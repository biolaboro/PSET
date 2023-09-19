rule assay:
    message:
        "output assay data to JSON"
    input:
        tsv=config["file"],
    output:
        jsn=expand(subroot / "assay.json", id=assays),
    params:
        flank=config["flank"],
        context=context,
        parent=root,
    threads: 1
    run:
        with open(input.tsv) as tsv:
            for assay in parse_assays(tsv, context=params.context, target_type=int):
                keys = assay.kflanks() if config["flank"] else ("amplicon",)
                camp = assay.camplicon()
                contexts = {}
                # calculate the actual amount of 5'/3'-context available per component query
                for key in keys:
                    coor = assay.acoors() if key == "amplicon" else assay.ccoors[key]
                    contexts[key] = (min(coor[0], params.context[0]), min(len(camp) - coor[1], params.context[1]))
                lcl = {key: str(next(assay.records(key, expand=1, context=params.context)).seq) for key in keys}
                qry = dict(lcl=lcl, glc=assay.components, ctx=contexts)
                # setup query for local and glocal alignment
                data = dict(zip(assay.key(), assay.val()), type=type(assay).__name__, coor=assay.ccoors, qry=qry)
                with params.parent.joinpath(assay.id).joinpath("assay.json").open("w") as jsn:
                    json.dump(data, jsn, indent=4)


rule local:
    message:
        "BLAST+: {wildcards.id}"
    input:
        jsn=subroot / "assay.json",
    output:
        tsv=subroot / "lcl.tsv",
    params:
        db=db,
        outfmt=f"7 {' '.join(fields_8CB)} staxids",
        flags=list(chain.from_iterable(confb)),
    log:
        log=subroot / "lcl.log",
    threads: 8
    shell:
        """
        {{
            jq -r '.qry.lcl | to_entries | map(">\(.key)\n\(.value)") | .[]' {input.jsn:q} | \
            blastn \
                -query - \
                -db {params.db:q} \
                -out {output.tsv:q} \
                -outfmt {params.outfmt:q} \
                -num_threads {threads} \
                {params.flags}
        }} 2> {log.log:q}
        """


rule filter:
    message:
        "filter by similarity scores: {wildcards.id}"
    input:
        jsn=subroot / "assay.json",
        tsv=rules.local.output.tsv,
    output:
        tsv=subroot / "flt.sim.tsv",
        txt=subroot / "flt.bat.txt",
    log:
        log=subroot / "flt.log",
    params:
        script=base / "scripts" / "pset" / "filter.py",
        sim=config["simlcl"],
    threads: 1
    shell:
        """
        python3 {params.script:q} \
            {input.jsn:q} \
            {input.tsv:q} \
            -tsv {output.tsv:q} \
            -txt {output.txt:q} \
            -sim {params.sim:q} \
            2> {log.log:q}
        """


rule library:
    message:
        "extract aligned subjects: {wildcards.id}"
    input:
        txt=rules.filter.output.txt,
    output:
        fna=subroot / "lib.seq.fasta",
        jsn=subroot / "lib.map.json",
        tsv=subroot / "lib.tax.tsv",
    params:
        db=db,
    threads: 1
    run:
        recs = {}
        maps = defaultdict(list)
        # output taxonomy mappings
        cmd = ("blastdbcmd", "-db", params.db, "-entry_batch", input.txt, "-outfmt", "%a\t%T", "-out", output.tsv)
        check_call(cmd)
        # output extracted sequences
        cmd = ("blastdbcmd", "-db", params.db, "-entry_batch", input.txt, "-outfmt", ">%a\n%s")
        with Popen(cmd, stdout=PIPE, universal_newlines=True) as proc:
            for rec in SeqIO.parse(proc.stdout, "fasta"):
                key = seguid(rec.seq)
                recs[key] = recs.get(key, rec)
                maps[key].append(rec.id)
        SeqIO.write(recs.values(), output.fna, "fasta")
        # output duplicate mappings
        with open(output.jsn, "w") as file:
            json.dump(maps, file, indent=4)
            print(file=file)


rule glocal:
    message:
        "glocal alignment: {wildcards.id}"
    input:
        jsn=subroot / "assay.json",
        lib=rules.library.output.fna,
    output:
        tsv=subroot / "glc.tsv",
    log:
        log=subroot / "glc.log",
    params:
        *chain.from_iterable(confg),
    threads: 8
    shell:
        """
        {{
            jq -r '.qry.glc | to_entries | map(">\(.key)\n\(.value)") | .[]' {input.jsn:q} | \
                glsearch36 -n -m 8CB -T {threads} {params:q} - {input.lib:q} | \
                awk -F $"\t" -v OFS=$"\t" '/^#/ {{ print; }} /^[^#]/ {{ $2=$2":"n++; print; }}'
        }} 1> {output.tsv:q} 2> {log.log:q}
        """


rule call:
    message:
        "output the assay calls: {wildcards.id}"
    input:
        assay=subroot / "assay.json",
        mapping=rules.library.output.jsn,
        taxa=rules.library.output.tsv,
        alignment=rules.glocal.output.tsv,
    output:
        jsn=subroot / "call.json",
    log:
        log=subroot / "call.log",
    params:
        **dkwargs,
        script=base / "scripts" / "pset" / "call.py",
        dburl=config["dburl"],
        sim=config["simglc"],
        xtaxa=config["xtaxa"].split(","),
    threads: 1
    shell:
        """
        {params.script:q} \
            -assay {input.assay:q} \
            -alignment {input.alignment:q} \
            -mapping {input.mapping:q} \
            -taxa {input.taxa:q} \
            -dburl {params.dburl:q} \
            -sim {params.sim:q} \
            -dFR {params.dFR:q} \
            -dF3F2 {params.dF3F2:q} \
            -dF2F1c {params.dF2F1c:q} \
            -dF1cB1c {params.dF1cB1c:q} \
            -xtaxa {params.xtaxa} \
            -out {output.jsn:q} \
            2> {log.log:q}
        """


rule target_pset:
    message:
        "target pset"
    input:
        expand(rules.call.output.jsn, id=assays),
