rule assay:
    message:
        "output assay data to JSON: {wildcards.id}"
    output:
        jsn=subroot / "assay.json",
    params:
        flank=config["flank"],
        context=context,
    threads: 1
    run:
        assay = assays[wildcards.id]
        keys = assay.flanks if config["flank"] else ("amplicon",)
        camp = assay.camplicon()
        contexts = OrderedDict()
        # calculate the actual amount of 5'/3'-context available per component query
        for key in keys:
            coor = assay.coordinates[key]
            contexts[key] = (min(coor[0], params.context[0]), min(len(camp) - coor[1], params.context[1]))
        lcl = {key: str(next(assay.records(key, expand=1, context=params.context)).seq) for key in keys}
        qry = dict(lcl=lcl, glc=assay.components, ctx=contexts)
        # setup query for local and glocal alignment
        data = dict(zip(assay.key(), assay.val()), coor=assay.coordinates, qry=qry, config=config)
        with open(output.jsn, "w") as file:
            json.dump(data, file, indent=4)


rule local:
    message:
        "BLAST+: {wildcards.id}"
    input:
        jsn=rules.assay.output.jsn,
    output:
        tsv=subroot / "lcl.tsv",
    params:
        db=db,
        outfmt=f"7 {' '.join(fields_8CB)} staxids",
        flags=list(chain.from_iterable(argsb)),
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
        jsn=rules.assay.output.jsn,
        tsv=rules.local.output.tsv,
    output:
        tsv=subroot / "flt.sim.tsv",
        txt=subroot / "flt.bat.txt",
    log:
        log=subroot / "flt.log",
    params:
        script=base / "scripts" / "filter.py",
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
        recs = OrderedDict()
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
        jsn=rules.assay.output.jsn,
        lib=rules.library.output.fna,
    output:
        tsv=subroot / "glc.tsv",
    log:
        log=subroot / "glc.log",
    params:
        *chain.from_iterable(argsg),
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
        assay=rules.assay.output.jsn,
        mapping=rules.library.output.jsn,
        taxa=rules.library.output.tsv,
        alignment=rules.glocal.output.tsv,
    output:
        jsn=subroot / "call.json",
    log:
        log=subroot / "call.log",
    params:
        script=base / "scripts" / "call.py",
        dburl=config["dburl"],
        sim=config["simglc"],
        dist=config["dist"],
        xtaxa=config["xtaxa"].replace(",", " "),
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
            -dist {params.dist:q} \
            -xtaxa {params.xtaxa} \
            -out {output.jsn:q} \
            2> {log.log:q}
        """


rule target_assay:
    message:
        "target assay"
    input:
        expand(rules.assay.output.jsn, id=assays),


rule target_pset:
    message:
        "target pset"
    input:
        expand(rules.call.output.jsn, id=assays),
