report: base.parent / "report" / "workflow.rst"


rule tab_alignments:
    message:
        "tabulate alignments: {wildcards.id}"
    input:
        jsn=rules.call.output.jsn,
    output:
        tsv=report(subroot / subroot / "aln.tsv", category="Assay", subcategory="{id}"),
    params:
        script=base / "scripts" / "aln.jq",
    threads: 1
    shell:
        """
        jq -r -f {params.script:q} {input.jsn:q} | tr '|' '\n' > {output.tsv:q}
        """


rule tab_hits:
    message:
        "tabulate hits: {wildcards.id}"
    input:
        jsn=rules.call.output.jsn,
    output:
        tsv=report(subroot / "hit.tsv", category="Assay", subcategory="{id}"),
    params:
        script=base / "scripts" / "hit.jq",
    threads: 1
    shell:
        """
        jq -r -f {params.script:q} {input.jsn:q} > {output.tsv:q}
        """


rule tab_confusion:
    message:
        "tabulate confusion matrix"
    input:
        jsn=expand(rules.call.output.jsn, id=assays),
    output:
        tsv=report(root / "con.tsv", category="Summary"),
    params:
        script=base / "scripts" / "con.jq",
    shell:
        """
        {{
            printf "db\tid\tTP\tFN\tFP\tTN\n";
            jq -r -f {params.script:q} {input.jsn:q};
        }} > {output.tsv:q}
        """


rule plot_mismatches:
    message:
        "plot assay component mismatches: {wildcards.id}"
    input:
        rules.call.output.jsn,
    output:
        png=report(subroot / "mis.png", caption=base.parent / "report" / "mismatch.rst", category="Assay", subcategory="{id}"),
    log:
        log=subroot / "mis.log",
    params:
        script=base / "scripts" / "plot_mismatches.py",
        dpi=config["dpi"],
    threads: 1
    shell:
        """
        {params.script:q} {input:q} -out {output.png:q} -dpi {params.dpi} 2> {log.log:q}
        """


rule plot_heat:
    message:
        "plot heatmap"
    input:
        expand(rules.call.output.jsn, id=assays),
    output:
        png=report(root / "heat.png", caption=base.parent / "report" / "heat.rst", category="Summary"),
    log:
        log=root / "mis.log",
    params:
        parents=lambda wildcards, input: [Path(ele).parent for ele in input],
        script=base / "scripts" / "plot_heat.py",
        dpi=config["dpi"],
    threads: 1
    shell:
        """
        {params.script:q} {params.parents:q} -out {output.png:q} -dpi {params.dpi} 2> {log.log:q}
        """


rule tab_config:
    input:
        jsn=expand(rules.assay.output.jsn, id=assays).pop(),
    output:
        tsv=report(root / "cfg.tsv", category="Meta"),
    params:
        cmd='(["key", "value"] | @tsv), (.config | to_entries[] | [.key, .value] | @tsv)',
    threads: 1
    shell:
        """
        jq -r {params.cmd:q} {input.jsn} > {output.tsv:q}
        """


rule tab_targets:
    input:
        jsn=expand(rules.assay.output.jsn, id=assays),
    output:
        tsv=report(root / "tar.tsv", category="Meta"),
    params:
        db=config["db"],
        taxdb=config["dburl"][len("sqlite:///") :],
    shell:
        """
        touch {output.tsv}
        {{
            jq '.targets[]' {input.jsn:q} | sort -u | \
                blastdbcmd -db {params.db:q} -taxidlist - -outfmt '%T' | sort -n | uniq -c | \
            awk -v OFS=$'\t' '{{ print($2, $1); }}' > {output.tsv}.tmp1
            cut -f 1 {output.tsv}.tmp1 | \
                xargs python -m taxa.taxa -database {params.taxdb:q} lineage | \
                awk '$1==$2 {{ print($3, $6); }}' > {output.tsv}.tmp2
            paste {output.tsv}.tmp1 {output.tsv}.tmp2 > {output.tsv}
            rm -f {output.tsv}.tmp[12]
        }} || true
        """


rule txt_blastdbinfo:
    output:
        txt=report(root / "bdb.tsv", category="Meta"),
    params:
        db=config["db"],
    shell:
        """
        blastdbcmd -db {params.db:q} -info > {output.txt}
        """


rule target_tsv:
    input:
        expand(rules.tab_hits.output.tsv, id=assays),
        expand(rules.tab_alignments.output.tsv, id=assays),
        rules.tab_confusion.output.tsv,
        rules.tab_config.output.tsv,
        rules.tab_targets.output.tsv,


rule target_report:
    input:
        expand(rules.tab_hits.output.tsv, id=assays),
        expand(rules.tab_alignments.output.tsv, id=assays),
        rules.tab_confusion.output.tsv,
        expand(rules.plot_mismatches.output.png, id=assays),
        rules.plot_heat.output.png,
        rules.tab_config.output.tsv,
        rules.tab_targets.output.tsv,
        rules.txt_blastdbinfo.output.txt,
