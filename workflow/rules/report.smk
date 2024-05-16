report: base.parent / "README.rst"


rule tab_alignments:
    message:
        "tabulate alignments: {wildcards.id}"
    input:
        jsn=rules.call.output.jsn,
    output:
        tsv=report(subroot / "aln.tsv", category="Assay", subcategory="{id}"),
    params:
        script=base / "scripts" / "report" / "aln.jq",
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
        script=base / "scripts" / "report" / "hit.jq",
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
        script=base / "scripts" / "report" / "con.jq",
    shell:
        """
        {{
            printf "id\tTP\tFN\tFP\tTN\n";
            jq -r -f {params.script:q} {input.jsn:q};
        }} > {output.tsv:q}
        """


rule plot_report:
    message:
        "generate plot report"
    input:
        json=[Path(ele).absolute() for ele in expand(rules.call.output.jsn, id=assays)],
    output:
        html=root.absolute() / "report.html",
    log:
        log=root / "report.log",
    params:
        script=base.absolute() / "scripts" / "report" / "plot.Rmd",
    threads: 8
    shell:
        """
        Rscript \
            -e "args = commandArgs(trailingOnly = T); rmarkdown::render(args[1], output_file=args[2], params=list(json=args[3:length(args)]));" \
            {params.script:q} \
            {output:q} \
            {input.json} \
            2> {log.log}
        """


rule tab_config:
    output:
        tsv=report(root / "cfg.tsv", category="Meta"),
    params:
        cmd='(["key", "value"] | @tsv), (. | to_entries[] | [.key, .value] | @tsv)',
        config=json.dumps(config),
    threads: 1
    shell:
        """
        echo {params.config:q} | jq -r {params.cmd:q} > {output.tsv:q}
        """


rule target_tsv:
    input:
        expand(rules.tab_hits.output.tsv, id=assays),
        expand(rules.tab_alignments.output.tsv, id=assays),
        rules.tab_confusion.output.tsv,
        rules.tab_config.output.tsv,


rule target_report:
    input:
        expand(rules.tab_hits.output.tsv, id=assays),
        expand(rules.tab_alignments.output.tsv, id=assays),
        rules.tab_confusion.output.tsv,
        rules.tab_config.output.tsv,
        rules.plot_report.output.html,
