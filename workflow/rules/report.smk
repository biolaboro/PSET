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
        script=base / "scripts" / "report" / "plot_mismatches.py",
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
        log=root / "heat.log",
    params:
        parents=lambda wildcards, input: [Path(ele).parent for ele in input],
        script=base / "scripts" / "report" / "plot_heat.py",
        dpi=config["dpi"],
    threads: 1
    shell:
        """
        {params.script:q} {params.parents:q} -out {output.png:q} -dpi {params.dpi} 2> {log.log:q}
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


rule report_offline:
    input:
        tsv=rules.tab_confusion.output.tsv,
        png=Path(rules.plot_heat.output.png).resolve(),
    output:
        html=root / "report_offline.html",
    threads: 1
    params:
        md=base.parent / "README.md",
        cmd='{printf("|"); for (i = 1; i <= NF; ++i) printf("%s|", $i); printf("\\n")}',
        css="<style>table, th, td { border: 1px solid; } </style>",
    shell:
        """
        {{
            cat {params.md:q};
            printf "\n## Confusion Matrix\n";
            NF="$(awk -F $'\t' 'NR == 1 {{ printf(NF); }}' {input.tsv:q})";
            awk -F $'\t' 'NR == 1 {params.cmd}' {input.tsv:q};
            printf "|"; seq 1 "$NF" | while read -r _; do printf -- "-|"; done; printf "\\n";
            awk -F $'\t' 'NR >= 2 {params.cmd}' {input.tsv:q};
            printf "\n## Heatmap\n";
            printf "\n![Heatmap](%s)\n" {input.png:q};
            printf "\n%s\n" {params.css:q}
        }} | pandoc --metadata title=PSET --embed-resources --standalone - -o {output.html:q}
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
        expand(rules.plot_mismatches.output.png, id=assays),
        rules.plot_heat.output.png,
        rules.tab_config.output.tsv,


rule target_report_offline:
    input:
        md=rules.report_offline.output.html,
