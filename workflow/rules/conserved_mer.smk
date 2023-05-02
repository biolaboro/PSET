from pathlib import Path

base = Path(workflow.basedir).parent
path = Path(config["path"])
root = Path(config.get("out", "results")) / config.get("lab", path.stem)


rule mash:
    message:
        "calculate mash dist"
    input:
        fna=path,
    output:
        tsv=root / "dst.tsv",
    log:
        log=root / "dst.log",
    threads: 8
    shell:
        """
        mash dist -p {threads} -C -i {input.fna:q} {input.fna:q} > {output.tsv:q} 2> {log.log:q}
        """


checkpoint cluster:
    message:
        "cluster distance matrix"
    input:
        tsv=rules.mash.output.tsv,
    output:
        directory(root / "cls"),
    log:
        log=root / "cls.log",
    params:
        root=base / "scripts",
    threads: 8
    shell:
        """
        awk -F $'\t' -v OFS=$'\t' '{{ sub(/:.*/, "", $1); sub(/:.*/, "", $2); print($1, $2, $3); }}' {input.tsv:q} | \
        {params.root:q}/dstcls.py -lab {output:q} -proc {threads} - 2> {log.log:q}
        """


rule library:
    input:
        fna=path,
        tsv=root / "cls" / "{cls}.tsv",
    output:
        fna=root / "cls" / "{cls}.lib.fna",
    log:
        log=root / "cls" / "{cls}.lib.log",
    shell:
        """
        awk 'NR > 2 {{ print $2 }}' {input.tsv:q} | xargs -L 8 samtools faidx {input.fna:q} > {output.fna:q} 2> {log.log:q}
        """


rule reference:
    input:
        fna=path,
        tsv=root / "cls" / "{cls}.tsv",
    output:
        fna=root / "cls" / "{cls}.ref.fna",
    log:
        log=root / "cls" / "{cls}.ref.log",
    shell:
        """
        awk 'NR == 2 {{ print $2 }}' {input.tsv:q} | xargs -L 1 samtools faidx {input.fna:q} > {output.fna:q} 2> {log.log:q}
        """


rule nucmer:
    message:
        """
        run nucmer
        """
    input:
        ref=rules.reference.output.fna,
        lib=rules.library.output.fna,
    output:
        dlt=root / "cls" / "{cls}.mer.delta",
    log:
        log=root / "cls" / "{cls}.mer.log",
    params:
        cls=lambda wildcards: root / "cls" / f"{wildcards.get('cls')}.mer",
        minmatch=config.get("minmatch", 20),
    threads: 8
    shell:
        """
        nucmer --maxmatch --minmatch {params.minmatch:q} --prefix {params.cls:q} --threads {threads} {input.ref:q} {input.lib:q}
        """


rule delta:
    input:
        dlt=rules.nucmer.output.dlt,
    output:
        dlt=root / "cls" / "{cls}.mer.100.delta",
    log:
        log=root / "cls" / "{cls}.mer.100.log",
    params:
        i=config.get("i", 100),
    threads: 1
    shell:
        """
        delta-filter -i {params.i:q} {input.dlt:q} > {output.dlt:q} 2> {log.log:q}
        """


rule coords:
    input:
        dlt=rules.delta.output.dlt,
    output:
        tsv=root / "cls" / "{cls}.rng.tsv",
    log:
        log=root / "cls" / "{cls}.rng.log",
    threads: 1
    shell:
        """
        show-coords -T {input.dlt:q} | awk -v OFS=$"\t" 'NR > 4 {{ print $1, $2, $9; }}' > {output.tsv:q} 2> {log.log:q}
        """


rule conserved:
    input:
        ref=rules.reference.output.fna,
        lib=rules.library.output.fna,
        tsv=rules.coords.output.tsv,
    output:
        jsn=root / "cls" / "{cls}.json",
    log:
        log=root / "cls" / "{cls}.log",
    params:
        script=base / "scripts" / "mum_to_conserved.py",
    threads: 8
    shell:
        """
        {params.script:q} {input.ref:q} {input.lib:q} {input.tsv:q} -proc {threads} > {output.jsn:q} 2> {log.log:q}
        """


rule tabulate:
    input:
        jsn=rules.conserved.output.jsn,
    output:
        tsv=root / "cls" / "{cls}.con.tsv",
    log:
        log=root / "cls" / "{cls}.con.log",
    params:
        script=base / "scripts" / "cluster.jq",
    shell:
        """
        jq -r -f {params.script:q} {input.jsn} > {output.tsv:q} 2> {log.log:q}
        """


def get_clusters(wildcards):
    pattern = Path(checkpoints.cluster.get(**wildcards).output[0]) / "{cls}.tsv"
    return glob_wildcards(pattern).cls


rule listing:
    input:
        tsv=lambda wildcards: expand(rules.tabulate.output.tsv, cls=get_clusters(wildcards)),
    output:
        txt=root / "lst.txt",
    log:
        log=root / "lst.log",
    shell:
        """
        for ele in {input.tsv:q}; do echo "$ele"; done > {output.txt:q} 2> {log.log:q}
        """


rule all:
    input:
        rules.listing.output.txt,
