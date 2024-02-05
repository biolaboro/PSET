import os
from pathlib import Path
from Bio import SeqIO

base = Path(workflow.basedir)
root = base.parent
path_out = Path(config.get("out", Path("resources") / "assay"))
os.makedirs(path_out, exist_ok=True)

path_conf = root / "conf" / "agen"

path_file = config["file"]

modes = config["mode"].split(",")
confs = {ele.stem: ele for ele in path_conf.glob("*.json")}
accs = [record.id for record in SeqIO.parse(path_file, "fasta")]


rule split:
    input:
        fna=path_file,
    output:
        fna=path_out / "{acc}" / "seq.fna",
    run:
        for record in SeqIO.parse(input.fna, "fasta"):
            SeqIO.write(record, output.fna, "fasta")


rule agen:
    message:
        "generate assay candidates"
    input:
        fna=path_out / "{acc}" / "seq.fna",
    output:
        jsn=path_out / "{acc}" / "{mode}.{conf}.json",
    log:
        log=path_out / "{acc}" / "{mode}.{conf}.log",
    params:
        script=root / "scripts" / "agen" / "agen.py",
        conf=lambda wildcards, output: confs[wildcards.conf],
        cstr=config.get("cstr", "GLOBAL:PRIMER_NUM_RETURN=1000"),
        limit=int(config.get("n_per_region", 1)),
        optional_loop=config.get("optional_loop", ""),
    threads: 8
    shell:
        """
        optional_loop={params.optional_loop:q}
        {params.script:q} \
            -mode {wildcards.mode:q} \
            -conf {params.conf:q} \
            ${{optional_loop:+--optional-loop}} \
            -cstr {params.cstr} \
            -limit {params.limit:q} \
            -proc {threads} \
            -out {output.jsn:q} \
            {input.fna:q} \
            2> {log.log:q}
        """


rule rank:
    input:
        jsn=expand(path_out / "{{acc}}" / "{{mode}}.{conf}.json", conf=confs),
    output:
        jsn=path_out / "{acc}" / "{mode}.json",
    log:
        log=path_out / "{acc}" / "{mode}.log",
    params:
        script=root / "scripts" / "agen" / "rank.jq",
        targets=config.get("targets", "1"),
        n_per_region=int(config.get("n_per_region", 1)),
        n_total=int(config.get("n_total", 10)),
    wildcard_constraints:
        mode="(LAMP|PCR)",
    shell:
        """
        jq \
            -r \
            -s \
            -f {params.script:q} \
            --arg targets {params.targets:q} \
            --argjson n_per_region {params.n_per_region:q} \
            --argjson n_total {params.n_total:q} \
            {input.jsn} \
            1> {output.jsn} \
            2> {log.log:q}
        """


rule tableize:
    input:
        jsn=rules.rank.output.jsn,
    output:
        tsv=path_out / "{acc}" / "{mode}.tab",
    params:
        script=root / "scripts" / "agen" / "tableize.jq",
        targets=config.get("targets", "1"),
    shell:
        """
        jq -r -f {params.script:q} --arg targets {params.targets:q} {input.jsn} > {output.tsv:q}
        """


rule assayfy:
    input:
        jsn=rules.rank.output.jsn,
    output:
        tsv=path_out / "{acc}" / "{mode}.tsv",
    params:
        script=root / "scripts" / "agen" / "assayfy.jq",
        targets=config.get("targets", "1"),
    shell:
        """
        jq -r -f {params.script:q} --arg targets {params.targets:q} {input.jsn} > {output.tsv:q}
        """


rule target:
    input:
        expand(rules.tableize.output.tsv, acc=accs, conf=confs, mode=modes),
        expand(rules.assayfy.output.tsv, acc=accs, conf=confs, mode=modes),
