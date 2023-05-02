from pathlib import Path

path = Path(config["path"])
base = Path(workflow.basedir).parent
out = config.get("out", base.parent / "results" / "conserved" / "msa")
lab = config.get("lab", path.stem)
root = Path(out)


rule mafft:
    input:
        config["path"],
    output:
        root / f"{lab}.msa.fna",
    log:
        root / f"{lab}.msa.log",
    threads: 8
    shell:
        """
        mafft --thread {threads} {input:q} 1> {output:q} 2> {log:q}
        """


rule conserved:
    input:
        rules.mafft.output,
    output:
        root / f"{lab}.con.fna",
    log:
        root / f"{lab}.con.log",
    threads: 1
    params:
        script=base / "scripts" / "msa_to_conserved.py",
        min=config.get("min", 100),
        thr=config.get("thr", 1),
    shell:
        """
        {params.script:q} -gapless -min {params.min:q} -thr {params.thr:q} {input:q} 1> {output:q} 2> {log:q}
        """


rule primer3:
    input:
        fna=rules.conserved.output,
        cfg=config.get("cfg", base.parent / "conf" / "primer3.json"),
    output:
        jsn=root / f"{lab}.primer3.json",
    log:
        log=root / f"{lab}.primer3.log",
    threads: 8
    params:
        script=base / "scripts",
    shell:
        """
        {params.script:q}/run_primer3.py {input.fna:q} {input.cfg:q} -proc {threads} 1> {output:q} 2> {log.log:q}
        """


rule to_assay:
    input:
        jsn=rules.primer3.output.jsn,
    output:
        tsv=root / f"{lab}.assay.tsv",
    log:
        log=root / f"{lab}.assay.log",
    threads: 1
    params:
        script=base / "scripts",
        max=config.get("max", 999),
        targets=config.get("targets", 1),
    shell:
        """
        {params.script:q}/primer3_to_definition.py {input.jsn:q} | \
            jq --arg max {params.max:q} --arg targets {params.targets:q} -r -f {params.script:q}/primer3_to_assay.jq \
                1> {output.tsv:q} \
                2> {log.log:q}
        """


rule all:
    input:
        rules.to_assay.output.tsv,
