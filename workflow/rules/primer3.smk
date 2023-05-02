# imports

from pathlib import Path

# globals

base = Path(workflow.basedir)
root = base.parent.parent
env = str(root / "workflow" / "envs" / "env.yml")
lab = config.get("lab", Path(config["fna"]).stem.split(".")[0])

# rules


rule primer3:
    input:
        fna=config["fna"],
        cfg=config["cfg"],
    output:
        jsn=root / "resources" / "assay" / f"{lab}.primer3.json",
    log:
        log=root / "resources" / "assay" / f"{lab}.primer3.log",
    threads: 1
    conda:
        env
    params:
        root=root / "workflow" / "scripts",
    shell:
        """
        {params.root:q}/run_primer3.py {input.fna:q} {input.cfg:q} -proc {threads} 1> {output:q} 2> {log.log:q}
        """


rule to_assay:
    input:
        jsn=rules.primer3.output.jsn,
    output:
        tsv=root / "resources" / "assay" / f"{lab}.assay.tsv",
    log:
        log=root / "resources" / "assay" / f"{lab}.assay.log",
    threads: 1
    conda:
        env
    params:
        root=root / "workflow" / "scripts",
        max=config.get("max", 999),
    shell:
        """
        {params.root:q}/primer3_to_definition.py {input.jsn:q} | \
            jq --arg max {params.max} -r -f {params.root:q}/primer3_to_assay.jq \
                1> {output.tsv:q} \
                2> {log.log:q}
        """
