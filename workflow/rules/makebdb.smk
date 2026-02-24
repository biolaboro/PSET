from pathlib import Path


root = Path(config["root"])
title = config["title"]


rule makeblastdb:
    message:
        """Create BLAST+ database with the FASTA/GenBank files and taxonomy info..."""
    input:
        seq=config["path"],
    output:
        log=root / f"{title}.log",
    log:
        log=root / f"{title}.err.log",
    params:
        out=root / title,
        taxid=config.get("taxid", 1),
        taxid_map=config.get("taxid_map", ""),
        title=title,
    threads: 1
    shell:
        """
        {{
            tax1={params.taxid_map:q}
            tax2={params.taxid:q}
            if [ -z "$tax1" ]; then tax="-taxid $tax2"; else tax="-taxid_map $tax1"; fi
            makeblastdb \
                -in {input.seq:q} \
                -input_type fasta \
                -dbtype nucl \
                -title {params.title:q} \
                -parse_seqids \
                -hash_index \
                -out {params.out} \
                -blastdb_version 5 \
                -logfile {output.log:q} \
                $tax
            rm -f {input.seq:q}
        }} 2> {log.log:q}
        """
