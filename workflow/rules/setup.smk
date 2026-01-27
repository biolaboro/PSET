# globals

root = Path(workflow.basedir).parent.parent / "resources"

# rules


rule taxa:
    message:
        """download NCBI Taxonomy files and build SQLite database..."""
    output:
        db=root / "taxa" / "taxa.db",
    log:
        log=root / "taxa" / "taxa.log",
    params:
        out=root / "taxa",
        url="ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy",
        tar="taxdump.tar.gz",
    threads: 1
    shell:
        """
        curl --output-dir {params.out} -O "{params.url}/{params.tar}" 2> {log.log:q} && \
        python3 -m taxa.taxa -database {output.db} create "{params.out}/{params.tar}" 2>> {log.log:q}
        """


rule ebov:
    message:
        """download EBOV nucleotide sequences and build BLAST+ database"""
    output:
        expand(
            root / "blast" / "EBOV" / "EBOV.{ext}",
            ext=("nhd", "nhi", "nhr", "nin", "nog", "nos", "not", "nsq", "ntf", "nto"),
        ),
    log:
        log=root / "blast" / "EBOV" / "EBOV.log",
    params:
        acc=root / "test" / "EBOV.ssv",
        out=root / "blast" / "EBOV" / "EBOV",
        db="nuccore",
    threads: 1
    shell:
        """
        cut -f 1 -d ' ' {params.acc:q} | \
        efetch -format fasta -db {params.db:q} | \
        makeblastdb \
            -parse_seqids \
            -hash_index \
            -blastdb_version 5 \
            -dbtype nucl \
            -title "$(basename {params.out:q})" \
            -out {params.out:q} \
            -logfile {log.log:q} \
            -taxid_map {params.acc:q}
        """


if config.get("db"):

    rule download:
        message:
            """download BLAST+ database from NCBI..."""
        output:
            log=root / "blast" / config["db"] / f"{config['db']}.stdout.log",
        log:
            log=root / "blast" / config["db"] / f"{config['db']}.stderr.log",
        params:
            db=config["db"],
            out=root / "blast" / config["db"],
            script=Path(workflow.basedir).parent / "scripts" / "util" / "dbdl.sh",
        threads: 4
        shell:
            """
            bash -x {params.script:q} {params.db:q} {params.out:q} {threads:q} > {output.log:q} 2> {log.log:q}
            """
