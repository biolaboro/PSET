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
            log=root / "blast" / config["db"] / f"{config['db']}.log",
        params:
            db=config["db"],
            decompress=config.get("decompress", ""),
            source=config.get("source", "ncbi"),
            blastdb_version=5,
            out=root / "blast" / config["db"],
        threads: 4
        shell:
            """
            decompress={params.decompress:q}
            mkdir -p {params.out:q} && \
            cd {params.out:q} || exit 1 && \
            update_blastdb.pl \
                --source {params.source:q} \
                "${{decompress:+--decompress}}" \
                --blastdb_version {params.blastdb_version:q} \
                --num_threads {threads} \
                {params.db:q} | \
                tee {output.log:q}
            # in case complete decompress fails for some reason or the flag wasn't added
            # this is known when the tar.gz files remain
            find . -name "*.tar.gz" -type f | xargs -L 1 -P {threads} -- tar --keep-newer-files -x -f
            find . -name "*.tar.gz" -type f -exec rm {{}} \;
            """
