import os
import shlex
from pathlib import Path

from Bio import SeqIO

from pset.assay import parse_assays
from pset.util import argify


base = Path(workflow.basedir).parent
resources = Path("resources")
results = Path("results")

lab = config.get("lab", Path(config.get("file", "")).stem)

confb = [ele.split("=") for ele in shlex.split(config.get("confb", ""))]

url = "https://iedb.org/downloader.php?file_name=doc/epitope_full_v3_json.zip"


rule download:
    message:
        "download the IEDB stuff"
    output:
        zip=resources / "iedb" / Path(url).name,
    params:
        dir=lambda wildcards, output: Path(output.zip).parent,
        name=Path(url).name,
        url=url,
    threads: 1
    shell:
        """
        curl --create-dirs --output-dir {params.dir:q} -o {params.name:q} {params.url:q}
        """


rule tableize:
    input:
        zip=rules.download.output.zip,
    output:
        tsv=resources / "iedb" / "iedb.tsv",
    params:
        script=base / "scripts" / "eset" / "tableize.jq",
    threads: 1
    shell:
        """
        unzip -p {input.zip:q} | jq -r -f {params.script:q} > {output.tsv:q}
        """


rule blastdb:
    input:
        tsv=rules.tableize.output.tsv,
    output:
        multiext(
            str(resources / "blast" / "iedb" / "iedb"),
            ".pdb",
            ".phd",
            ".phi",
            ".phr",
            ".pin",
            ".pog",
            ".pos",
            ".pot",
            ".psq",
            ".ptf",
            ".pto",
        ),
        txt=resources / "blast" / "iedb" / "iedb.txt",
    log:
        log=resources / "blast" / "iedb" / "iedb.log",
    params:
        out=lambda wildcards, output: Path(output.txt).parent / "iedb",
    threads: 1
    shell:
        """
        awk -F $'\t' 'NR > 1 {{ print($1, $NF); }}' {input.tsv:q} > {output.txt:q} && \
        awk -F $'\t' 'NR > 1 {{ printf(">%d [molecule=%s] [organism=%s]\\n%s\\n", $1, $2, $3, $4); }}' {input.tsv:q} | \
            makeblastdb \
                -in - \
                -input_type fasta \
                -dbtype prot \
                -title iedb \
                -parse_seqids \
                -hash_index \
                -out {params.out:q} \
                -blastdb_version 5 \
                -logfile {log.log:q} \
                -taxid_map {output.txt:q}
        """


rule query:
    input:
        tsv=config.get("file", ""),
    output:
        fna=results / lab / "query.fna",
    run:
        os.makedirs(Path(output.fna).parent, exist_ok=True)
        with open(input.tsv) as file:
            records = []
            for assay in parse_assays(file):
                for record in assay.records():
                    record.id = assay.id + "-" + record.id
                    records.append(record)
            SeqIO.write(records, output.fna, "fasta")


rule blastx:
    input:
        fna=rules.query.output.fna,
    output:
        tsv=results / lab / "blastx.tsv",
    log:
        log=results / lab / "blastx.log",
    threads: 8
    params:
        db=config.get("db", resources / "blast" / "iedb" / "iedb"),
        outfmt="7 std slen staxids stitle btop",
        flags=list(chain.from_iterable(confb)),
    shell:
        """
            blastx \
                -query {input.fna:q} \
                -db {params.db:q} \
                -outfmt {params.outfmt:q} \
                -num_threads {threads} \
                -out {output.tsv:q} \
                {params.flags} \
                2> {log.log:q}
        """
