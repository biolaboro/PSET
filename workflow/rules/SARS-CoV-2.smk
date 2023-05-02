import json
import lzma
from pathlib import Path
from subprocess import PIPE, Popen
from getpass import getpass

# setup

path_provision = Path("resources") / "provision.json.xz"
path_blast = Path("resources") / "blast" / "SARS-CoV-2"
path_meta = Path("resources") / "meta" / "SARS-CoV-2.tsv"


# rules


rule download:
    message:
        "download GISAID sequence data"
    output:
        xz=path_provision,
    params:
        url="https://www.epicov.org/epi3/3p/resseq02/export/",
        fname=path_provision.name,
        uname=config["user"],
    shell:
        "curl -su {params.uname:q} {params.url:q}/{params.fname:q} > {output.xz:q}"


rule makeblastdb:
    message:
        "make BLAST+ database"
    input:
        xz=rules.download.output,
    output:
        directory(path_blast),
    params:
        title="SARS-CoV-2",
        out=path_blast / "SARS-CoV-2",
        taxid="2697049",
    log:
        log=path_blast / "SARS-CoV-2.log",
    run:
        cmd = (
            "makeblastdb",
            "-in",
            "-",
            "-input_type",
            "fasta",
            "-dbtype",
            "nucl",
            "-title",
            params.title,
            "-parse_seqids",
            "-hash_index",
            "-out",
            params.out,
            "-blastdb_version",
            "5",
            "-logfile",
            log.log,
            "-taxid",
            params.taxid,
        )
        with Popen(cmd, stdin=PIPE, universal_newlines=True) as pipe:
            with pipe.stdin as stdin:
                with lzma.open(input.xz[0], "rt") as file:
                    for record in map(json.loads, file):
                        if record["covv_host"] == "Human" and record["is_high_coverage"]:
                            print(">", record["covv_accession_id"], "\n", record["sequence"], "\n", sep="", file=stdin)


rule metadata:
    message:
        "xz.json -> tsv"
    input:
        xz=rules.download.output,
    output:
        tsv=path_meta,
    run:
        with lzma.open(input.xz[0], "rt") as file1, open(output.tsv, "w") as file2:
            record = json.loads(next(file1))
            fields = [key for key in record if key != "sequence"]
            print(*fields, sep="\t", file=file2)
            print(*(record[key] for key in fields), sep="\t", file=file2)
            for record in map(json.loads, file1):
                print(*(record[key] for key in fields), sep="\t", file=file2)


rule update:
    input:
        path_provision,
        path_blast,
        path_meta,
