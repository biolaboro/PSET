# imports

import json

# globals

root = Path(workflow.basedir).parent.parent
env = str(root / "workflow" / "envs" / "env.yml")

# rules


def listify(path):
    with open(path) as file:
        yield from map(str.strip, file)


def db_files(path, db):
    with open(path) as file:
        obj = json.load(file)
    obj = next(ele for ele in obj if ele["dbname"] == db)
    yield from (Path(ele).name for ele in obj["files"])


rule taxa:
    output:
        db=root / "resources" / "taxa.db",
    log:
        log=root / "resources" / "log_taxa.txt",
    params:
        out=root / "resources",
        url="ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy",
        tar="taxdump.tar.gz",
    threads: 1
    conda:
        env
    shell:
        """
        curl --output-dir {params.out} -O "{params.url}/{params.tar}" 2> {log.log:q} && \
        python3 -m taxa.taxa -database {output.db} create "{params.out}/{params.tar}" 2>> {log.log:q}
        """


rule data:
    output:
        expand(root / "resources" / "blast" / "EBOV" / "EBOV.{ext}", ext=("nhd", "nhi", "nhr", "nin", "nog", "nos", "not", "nsq", "ntf", "nto")),
        ssv=root / "resources" / "blast" / "EBOV" / "EBOV.ssv",
        log=root / "resources" / "blast" / "EBOV" / "EBOV.log",
    log:
        log=root / "resources" / "log_data.txt",
    params:
        out=root / "resources" / "blast" / "EBOV" / "EBOV",
        db="nuccore",
        query="txid11266[Organism:exp] AND biomol_genomic[PROP] NOT gbdiv_pat[PROP] NOT gbdiv_syn[PROP]",
        n=500,
    threads: 1
    conda:
        env
    shell:
        """
        esearch -db {params.db:q} -query {params.query:q} | \
        efetch -format uid | \
        head -n {params.n} | \
        esummary -db nuccore | \
        xtract -pattern DocumentSummary -element AccessionVersion TaxId > {output.ssv:q}
        cut -f 1 -d " " {output.ssv:q} | \
        epost -db {params.db:q} | \
        efetch -format fasta | \
        makeblastdb \
            -parse_seqids \
            -hash_index \
            -blastdb_version 5 \
            -dbtype nucl \
            -title "$(basename {params.out:q})" \
            -out {params.out:q} \
            -logfile {output.log} \
            -taxid_map {output.ssv:q}
        """


rule json:
    output:
        jsn=root / "resources" / "blast" / "db.json",
    params:
        url="ftp://ftp.ncbi.nlm.nih.gov/blast/db/v5/",
    shell:
        """
        curl -s --list-only {params.url:q} | \
        grep -E "\-metadata\.json$" | \
        xargs -I % curl https://ftp.ncbi.nlm.nih.gov/blast/db/v5/% | \
        jq -s . > {output.jsn:q}
        """


if config.get("db"):

    rule download:
        input:
            jsn=rules.json.output.jsn,
        output:
            txt=root / "resources" / "blast" / config["db"] / "lst.txt",
        params:
            db=config["db"],
            out=root / "resources" / "blast" / config["db"],
        threads: 4
        shell:
            """
            for ext in "" ".md5"; do
                jq --arg db {params.db:q} --arg ext "$ext" 'map(select(.dbname == $db))[].files[] | "\(.)\($ext)"' {input.jsn:q} | \
                xargs -n 1 -P {threads} curl --output-dir {params.out} -O
            done
            find {params.out:q} -name "{params.db}*.tar.gz" | \
            while read -r ele; do
                x="$(openssl dgst -md5 "$ele" | cut -f 2 -d " ")"
                y="$(cat "$ele.md5" | cut -f 1 -d " ")"
                [[ "$x" == "$y" ]] || exit 1
                python -m tarfile --list "$ele"
                python -m tarfile --extract "$ele" {params.out:q}
            done | sort -u > {output.txt:q}
            """
