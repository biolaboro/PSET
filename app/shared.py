import os
import json
from math import ceil
from string import ascii_uppercase
from itertools import product
from pathlib import Path
from collections import Counter, defaultdict
from subprocess import PIPE, Popen, CalledProcessError

import pandas as pd
import networkx as nx
from taxa.taxa import ancestors, descendants
from shiny import ui

CALLS = ("TP", "TN", "FP", "FN", "XX")
COMPONENTS = ("F", "P", "R", "F3", "F2", "LF", "F1c", "B1c", "LB", "B2", "B3")
BLAST_DIR = Path("resources") / "blast"
CPU_COUNT = os.cpu_count()
DBURL = "sqlite:///resources/taxa/taxa.db"
PATH_RESULTS = Path("results")
CONFUSION_AGGREGATE_KEYS = ("batch", "db", "tax", "sci", "rank", "genus", "species", "count")
DPI_APP = 96
DPI_PLOT = 300
USER_KEY = "user"  # this holds the header key to retrieve the user id
USER_DEFAULT = "user"  # this is the default user id value if not in the header


def add_key(df, col="id", key_name="key", key_values=ascii_uppercase):
    """Calculate a unique key column.

    Args:
        df (pandas.DataFrame): the data frame
        col (str): the column name to generate a key off of
        key_name (str): the name of the generated key
        key_values (list): the list of strings to derive key value from 

    Returns:
        pandas.DataFrame: the data frame with the new key column
    """
    values = df[col].unique()
    keys = dict(
        zip(
            values,
            map("".join, product(key_values, repeat=ceil(len(values) / len(key_values))))
        )
    )
    df[key_name] = [keys[ele] for ele in df[col]]
    df.insert(0, key_name, df.pop(key_name))
    return df


def read_func(path, read_func, *args, **kwargs):
    """Wrapper function for pandas read_* methods...

    Args:
        path (str): the file path
        read_func (function): the read_* function
        args: forwarded to the read_* function
        kwargs: forwarded to the read_* function

    Returns:
        pandas.DataFrame: the data frame
    """
    try:
        return read_func(path, *args, **kwargs)
    except pd.errors.EmptyDataError:
        return pd.DataFrame()


def read_hits(path):
    """Read hits from PSET call.json file as dictionaries.

    Args:
        path (str): the file

    Yiels:
        dict: the next hit result
    """
    keys = ("id", "com", "psim", "astr", "call", "acc")
    with path.open() as file:
        obj = json.load(file)
        for nth, hit in enumerate(obj["hits"], start=1):
            for com, val in hit["evals"].items():
                for call, accs in hit["calls"].items():
                    for acc in accs:
                        yield dict(
                            zip(keys, (obj["assay"]["id"], com, val["psim"], val["astr"], call, acc)),
                            nth=nth,
                            batch=path.parts[-4],
                            db=path.parts[-3]
                        )


def blastdb_taxidlist(db, taxa, outfmt="%T"):
    """Call blastdbcmd to list taxonomy metadata.

    Args:
        db (str): the BLAST+ database
        taxa (list): the list of taxonomy identifiers
        outfmt (str): the blastdbcmd outfmt parameter

    Yields:
        str: the next output line
    """
    try:
        cmd = ("blastdbcmd", "-db", db, "-outfmt", outfmt, "-no_taxid_expansion", "-taxidlist", "-")
        with Popen(cmd, stdin=PIPE, stdout=PIPE, universal_newlines=True) as proc:
            with proc.stdin as file:
                print(*taxa, file=file, sep="\n")
            with proc.stdout as file:
                yield from file
    except CalledProcessError:
        yield from ()


def count_nntaxa_in_blastdb(curs, db, taxon, near_neighbors=True):
    """Count the taxon and its relatives in the BLAST+ database.

    Args:
        curs: the open database cursor
        db: the BLAST+ database
        taxon: the taxonomy identifier
        near_neighbors: the flag to include near neibors
    
    Return:
        Counter: the taxon counts
    """
    parent_tax_id = ancestors(curs, taxon)[-1]["parent_" * near_neighbors + "tax_id"]
    lineage = ancestors(curs, parent_tax_id)
    infos = [
        (
            *curs.execute(
                "SELECT id, name_class FROM tax_name WHERE tax_id == :tax_id AND name_class == 'scientific name';",
                dict(tax_id=ele["tax_id"])
            ).fetchone(),
            curs.execute("SELECT rank FROM tax_node WHERE tax_id == :tax_id;", dict(tax_id=ele["tax_id"])).fetchone()[0]
        )
        for ele in lineage
    ]
    meta = {
        **{ele["tax_id"]: dict(ele) for ele in descendants(curs, parent_tax_id)},
        **{ele["tax_id"]: dict(ele, id=info[0], name_class=info[1], rank=info[2]) for ele, info in zip(lineage, infos)}
    }

    T = nx.DiGraph((val["parent_tax_id"], key) for key, val in meta.items())

    if C := Counter(blastdb_taxidlist(db, (target for target, degree in T.out_degree() if not degree))):
        D = defaultdict(int)

        for key, val in sorted(C.items(), key=lambda e: e[1]):
            for ele in nx.shortest_path(T, source=1, target=int(key)):
                D[ele] += val

        for key, val in meta.items():
            val["count"] = D.get(key, 0)
            val["target"] = nx.has_path(T, source=int(taxon), target=int(key))

        return meta

    return C


def blastdbcmd_info():
    """List available BLAST+ databases and metadata.
    
    Returns:
        pd.DataFrame: the info
    """
    fields = "\t".join(("%f", "%p", "%t", "%d", "%l", "%n", "%U", "%v"))
    cmd = ("blastdbcmd", "-recursive", "-remove_redundant_dbs", "-list", BLAST_DIR, "-list_outfmt", fields)
    with Popen(cmd, stdout=PIPE, universal_newlines=True) as proc:
        with proc.stdout as file:
            data = pd.read_table(file, names=("path", "type", "title", "date", "bases", "sequences", "bytes", "version"))
    return data


def nucl_db_v5_choices():
    """List available nucleotide BLAST+ v5 databases and metadata.
    
    Returns:
        dict: the info
    """
    data = blastdbcmd_info()
    data = data[(data["type"] == "Nucleotide") & (data["version"] == 5)]
    return dict(zip(data["path"], data["title"]))


def monitor_snakemake(session, cmd, msg=None):
    """Run snakemake and monitor progress.
    
    Args:
        session: the Shiny session
        cmd: the command tuple to run
        msg: the message to display
    """
    with ui.Progress(min=0, max=100, session=session) as prog:
        prog.set(message=msg, detail=" ".join(cmd))
        with Popen(cmd, universal_newlines=True, bufsize=1, stderr=PIPE) as proc:
            with proc.stderr as file:
                for line in file:
                    line = line.strip()
                    print(line)
                    if line.endswith("%) done"):
                        prog.set(message=msg, detail=line, value=float(line.split(" ")[-2][1:-2]))
