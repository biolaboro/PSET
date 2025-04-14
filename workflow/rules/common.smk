# imports

import json
import shlex
from collections import defaultdict
from itertools import chain
from pathlib import Path
from subprocess import PIPE, Popen, check_call

from Bio import SearchIO, SeqIO
from Bio.SeqUtils.CheckSum import seguid

from pset.assay import Assay, parse_assays
from pset.util import fields_8CB


# parse args
db = Path(config["db"])
context = tuple(map(int, config["context"].split(",")))
confb = [ele.split("=") for ele in shlex.split(config.get("confb", ""))]
confg = [ele.split("=") for ele in shlex.split(config.get("confg", ""))]
dkwargs = {key: tuple(map(int, config[key].split(",", maxsplit=1))) for key in ("dFR", "dF3F2", "dF2F1c", "dF1cB1c")}

# globals
root = Path(config["out"]) / db.name
subroot = Path(config["out"]) / db.name / "{id}"
base = Path(workflow.basedir)

# load assays
ids = set(config.get("assays", "").split(",")) - {""}
with open(config["file"]) as file:
    # assays = {ele.id: ele for ele in parse_assays(file, context=context, target_type=int) if not ids or ele.id in ids}
    assays, problems = [], []
    for idx, ele in enumerate(parse_assays(file, context=context, target_type=int), start=1):
        rec, row = ele
        if rec is None:
            problems.append((idx, row))
        elif not ids or rec.id in ids:
            print(rec.id)
            assays.append(rec)
    nprob = len(problems)
    w = len(str(len(assays) + nprob))
    nprob and print(f"problems with {nprob} assay{'s' * (nprob - 1)}:", file=sys.stderr)
    print(*(f"\t@ row = {idx:0{w}}, id = {row.get('id')}" for idx, row in problems), sep="\n", file=sys.stderr)
    assays = { ele.id: ele for ele in assays }
