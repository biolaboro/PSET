# imports

import json
from collections import OrderedDict, defaultdict
from itertools import chain
from pathlib import Path
from subprocess import PIPE, Popen, check_call

from Bio import SearchIO, SeqIO
from Bio.SeqUtils.CheckSum import seguid


from pset.assay import Assay, parse_assays
from pset.util import argify, fields_8CB


# parse args
db = Path(config["db"])
context = tuple(map(int, config["context"].split(",")))
argsb = list(argify(*filter(len, config.get("confb", "").split(","))))
confb = {ele[0]: ele[1] if len(ele) > 1 else None for ele in argsb}
argsg = list(argify(*filter(len, config.get("confg", "").split(","))))
confg = {ele[0]: ele[1] if len(ele) > 1 else None for ele in argsg}

# globals
root = Path(config["out"]) / db.name
subroot = Path(config["out"]) / db.name / "{id}"
base = Path(workflow.basedir)

# load assays
ids = config.get("assays")
ids = set(ids.split(",")) if ids else ids
with open(config["file"]) as file:
    assays = OrderedDict((ele.id, ele) for ele in parse_assays(file, context=context) if not ids or ele.id in ids)
