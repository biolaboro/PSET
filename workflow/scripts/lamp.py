#!/usr/bin/env python3

import json
import sys
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, FileType
from collections import Counter, defaultdict
from csv import DictReader
from pathlib import Path

import pandas as pd


def iter_calls(path, id_to_lamp):
    for path in Path(path).rglob("call.json"):
        with path.open() as file:
            obj = json.load(file)
            lamp = id_to_lamp[obj["assay"]["id"]]
            yield from ((lamp, key, ele) for hit in obj["hits"] for key, val in hit["calls"].items() for ele in val)


def parse_args(argv):
    parser = ArgumentParser(description="filter BLAST+", formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument("file", help="the assay file")
    parser.add_argument("call")
    parser.add_argument("-out", type=FileType("w"), default="-")
    return parser.parse_args(argv)


def main(argv):
    args = parse_args(argv[1:])

    with open(args.file) as file:
        id_to_lamp = {row["id"]: row["lamp_id"] for row in DictReader(file, delimiter="\t")}

    data = defaultdict(list)
    for path in Path(args.call).rglob("call.json"):
        with path.open() as file:
            obj = json.load(file)
            for hit in obj["hits"]:
                for key, val in hit["calls"].items():
                    for ele in val:
                        id = obj["assay"]["id"]
                        data[(id_to_lamp[id], ele)].append(key)

    data = {key: Counter(val) for key, val in data.items()}

    counts = defaultdict(int)
    for key, val in data.items():
        if "TP" in val:
            counts[(key[0], ("TP" if val["TP"] == 3 else "FN"))] += 1
        elif "FP" in val:
            counts[(key[0], ("FP" if val["FP"] == 3 else "TN"))] += 1
        else:
            counts[(key[0], key[1])] += 1

    df_data = pd.DataFrame(dict(id=key[0], call=key[1], n=val) for key, val in counts.items())
    calls = ("TP", "TN", "FP", "FN")
    df_data = df_data.assign(call=pd.Categorical(df_data["call"], categories=calls, ordered=True))
    columns = ["id", *calls]
    df_data = df_data.pivot_table("n", "id", "call", fill_value=0).reset_index().reindex(columns=columns, fill_value=0)

    df_data.to_csv(args.out, sep="\t", index=None)

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
