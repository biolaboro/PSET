#!/usr/bin/env python3
import json
import sys
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser
from multiprocessing import Pool

import pandas as pd

from pset.assay import align_type


def read_json(path):
    def process_data_json(obj):
        db = obj["db"]
        id = obj["assay"]["id"]
        for hit in obj.get("hits", []):
            evals = hit.get("evals", {})
            calls = hit.get("calls", {})
            # { "calls": { call: [] } }
            n = sum(len(calls.get(key, [])) for key in ("TP", "FN", "TPN", "FNN"))
            if n > 0:
                for key, val in evals.items():
                    if key != "amplicon":
                        for idx, ele in enumerate(zip(val["qaln"], val["saln"]), start=1):
                            result = align_type(*ele)
                            yield dict(id=id, db=db, com=key, pos=idx, mut=result.name, n=n, qlen=len(evals[key]["qaln"]))

    with open(path) as file:
        return pd.DataFrame(process_data_json(json.load(file)))


def func(df, n):
    group_total = df["ncall"].sum()
    return ((df["nmut"] * df["ncall"]).sum() / group_total) * (group_total / n)


def parse_args(argv):
    parser = ArgumentParser(description="", formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument("-json", help="the call JSON file", nargs="+")
    parser.add_argument("-out", help="the output file prefix", default="out.png")
    parser.add_argument("-proc", type=int, default=1)
    return parser.parse_args(argv)


def main(argv):
    args = parse_args(argv[1:])

    with Pool(args.proc) as pool:
        results = pool.map(read_json, args.json)
        df = pd.concat(results)

    mpos = dict(df[["com", "pos"]].groupby(["com"]).max().reset_index().to_dict(orient="split")["data"])
    df.loc[df["mut"].isin(("SIM", "UNK")), "n"] = 0
    df = (
        pd.pivot(df.groupby(["db", "id", "com", "pos"])["n"].sum().reset_index(name="n"), index=["db", "id"], columns=["com", "pos"], values="n")
        .reset_index()
        .fillna(0)
    )
    df.to_csv("out.tsv", sep="\t")

    # Create a Pandas Excel writer using XlsxWriter as the engine.
    from xlsxwriter import Workbook
    from xlsxwriter.utility import xl_rowcol_to_cell

    workbook = Workbook(args.out)
    worksheet_data = workbook.add_worksheet("data")
    worksheet_plot = workbook.add_worksheet("plot")

    for idx, row in df.iterrows():
        data = [row[(key, i)] for key, val in mpos.items() for i in range(1, val + 1)]
        worksheet_data.write_row(f"A{idx + 1}", data=data)

    keys = ("primer1", "probe", "primer2")
    worksheet_plot.write_row("A1", ("id", *keys))
    for idx in range(len(df)):
        pcol = 1
        worksheet_plot.write_row(row=idx + 1, col=0, data=[df["id"][idx]])
        for col, key in enumerate(keys):
            val = mpos[key]
            col1 = xl_rowcol_to_cell(idx, pcol)
            col2 = xl_rowcol_to_cell(idx, pcol + val - 1)
            worksheet_plot.add_sparkline(row=idx + 1, col=col + 1, options=dict(range=f"data!{col1}:{col2}", type="column"))
            pcol += val

    workbook.close()

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
