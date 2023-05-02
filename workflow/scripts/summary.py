#!/usr/bin/env python3

import json
import sys
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser
from functools import partial
from multiprocessing import Pool

import pandas as pd

from pset.assay import Assay


def read_data_json_1(path):
    def process_data_json(obj):
        db = obj["db"]
        id = obj["assay"]["id"]
        assay = Assay.factory(obj["assay"]["id"], obj["assay"]["definition", obj["assay"]["targets"]])
        for hit in obj.get("hits", []):
            evals = hit.get("evals", {})
            calls = hit.get("calls", {})
            for com, eval in evals.items():
                for mut, nmut in eval["atypes"].items():
                    for call, accs in calls.items():
                        yield dict(db=db, id=id, com=com, mut=mut, nmut=nmut, qsim=eval["qsim"], call=call[:2], ncall=len(accs))

    with open(path) as file:
        return pd.DataFrame(process_data_json(json.load(file)))


def read_data_json_2(path):
    def process_data_json(obj):
        db = obj["db"]
        id = obj["assay"]["id"]
        nn = set(obj["near"])
        for hit in obj["hits"]:
            for key, val in hit["calls"].items():
                for ele in val:
                    yield dict(db=db, id=id, acc=ele, nn=ele in nn, call=key[:2])

    with open(path) as file:
        return pd.DataFrame(process_data_json(json.load(file)))


def read_data_json_3(path):
    def process_data_json(obj):
        db = obj["db"]
        id = obj["assay"]["id"]
        for hit in obj["hits"]:
            evals = hit.get("evals", {})
            calls = hit.get("calls", {})
            # ignore amplicon entry
            del evals["amplicon"]
            # ignore XX
            if "XX" in calls:
                del calls["XX"]
            # the number of non XX-calls
            n = sum(map(len, calls.values()))
            # the total length of alignment areas
            l = sum(len(val["saln"]) for val in evals.values())
            # the number of matches across alignment areas
            m = sum(val["atypes"].get("SIM", 0) for val in evals.values())
            yield dict(db=db, id=id, n=n, asim=m / l)

    with open(path) as file:
        df = pd.DataFrame(process_data_json(json.load(file)))
        df = df.groupby(["db", "id"]).apply(lambda x: (x["asim"] * x["n"]).sum() / x["n"].sum())
        return df


def func(df, n):
    group_total = df["ncall"].sum()
    return ((df["nmut"] * df["ncall"]).sum() / group_total) * (group_total / n)


def parse_args(argv):
    parser = ArgumentParser(description="", formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument("-json", help="the JSON files", nargs="+")
    parser.add_argument("-proc", type=int, default=1)
    parser.add_argument("-out")
    return parser.parse_args(argv)


def main(argv):
    args = parse_args(argv[1:])

    with Pool(args.proc) as pool:
        results = pool.map(read_data_json_1, args.json)
        df_data_1 = pd.concat(results)
        results = pool.map(read_data_json_2, args.json)
        df_data_2 = pd.concat(results)
        results = pool.map(read_data_json_3, args.json)
        df_data_3 = pd.concat(results).reset_index(name="asim")

    trues = {"TP", "FN"}

    df_nn = df_data_2.loc[df_data_2["nn"]].groupby(["db", "id", "nn", "call"]).size().reset_index(name="n")
    df_nn["call"] = "NN." + df_nn["call"]
    df_nn = pd.pivot_table(df_nn, values="n", index=["db", "id"], columns="call").reset_index()

    categories = ("TP", "TN", "FP", "FN", "P")
    df_data_2 = df_data_2.assign(call=pd.Categorical(df_data_2["call"], categories=categories, ordered=True))
    df_data_2.reset_index(inplace=True)
    df_confusion = pd.pivot_table(
        df_data_2.groupby(["db", "id", "call"]).size().reset_index(name="n"), values="n", index=["db", "id"], columns="call"
    )

    # the number of true accessions
    df_data_2 = df_data_2.loc[df_data_2["call"].isin(trues)]
    n = df_data_2["acc"].nunique()
    # the number of true accessions by assay
    df_counts = df_data_2.groupby(["db", "id", "call"]).size().reset_index(name="count")
    df_counts["cov"] = df_counts["count"] / n
    df_confusion["P"] = n

    print(df_data_1)
    df_data_1 = df_data_1.loc[df_data_1["call"].isin(trues)]
    df_data_1["mut"] = df_data_1["mut"].map({"TRS": "SUB", "TRV": "SUB", "DIS": "SUB", "INS": "IND", "DEL": "IND"})
    df_data_1 = df_data_1.assign(call=pd.Categorical(df_data_1["mut"], categories=["IND", "SUB"], ordered=True))
    df_data_1 = df_data_1.groupby(["db", "id", "com", "mut"], group_keys=False).apply(partial(func, n=n)).reset_index(name="n")
    df_data_1.reset_index(inplace=True)
    df_data_1 = pd.pivot_table(df_data_1, values="n", index=["db", "id"], columns=["com", "mut"], dropna=False)

    df_confusion = df_confusion.reset_index().set_index(["db", "id"])
    df_nn = df_nn.set_index(["db", "id"])
    df_data_3 = df_data_3.set_index(["db", "id"])
    df_data_1 = df_data_1.reset_index().set_index(["db", "id"])

    df_result = pd.merge(pd.merge(df_confusion, df_nn, on=["db", "id"], how="left"), df_data_3, on=["db", "id"], how="left")
    df_result = pd.merge(pd.concat([df_result], keys=[""], axis=1), df_data_1, on=["db", "id"], how="left")
    df_result.fillna(0).to_excel(args.out, sheet_name="summary")

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
