#!/usr/bin/env python3
import json
import sys
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser

import pandas as pd
from plotnine import (
    aes,
    element_blank,
    element_rect,
    facet_wrap,
    geom_bar,
    geom_blank,
    ggplot,
    ggsave,
    scale_fill_discrete,
    scale_y_continuous,
    theme,
)

from pset.assay import AlignType, Assay, align_type


def process_data_json(obj):
    db = obj["db"]
    id = obj["assay"]["id"]
    for hit in obj.get("hits", []):
        evals = hit.get("evals", {})
        calls = hit.get("calls", {})
        # { "calls": { call: [] } }
        n = len(calls.get("TP", [])) + len(calls.get("FN", []))
        if n > 0:
            for key, val in evals.items():
                for idx, ele in enumerate(zip(val["qaln"], val["saln"]), start=1):
                    result = align_type(*ele)
                    if result != AlignType.SIM:
                        yield dict(id=id, db=db, com=key, pos=idx, mut=result.name, n=n)


def parse_args(argv):
    parser = ArgumentParser(description="", formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument("data", help="the call JSON file")
    parser.add_argument("-out", help="the output file prefix", default="out.png")
    parser.add_argument("-dpi", help="the image DPI", type=int, default=300)
    return parser.parse_args(argv)


def main(argv):
    args = parse_args(argv[1:])

    with open(args.data) as file:
        obj = json.load(file)

    assay = obj["assay"]
    assay = Assay.factory(assay["id"], assay["definition"], set(assay.get("targets", [1])))
    df_data = pd.DataFrame(process_data_json(obj))

    aln_type_exclude = ("SIM", "UNK", "DIS")
    categories = tuple(ele.name for ele in AlignType if not ele.name in set(aln_type_exclude))

    if not df_data.empty:
        df_data = df_data[df_data["com"] != "amplicon"]
        df_data = df_data[~df_data["mut"].isin(aln_type_exclude)]
        df_data = df_data.assign(call=pd.Categorical(df_data["mut"], categories=categories))
        df_data.reset_index(inplace=True)

    records = tuple(assay.records(*assay.primers))
    df_dummy = pd.concat(
        (
            pd.DataFrame((dict(com=ele.id, pos=1, n=0, mut="GAP") for ele in records)),
            pd.DataFrame((dict(com=ele.id, pos=len(ele) + 1, n=0, mut="GAP") for ele in records)),
        )
    )

    df_data = df_dummy if df_data.empty else df_data

    plot = (
        ggplot(df_data, aes("pos", "n", fill="mut"))
        + geom_blank(data=df_dummy)
        + geom_bar(stat="identity")
        + facet_wrap("~ com", ncol=1, scales="free_x")
        + scale_y_continuous()
        + scale_fill_discrete(limits=categories)
        + theme(panel_background=element_blank(), panel_border=element_rect(colour="black"), subplots_adjust=dict(hspace=0.5))
    )
    ggsave(plot=plot, width=6.5, height=6.5, units="in", filename=args.out, dpi=args.dpi)

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
