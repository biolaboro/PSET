#!/usr/bin/env python3
import json
import sys
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser
from multiprocessing import Pool
from pathlib import Path

import pandas as pd
from plotnine import (
    aes,
    element_blank,
    element_rect,
    element_text,
    facet_grid,
    geom_tile,
    ggplot,
    ggsave,
    theme,
)

path = sys.argv[1]


def read_data(path):
    with path.joinpath("lib.map.json").open() as file:
        df_meta = pd.DataFrame(dict(key=key, acc=ele) for key, val in json.load(file).items() for ele in val)
    df_data = pd.read_csv(path / "hit.tsv", delimiter="\t")
    if not df_data.empty:
        df_data["id"] = path.stem
        df_data = df_data.merge(df_meta, on="acc", how="left")
        df_data["n"] = df_data.groupby("key")["key"].transform("count")
        df_data = df_data[["id", "key", "call", "n", "heat"]].drop_duplicates()
    return df_data


def parse_args(argv):
    parser = ArgumentParser(description="", formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument("path", help="the paths to the PSET results", nargs="+")
    parser.add_argument("-out", help="the output file", default="heat.png")
    parser.add_argument("-scales", help="the figure facet scales parameter", default="free_y")
    parser.add_argument("-dpi", help="the image DPI", type=int, default=300)
    parser.add_argument("-proc", type=int, default=1)
    return parser.parse_args(argv)


def main(argv):
    args = parse_args(argv[1:])

    with Pool(args.proc) as pool:
        results = pool.map(read_data, map(Path, args.path))
        df_data = pd.concat(results)

    calls = ("TP", "TN", "FP", "FN", "XX")
    df_data = df_data.assign(call=pd.Categorical(df_data["call"], categories=calls)).reset_index()

    # tile
    plot = (
        ggplot(df_data, aes("id", "key", fill="heat"))
        + geom_tile()
        + facet_grid("call ~ .", scales=args.scales)
        + theme(
            axis_text_x=element_text(rotation=90, size=4),
            axis_text_y=element_blank(),
            panel_background=element_blank(),
            panel_border=element_rect(colour="black"),
            strip_background=element_rect(color="black"),
        )
    )
    ggsave(plot=plot, width=6.5, height=6.5, units="in", filename=args.out, dpi=args.dpi)

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
