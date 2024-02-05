#!/usr/bin/env python3

import json
import sys
import textwrap
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser
from functools import partial
from multiprocessing import Pool
from pathlib import Path

import jq
import pandas as pd
from plotnine import (
    aes,
    element_blank,
    element_rect,
    element_text,
    facet_grid,
    geom_bar,
    geom_tile,
    ggplot,
    guide_legend,
    guides,
    scale_fill_gradientn,
    theme,
)

jq_prog_load_calls = """
    .assay.id as $id |
    .hits | [
        map(.evals | to_entries | map({id: $id, com: .key, heat: .value.psim, qaln: .value.qaln, saln: .value.saln})),
        map(.calls | to_entries | map(.key as $key | .value | map({call: $key, acc: .})[]))
    ] |
        transpose |
        map(.[0][] as $x | .[1][] as $y | ($x + $y))[]
    """

jq_prog_load_muts = """
    def qpos(str):
        # "GATTACA"    -> [1,2,3,4,5,6,7]
        # "GA--TT-ACA" -> [1,2,2,2,3,4,4,5,6,7]
        str | split("") | map(if . == "-" then 0 else 1 end) | to_entries as $arr |
          $arr | map($arr[:(.key + 1)] | map(.value) | add);

    .assay.id as $id |
    .atype as $atype |
    .hits | [
        map(.evals | to_entries | map({com: .key, psim: .value.psim, astr: .value.astr, qaln: .value.qaln})),
        # select TP/FN only
        map(.calls | to_entries | map(select(.key == "TP" or .key == "FN") | {call: .key, acc: .value[0], nseq: (.value | length)}))
    ] |
        transpose |
        # select non-perfect alignments
        map(.[0][] as $x | .[1][] as $y | select($x.psim < 1) | ({$id} + $x + $y) | del(.psim)) |
        map(
            qpos(.qaln) as $pos |
            . +
                (
                  .astr | split("") | map(tonumber) | to_entries |
                      map(
                          $atype[(.value - 1)] as $mut |
                          select($mut != "IDN" and $mut != "SIM") |
                          { pos: $pos[.key], mut: $mut }
                      )[]
                ) |
                del(.astr, .qaln)
        )[]
    """

jq_prog_load_nmut = """
    .assay.id as $id |
    .atype as $atype |
    .hits | [
        map(.evals | to_entries | map({com: .key, nmut: .value.qaln | split("") | map(select(test("[a-z]"))) | length})),
        # select TP/FN only
        map(.calls | to_entries | map(select(.key == "TP" or .key == "FN") | {call: .key, acc: .value[0], nseq: (.value | length)}))
    ] |
        transpose |
        map(.[0][] as $x | .[1][] as $y | ({$id} + $x + $y))[]
"""


def run_jq(root, prog):
    compiled = jq.compile(prog)
    with Path(root).joinpath("call.json").open() as file:
        return pd.DataFrame(compiled.input_value(json.load(file)))


def plot_heat(data):
    # data["call"] = data["call"].astype(pd.CategoricalDtype(categories=["TP", "FN", "FP", "TN", "XX"], ordered=True))
    # print(data["call"])
    return (
        ggplot(data, aes("acc", "com"))
        + geom_tile(aes(fill="heat"), show_legend=dict(fill=True))
        + facet_grid(
            "id ~ call",
            space="free_y",
            scales="free_y",
            labeller=lambda x: textwrap.fill(x, width=10),
        )
        + scale_fill_gradientn(
            colors=("red", "red", "orange", "green", "darkgreen"),
            breaks=(0.0, 0.75, 0.8, 0.99, 1.0),
            values=(0.0, 0.75, 0.8, 0.99, 1.0),
            limits=(0, 1),
        )
        + guides(fill=guide_legend(reverse=True))
        + theme(
            axis_text_x=element_blank(),
            axis_text_y=element_text(size=6),
            panel_background=element_rect(fill="white"),
            panel_border=element_rect(colour="black", linetype="dotted"),
            strip_background=element_rect(fill="white"),
            strip_text_x=element_text(size=6),
            strip_text_y=element_text(size=6, rotation=0),
        )
    )


def plot_muts(data):
    return (
        ggplot(data, aes(x="pos"))
        + geom_bar(aes(fill="mut"))
        + facet_grid("id ~ com")
        + guides(fill=guide_legend(nrow=1))
        + theme(
            axis_text_y=element_text(size=6),
            panel_border=element_rect(linetype="dotted"),
            strip_background=element_rect(fill="white"),
            strip_text_y=element_text(hjust=0, angle=0),
        )
    )


def plot_nmut_com(data):
    return (
        ggplot(data.groupby(["id", "com", "nmut"], as_index=False).sum("nseq"), aes("nmut", "nseq"))
        + geom_bar(stat="identity")
        + facet_grid("id ~ com")
        + theme(
            axis_text_y=element_text(size=6),
            panel_border=element_rect(linetype="dotted"),
            strip_background=element_rect(fill="white"),
            strip_text_y=element_text(hjust=0, angle=0),
        )
    )


def plot_nmut_acc(data):
    return (
        ggplot(data.groupby(["id", "com", "nmut"], as_index=False).sum("nseq"), aes("nmut", "nseq"))
        + geom_bar(aes(fill="com"), stat="identity")
        + facet_grid("id ~ .")
        + theme(
            axis_text_y=element_text(size=6),
            panel_border=element_rect(linetype="dotted"),
            strip_background=element_rect(fill="white"),
            strip_text_y=element_text(hjust=0, angle=0),
            text=element_text(family="mono"),
        )
    )


def plot_cover(data):
    data = data.loc[data["call"] == "TP"]
    data["n"] = data.groupby("id").transform("size")
    data = data.sort_values(by=["n"], ascending=False).drop_duplicates(subset=["acc"])
    return (
        ggplot(data, aes("id")) +
        geom_bar() +
        theme(
            panel_border=element_rect(linetype="dotted")
        )
    )


def plot_data(path, mode, proc=1):
    work = dict(
        cover=(jq_prog_load_calls, plot_cover),
        heat=(jq_prog_load_calls, plot_heat),
        muts=(jq_prog_load_muts, plot_muts),
        nmut_com=(jq_prog_load_nmut, plot_nmut_com),
        nmut_acc=(jq_prog_load_nmut, plot_nmut_acc),
    )
    prog, plotter = work[mode]
    with Pool(proc) as pool:
        data = pd.concat((ele for ele in pool.map(partial(run_jq, prog=prog), path) if not ele.empty))
    return plotter(data)


def parse_args(argv):
    parser = ArgumentParser(description="", formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument("path", help="the paths to the PSET results", nargs="+")
    choices = ("cover", "heat", "muts", "nmut_com", "nmut_acc")
    parser.add_argument("-mode", choices=choices, default=list(choices)[0])
    parser.add_argument("-variant", choices=choices, default=list(choices)[0])
    parser.add_argument("-out", help="the output file", default="heat.png")
    parser.add_argument("-width", help="the image width in inches", type=float, default=6.5)
    parser.add_argument("-height", help="the image height in inches", type=float, default=6.5)
    parser.add_argument("-dpi", help="the image DPI", type=int, default=300)
    parser.add_argument("-proc", type=int, default=1)
    return parser.parse_args(argv)


def main(argv):
    args = parse_args(argv[1:])
    plot = plot_data(args.path, args.mode, args.proc)
    plot.save(width=args.width, height=args.height, units="in", filename=args.out, dpi=args.dpi)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
