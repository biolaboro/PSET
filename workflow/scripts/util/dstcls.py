#!/usr/bin/env python3

import sys
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, FileType
from os import makedirs
from pathlib import Path

import pandas as pd
from sklearn import manifold
from sklearn.mixture import GaussianMixture
from sklearn.model_selection import GridSearchCV


def parse_args(argv):
    parser = ArgumentParser(description="", formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument("file", help="the FASTA file", type=FileType())
    parser.add_argument("-lab", default="lab", type=Path)
    parser.add_argument("-comp", type=int, default=8)
    parser.add_argument("-eps", type=int, default=1e-9)
    parser.add_argument("-init", type=int, default=10)
    parser.add_argument("-iter", type=int, default=10000)
    parser.add_argument("-seed", type=int, default=7080)
    parser.add_argument("-proc", type=int, default=1)
    return parser.parse_args(argv)


def main(argv):
    args = parse_args(argv[1:])
    makedirs(args.lab, exist_ok=True)

    df = pd.read_csv(args.file, sep="\t", names=("qry", "ref", "dst"))
    X = df.pivot_table("dst", index="qry", columns="ref")
    index = list(X.index)

    mds_args = dict(
        n_components=2,
        metric=True,
        n_init=args.init,
        max_iter=args.iter,
        eps=args.eps,
        n_jobs=args.proc,
        random_state=args.seed,
        dissimilarity="precomputed",
        normalized_stress="auto",
    )
    mds = manifold.MDS(**mds_args)
    mds_args.update(dict(metric=False, n_init=1))
    init = mds.fit(X).embedding_
    mds = manifold.MDS(**mds_args)
    X = mds.fit_transform(X, init=init)

    param_grid = {
        "n_components": range(1, args.comp + 1),
        "covariance_type": ["spherical", "tied", "diag", "full"],
    }
    grid_search = GridSearchCV(GaussianMixture(), param_grid=param_grid, scoring=lambda est, X: -est.bic(X), n_jobs=args.proc)
    grid_search.fit(X)

    df = pd.DataFrame(grid_search.cv_results_)[["param_n_components", "param_covariance_type", "mean_test_score"]]
    df["mean_test_score"] = -df["mean_test_score"]
    df = df.rename(columns=dict(mean_test_score="BIC"))
    print(df.sort_values(by="BIC"), file=sys.stderr)

    preds = grid_search.best_estimator_.predict(X)
    df = pd.DataFrame(zip(index, preds), columns=("acc", "cls")).groupby("cls")
    for cls, grp in df:
        path = args.lab / f"{cls}.tsv"
        grp.to_csv(path, sep="\t", index_label="idx")

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
