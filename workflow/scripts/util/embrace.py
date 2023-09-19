#!/usr/bin/env python3

import sys
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, FileType
from csv import DictReader
from functools import partial, reduce
from itertools import chain
from multiprocessing import Pool
from subprocess import DEVNULL, PIPE, Popen, run
from tempfile import NamedTemporaryFile

import pandas as pd
from Bio import SearchIO, Seq, SeqIO, pairwise2

from pset.assay import LAMP, PCR, Oligo
from pset.util import argify, fields_8CB, iter_hsps

key_combos = (
    (("A",), Oligo),
    (("F", "P", "R"), PCR),
    (("F", "R"), PCR),
    (("F3", "F2", "LF", "F1c", "B1c", "LB", "B2", "B3"), LAMP),
    (("F3", "F2", "LF", "F1c", "B1c", "B2", "B3"), LAMP),
    (("F3", "F2", "F1c", "B1c", "LB", "B2", "B3"), LAMP),
    (("F3", "F2", "F1c", "B1c", "B2", "B3"), LAMP),
)

key_brackets = dict(
    zip(
        ("A", "F", "P", "R", "F3", "F2", "LF", "F1c", "B1c", "LB", "B2", "B3"),
        ("[]", "[]", "()", "[]", "[]", "[]", "()", "[]", "[]", "()", "[]", "[]"),
    )
)


def process(row, **kwargs):
    # infer assay type from fields
    keys = ()
    for combo, fn in key_combos:
        if sum(bool(row.get(key)) for key in combo) == len(combo):
            keys = combo
            definition = []
            for key in combo:
                b = key_brackets[key]
                definition.append(b[0] + row[key] + b[1])
            assay = fn("".join(definition))
            break

    if not keys:
        raise ValueError("incomplete assay key set")

    # find definition subject
    names = "qaccver saccver sstart send sstrand pident qcovhsp"
    types = (str, str, int, int, str, float, float)
    dtypes = dict(zip(names, types))
    targets = row.get("targets", "").replace(";", ",")
    cmd = (
        "blastn",
        "-subject_besthit",
        *("-db", kwargs["db"]),
        *("-task", "blastn"),
        *("-outfmt", f"6 " + names),
        *(("-taxids", targets) if targets else ()),
        *chain.from_iterable(kwargs["confb"]),
    )
    # print(cmd, file=sys.stderr)
    results = []
    for key in keys:
        with NamedTemporaryFile() as temp1, NamedTemporaryFile() as temp2:
            temp1.write(f">{key}\n{row[key]}\n".encode())
            temp1.flush()
            run((*cmd, "-query", temp1.name, "-out", temp2.name), universal_newlines=True)
            results.append(pd.read_csv(temp2, sep="\t", names=names.split(" "), dtype=dtypes))

    df = reduce(lambda x, y: pd.merge(x, y, on="saccver"), results)

    result = {}
    if not df.empty:
        with NamedTemporaryFile() as temp:
            cmd = ("blastdbcmd", "-db", kwargs["db"], "-entry", df.head(1).iloc[0]["saccver"], "-out", temp.name)
            # print(cmd[:-1], "temp", file=sys.stderr)
            run(cmd)
            record = SeqIO.read(temp.name, "fasta")
            cmd = ("glsearch36", "-n", "-m", "8CB", *chain.from_iterable(kwargs["confg"]), "-", temp.name)
            # print(cmd[:-1], "temp", file=sys.stderr)
            temp.seek(0)
            with Popen(cmd, stdin=PIPE, stdout=PIPE, stderr=DEVNULL, universal_newlines=True) as proc:
                with proc.stdin as file:
                    for key in keys:
                        print(f">{key}", row[key], sep="\n", file=file)
                with proc.stdout as file:
                    hsps = list(iter_hsps(SearchIO.parse(file, "blast-tab", comments=True, fields=fields_8CB)))
                    hit = next(assay.hits(hsps, **kwargs["dkwargs"]), ())

        if hit:
            p1, p2 = hit[0].hit_start - kwargs["context"], hit[-1].hit_end + kwargs["context"]
            p1 = 0 if p1 < 0 else p1
            definition = str(record.seq[p1:p2])
            sseqs = [record.seq[hsp.hit_start : hsp.hit_end] for hsp in hit]
            if combo[0] != hit[0].query_id:
                definition = Seq.reverse_complement(definition)
                sseqs = [ele.reverse_complement() for ele in sseqs][::-1]
                hit = hit[::-1]

            if "P" in assay:
                P1, P2 = row["P"], Seq.reverse_complement(row["P"])
                S1, S2 = pairwise2.align.globalxx(P1, sseqs[1])[0].score, pairwise2.align.globalxx(P2, sseqs[1])[0].score
                row["P"] = P2 if S1 < S2 else P1

            for key, hsp, sseq in zip(combo, hit, sseqs):
                b = key_brackets[key]
                sseq = str(sseq)
                definition = definition.replace(sseq, b[0] + row[key] + b[1])

            assay = assay.factory(definition)
            result = dict(defintion=definition)
            result["hit"] = f"{record.id}:{p1}-{p2}"
            for hsp in hit:
                result[f"{hsp.query_id}.pid"] = hsp.ident_pct
                result[f"{hsp.query_id}.gap"] = hsp.gapopen_num

    return result


def main(argv):
    args = parse_args(argv[1:])

    argsb = list(argify(*filter(len, args.confb.split(","))))
    confb = tuple((ele[0], ele[1]) if len(ele) > 1 else None for ele in argsb)
    argsg = list(argify(*filter(len, args.confg.split(","))))
    confg = tuple((ele[0], ele[1]) if len(ele) > 1 else None for ele in argsg)
    dkwargs = dict(dFR=args.dFR, dF3F2=args.dF3F2, dF2F1c=args.dF2F1c, dF1cB1c=args.dF1cB1c)
    kwargs = dict(db=args.db, context=args.context, confb=confb, confg=confg, dkwargs=dkwargs)

    with args.file as file:
        rows = list(DictReader(file, delimiter="\t"))

    with Pool(processes=args.nproc) as pool:
        df = pd.DataFrame(pool.map(partial(process, **kwargs), rows))

    df.to_csv(args.out, sep="\t", index=False)

    return 0


def parse_args(argv):
    parser = ArgumentParser(description="", formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument("file", help="the input file", type=FileType(), default="-")
    parser.add_argument("db")
    parser.add_argument("-context", help="the amount of 5'/3' definition context", type=int, default=6)
    parser.add_argument("-confb", help="the list of additional flags to configure blastn", default="qcov_hsp_perc=100")
    parser.add_argument("-confg", help="the list of additional flags to configure glsearch36", default="E=10000")
    parser.add_argument("-dFR", help="the min/max distance range (inclusive) for F/R primers", default=(1, 1000), type=int, nargs=2)
    parser.add_argument("-dF3F2", help="the min/max distance range (inclusive) for F3/F2 primers", default=(20, 80), type=int, nargs=2)
    parser.add_argument("-dF2F1c", help="the min/max distance range (inclusive) for F2/F1c primers", default=(20, 80), type=int, nargs=2)
    parser.add_argument("-dF1cB1c", help="the min/max distance range (inclusive) for F1c/B1c primers", default=(1, 100), type=int, nargs=2)
    parser.add_argument("-nproc", help="the number of processes", type=int, default=1)
    parser.add_argument("-out", help="the output file", type=FileType("w"), default="-")
    return parser.parse_args(argv)


if __name__ == "__main__":
    sys.exit(main(sys.argv))
