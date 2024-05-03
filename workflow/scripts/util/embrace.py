#!/usr/bin/env python3

import sys
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, FileType
from csv import DictReader, DictWriter
from functools import reduce
from itertools import chain, combinations, product
from signal import SIGINT, signal
from subprocess import PIPE, Popen, check_output, run
from tempfile import NamedTemporaryFile

import pandas as pd
from Bio import Seq

from pset.assay import Assay
from pset.util import argify

key_combos = (
    ("F3", "F2", "LF", "F1c", "B1c", "LB", "B2", "B3"),
    ("F3", "F2", "LF", "F1c", "B1c", "B2", "B3"),
    ("F3", "F2", "F1c", "B1c", "LB", "B2", "B3"),
    ("F3", "F2", "F1c", "B1c", "B2", "B3"),
    ("F3", "FIP", "BIP", "B3"),
    ("F", "P", "R"),
    ("F", "R"),
    ("O",),
)

key_brackets = dict(
    zip(
        ("O", "F", "P", "R", "F3", "F2", "LF", "F1c", "B1c", "LB", "B2", "B3"),
        ("[]", "[]", "()", "[]", "[]", "[]", "()", "[]", "[]", "()", "[]", "[]"),
    )
)


def blast_keys(row, keys, **kwargs):
    # find definition subject
    names = ("qaccver", "saccver", "sstart", "send", "sstrand", "pident", "qcovhsp")
    types = (str, str, int, int, str, float, float)
    dtypes = dict(zip(names, types))
    cmd = (
        "blastn",
        "-subject_besthit",
        *("-db", kwargs["db"]),
        *("-task", "blastn"),
        *("-qcov_hsp_perc", "100"),
        *("-num_threads", str(kwargs["nproc"])),
        *("-outfmt", f"6 " + " ".join(names)),
    )
    for key in keys:
        with NamedTemporaryFile() as temp1, NamedTemporaryFile() as temp2:
            temp1.write(f">{key}\n{row[key]}\n".encode())
            temp1.flush()
            run((*cmd, "-query", temp1.name, "-out", temp2.name), universal_newlines=True)
            result = pd.read_csv(temp2, sep="\t", names=names, dtype=dtypes)
            result.rename(columns={ele: f"{ele}_{key}" for ele in set(names) - {"saccver"}}, inplace=True)
            yield result


def deconstruct(row, **kwargs):
    # print("deconstruct")
    keys = ("F3", "LF" if row.get("LF") else "", "LB" if row.get("LB") else "", "B3")
    keys = tuple(ele for ele in keys if ele)
    results = list(blast_keys(row, keys, **kwargs))
    df = reduce(lambda x, y: pd.merge(x, y, on="saccver"), results)

    update = {}

    if df.shape[0] > 0:
        saccver = df[f"saccver"].iloc[0]
        coors = df[[f"sstart_{keys[0]}", f"send_{keys[0]}", f"sstart_{keys[-1]}", f"send_{keys[-1]}"]].iloc[0]
        sstart, send = min(coors), max(coors)
        cmd = ("blastdbcmd", "-db", kwargs["db"], "-entry", saccver, "-outfmt", "%s", "-range", f"{sstart}-{send}")
        record = check_output(cmd, universal_newlines=True)
        outfmt = "qseqid sstart send sstrand qseq"
        with NamedTemporaryFile() as subject, NamedTemporaryFile() as query:
            subject.write(f">sbj\n{record}\n".encode())
            subject.flush()
            for key in keys:
                query.write(f">{key}\n{row[key]}\n".encode())
            for key in ("FIP", "BIP"):
                if row.get(key):
                    query.write(f">IP\n{row[key]}\n".encode())
            query.flush()
            cmd = (
                'blastn',
                '-subject', subject.name,
                '-query', query.name,
                '-task', "blastn-short",
                '-evalue', '0.01',
                '-perc_identity', '90',
                '-outfmt', f"6 {outfmt}"
            )
            with Popen(cmd, universal_newlines=True, stdout=PIPE) as proc:
                with proc.stdout as file:
                    df = pd.read_csv(file, sep="\t", names=outfmt.split(" "), dtype=dict(sstart=int, send=int))

        df.loc[df['sstrand'] == "minus", ['sstart', 'send']] = df.loc[df['sstrand'] == "minus", ['send', 'sstart']].values
        df = df.groupby("qseqid")
        if "IP" in df.groups:
            for F3, B3 in product(df.get_group("F3").to_dict('records'), df.get_group("B3").to_dict('records')):
                F3, B3 = (F3, B3) if int(F3["sstart"]) < int(B3["sstart"]) else (B3, F3)
                for item in combinations(df.get_group("IP").to_dict("records"), 4):
                    item = list(sorted(item, key=lambda x: int(x["sstart"])))
                    # print(*zip(
                    #     (
                    #         f'{len(item[0]["qseq"])} + {len(item[1]["qseq"])} == {len(row["FIP"])}',
                    #         f'{len(item[2]["qseq"])} + {len(item[3]["qseq"])} == {len(row["BIP"])}',
                    #         f'({F3["sstrand"]} == {item[0]["sstrand"]} == {item[2]["sstrand"]})',
                    #         f'({B3["sstrand"]} == {item[3]["sstrand"]} == {item[1]["sstrand"]})',
                    #         f'{F3["send"]} <= {item[0]["sstart"]} <= {item[0]["send"]} <= {item[1]["sstart"]} <= {item[1]["send"]} <= {item[2]["sstart"]} <= {item[2]["send"]} <= {item[3]["sstart"]} <= {item[3]["send"]} <= {B3["sstart"]}'
                    #     ),
                    #     (
                    #         len(item[0]["qseq"]) + len(item[1]["qseq"]) == len(row["FIP"]),
                    #         len(item[2]["qseq"]) + len(item[3]["qseq"]) == len(row["BIP"]),
                    #         (F3["sstrand"] == item[0]["sstrand"] == item[2]["sstrand"]),
                    #         (B3["sstrand"] == item[3]["sstrand"] == item[1]["sstrand"]),
                    #         F3["send"] <= item[0]["sstart"] <=
                    #         item[0]["send"] <= item[1]["sstart"] <=
                    #         item[1]["send"] <= item[2]["sstart"] <=
                    #         item[2]["send"] <= item[3]["sstart"] <=
                    #         item[3]["send"] <= B3["sstart"]
                    #     )
                    # ), sep="\n")
                    # print()
                    if (
                        len(item[0]["qseq"]) + len(item[1]["qseq"]) == len(row["FIP"]) and
                        len(item[2]["qseq"]) + len(item[3]["qseq"]) == len(row["BIP"]) and
                        (F3["sstrand"] == item[0]["sstrand"] == item[2]["sstrand"]) and
                        # (B3["sstrand"] == item[3]["sstrand"] == item[1]["sstrand"]) and
                        F3["send"] <= item[0]["sstart"] <=
                        item[0]["send"] <= item[1]["sstart"] <=
                        item[1]["send"] <= item[2]["sstart"] <=
                        item[2]["send"] <= item[3]["sstart"] <=
                        item[3]["send"] <= B3["sstart"]
                    ):
                        # print("F3>", F3["sstrand"], F3["sstart"], F3["send"])
                        # print("F2>", item[0]["sstrand"], item[0]["sstart"], item[0]["send"])
                        # print("F1c>", item[1]["sstrand"], item[1]["sstart"], item[1]["send"])
                        # print("B1c>", item[2]["sstrand"], item[2]["sstart"], item[2]["send"])
                        # print("B2>", item[3]["sstrand"], item[3]["sstart"], item[3]["send"])
                        # print("B3>", B3["sstrand"], B3["sstart"], B3["send"])
                        update = dict(F2=item[0]["qseq"], F1c=item[1]["qseq"], B1c=item[2]["qseq"], B2=item[3]["qseq"])

    return update


def process(row, **kwargs):
    try:
        Assay.factory(row["definition"])
        return row
    except:
        pass

    # print()
    # print()
    # print(row)
    if row.get("FIP") or row.get("BIP"):
        row.update(deconstruct(row, **kwargs))

    # check and get complete set of assay component keys
    keys = next(ele for ele in key_combos if sum(bool(row.get(key)) for key in ele) == len(ele))
    try:
        results = list(blast_keys(row, keys, **kwargs))
    except StopIteration:
        results = []

    df = reduce(lambda x, y: pd.merge(x, y, on="saccver"), results)
    if df.shape[0] > 0 and len(keys) == len(results):
        saccver = df[f"saccver"].iloc[0]
        coors = df[[f"sstart_{keys[0]}", f"send_{keys[0]}", f"sstart_{keys[-1]}", f"send_{keys[-1]}"]].iloc[0]
        sstart, send = min(coors) - kwargs["context"], max(coors) + kwargs["context"]
        # get the sequence from blastdbcmd
        cmd = ("blastdbcmd", "-db", kwargs["db"], "-entry", saccver, "-outfmt", "%s")
        record = check_output(cmd, universal_newlines=True)
        for key in keys:
            start, end = df[[f"sstart_{key}", f"send_{key}"]].iloc[0]
            start, end = min(start, end), max(start, end)
            # print(key, row[key])
            row[key] = row[key] if df[f"sstrand_{key}"].iloc[0] == "plus" else Seq.reverse_complement(row[key])
            record = record[:start - 1] + row[key] + record[end:]
        for key in keys:
            b = key_brackets[key]
            # print(key, row[key], b)
            record = record.replace(row[key], b[0] + row[key] + b[1])
        try:
            assay = Assay.factory(record[sstart - 1:send + 2 * len(keys)], id=row.get("id", ""))
            # print(assay)
            row["definition"] = assay.definition
            row["reference"] = f"{saccver}:{start}-{end}"
        except:
            row["definition"] = ""
            row["reference"] = ""

    # print()
    # print()

    return row


def main(argv):
    args = parse_args(argv[1:])

    kwargs = dict(db=args.db, context=args.context, nproc=args.nproc)

    with args.file as file:
        rows = list(DictReader(file, delimiter="\t"))

    with args.out as file:
        fieldnames = list(rows[0])
        if "definition" not in fieldnames:
            fieldnames.append("definition")
        if "reference" not in fieldnames:
            fieldnames.append("reference")
        writer = DictWriter(file, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            result = process(row, **kwargs)
            writer.writerow(result)

    return 0


def parse_args(argv):
    parser = ArgumentParser(description="", formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument("file", help="the input file", type=FileType(), default="-")
    parser.add_argument("db", help="the BLAST+ database")
    parser.add_argument("-context", help="the amount of 5'/3' definition context", type=int, default=6)
    parser.add_argument("-nproc", help="the number of processes", type=int, default=1)
    parser.add_argument("-out", help="the output file", type=FileType("w"), default="-")
    return parser.parse_args(argv)


if __name__ == "__main__":
    signal(SIGINT, lambda signum, _: sys.exit(signum))
    sys.exit(main(sys.argv))
