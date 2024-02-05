#!/usr/bin/env python3

import json
import sys
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, FileType
from collections import defaultdict
from itertools import product
from multiprocessing import Pool
from operator import itemgetter

import primer3
from Bio import SeqIO
from Bio.Seq import reverse_complement as rc
from Bio.SeqRecord import SeqRecord

from pset.assay import LAMP, PCR, Assay, amb


def to_seq_args(records, l=None, r=None):
    # left/right -> (start-distance, size)
    for record in records:
        seq_args = dict(SEQUENCE_TEMPLATE=str(record.seq))
        if l and r:
            seq_args["SEQUENCE_PRIMER_PAIR_OK_REGION_LIST"] = (l[0], l[1] - l[0], len(record) - r[1], r[1] - r[0])
        elif l and not r:
            seq_args["SEQUENCE_INCLUDED_REGION"] = (l[0], l[1] - l[0])
        elif not l and r:
            seq_args["SEQUENCE_INCLUDED_REGION"] = (len(record) - r[1], r[1] - r[0])
        yield seq_args


def iter_coors(results, key="PRIMER_PAIR_NUM_RETURNED"):
    yield from ((i, j) for i in range(len(results)) for j in range(results[i][key]))


def primer_right_normalize(idx_len):
    i, l = idx_len
    return i - l + 1, l


def subamplicon(coor, results, records):
    i, j = coor
    # The selected left primer (the primer to the left in the input sequence).
    # i is the 0-based index of the start base of the primer, and n is its length.
    idx1, len1 = results[i][f"PRIMER_LEFT_{j}"]
    # The selected right primer (the primer to the right in the input sequence).
    # i is the 0-based index of the last base of the primer, and n is its length.
    idx2, _ = primer_right_normalize(results[i][f"PRIMER_RIGHT_{j}"])
    # substring range
    range = (idx1 + len1, idx2)
    # record
    id = f"{coor[0]}-{coor[1]},{range[0]}-{range[1]}"
    annotations = dict(range=range, coor=coor)
    record = SeqRecord(records[i][slice(*range)].seq, id=id, annotations=annotations)
    return record


def main_pcr(args, conf, records):
    # Primer3
    with Pool(args.proc) as pool:
        global_args = conf["GLOBAL"]
        global_args["PRIMER_PICK_LEFT_PRIMER"] = 1
        global_args["PRIMER_PICK_RIGHT_PRIMER"] = 1
        iterable = ((dict(SEQUENCE_TEMPLATE=str(record.seq)), global_args) for record in records)
        results = list(pool.starmap(primer3.bindings.design_primers, iterable))

    assays = []
    counter = defaultdict(int)
    context = conf.get("PSET", {}).get("context", (0, 0))
    for i, e in enumerate(zip(records, results)):
        record, result = e
        for j in range(result["PRIMER_PAIR_NUM_RETURNED"]):
            # context
            definition = str(record.seq)
            pos1 = result[f"PRIMER_LEFT_{j}"][0]
            pos2 = result[f"PRIMER_RIGHT_{j}"][0] + 1
            pos1 -= context[0]
            pos1 = 0 if pos1 < 0 else pos1
            pos2 += context[1]
            # define
            definition = definition[pos1:pos2]
            seq = result[f"PRIMER_LEFT_{j}_SEQUENCE"]
            definition = definition.replace(seq, f"[{seq}]")
            seq = result.get(f"PRIMER_INTERNAL_{j}_SEQUENCE")
            definition = definition.replace(seq, f"({seq})") if seq else definition
            seq = rc(result[f"PRIMER_RIGHT_{j}_SEQUENCE"])
            definition = definition.replace(seq, f"[{seq}]")
            # score
            penalty = result[f"PRIMER_PAIR_{j}_PENALTY"]
            # output
            counter[record.id] += 1
            assay = Assay.factory(definition)
            assays.append(((i, counter[record.id], pos1, pos2), assay, penalty))

    assays.sort(key=lambda x: x[-1])
    for coor, assay, penalty in assays:
        i, n, pos1, pos2 = coor
        width1 = len(str(len(records[i])))
        width2 = len(str(counter[records[i].id]))
        assay.id = f"{records[i].id}_{pos1:0{width1}d}-{pos2:0{width1}d}_{n:0{width2}d}"
        yield (assay, penalty)


def main_lamp(args, conf, records):
    stats = defaultdict(int)

    # Primer3
    with Pool(args.proc) as pool:
        # Calculate F3/B3
        conf["F3B3"]["PRIMER_PICK_LEFT_PRIMER"] = 1
        conf["F3B3"]["PRIMER_PICK_RIGHT_PRIMER"] = 1
        iterable = ((ele, conf["F3B3"]) for ele in to_seq_args(records))
        results1 = pool.starmap(primer3.bindings.design_primers, iterable)
        records1 = [subamplicon(coor, results1, records) for coor in iter_coors(results1)]
        print("F3/B3", len(records1), file=sys.stderr)
        # Calculate F2/B2
        conf["F2B2"]["PRIMER_PICK_LEFT_PRIMER"] = 1
        conf["F2B2"]["PRIMER_PICK_RIGHT_PRIMER"] = 1
        dist_range = conf["LAMP"]["F3F2_DIST_RANGE"]
        iterable = ((ele, conf["F2B2"]) for ele in to_seq_args(records1, l=dist_range, r=dist_range))
        results2 = pool.starmap(primer3.bindings.design_primers, iterable)
        records2 = [subamplicon(coor, results2, records1) for coor in iter_coors(results2)]
        print("F2/B2", len(records2), file=sys.stderr)
        # Calculate LF
        conf["LFLB"]["PRIMER_PICK_LEFT_PRIMER"] = 0
        conf["LFLB"]["PRIMER_PICK_RIGHT_PRIMER"] = 1
        iterable = ((ele, conf["LFLB"]) for ele in to_seq_args(records2))
        results3LF = list(pool.starmap(primer3.bindings.design_primers, iterable))
        print("LF", sum(result["PRIMER_RIGHT_NUM_RETURNED"] for result in results3LF), file=sys.stderr)
        # Calculate F1c
        conf["F1cB1c"]["PRIMER_PICK_LEFT_PRIMER"] = 0
        conf["F1cB1c"]["PRIMER_PICK_RIGHT_PRIMER"] = 1
        dist_range = conf["LAMP"]["F2F1c_DIST_RANGE"]
        iterable = ((ele, conf["F1cB1c"]) for ele in to_seq_args(records2, l=dist_range))
        results3F1c = list(pool.starmap(primer3.bindings.design_primers, iterable))
        print("F1c", sum(result["PRIMER_RIGHT_NUM_RETURNED"] for result in results3F1c), file=sys.stderr)
        # Calculate LB
        conf["LFLB"]["PRIMER_PICK_LEFT_PRIMER"] = 1
        conf["LFLB"]["PRIMER_PICK_RIGHT_PRIMER"] = 0
        iterable = ((ele, conf["LFLB"]) for ele in to_seq_args(records2))
        results3LB = list(pool.starmap(primer3.bindings.design_primers, iterable))
        print("LB", sum(result["PRIMER_LEFT_NUM_RETURNED"] for result in results3LB), file=sys.stderr)
        # Calculate B1c
        conf["F1cB1c"]["PRIMER_PICK_LEFT_PRIMER"] = 1
        conf["F1cB1c"]["PRIMER_PICK_RIGHT_PRIMER"] = 0
        dist_range = conf["LAMP"]["F2F1c_DIST_RANGE"]
        iterable = ((ele, conf["F1cB1c"]) for ele in to_seq_args(records2, r=dist_range))
        results3B1c = list(pool.starmap(primer3.bindings.design_primers, iterable))
        print("B1c", sum(result["PRIMER_LEFT_NUM_RETURNED"] for result in results3B1c), file=sys.stderr)

    # process combos
    sublamps = []

    oligo_calc = primer3.thermoanalysis.ThermoAnalysis(**conf["THERMO"])
    counter = defaultdict(int)
    for i in range(len(records2)):
        iP2, jP2 = records2[i].annotations["coor"]
        iP3, jP3 = records1[iP2].annotations["coor"]
        pos1 = results1[iP3][f"PRIMER_LEFT_{jP3}"][0]
        pos2 = results1[iP3][f"PRIMER_RIGHT_{jP3}"][0] + 1

        # check F2/B2 stability
        seq = results2[iP2][f"PRIMER_LEFT_{jP2}_SEQUENCE"]
        if oligo_calc.calc_end_stability(seq, seq).dh >= 0:
            stats["F2:calc_end_stability"] += 1
            continue
        seq = rc(results2[iP2][f"PRIMER_RIGHT_{jP2}_SEQUENCE"])
        if oligo_calc.calc_end_stability(seq, seq).dh >= 0:
            stats["B2:calc_end_stability"] += 1
            continue

        # sort primers
        # check LF-F1c overlap/order
        # check F1c stability
        result3LF, result3F1c = results3LF[i], results3F1c[i]
        n1, n2 = result3LF["PRIMER_RIGHT_NUM_RETURNED"], result3F1c["PRIMER_RIGHT_NUM_RETURNED"]
        F = []

        if args.optional_loop:
            for j2 in range(n2):
                F1c = primer_right_normalize(result3F1c[f"PRIMER_RIGHT_{j2}"])
                seq = result3F1c[f"PRIMER_RIGHT_{j2}_SEQUENCE"]
                is_stable = oligo_calc.calc_end_stability(seq, seq).dh < 0
                if is_stable:
                    F.append((-1, j2))
                else:
                    stats["LF-F1c:unstable"] += not is_stable

        for j1, j2 in product(range(n1), range(n2)):
            LF = primer_right_normalize(result3LF[f"PRIMER_RIGHT_{j1}"])
            F1c = primer_right_normalize(result3F1c[f"PRIMER_RIGHT_{j2}"])
            seq = result3F1c[f"PRIMER_RIGHT_{j2}_SEQUENCE"]
            is_ordered = sum(LF) < F1c[0]
            is_stable = oligo_calc.calc_end_stability(seq, seq).dh < 0
            if is_ordered and is_stable:
                F.append((j1, j2))
            else:
                stats["LF-F1c:overlap"] += not is_ordered
                stats["LF-F1c:unstable"] += not is_stable

        if not F:
            stats["LF-F1c:none"] += 1
            continue

        # sort primers
        # check B1c-LB overlap/order
        # check B1c stability
        result3B1c, result3LB = results3B1c[i], results3LB[i]
        n1, n2 = result3B1c["PRIMER_LEFT_NUM_RETURNED"], result3LB["PRIMER_LEFT_NUM_RETURNED"]
        B = []

        if args.optional_loop:
            for j1 in range(n1):
                B1c = result3B1c[f"PRIMER_LEFT_{j1}"]
                seq = rc(result3B1c[f"PRIMER_LEFT_{j1}_SEQUENCE"])
                is_stable = oligo_calc.calc_end_stability(seq, seq).dh < 0
                if is_ordered and is_stable:
                    B.append((j1, -1))
                else:
                    stats["B1c-LB:unstable"] += not is_stable

        for j1, j2 in product(range(n1), range(n2)):
            B1c = result3B1c[f"PRIMER_LEFT_{j1}"]
            LB = result3LB[f"PRIMER_LEFT_{j2}"]
            seq = rc(result3B1c[f"PRIMER_LEFT_{j1}_SEQUENCE"])
            is_ordered = sum(B1c) < LB[0]
            is_stable = oligo_calc.calc_end_stability(seq, seq).dh < 0
            if is_ordered and is_stable:
                B.append((j1, j2))
            else:
                stats["B1c-LB:overlap"] += not is_ordered
                stats["B1c-LB:unstable"] += not is_stable

        if not B:
            stats["B1c-LB:none"] += 1

        # check F1c-B1c distance
        key = (records[iP3].id, pos1, pos2)
        for f, b in product(F, B):
            _, f2, b1, __ = *f, *b
            F1c = result3F1c[f"PRIMER_RIGHT_{f2}"]
            B1c = result3B1c[f"PRIMER_LEFT_{b1}"]
            if conf["LAMP"]["F1cB1c_DIST_RANGE"][0] <= B1c[0] - sum(F1c) + 1 <= conf["LAMP"]["F1cB1c_DIST_RANGE"][1] and counter[key] < args.limit:
                sublamps.append((i, f, b))
                counter[key] += 1
            else:
                stats["F1c-B1c:distance"] += 1

    print(*(ele for ele in stats.items() if ele[1]), file=sys.stderr, sep="\n")

    print("output", file=sys.stderr)
    assays = []
    context = conf.get("PSET", {}).get("context", (0, 0))
    brackets = "[]", "[]", "()", "[]", "[]", "()", "[]", "[]"
    counter.clear()
    for i, coor1, coor2 in sublamps:
        iP2, jP2 = records2[i].annotations["coor"]
        iP3, jP3 = records1[iP2].annotations["coor"]
        pos1 = results1[iP3][f"PRIMER_LEFT_{jP3}"][0]
        pos2 = results1[iP3][f"PRIMER_RIGHT_{jP3}"][0] + 1

        sequences = (
            results1[iP3][f"PRIMER_LEFT_{jP3}_SEQUENCE"],  # F3
            results2[iP2][f"PRIMER_LEFT_{jP2}_SEQUENCE"],  # F2
            " " if coor1[0] == -1 else rc(results3LF[i][f"PRIMER_RIGHT_{coor1[0]}_SEQUENCE"]),  # LF
            rc(results3F1c[i][f"PRIMER_RIGHT_{coor1[1]}_SEQUENCE"]),  # F1c
            results3B1c[i][f"PRIMER_LEFT_{coor2[0]}_SEQUENCE"],  # B1c
            " " if coor2[1] == -1 else results3LB[i][f"PRIMER_LEFT_{coor2[1]}_SEQUENCE"],  # LB
            rc(results2[iP2][f"PRIMER_RIGHT_{jP2}_SEQUENCE"]),  # B2
            rc(results1[iP3][f"PRIMER_RIGHT_{jP3}_SEQUENCE"]),  # B3
        )
        penalty = sum(
            (
                results1[iP3][f"PRIMER_PAIR_{jP3}_PENALTY"],
                results2[iP2][f"PRIMER_PAIR_{jP2}_PENALTY"],
                0 if coor1[0] == -1 else results3LF[i][f"PRIMER_RIGHT_{coor1[0]}_PENALTY"],
                results3F1c[i][f"PRIMER_RIGHT_{coor1[1]}_PENALTY"],
                results3B1c[i][f"PRIMER_LEFT_{coor2[0]}_PENALTY"],
                0 if coor2[1] == -1 else results3LB[i][f"PRIMER_LEFT_{coor2[1]}_PENALTY"],
            )
        )
        # define
        cpos1 = pos1 - context[0]
        cpos1 = 0 if cpos1 < 0 else cpos1
        cpos2 = pos2 + context[1]
        definition = str(records[iP3][cpos1:cpos2].seq)
        for delimiter, sequence in zip(brackets, sequences):
            definition = definition.replace(sequence, f"{delimiter[0]}{sequence}{delimiter[1]}")
        # counter
        counter[records[iP3].id] += 1
        # append
        assay = Assay.factory(definition)
        assays.append(((iP3, counter[records[iP3].id], pos1, pos2), assay, penalty))

    assays.sort(key=lambda x: x[-1])
    for coor, assay, penalty in assays:
        i, n, pos1, pos2 = coor
        width1 = len(str(len(records[i])))
        width2 = len(str(counter[records[i].id]))
        assay.id = f"{records[i].id}_{pos1:0{width1}d}-{pos2:0{width1}d}_{n:0{width2}d}"
        yield (assay, penalty)


def parse_config(path, cstr):
    cobj = defaultdict(dict)
    for ele in cstr:
        tokens = ele.split(":")
        key, val = tokens[1].split("=", maxsplit=1)
        val = val.split(",")
        val = list(map(int, val)) if val[0].isdigit() else val
        val = val[0] if len(val) == 1 else val
        cobj[tokens[0]][key] = val
    with open(path) as file:
        config = json.load(file)
        config["GLOBAL"].update(cobj["GLOBAL"])
        for key in ("F3B3", "F2B2", "F1cB1c", "LFLB"):
            config[key] = {**config["GLOBAL"], **config[key], **cobj[key]}
        for key in ("THERMO", "PSET", "LAMP"):
            config[key].update(cobj[key])
    return config


def parse_args(argv):
    parser = ArgumentParser(description="FASTA file -> primers", formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument("seqs", help="the FASTA file", type=FileType())
    choices = ("PCR", "LAMP")
    parser.add_argument("-mode", help="the assay mode", choices=choices, default=choices[0])
    parser.add_argument("-conf", help="the JSON file of Primer3 global input tags (https://primer3.org/manual.html#globalTags)")
    parser.add_argument("--optional-loop", action="store_true", help="the flag to make LF/LB component generation optional")
    parser.add_argument(
        "-cstr",
        help="the configuration override string(s), format='section:key=val,...', example: 'GLOBAL:PRIMER_NUM_RETURN=500'",
        nargs="+",
        default=[],
    )
    parser.add_argument("-limit", help="the limit for the number of primer combinations", default=10, type=int)
    parser.add_argument("-proc", help="the number of processes", default=1, type=int)
    parser.add_argument("-out", type=FileType("w"), default="-", help="the output JSON file")
    return parser.parse_args(argv)


def main(argv):
    # parse args/conf
    args = parse_args(argv[1:])
    conf = parse_config(args.conf, args.cstr)

    # load FASTA
    with args.seqs as file:
        records = list(SeqIO.parse(file, "fasta"))

    if args.mode == "PCR":
        results = main_pcr(args, conf, records)
    elif args.mode == "LAMP":
        results = main_lamp(args, conf, records)

    # output assays
    oligo_calc = primer3.thermoanalysis.ThermoAnalysis(**conf["THERMO"])
    rc_keys = {PCR.parse_regex(PCR.regex)[-1][1], *itemgetter(2, 3, 6, 7)([ele[1] for ele in LAMP.parse_regex(LAMP.regex)])}
    with args.out as file:
        json.dump(
            dict(
                args={key: getattr(val, "name", val) for key, val in vars(args).items()},
                conf=conf,
                results=[
                    dict(
                        id=assay.id,
                        definition=assay.definition,
                        penalty=penalty,
                        sequences={
                            key: dict(
                                pos=assay.ccoors[key],
                                len=len(val),
                                Tm=oligo_calc.calc_tm(val),
                                dG5=oligo_calc.calc_end_stability(rc(val), rc(val)).dh,
                                dG3=oligo_calc.calc_end_stability(val, val).dh,
                                gc=sum((amb[ele].count("C") + amb[ele].count("G")) / len(amb[ele]) for ele in val) / len(val),
                                seq=rc(val) if key in rc_keys else val,
                            )
                            for key, val in assay.components.items()
                        },
                        **(dict(fip=assay.fip(), bip=assay.bip()) if args.mode == "LAMP" else dict()),
                    )
                    for assay, penalty in results
                ],
            ),
            fp=file,
            indent=True,
        )

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
