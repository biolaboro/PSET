#!/usr/bin/env python3

"""This module contains assay processing functions.
"""

import sys
from collections import OrderedDict
from csv import DictReader
from enum import Enum, auto
from functools import reduce
from itertools import product
from operator import mul
from re import findall, match, search, split, sub

from Bio.Data.IUPACData import ambiguous_dna_letters as dna_key
from Bio.Data.IUPACData import ambiguous_dna_values as dna_val
from Bio.Seq import Seq, reverse_complement
from Bio.SeqRecord import SeqRecord
from pset.util import iter_limit, iter_lines, sniff_lines

amb = {
    **{key: "".join(sorted(val)) for key, val in dna_val.items()},
    **{key.lower(): "".join(sorted(val)) for key, val in dna_val.items()}
}

dna = set("ACGTacgt")

trs = set(zip("ACGTACGTacgtacgt", "GTACgtacGTACgtac"))

trv = set(zip("AACGCTGTAACGCTGTaacgctgtaacgctgt", "CTGTAACGctgtaacgCTGTAACGctgtaacg"))

unk = set("NXnx")


class Assay(object):
    """The base assay object. Each subclass perfoms additional definition parsing.
    """

    def __init__(self, id, definition, targets):
        super().__init__()
        self.id = id
        self.definition = definition
        self.targets = targets
        self.components = OrderedDict()
        self.coordinates = OrderedDict()
        self.regex = ""
        self.delimiters = ()

    def __eq__(self, y):
        return self.val() == y.val()

    def __hash__(self):
        return hash(self.val())

    def __repr__(self):
        return str(self.val())

    def __getitem__(self, key):
        if key == "amplicon":
            result = self.amplicon()
        elif key == "camplicon":
            result = self.camplicon()
        else:
            result = self.components[key]
        return result

    @staticmethod
    def key():
        return "id", "definition", "targets", "type"

    def val(self):
        return self.id, self.definition, tuple(sorted(self.targets)), type(self).__name__

    def description(self):
        return f"[id={self.id}] [targets={self.targets}] [type={type(self).__name__}]"

    def adefinition(self):
        return search(r"(\[.+\])", self.definition).group(1)

    def amplicon(self):
        return "".join(ele for ele in search(r"\[(.+)\]", self.definition).group(1) if ele.isalpha())

    def camplicon(self):
        return "".join(ele for ele in self.definition if ele.isalpha())

    def context5(self):
        match = search(r"^([^[]+)\[", self.definition)
        return match.group(1) if match else ""

    def context3(self):
        match = search(r"\]([^[]+)$", self.definition)
        return match.group(1) if match else ""

    def conponent(self, key, ctx5=0, ctx3=0):
        start, end = self.coordinates[key]
        start, end = start - ctx5, end + ctx3
        start = 0 if start < 0 else start
        return self.camplicon()[start:end]

    def delim_coor(self):
        i = -1
        delims = OrderedDict()
        for e in self.definition:
            if e.isalpha():
                i += 1
            else:
                delims[i] = delims.get(i, "") + e
        return delims

    def delim_info(self):
        yield from (
            dict(zip(("open", "component", "close"), match))
            for match in findall(r"(.)\(\?P<([^<>]+)>[^)]+\)(.)", sub(r"\\(.)", r"\1", self.regex))
        )

    def delimify(self, qaln, saln):
        o = {ele["open"] for ele in self.delimiters}
        c = {ele["close"] for ele in self.delimiters}
        coors = self.delim_coor()
        p = -1
        aln = list(zip(qaln, saln))
        aln += [aln[-1]]
        if not self.context5():
            yield self.delimiters[0]["open"]
        for i in range(len(aln) - 1):
            p += aln[i][0] != "-"
            yield aln[i][1]
            for d in coors.get(p, ""):
                if (d in c and aln[i][0] != "-") or (d in o and aln[i + 1][0] != "-"):
                    yield d

    def records(self, *keys, expand=0):
        for key in keys:
            val = self[key]
            if expand:
                nexp = num_expansions(val)
                nexp = nexp if expand < 0 else min(expand, nexp)
                width = len(str(nexp))
                for idx, ele in enumerate(gen_expansions(val, limit=expand), start=1):
                    name = "{}-{:0{width}}".format(key, idx, width=width)
                    yield SeqRecord(Seq(ele), id=name, name=name, description=self.description())
            else:
                yield SeqRecord(Seq(val), id=key, name=key, description=self.description())

    def update(self):
        pos, off = 0, {}
        for idx, ele in enumerate(self.definition, start=1):
            pos += ele.isalpha()
            off[idx] = pos

        match = search(self.regex, self.definition)

        for key, val in match.groupdict().items():
            self.components[key] = "".join(ele for ele in val if ele.isalpha())
            self.coordinates[key] = tuple(map(off.get, match.span(key)))

        self.components = OrderedDict((ele["component"], self.components[ele["component"]]) for ele in self.delimiters)
        self.coordinates = OrderedDict((ele["component"], self.coordinates[ele["component"]]) for ele in self.delimiters)

    def is_range_valid(self, r):
        return 0 <= r[0] < r[1]

    def is_range_disjoint(self, r1, r2):
        return r1[0] < r1[1] < r2[0] < r2[1] or r2[0] < r2[1] < r1[0] < r1[1]

    def is_primer_hit(self, p1, p2, dist):
        rng1 = (p1.hit_start, p1.hit_end)
        rng2 = (p2.hit_start, p2.hit_end)
        return all((
            p1.query_strand == p2.query_strand,
            self.is_range_valid(rng1),
            self.is_range_valid(rng2),
            self.is_range_disjoint(rng1, rng2),
            abs(p2.hit_start - p1.hit_end + 1) <= dist
        ))

    def primer_hits(self, primer1, primer2, dist):
        yield from (ele for ele in product(primer1, primer2) if self.is_primer_hit(*ele, dist))

    @staticmethod
    def factory(id, definition, targets={1}):
        definition = definition.upper()
        assay = (
            Amplicon if match(Amplicon.regex, definition) else
            next((ele for ele in (Lamp, *Outer.__subclasses__()) if match(ele.regex, definition)), Outer)
        )
        return assay(id, definition, targets)


class Amplicon(Assay):
    regex = rf"^[{dna_key}]*\[(?P<amplicon>[{dna_key}]+)\][{dna_key}]*$"  # Mark Sanders pro-tip

    def __init__(self, id, definition, targets):
        super().__init__(id, definition, targets)
        self.regex = Amplicon.regex
        self.delimiters = list(self.delim_info())
        self.update()

    def hits(self, rows, distance):
        yield from ((row, ) for row in rows)


class Outer(Assay):
    regex = rf"^[{dna_key}]*\[(?P<outer1>[{dna_key}]+)\][{dna_key}]*\[(?P<outer2>[{dna_key}]+)\][{dna_key}]*$"

    def __init__(self, id, definition, targets):
        super().__init__(id, definition, targets)
        self.regex = Outer.regex
        self.delimiters = list(self.delim_info())
        self.update()

    def hits(self, rows, dist):
        outer1 = (row for row in rows if row.query_id == "outer1")
        outer2 = (row for row in rows if row.query_id == "outer2")
        yield from self.primer_hits(outer1, outer2, dist)


class Inner(Outer):
    regex = rf"[{dna_key}\[\]]*{{(?P<inner1>[{dna_key}\]]*)}}[{dna_key}]*{{(?P<inner2>[{dna_key}\[]*)}}[{dna_key}\[\]]*"

    def __init__(self, id, definition, targets):
        _definition = "".join([ele for ele in definition if ele not in "{}"])
        super().__init__(id, _definition, targets)
        self.regex = Inner.regex
        self.delimiters = (self.delimiters[0], *self.delim_info(), self.delimiters[-1])
        self.definition = definition
        self.update()

    def hits(self, rows, dist):
        inner1 = (row for row in rows if row.query_id == "inner1")
        inner2 = (row for row in rows if row.query_id == "inner2")
        for outer1, outer2 in super().hits(rows, dist):
            bounds1 = sorted((outer1.hit_start, outer1.hit_end, outer2.hit_start, outer2.hit_end))
            for inner1, inner2 in super().primer_hits(inner1, inner2, dist):
                bounds2 = sorted((inner1.hit_start, inner1.hit_end, inner2.hit_start, inner2.hit_end))
                if bounds1[0] <= bounds2[0] < bounds2[-1] <= bounds1[-1]:
                    yield outer1, inner1, inner2, outer2


class Lamp(Assay):
    regex = (
        rf"^[{dna_key}]*\[(?P<F3>[{dna_key}]+)\][{dna_key}]*"
        rf"\[(?P<F2>[{dna_key}]+)\][{dna_key}]*\<(?P<LF>[{dna_key}]+)\>[{dna_key}]*\[(?P<F1>[{dna_key}]+)\]"
        rf"[{dna_key}]*"
        rf"\[(?P<B1c>[{dna_key}]+)\][{dna_key}]*\<(?P<LB>[{dna_key}]+)\>[{dna_key}]*\[(?P<B2c>[{dna_key}]+)\]"
        rf"[{dna_key}]*\[(?P<B3c>[{dna_key}]+)\][{dna_key}]*$"
    )

    def __init__(self, id, definition, targets):
        super().__init__(id, definition, targets)
        self.regex = Lamp.regex
        self.delimiters = list(self.delim_info())
        self.update()

    def hits(self, rows, dist):
        F3 = (row for row in rows if row.query_id == "F3")
        B3c = (row for row in rows if row.query_id == "B3c")
        F2 = (row for row in rows if row.query_id == "F2")
        B2c = (row for row in rows if row.query_id == "B2c")
        LF = (row for row in rows if row.query_id == "LF")
        LB = (row for row in rows if row.query_id == "LB")
        F1 = (row for row in rows if row.query_id == "F1")
        B1c = (row for row in rows if row.query_id == "B1c")
        for h3, h3c in self.primer_hits(F3, B3c, dist):
            b3 = sorted((h3.hit_start, h3.hit_end, h3c.hit_start, h3c.hit_end))
            for h2, h2c in self.primer_hits(F2, B2c, dist):
                b2 = sorted((h2.hit_start, h2.hit_end, h2c.hit_start, h2c.hit_end))
                for hf, hb in self.primer_hits(LF, LB, dist):
                    bl = sorted((hf.hit_start, hf.hit_end, hb.hit_start, hb.hit_end))
                    for h1, h1c in self.primer_hits(F1, B1c, dist):
                        b1 = sorted((h1.hit_start, h1.hit_end, h1c.hit_start, h1c.hit_end))
                        if (
                            b3[1] <= b2[0] <=
                            b2[1] <= bl[0] <=
                            bl[1] <= b1[0] <=
                            b1[3] <= bl[2] <=
                            bl[3] <= b2[2] <=
                            b2[2] <= b2[3] <=
                            b2[3] <= b3[2]
                        ):
                            yield h3, h2, hf, h1, h1c, hb, h2c, h3c


class Probe(Outer):
    regex = rf"[{dna_key}]*\[[{dna_key}\]]*\((?P<probe>[{dna_key}\[\]]+)\)[{dna_key}\[]*\][{dna_key}]*"

    def __init__(self, id, definition, targets):
        _definition = "".join([ele for ele in definition if ele not in "()"])
        super().__init__(id, _definition, targets)
        self.regex = Probe.regex
        self.delimiters = (self.delimiters[0], *self.delim_info(), self.delimiters[-1])
        self.definition = definition
        self.update()

    def hits(self, rows, dist):
        probes = (row for row in rows if row.query_id == "probe")
        for outer1, outer2 in super().hits(rows, dist):
            bounds = sorted((outer1.hit_start, outer1.hit_end, outer2.hit_start, outer2.hit_end))
            for probe in probes:
                if bounds[0] <= probe.hit_start < probe.hit_end <= bounds[-1]:
                    yield outer1, probe, outer2


class AlignType(Enum):
    TRS = auto()
    TRV = auto()
    INS = auto()
    DEL = auto()
    DIS = auto()
    SIM = auto()
    GAP = auto()
    UNK = auto()


def align_type(ref, alt):
    if ref == "-" and alt == "-":
        return AlignType.GAP
    elif ref == "-":
        return AlignType.INS
    elif alt == "-":
        return AlignType.DEL
    elif alt in unk:
        return AlignType.UNK
    elif (ref, alt) in trs:
        return AlignType.TRS
    elif (ref, alt) in trv:
        return AlignType.TRV
    else:
        return AlignType.SIM if is_similar(ref, alt) else AlignType.DIS


def decode_btop(seq, btop, sbj=True):
    i = 0
    for ele in filter(len, split("([^0-9]{2})", btop)):
        if ele.isdigit():
            n = int(ele)
            yield seq[i:i + n]
            i += n
        elif ele[0] == "-":
            yield ele[sbj].lower()
        elif is_similar(*ele):
            yield ele[sbj]
            i += 1
        else:
            yield ele[sbj].lower()
            i += 1


def is_similar(ref, alt):
    return alt not in unk and bool({*amb.get(ref, ref)} & {*amb.get(alt, alt)})


def gen_expansions(seq, limit=None):
    indexes = [idx for idx, ele in enumerate(seq) if ele in amb]
    codes = [amb[seq[idx]] for idx in indexes]
    for combo in iter_limit(product(*codes), limit=limit):
        edits = dict(zip(indexes, combo))
        yield "".join(edits.get(idx, ele) for idx, ele in enumerate(seq))


def num_expansions(seq):
    return reduce(mul, (len(amb.get(ch, ch)) for ch in seq), 1)


def parse_assays(file, comment="#", sniff=True, **kwargs):
    kwargs["dialect"] = sniff_lines(file, comment=comment) if sniff else None
    reader = DictReader(iter_lines(file, comment=comment), **kwargs)
    reader.fieldnames = [field.lower().strip() for field in reader.fieldnames]
    for row in reader:
        id = row["id"].strip()
        definition = row["definition"].strip()
        targets = set(map(int, row.get("targets", "1").split(";")))
        yield Assay.factory(id, definition, targets)


def write_assays(records, file=sys.stdout, sep="\t"):
    print("id", "definition", "targets", sep=sep, file=file)
    for record in records:
        targets = ";".join(sorted(map(str, record.targets)))
        print(record.id, record.definition, targets, sep=sep, file=file)
