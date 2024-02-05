#!/usr/bin/env python3

"""This module contains assay processing functions.
"""

import json
import re
import sys
from csv import DictReader
from enum import Enum, auto
from functools import reduce
from itertools import chain, product
from operator import mul

from Bio.Data.IUPACData import ambiguous_dna_letters as dna_key
from Bio.Data.IUPACData import ambiguous_dna_values as dna_val
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from pset.util import iter_limit, iter_lines, sniff_lines

# dict of DNA codes to unambiguous representations (case-insensitive): "A" -> "A", "N" -> "ACGT", ...
amb = {**{key: "".join(sorted(val)) for key, val in dna_val.items()}, **{key.lower(): "".join(sorted(val)) for key, val in dna_val.items()}}

# set of canonical DNA codes
dna = set("ACGTacgt")

# set of transition mutation pairs
trs = set(zip("ACGTACGTacgtacgt", "GTACgtacGTACgtac"))

# set of transversion mutation pairs
trv = set(zip("AACGCTGTAACGCTGTaacgctgtaacgctgt", "CTGTAACGctgtaacgCTGTAACGctgtaacg"))

# set of unknown DNA codes
unk = set("NXnx")


class Assay(object):
    """The base assay object.

    Use the factory method.

    Attributes:
        definition (str): the string that delimits primers within an amplicon and implicitly sets assay type
        targets (set[Any]): the set of targeted taxonomic identifiers
        id (str): the identifier
        dcoors (dict[str, tuple[int, int]]): the dictionary mapping each assay component key to ranges on the definition
        ccoors (dict[str, tuple[int, int]]): the dictionary mapping each assay component key to ranges on the camplicon
        components (dict[str, str]): the dictionary mapping each assay component key to sequence
    """

    # the regex for extracting primers according to named groups from the definition
    regex = ""

    def __init__(self, definition, targets=None, id=None):
        """Initialized the assay.

        Args:
            definition (str): the assay definition
            targets (set[Any] or None): the set of targeted taxonomic identifiers
            id (str or None): the assay identifier
        """
        super().__init__()

        self.definition = definition
        self.targets = set(targets or ())
        self.id = id or ""

        match = re.search(self.regex, self.definition)
        groups = {key: val for key, val in match.groupdict().items() if val}

        pos, off = 0, {}
        for idx, ele in enumerate(self.definition, start=1):
            pos += ele.isalpha()
            off[idx] = pos

        self.dcoors = {key: match.span(key) for key, val in groups.items() if val}
        self.ccoors = {key: tuple(map(off.get, match.span(key))) for key, val in groups.items() if val}
        self.components = {key: "".join(filter(str.isalpha, val)) for key, val in groups.items() if val}

    def __eq__(self, y):
        return self.val() == y.val()

    def __hash__(self):
        return hash(self.val())

    def __repr__(self):
        return str(self.val())

    def __contains__(self, key):
        return key in self.components

    def __getitem__(self, key):
        return self.amplicon() if key == "amplicon" else self.components[key]

    def description(self):
        """Calculate a header for FASTA output.

        Returns:
            str: the description
        """
        return f"[id={self.id}] [targets={','.join(map(str, self.targets))}] [type={type(self).__name__}]"

    def acoors(self):
        coordinates = tuple(self.ccoors.values())
        return (coordinates[0][0], coordinates[-1][-1])

    def bcoors(self):
        i = -1
        delims = {}
        for e in self.definition:
            if e.isalpha():
                i += 1
            else:
                delims[i] = delims.get(i, "") + e
        return delims

    @staticmethod
    def parse_regex(regex):
        return re.findall(r"(.)\(\?P<([^<>]+)>[^)]+\)(.)", re.sub(r"\\(.)", r"\1", regex.replace(")?", ")")))

    def brackets(self):
        """Calculate the set of open and close brackets for each component based on the definition regex."""
        return [dict(zip(("open", "component", "close"), match)) for match in self.parse_regex(self.regex)]

    def embrace(self, qaln, saln):
        brackets = [ele for ele in self.brackets() if ele["component"] in self]
        o = {ele["open"] for ele in brackets}
        c = {ele["close"] for ele in brackets}
        coors = self.bcoors()
        p = -1
        aln = list(zip(qaln, saln))
        aln += [aln[-1]]
        if not self.definition[: self.acoors()[0]]:  # 5'-context
            yield brackets[0]["open"]
        for i in range(len(aln) - 1):
            p += aln[i][0] != "-"
            yield aln[i][1]
            for d in coors.get(p, ""):
                if (d in c and aln[i][0] != "-") or (d in o and aln[i + 1][0] != "-"):
                    yield d

    def amplicon(self):
        return self.camplicon()[slice(*self.acoors())]

    def camplicon(self):
        """Calculate the definition sequence without primer brackets.

        The camplicon is the contextualized amplicon, which is the amplicon with 5'/3'-context.

        Returns:
            str: the contextualized amplicon
        """
        return "".join(filter(str.isalpha, self.definition))

    def records(self, *keys, expand=0, context=(0, 0)):
        """Generate the sequence for the key specified.

        Args:
            keys: the component keys
            expand (int, optional): the number of unambiguous sequence expansions to generate, defaults to 0
            context (tuple[int], optional): The amount of 5'/3'-context to include. Defaults to (0, 0).

        Yields:
            (SeqRecord): the next sequence (expansion) of the corresponding key
        """
        keys = keys or self.components
        for key in keys:
            val = self[key]
            coor = self.acoors() if key == "amplicon" else self.ccoors[key]
            idx5 = coor[0] - context[0]
            idx5 = 0 if idx5 < 0 else idx5
            val = self.camplicon()[idx5 : coor[1] + context[1]]
            if expand > 0:
                nexp = num_expansions(val)
                width = len(str(nexp))
                for idx, ele in enumerate(gen_expansions(val, limit=expand), start=1):
                    name = "{}-{:0{width}}".format(key, idx, width=width)
                    yield SeqRecord(Seq(ele), id=name, name=name, description=self.description())
            else:
                yield SeqRecord(Seq(val), id=key, name=key, description=self.description())

    def val(self):
        """Return a tuple of values for calculating the hash of this object.

        Returns:
            tuple: the values
        """
        return self.id, tuple(sorted(self.targets)), self.definition

    def kprobes(self):
        return tuple(ele["component"] for ele in self.brackets() if ele["open"] == "(" and ele["close"] == ")")

    def kprimers(self):
        return tuple(ele["component"] for ele in self.brackets() if ele["open"] == "[" and ele["close"] == "]")

    def kflanks(self):
        keys = self.kprimers()
        return keys[0], keys[-1]

    @staticmethod
    def key():
        """Return a tuple of keys corresponding to the computed value tuple.

        Returns:
            tuple: the key names
        """
        return "id", "targets", "definition"

    @staticmethod
    def filter_hsps(hsps, query_id):
        yield from (hsp for hsp in hsps if hsp.query_id == query_id)

    @staticmethod
    def primer_hits(primers1, primers2, dist):
        """Generate all possible primer hits.

        Args:
            primers1 (list[HSP]): the list of primer HSP objects
            primers2 (list[HSP]): the list of primer HSP objects
            dist (Tuple[int,int]): the min/max distance range (inclusive)

        Yields:
            tuple[HSP]: the next combination of primers forming a hit
        """
        for p1, p2 in product(primers1, primers2):
            p1, p2 = (p1, p2) if p1.hit_end < p2.hit_start else (p2, p1)
            if p1.query_strand == p2.query_strand and dist[0] <= p2.hit_end - p1.hit_start <= dist[1]:
                yield p1, p2

    @staticmethod
    def from_json(path):
        """Load an assay from JSON.

        Args:
            path (str): the path to the JSON file

        Returns:
            Assay: the assay object
        """
        with open(path) as file:
            obj = json.load(file)
            return Assay.factory(obj["definition"], obj.get("targets"), obj.get("id"))

    @staticmethod
    def factory(definition, targets=None, id=None):
        """Create an assay object.

        Args:
            definition (str): the assay definition
            targets (set[Any] or None): the set of targeted taxonomic identifiers
            id (str): the assay identifier

        Returns:
            _type_: _description_
        """
        definition = definition.upper()
        assay = next((ele for ele in subclasses(Assay) if re.match(ele.regex, definition)))
        return assay(definition, targets, id)


class Oligo(Assay):
    """The oligo assay for testing a single primer/probe."""

    regex = rf"^[{dna_key}]*\[(?P<O>[{dna_key}]+)\][{dna_key}]*$"  # Mark Sanders pro-tip

    def __init__(self, definition, targets=None, id=None):
        super().__init__(definition, targets, id)

    def hits(self, hsps, **kwargs):
        # no need to check for arrangement, so just return each
        yield from ((hsp,) for hsp in hsps)


class PCR(Assay):
    """The assay for testing a PCR primer set."""

    regex = rf"^[{dna_key}]*" rf"\[(?P<F>[{dna_key}]+)\][{dna_key}]*(?:\((?P<P>[{dna_key}]+)\))?[{dna_key}]*\[(?P<R>[{dna_key}]+)\]" rf"[{dna_key}]*$"

    def __init__(self, definition, targets=None, id=None):
        super().__init__(definition, targets, id)

    def has_probe(self):
        return self.kprobes()[0] in self

    def hits(self, hsps, **kwargs):
        k1, k2 = self.kflanks()
        hits = self.primer_hits(self.filter_hsps(hsps, k1), self.filter_hsps(hsps, k2), dist=kwargs["dFR"])
        if self.has_probe():
            probes = list(self.filter_hsps(hsps, self.kprobes()[0]))
            hits = ((p1, pr, p2) for p1, p2 in hits for pr in probes if p1.hit_end <= pr.hit_start and pr.hit_end <= p2.hit_start)
        yield from hits


class LAMP(Assay):
    """The assay for testing a LAMP primer set."""

    regex = (
        rf"^[{dna_key}]*"
        rf"\[(?P<F3>[{dna_key}]+)\]"
        rf"[{dna_key}]*"
        rf"\[(?P<F2>[{dna_key}]+)\][{dna_key}]*(?:\((?P<LF>[{dna_key}]+)\))?[{dna_key}]*\[(?P<F1c>[{dna_key}]+)\]"
        rf"[{dna_key}]*"
        rf"\[(?P<B1c>[{dna_key}]+)\][{dna_key}]*(?:\((?P<LB>[{dna_key}]+)\))?[{dna_key}]*\[(?P<B2>[{dna_key}]+)\]"
        rf"[{dna_key}]*"
        rf"\[(?P<B3>[{dna_key}]+)\]"
        rf"[{dna_key}]*$"
    )

    def __init__(self, definition, targets=None, id=None):
        super().__init__(definition, targets, id)
        # decompose definition into subassays
        # [0] -> outer flanking primer set
        # [1] -> inner primer set near 5'-end
        # [2] -> inner primer set near 3'-end
        keys = self.kprimers()
        # definition indexes of each subassay primer set
        idxs = ((0, 5), (1, 2), (3, 4))
        # characters to strip from each subassay definition interprimer region
        excl = ("[]()", "[]", "[]")
        coor = self.dcoors
        self.subassays = []
        for i in range(3):
            i1, i2 = idxs[i]
            k1, k2 = keys[i1], keys[i2]
            subdefinition = (
                ""
                + definition[coor[k1][0] - 1 : coor[k1][1] + 1]
                + "".join(ele for ele in definition[coor[k1][1] + 1 : coor[k2][0] - 1] if ele not in excl[i])
                + definition[coor[k2][0] - 1 : coor[k2][1] + 1]
            )
            self.subassays.append(Assay.factory(subdefinition, targets, f"{self.id}-{k1}/{k2}"))

    def fip(self):
        _, k2, k3, __, ___, ____ = self.kprimers()
        return str(Seq(self[k3]).reverse_complement()) + self[k2]

    def bip(self):
        _, __, ___, k4, k5, ____ = self.kprimers()
        return self[k4] + str(Seq(self[k5]).reverse_complement())

    def ksubcomponents(self):
        keys2 = tuple(self.subassays[0].components)
        keys2 = (keys2[0], *self.subassays[1].components, *self.subassays[2].components, keys2[-1])
        return keys2

    def translate_hsps(self, hsps):
        keys1 = tuple(self.components)
        keys2 = self.ksubcomponents()
        for hsp in hsps:
            hsp.query_id = keys2[keys1.index(hsp.query_id)]
            yield hsp

    def hits(self, hsps, **kwargs):
        k1, k2, k3, k4, k5, k6 = self.kprimers()
        kp1, kp2 = self.kprobes()

        # [k1]-[k2]-[k3]-[k4]-[k5]-[k6]
        #      |-sub2--| |--sub3-|
        # |-----------sub1------------|

        # note: subassays [1] and [2] store optional loop as a probe
        # note: need to rename query_id from LAMP keys to subassay keys

        # procedure:
        # filter subhsps corresponding to the subassay
        # rename of each hsp query_id to subassay expectation
        # calculate subassay hits
        # undo rename

        # calculate hits to subassay 1
        keys = (k1, k6)
        subhsps = list(self.translate_hsps(chain.from_iterable((self.filter_hsps(hsps, key) for key in keys))))
        subhits1 = list(self.subassays[0].hits(subhsps, **kwargs))

        # for [1] and [2]:
        # len(ele) == 2 => no loop, otherwise check strand consistent with proper LAMP design

        # calculate hits to subassay 2
        keys = (k2, k3, kp1) if kp1 in self else (k2, k3)
        subhsps = list(self.translate_hsps(chain.from_iterable((self.filter_hsps(hsps, key) for key in keys))))
        subhits2 = list(self.subassays[1].hits(subhsps, **kwargs))
        subhits2 = [ele for ele in subhits2 if len(ele) == 2 or ele[1].query_strand == ele[2].query_strand]
        # calculate hits to subassay 3
        keys = (k4, k5, kp2) if kp2 in self else (k4, k5)
        subhsps = list(self.translate_hsps(chain.from_iterable((self.filter_hsps(hsps, key) for key in keys))))
        subhits3 = list(self.subassays[2].hits(subhsps, **kwargs))
        subhits3 = [ele for ele in subhits3 if len(ele) == 2 or ele[1].query_strand == ele[0].query_strand]

        # A     -> F3
        # B[0]  -> F2
        # B[-1] -> F1c
        # C[0]  -> B1c
        # C[-1] -> B2
        # D     -> B3
        dF3F2 = kwargs["dF3F2"]
        dF2F1c = kwargs["dF2F1c"]
        dF1cB1c = kwargs["dF1cB1c"]
        for A, D in subhits1:
            for B, C in product(subhits2, subhits3):
                B, C = (B, C) if B[-1].hit_end < C[0].hit_start else (C, B)
                if (
                    # check arrangement
                    A.hit_end < B[0].hit_start
                    and B[-1].hit_end < C[0].hit_start
                    and C[-1].hit_end < D.hit_start
                    # check distances
                    and dF3F2[0] <= B[0].hit_start - A.hit_end <= dF3F2[1]
                    and dF3F2[0] <= D.hit_start - C[-1].hit_end <= dF3F2[1]
                    and dF2F1c[0] <= B[-1].hit_start - B[0].hit_end <= dF2F1c[1]
                    and dF2F1c[0] <= C[-1].hit_start - C[0].hit_end <= dF2F1c[1]
                    and dF1cB1c[0] <= C[0].hit_start - B[-1].hit_end <= dF1cB1c[1]
                ):
                    hit = (*A, *B, *C, *D)
                    # translate back from subassay component key to assay component key
                    for key, ele in zip(self.components, hit):
                        ele.query_id = key
                    yield hit


class AlignType(Enum):
    """Enumerate alignment event types."""

    TRS = auto()  # transition
    TRV = auto()  # transversion
    INS = auto()  # insertion
    DEL = auto()  # deletion
    DIS = auto()  # insertion
    IDN = auto()  # identity
    SIM = auto()  # similar
    GAP = auto()  # gap
    UNK = auto()  # unknown


def subclasses(cls):
    """https://stackoverflow.com/a/3862957"""
    return set(cls.__subclasses__()).union((s for c in cls.__subclasses__() for s in subclasses(c)))


def align_type(ref, alt):
    """Calculate the alignment event type between two DNA codes.

    Args:
        ref (str): the reference allele
        alt (str): the alternate allele

    Returns:
        AlignType: the alignment event type
    """
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
    elif ref == alt:
        return AlignType.IDN
    else:
        return AlignType.SIM if is_similar(ref, alt) else AlignType.DIS


def decode_btop(seq, btop, sbj=True):
    """Decode BLAST trace-back operations alignment representation.

    https://www.ncbi.nlm.nih.gov/books/NBK569862/#_ckbk_Dispsrchcustom_Traceback_operations_

    Lower case letters on the subject indicate dissimilarity to the query.

    Args:
        seq (str): the query string
        btop (str): the BTOP string
        sbj (bool): the flag to recover the subject or query alignment

    Yields:
        str: the next aligned DNA letter
    """
    i = 0
    for ele in filter(len, re.split("([^0-9]{2})", btop)):
        if ele.isdigit():  # identity of length n
            n = int(ele)
            yield seq[i : i + n]
            i += n
        elif ele[0] == "-":  # gap
            yield ele[sbj].lower()
        elif is_similar(*ele):  # similarity
            yield ele[sbj]
            i += 1
        else:  # dissimilarity
            yield ele[sbj].lower()
            i += 1


def is_similar(ref, alt):
    """Calculate wheteher two aligned DNA codes are similar.

    Args:
        ref (str): the reference allele
        alt (str): the alternate allele

    Returns:
        bool: the result
    """
    return alt not in unk and bool({*amb.get(ref, ref)} & {*amb.get(alt, alt)})


def gen_expansions(seq, limit=-1):
    """Generate all sequences represented by one that contains ambiguous DNA codes.

    Args:
        seq (str): the DNA sequence
        limit (int, optional): The limit on the number of expansions, defaults to -1

    Yields:
        str: the next sequence expansion
    """
    indexes = [idx for idx, ele in enumerate(seq) if ele in amb]
    codes = [amb[seq[idx]] for idx in indexes]
    for combo in iter_limit(product(*codes), limit=limit):
        edits = dict(zip(indexes, combo))
        yield "".join(edits.get(idx, ele) for idx, ele in enumerate(seq))


def num_expansions(seq):
    """Calculate the number of sequences represented by one that contains ambiguous DNA codes.

    Args:
        seq (str): the DNA sequence

    Returns:
        int: the number of expansions
    """
    return reduce(mul, (len(amb.get(ch, ch)) for ch in seq), 1)


def parse_assays(file, context=(0, 0), comment="#", sniff=True, target_type=str, **kwargs):
    """Parse assays from a tabular file.

    Args:
        file (file): the file handle
        context (tuple[int], optional): the amount of 5'/3'-context to include. Defaults to (0, 0).
        comment (str, optional): the comment line character. Defaults to "#".
        sniff (bool, optional): the flag to sniff the tabular file format dialect. Defaults to True.
        target_type (func): the function that converts each target value. Defaults to int(val).

    Yields:
        Assay: the next assay
    """
    kwargs["dialect"] = sniff_lines(file, comment=comment) if sniff else None
    reader = DictReader(iter_lines(file, comment=comment), **kwargs)
    reader.fieldnames = [field.lower().strip() for field in reader.fieldnames]
    for row in reader:
        definition = row["definition"].strip()
        tokens = re.split(r"[\[\]]", definition)
        idx5 = len(tokens[0]) - context[0]
        idx5 = 0 if idx5 < 0 else idx5
        definition = tokens[0][idx5:] + definition[len(tokens[0]) : len(definition) - len(tokens[-1])] + tokens[-1][: context[1]]
        targets = set(map(target_type, row["targets"].split(";"))) if "targets" in row else ()
        id = row.get("id")
        yield Assay.factory(definition, targets, id)


def write_assays(records, file=sys.stdout, sep="\t"):
    """Write the delimited assay text file that contains "id", "definition", and "target" columns.

    Args:
        records (iterable[Assay]): the assay objects
        file (file, optional): the file handle, defaults to sys.stdout
        sep (str, optional): the tabular file separator, defaults to "\t"
    """
    print(*Assay.key(), sep=sep, file=file)
    for record in records:
        val = record.val()
        print(val[0], ";".join(val[1]), val[2], sep=sep, file=file)


def parse_args(argv):
    from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, FileType

    parser = ArgumentParser(description="assay util...", formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument("file", help="the assay file", type=FileType())
    choices = ("records", "coordinates")
    parser.add_argument("-command", choices=choices, default=choices[0])
    parser.add_argument("-ids", help="the assay identifiers to process, default includes all", nargs="+", default=[])
    parser.add_argument("-keys", help="the assay component keys, defaults to all", nargs="+", default=())
    parser.add_argument("-context", help="the amount of 5'/3'-context to include", type=str, default="0,0")
    parser.add_argument("-expand", help="the number of ambiguous sequence expansions", type=int, default=-1)
    return parser.parse_args(argv)


def main(argv):
    from Bio import SeqIO

    args = parse_args(argv[1:])
    context = tuple(map(int, args.context.split(",")))

    ids = set(args.ids)
    with args.file as file:
        assays = list(ele for ele in parse_assays(file, context=context, sniff=False, delimiter="\t") if not ids or ele.id in ids)

    if args.command == "records":
        for assay in assays:
            SeqIO.write(assay.records(*args.keys, expand=args.expand, context=context), sys.stdout, "fasta")
    elif args.command == "coordinates":
        for assay in assays:
            keys = (key for key in args.keys if key in assay.components) if args.keys else assay.components
            for key in keys:
                print(key, *assay.ccoors[key], sep="\t")


if __name__ == "__main__":
    sys.exit(main(sys.argv))
