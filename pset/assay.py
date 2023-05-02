#!/usr/bin/env python3

"""This module contains assay processing functions.
"""

import json
import re
import sys
from collections import OrderedDict
from csv import DictReader
from enum import Enum, auto
from functools import reduce
from itertools import product
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

    Each subclass perfoms additional definition parsing.
    Use the factory method...

    Attributes:
        id (str): the identifier
        definition (str): the string that delimits primers within an amplicon and implicitly sets assay type
        targets (set[int]): the set of integers for corresponding to the identifiers of taxonomic targets
        components (dict[str, str]): the dictionary mapping each assay component key to sequence
        coordinates (dict[str, tuple[int, int]]): the dictionary mapping each assay component key to ranges on the camplicon
        primers (list[str]): the list of primers keys
        flanks (list[str]): the primers flanking the amplicon
        regex (str): the regex for extracting primers according to named groups from the definition
        delimiters (list[dict[str, str]]):
            the list of dictionaries where each one has a "component", "open", and "close" entry
            where the first specifies the component name and the latter two specify the delimiter characters
    """

    def __init__(self, id, definition, targets):
        """Initialized the assay.

        Call update to set internal attributes.
        Preferred initialization method is to call factory.

        Args:
            id (str): the assay identifier
            definition (str): the assay definition
            targets (set[int]): the assay targets as the set of integer taxonomic identifiers
        """
        super().__init__()
        self.id = id
        self.definition = definition
        self.targets = targets
        self.components = OrderedDict()
        self.coordinates = OrderedDict()
        self.primers = ()
        self.flanks = ()
        self.regex = ""
        self.delimiters = ()

    def __eq__(self, y):
        return self.val() == y.val()

    def __hash__(self):
        return hash(self.val())

    def __repr__(self):
        return str(self.val())

    def __getitem__(self, key):
        return self.components[key]

    @staticmethod
    def key():
        """Return a tuple of keys corresponding to the computed value tuple.

        Returns:
            tuple: the key names
        """
        return "id", "definition", "targets", "type"

    def val(self):
        """Return a tuple of values for calculating the hash of this object.

        Returns:
            tuple: the values
        """
        return self.id, self.definition, tuple(sorted(self.targets)), type(self).__name__

    def description(self):
        """Calculate a header for FASTA output.

        Returns:
            str: the description
        """
        return f"[id={self.id}] [targets={self.targets}] [type={type(self).__name__}]"

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
        """Calculate the set of open and close delimiters for each component based on the definition regex.

        Yields:
            dict[tuple]: the next entry that includes the "component" key together with "open" and "close" delimiters
        """
        yield from (
            dict(zip(("open", "component", "close"), match))
            for match in re.findall(r"(.)\(\?P<([^<>]+)>[^)]+\)(.)", re.sub(r"\\(.)", r"\1", self.regex))
        )

    def delimify(self, qaln, saln):
        o = {ele["open"] for ele in self.delimiters}
        c = {ele["close"] for ele in self.delimiters}
        coors = self.delim_coor()
        p = -1
        aln = list(zip(qaln, saln))
        aln += [aln[-1]]
        if not self.definition[: self.coordinates["amplicon"][0]]:  # 5'-context
            yield self.delimiters[0]["open"]
        for i in range(len(aln) - 1):
            p += aln[i][0] != "-"
            yield aln[i][1]
            for d in coors.get(p, ""):
                if (d in c and aln[i][0] != "-") or (d in o and aln[i + 1][0] != "-"):
                    yield d

    def camplicon(self):
        """Calculate the definition sequence without primer delimiters.

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
        for key in keys:
            val = self[key]
            coor = self.coordinates[key]
            idx5 = coor[0] - context[0]
            idx5 = 0 if idx5 < 0 else idx5
            val = self.camplicon()[idx5 : coor[1] + context[1]]
            if expand:
                nexp = num_expansions(val)
                width = len(str(nexp))
                for idx, ele in enumerate(gen_expansions(val, limit=expand), start=1):
                    name = "{}-{:0{width}}".format(key, idx, width=width)
                    yield SeqRecord(Seq(ele), id=name, name=name, description=self.description())
            else:
                yield SeqRecord(Seq(val), id=key, name=key, description=self.description())

    def update(self):
        """Set components and coordinates based on new definition.

        Parses the definition based on the regex to update the components and coordinates attributes.
        """
        pos, off = 0, OrderedDict()
        for idx, ele in enumerate(self.definition, start=1):
            pos += ele.isalpha()
            off[idx] = pos

        match = re.search(self.regex, self.definition)
        for key, val in match.groupdict().items():
            self.coordinates[key] = tuple(map(off.get, match.span(key)))
            self.components[key] = "".join(filter(str.isalpha, val))

        self.components = OrderedDict((ele["component"], self.components[ele["component"]]) for ele in self.delimiters)
        self.coordinates = OrderedDict((ele["component"], self.coordinates[ele["component"]]) for ele in self.delimiters)

        self.primers = tuple(self.coordinates)
        self.flanks = (self.primers[0], self.primers[-1]) if len(self.primers) > 1 else ()

        coordinates = tuple(self.coordinates.values())
        self.coordinates["amplicon"] = (coordinates[0][0], coordinates[-1][-1])
        self.components["amplicon"] = self.camplicon()[slice(*self.coordinates["amplicon"])]

    def is_range_valid(self, r):
        """Test whether the range is valid.

        Args:
            r (int): the range

        Returns:
            bool: the result
        """
        return 0 <= r[0] < r[1]

    def is_range_disjoint(self, r1, r2):
        """Test whether the ranges are disjoint.

        Args:
            r1 (int): the range
            r2 (int): the range

        Returns:
            bool: the result
        """
        return r1[0] < r1[1] < r2[0] < r2[1] or r2[0] < r2[1] < r1[0] < r1[1]

    def is_primer_hit(self, p1, p2, dist):
        """Test whether the primer alignment ranges are correctly arranged.

        Note: this is written with respect to results from the FASTA alignment suite.

        Args:
            p1 (HSP): the primer HSP
            p2 (HSP): the prtime HSP
            dist (int): the maximum allowed distance between the primers

        Returns:
            bool: the result
        """
        rng1 = (p1.hit_start, p1.hit_end)
        rng2 = (p2.hit_start, p2.hit_end)
        return all(
            (
                p1.query_strand == p2.query_strand,
                self.is_range_valid(rng1),
                self.is_range_valid(rng2),
                self.is_range_disjoint(rng1, rng2),
                abs(p2.hit_start - p1.hit_end + 1) <= dist,
            )
        )

    def primer_hits(self, primer1, primer2, dist):
        """Generate all possible primer hits.

        Args:
            primer1 (list[HSP]): the list of primer HSP objects
            primer2 (list[HSP]): the list of primer HSP objects
            dist (int): the maximum allowed distance between the primers

        Yields:
            tuple[HSP, ...]: the next combination of primers forming a hit
        """
        yield from (ele for ele in product(primer1, primer2) if self.is_primer_hit(*ele, dist))

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
            return Assay.factory(obj["id"], obj["definition"], set(obj.get("targets", [1])))

    @staticmethod
    def factory(id, definition, targets={}):
        """Create an assay object.

        Args:
            id (str): the assay identifier
            definition (str): the assay definition
            targets (set[int]): the assay targets as the set of integer taxonomic identifiers, defaults to {1}

        Returns:
            _type_: _description_
        """
        definition = definition.upper()
        assay = (
            Amplicon
            if re.match(Amplicon.regex, definition)
            else next((ele for ele in Primer.__subclasses__() if re.match(ele.regex, definition)), Primer)
        )
        return assay(id, definition, targets)


class Amplicon(Assay):
    """The amplicon assay for testing alignments to the amplicon."""

    regex = rf"^[{dna_key}]*\[(?P<amplicon>[{dna_key}]+)\][{dna_key}]*$"  # Mark Sanders pro-tip

    def __init__(self, id, definition, targets):
        super().__init__(id, definition, targets)
        self.regex = Amplicon.regex
        self.delimiters = list(self.delim_info())
        self.update()

    def hits(self, rows, **kwargs):
        # no need to check for arrangement, so just return each
        yield from ((row,) for row in rows)


class Primer(Assay):
    """The assay for testing a primer set."""

    regex = rf"^[{dna_key}]*\[(?P<primer1>[{dna_key}]+)\][{dna_key}]*\[(?P<primer2>[{dna_key}]+)\][{dna_key}]*$"

    def __init__(self, id, definition, targets):
        super().__init__(id, definition, targets)
        self.regex = Primer.regex
        self.delimiters = list(self.delim_info())
        self.update()

    def hits(self, rows, **kwargs):
        primer1 = (row for row in rows if row.query_id == "primer1")
        primer2 = (row for row in rows if row.query_id == "primer2")
        yield from self.primer_hits(primer1, primer2, kwargs["dist"])


class Probe(Primer):
    """The assay for testing a primer set with a probe"""

    regex = rf"[{dna_key}]*\[[{dna_key}\]]*\((?P<probe>[{dna_key}\[\]]+)\)[{dna_key}\[]*\][{dna_key}]*"

    def __init__(self, id, definition, targets):
        _definition = "".join([ele for ele in definition if ele not in "()"])
        super().__init__(id, _definition, targets)
        self.regex = Probe.regex
        self.delimiters = (self.delimiters[0], *self.delim_info(), self.delimiters[-1])
        self.definition = definition
        self.update()

    def hits(self, rows, **kwargs):
        probes = (row for row in rows if row.query_id == "probe")
        # get hits based on primer hits
        for primer1, primer2 in super().hits(rows, **kwargs):
            # check whether each probe falls within the primers
            bounds = sorted((primer1.hit_start, primer1.hit_end, primer2.hit_start, primer2.hit_end))
            for probe in probes:
                if bounds[0] <= probe.hit_start < probe.hit_end <= bounds[-1]:
                    yield primer1, probe, primer2


class Nested(Primer):
    """The assay for testing nested primer sets."""

    regex = rf"[{dna_key}\[\]]*{{(?P<nested1>[{dna_key}\]]*)}}[{dna_key}]*{{(?P<nested2>[{dna_key}\[]*)}}[{dna_key}\[\]]*"

    def __init__(self, id, definition, targets):
        _definition = "".join([ele for ele in definition if ele not in "{}"])
        super().__init__(id, _definition, targets)
        self.regex = Nested.regex
        self.delimiters = (self.delimiters[0], *self.delim_info(), self.delimiters[-1])
        self.definition = definition
        self.update()

    def hits(self, rows, **kwargs):
        nested1 = (row for row in rows if row.query_id == "nested1")
        nested2 = (row for row in rows if row.query_id == "nested2")
        # get hits based on primer primer hits
        for primer1, primer2 in super().hits(rows, **kwargs):
            bounds1 = sorted((primer1.hit_start, primer1.hit_end, primer2.hit_start, primer2.hit_end))
            # get hits based on nested primer hits
            for nested1, nested2 in super().primer_hits(nested1, nested2, **kwargs):
                # check whether nested primers fall within the primer primers
                bounds2 = sorted((nested1.hit_start, nested1.hit_end, nested2.hit_start, nested2.hit_end))
                if bounds1[0] <= bounds2[0] < bounds2[-1] <= bounds1[-1]:
                    yield primer1, nested1, nested2, primer2


class AlignType(Enum):
    """Enumerate alignment event types."""

    TRS = auto()  # transition
    TRV = auto()  # transversion
    INS = auto()  # insertion
    DEL = auto()  # deletion
    DIS = auto()  # insertion
    SIM = auto()  # similar
    GAP = auto()  # gap
    UNK = auto()  # unknown


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


def parse_assays(file, context=(0, 0), comment="#", sniff=True, **kwargs):
    """_summary_

    Args:
        file (file): the file handle
        context (tuple[int], optional): The amount of 5'/3'-context to include. Defaults to (0, 0).
        comment (str, optional): The comment line character. Defaults to "#".
        sniff (bool, optional): The flag to sniff the tabular file format dialect. Defaults to True.

    Yields:
        Assay: the next assay
    """
    kwargs["dialect"] = sniff_lines(file, comment=comment) if sniff else None
    reader = DictReader(iter_lines(file, comment=comment), **kwargs)
    reader.fieldnames = [field.lower().strip() for field in reader.fieldnames]
    for row in reader:
        id = row["id"].strip()
        definition = row["definition"].strip()
        tokens = re.split(r"[\[\]]", definition)
        idx5 = len(tokens[0]) - context[0]
        idx5 = 0 if idx5 < 0 else idx5
        definition = tokens[0][idx5:] + definition[len(tokens[0]) : len(definition) - len(tokens[-1])] + tokens[-1][: context[1]]
        targets = row.get("targets")
        targets = set(map(int, row.get("targets", "").split(";"))) if targets else {}
        yield Assay.factory(id, definition, targets)


def write_assays(records, file=sys.stdout, sep="\t"):
    """Write the delimited assay text file that contains "id", "definition", and "target" columns.

    Args:
        records (iterable[Assay]): the assay objects
        file (file, optional): the file handle, defaults to sys.stdout
        sep (str, optional): the tabular file separator, defaults to "\t"
    """
    print("id", "definition", "targets", sep=sep, file=file)
    for record in records:
        targets = ";".join(sorted(map(str, record.targets)))
        print(record.id, record.definition, targets, sep=sep, file=file)


def parse_args(argv):
    from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, FileType

    parser = ArgumentParser(description="assay util...", formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument("file", help="the assay file", type=FileType())
    parser.add_argument("-ids", help="the assay identifiers to process, default includes all", nargs="+", default=[])
    parser.add_argument("-keys", help="the assay component keys", nargs="+")
    parser.add_argument("-context", help="the amount of 5'/3'-context to include", type=str, default="0,0")
    subparsers = parser.add_subparsers(help="sub-program", dest="command")
    subparser = subparsers.add_parser("records", help="the program to get assay sequences")
    subparser.add_argument("-expand", help="the number of ambiguous sequence expansions", type=int, default=-1)
    subparser = subparsers.add_parser("coordinates", help="the program to output assay component coordinates")
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
            keys = (key for key in args.keys if key in assay.components) if args.keys else assay.components
            SeqIO.write(assay.records(*keys, expand=args.expand, context=context), sys.stdout, "fasta")
    elif args.command == "coordinates":
        for assay in assays:
            keys = (key for key in args.keys if key in assay.components) if args.keys else assay.components
            for key in keys:
                print(key, *assay.coordinates[key], sep="\t")


if __name__ == "__main__":
    sys.exit(main(sys.argv))
