#!/usr/bin/env python3

import os
import unittest
from collections import defaultdict
from operator import attrgetter
from pathlib import Path
from subprocess import PIPE, Popen, check_call
from tempfile import NamedTemporaryFile

from Bio import SearchIO, SeqIO

from pset.assay import (
    Assay,
    decode_btop,
    gen_expansions,
    is_similar,
    num_expansions,
    parse_assays,
)
from pset.util import fields_8CB

root = Path(__file__).parent / "data"


def perfect_hits(assay, subject):
    hsps = defaultdict(list)
    with NamedTemporaryFile("w") as file:
        print(subject, file=file, flush=True)
        cmd = ("glsearch36", "-n", "-m", "8CB", "@", file.name)
        with Popen(cmd, stdout=PIPE, stdin=PIPE, universal_newlines=True) as proc:
            with proc.stdin as stdin:
                SeqIO.write(assay.records(*assay.components), stdin, "fasta")
            with proc.stdout as stdout:
                for result in SearchIO.parse(stdout, "blast-tab", fields=fields_8CB, comments=True):
                    for hit in result.hits:
                        for hsp in hit.hsps:
                            if hsp.ident_pct == 100:
                                hsps[hsp.hit_id].append(hsp)

    return next(assay.hits(next(iter(hsps.values()), []), dist=1000), None)


class Test(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        os.chdir(root)
        check_call(root / "setup.sh")

    def test_factory(self):
        print("\n")
        print(Assay.factory("assay", "[A]", [1]).components, sep="\n")
        print(Assay.factory("assay", "N[A]", [1]).components, sep="\n")
        print(Assay.factory("assay", "[A]N", [1]).components, sep="\n")
        print(Assay.factory("assay", "N[A]N", [1]).components, sep="\n")
        print(Assay.factory("assay", "[A][T]", [1]).components, sep="\n")
        print(Assay.factory("assay", "N[A][T]", [1]).components, sep="\n")
        print(Assay.factory("assay", "[A]N[T]", [1]).components, sep="\n")
        print(Assay.factory("assay", "[A][T]N", [1]).components, sep="\n")
        print(Assay.factory("assay", "N[A]N[T]N", [1]).components, sep="\n")
        print(Assay.factory("assay", "[A(T]CG[T)A]", [1]).components, sep="\n")
        print(Assay.factory("assay", "N[A(T]CG[T)A]", [1]).components, sep="\n")
        print(Assay.factory("assay", "N[A(T]CG[T)A]N", [1]).components, sep="\n")
        print(Assay.factory("assay", "[A{AA]C}{G[GGG}T]", [1]).components, sep="\n")
        print(Assay.factory("assay", "N[A{AA]C}{G[GGG}T]", [1]).components, sep="\n")
        print(Assay.factory("assay", "[A{AA]C}{G[GGG}T]N", [1]).components, sep="\n")
        print(Assay.factory("assay", "N[A{AA]C}{G[GGG}T]N", [1]).components, sep="\n")

    def test_amplicon(self):
        self.assertEqual(Assay.factory("assay", "[A{AA]C}{G[GGG}T]").components["amplicon"], "AAACGGGGT")
        self.assertEqual(Assay.factory("assay", "N[A{AA]C}{G[GGG}T]N").components["amplicon"], "AAACGGGGT")

    def test_num_expansions(self):
        """Test expansion number calculation..."""
        seq = "ACGTWSMKRYBDHVN"
        num = num_expansions(seq)
        ans = 1 * 1 * 1 * 1 * 2 * 2 * 2 * 2 * 2 * 2 * 3 * 3 * 3 * 3 * 4

        self.assertEqual(num, ans)
        self.assertEqual(num, len(list(gen_expansions(seq))))

    def test_gen_expansions(self):
        seq = "AWBN"
        ans = [
            "AACA",
            "AACC",
            "AACG",
            "AACT",
            "AAGA",
            "AAGC",
            "AAGG",
            "AAGT",
            "AATA",
            "AATC",
            "AATG",
            "AATT",
            "ATCA",
            "ATCC",
            "ATCG",
            "ATCT",
            "ATGA",
            "ATGC",
            "ATGG",
            "ATGT",
            "ATTA",
            "ATTC",
            "ATTG",
            "ATTT",
        ]
        self.assertListEqual(list(gen_expansions(seq)), ans)

    def test_is_similar(self):
        self.assertTrue(all(is_similar(*ele) for ele in zip("ACGT", "ACGT")))
        self.assertTrue(all(is_similar(*ele) for ele in zip("aaAA", "aAaA")))
        self.assertTrue(all(is_similar(*ele) for ele in zip("xxXXnnNN", "AaAaAaA")))
        self.assertTrue(all(not is_similar(*ele) for ele in zip("XNXNxnxn", "XNxnXNxn")))
        self.assertTrue(all(not is_similar(*ele) for ele in zip("AaAaAaA", "xxXXnnNN")))
        self.assertTrue(all(not is_similar(*ele) for ele in zip("acgtACGT", "TGCAtgca")))

    def test_invalid_definition(self):
        val = [
            "[][]",
            "[]()[]",
            "[]{}{}[]",
            "[GATTACA",
            "GATTACA]",
            "[GATTACA[]",
            "[GATTACA]]",
            "[GATTACA]CAT[0123456789]",
            "[GATTACA[CAT]GATTACA]",
            "(GATTACA)C[ATG]A[TTACA]",
            "[GATTACA]C[AT(G]ATTACA)",
            "{GAT}T[ACA]CATG[AT]T{ACA}",
            "{GATT[A}CA]CATG[A{T]TACA}",
        ]
        for idx, ele in enumerate(val, start=1):
            with self.assertRaises(Exception):
                print(ele, Assay.factory(f"assay-{idx}", ele, [1]))

    def test_hits(self):
        #                    5'<-> 3'
        # CCCCGCCGCGGGCCGGCGCG <-> CGCGCCGGCCCGCGGCGGGG
        # TTTGGGTGTGGTGTGGTGGT <-> ACCACCACACCACACCCAAA
        assay = Assay.factory("assay", "[CCCCGCCGCGGGCCGGCGCG]AAAAAAAAAATTATATTTATATATTATATTAAAAAAAAAA[TTTGGGTGTGGTGTGGTGGT]")
        self.assertIsNotNone(perfect_hits(assay, ">sub\nCCCCGCCGCGGGCCGGCGCGAAAAAAAAAATTATATTTATATATTATATTAAAAAAAAAATTTGGGTGTGGTGTGGTGGT\n"))
        self.assertIsNotNone(perfect_hits(assay, ">sub\nACCACCACACCACACCCAAATTTTTTTTTTAATATAATATATAAATATAATTTTTTTTTTCGCGCCGGCCCGCGGCGGGG\n"))
        self.assertIsNone(perfect_hits(assay, ">sub\nACCACCACACCACACCCAAAAAAAAAAAAATTATATTTATATATTATATTAAAAAAAAAACCCCGCCGCGGGCCGGCGCG\n"))
        self.assertIsNone(perfect_hits(assay, ">sub\nCGCGCCGGCCCGCGGCGGGGTTTTTTTTTTAATATAATATATAAATATAATTTTTTTTTTTTTGGGTGTGGTGTGGTGGT\n"))
        assay = Assay.factory("assay", "[CCCCGCCGCGGGCCGGCGCG]AAAAAAAAAA(TTATATTTATATATTATATT)AAAAAAAAAA[TTTGGGTGTGGTGTGGTGGT]")
        self.assertIsNotNone(perfect_hits(assay, ">sub\nCCCCGCCGCGGGCCGGCGCGAAAAAAAAAATTATATTTATATATTATATTAAAAAAAAAATTTGGGTGTGGTGTGGTGGT\n"))
        self.assertIsNotNone(perfect_hits(assay, ">sub\nACCACCACACCACACCCAAATTTTTTTTTTAATATAATATATAAATATAATTTTTTTTTTCGCGCCGGCCCGCGGCGGGG\n"))
        self.assertIsNone(perfect_hits(assay, ">sub\nTTTGGGTGTGGTGTGGTGGTAAAAAAAAAATTATATTTATATATTATATTAAAAAAAAAACGCGCCGGCCCGCGGCGGGG\n"))
        self.assertIsNone(perfect_hits(assay, ">sub\nCGCGCCGGCCCGCGGCGGGGTTTTTTTTTTAATATAATATATAAATATAATTTTTTTTTTTTTGGGTGTGGTGTGGTGGT\n"))
        self.assertIsNotNone(perfect_hits(assay, ">sub\nCCCCGCCGCGGGCCGGCGCGAAAAAAAAAAAATATAATATATAAATATAAAAAAAAAAAATTTGGGTGTGGTGTGGTGGT\n"))
        self.assertIsNotNone(perfect_hits(assay, ">sub\nACCACCACACCACACCCAAATTTTTTTTTTTTATATTTATATATTATATTTTTTTTTTTTCGCGCCGGCCCGCGGCGGGG\n"))
        self.assertIsNone(perfect_hits(assay, ">sub\nACCACCACACCACACCCAAAAAAAAAAAAAAATATAATATATAAATATAAAAAAAAAAAACCCCGCCGCGGGCCGGCGCG\n"))
        self.assertIsNone(perfect_hits(assay, ">sub\nCCCCGCCGCGGGCCGGCGCGTTTTTTTTTTTTATATTTATATATTATATTTTTTTTTTTTACCACCACACCACACCCAAA\n"))
        # from brad
        assay = Assay.factory(
            "assay",
            "CAGCTG[TGCGTCGGCAAACCAATGCT]ATTGAATCACTAGAAGGTCGAGTAACAACTCTTGA(GGCCAGCTTAAAACCCGTTC)AAGACATGGCAAAGACCATATCATCCCTGAATCGCAGCTGT[GCCGAAATGGTTGCAAAATACG]ACCTAC",
        )
        self.assertIsNotNone(
            perfect_hits(
                assay,
                ">sub\nTGCGTCGGCAAACCAATGCTATTGCAGCAGGATAGGACTTATAGACATCATGGACCCGTGAGGCCAGCTTAAAACCCGTTCAAGACATGGCAAAGACCATATCATCCCTGAATCGCAGCTGTGCCGAAATGGTTGCAAAATACGACCTAC\n",
            )
        )
        self.assertIsNotNone(
            perfect_hits(
                assay,
                ">sub\nCAGCTGTGCGTCGGCAAACCAATGCTATTCAAGACAACACGTAAAAGTGATACAACTCTTGAGGCCAGCTTAAAACCCGTTCAAGACCACTGGGCGAGCAACTGCCACGCCGAAATGGTTGCAAAATACGACCTAC\n",
            )
        )
        self.assertIsNone(
            perfect_hits(
                assay,
                ">sub\nCAGCTGTGCGTCGGGGCAATGCTATTGAATCACTAGAAGGTCGAGTAACAACTCTTGAGGCCAGCTTAAAACACTTTCAAGACATGGCAAAGACCATATCATCCCTGAATCGCAGCTGTGCCTTTCAAAATACGACCTAC\n",
            )
        )
        self.assertIsNone(
            perfect_hits(
                assay,
                ">sub\n CAGCTGTGCGTCGGCAGGGATGCTATTGAATCACTAGAAGGTCGAGTAACAACTCTTGAGGCCAGCTTAAAACGGTTCAAGACATGGCAAAGACCATATCATCCCTGAATCGCAGCTGTGCCGAAATGCTCCAAAATACGACCTAC\n",
            )
        )

    def test_decode_btop(self):
        self.assertEqual("".join(decode_btop("G", "1")), "G")
        self.assertEqual("".join(decode_btop("GATTACA", "7")), "GATTACA")
        self.assertEqual("".join(decode_btop("GATTACA", "GC6")), "cATTACA")
        self.assertEqual("".join(decode_btop("GATTACA", "2TATA3")), "GAaaACA")
        self.assertEqual("".join(decode_btop("GATTACA", "6AC")), "GATTACc")
        self.assertEqual("".join(decode_btop("GAACA", "2-T-T3")), "GAttACA")
        self.assertEqual("".join(decode_btop("GATTACA", "2T-T-3")), "GA--ACA")
        self.assertEqual("".join(decode_btop("GATACA", "2T--T3")), "GA-tACA")
        self.assertEqual("".join(decode_btop("GATTANA", "5NC1")), "GATTACA")
        self.assertEqual("".join(decode_btop("GATTACA", "5CN1")), "GATTAnA")

    def test_decode_btop_blastn(self):
        queries = SeqIO.to_dict(SeqIO.parse(root.joinpath("qry.fasta"), "fasta"))
        fields = ["qaccver", "saccver", "qstart", "qend", "qseq", "sseq", "length", "btop"]
        print()
        with root.joinpath("blastn.tsv").open() as file:
            for rec in SearchIO.parse(file, "blast-tab", fields=fields, comments=True):
                for hit in rec.hits:
                    for hsp in hit.hsps:
                        qry = str(queries[hsp.query_id].seq)[hsp.query_start : hsp.query_end]
                        qaln = str(hsp.query.seq)
                        saln = str(hsp.hit.seq)
                        print(qaln)
                        print("".join(decode_btop(qry, hsp.btop, sbj=False)).upper())
                        print("".join(decode_btop(qry, hsp.btop)).upper())
                        print(saln)
                        print()
                        self.assertEqual(qry, qaln.replace("-", ""))
                        self.assertEqual("".join(decode_btop(qry, hsp.btop, sbj=False)).upper(), qaln)
                        self.assertEqual("".join(decode_btop(qry, hsp.btop)).upper(), saln)

    def test_decode_btop_fasta36(self):
        queries = SeqIO.to_dict(SeqIO.parse(root.joinpath("qry.fasta"), "fasta"))
        hsp_getter = attrgetter("query_id", "hit_id", "query_start", "query_end", "hit_start", "hit_end")
        print()
        for prog in ("fasta", "glsearch"):
            with root.joinpath(f"{prog}36.tsv").open() as file1, root.joinpath(f"{prog}36.txt").open() as file2:
                tab = SearchIO.parse(file1, "blast-tab", fields=fields_8CB, comments=True)
                m10 = SearchIO.parse(file2, "fasta-m10")
                for rec1, rec2 in zip(tab, m10):
                    for hit1, hit2 in zip(rec1.hits, rec2.hits):
                        for hsp1, hsp2 in zip(hit1.hsps, hit2.hsps):
                            self.assertEqual(hsp_getter(hsp1), hsp_getter(hsp2))
                            qry = queries[hsp1.query_id].seq[hsp1.query_start : hsp1.query_end]
                            qry = str(qry if hsp1.query_strand >= 0 else qry.reverse_complement())
                            qaln = str(hsp2.query.seq)
                            saln = str(hsp2.hit.seq)
                            print(qaln)
                            print("".join(decode_btop(qry, hsp1.btop, sbj=False)).upper())
                            print("".join(decode_btop(qry, hsp1.btop)).upper())
                            print(saln)
                            print()
                            self.assertEqual(qry, qaln.replace("-", ""))
                            self.assertEqual("".join(decode_btop(qry, hsp1.btop, sbj=False)).upper(), qaln)
                            self.assertEqual("".join(decode_btop(qry, hsp1.btop)).upper(), saln)

    def test_delimify(self):
        print()

        # 001233345567789AABC
        # C[ATG(]AT)TA[CAC]AT
        assay = Assay.factory("assay", "C[ATG(]AT)TA[CAC]AT")

        # 0123456789ABC    0 123  45 67 89A BC
        # CATGATTACACAT    C[ATG(]AT)TA[CAC]AT
        # CAT-ATTACA-AT -> C[AT-(]AT)TA[CA-]AT
        qaln, saln = "CATGATTACACAT", "CAT-ATTACA-AT"
        self.assertEqual("".join(assay.delimify(qaln, qaln)), "C[ATG(]AT)TA[CAC]AT")
        self.assertEqual("".join(assay.delimify(qaln, saln)), "C[AT-(]AT)TA[CA-]AT")

        # 01234566789ABC    0 123  45 667 89A BC
        # CATGATT-ACACAT    C[ATG(]AT)T-A[CAC]AT
        # CAT-ATTTACA-AT -> C[AT-(]AT)TTA[CA-]AT
        qaln, saln = "CATGATT-ACACAT", "CAT-ATTTACA-AT"
        self.assertEqual("".join(assay.delimify(qaln, qaln)), "C[ATG(]AT)T-A[CAC]AT")
        self.assertEqual("".join(assay.delimify(qaln, saln)), "C[AT-(]AT)TTA[CA-]AT")

        # 0001234567798ABC    000 123  45 677798A BC
        # --CATGATTA-CACAT    --C[ATG(]AT)TA-[CAC]AT
        # CCCATGA--ATC-C-T -> CCC[ATG(]A-)-AT[C-C]-T
        qaln, saln = "--CATGATTA-CACAT", "CCCATGA--ATC-C-T"
        self.assertEqual("".join(assay.delimify(qaln, qaln)), "--C[ATG(]AT)TA-[CAC]AT")
        self.assertEqual("".join(assay.delimify(qaln, saln)), "CCC[ATG(]A-)-AT[C-C]-T")

        # 000122222333456789ABC    000 1111123  3345 67 89A BC
        # C--A----TG--ATTACACAT    C--[A----TG]--(AT)TA[CAC]AT
        # CAAATTTTTGGGATTACACAT -> CAA[ATTTTTG]GG(AT)TA[CAC]AT
        qaln, saln = "C--A----TG--ATTACACAT", "CAAATTTTTGGGATTACACAT"
        self.assertEqual("".join(assay.delimify(qaln, qaln)), "C--[A----TG]--(AT)TA[CAC]AT")
        self.assertEqual("".join(assay.delimify(qaln, saln)), "CAA[ATTTTTG]GG(AT)TA[CAC]AT")

        # 0122222333456789ABC    0 1222223 33 45 67 89A BC
        # CAT----G--ATTACACAT    C[AT----G]--(AT)TA[CAC]AT
        # CATTTTTGGGATTACACAT -> C[ATTTTTG]GG(AT)TA[CAC]AT
        qaln, saln = "CAT----G--ATTACACAT", "CATTTTTGGGATTACACAT"
        self.assertEqual("".join(assay.delimify(qaln, qaln)), "C[AT----G]--(AT)TA[CAC]AT")
        self.assertEqual("".join(assay.delimify(qaln, saln)), "C[ATTTTTG]GG(AT)TA[CAC]AT")

        assay = Assay.factory("assay", "[ATG(]AT)TA[CAC]")
        qaln, saln = "AT----G--ATTACAC", "ATTTTTGGGATTACAC"
        self.assertEqual("".join(assay.delimify(qaln, qaln)), "[AT----G]--(AT)TA[CAC]")
        self.assertEqual("".join(assay.delimify(qaln, saln)), "[ATTTTTG]GG(AT)TA[CAC]")

    def test_parse_assay(self):
        path = root / "assay.tsv"
        with path.open() as file:
            self.assertListEqual(
                list(parse_assays(file)),
                [
                    Assay.factory("assay-1", "[GAT]T[ACA]", {666}),
                    Assay.factory("assay-2", "[GAT](T)[ACA]", {666, 662}),
                    Assay.factory("assay-3", "[GAT](T)[ACA]", {666, 662}),
                    Assay.factory("assay-4", "[GATTACA]", {666, 662}),
                    Assay.factory("assay-5", "[GA{T]T}A{[CA}]", {666, 662}),
                ],
            )
            file.seek(0)
            self.assertListEqual(
                list(parse_assays(file, context=(0, 1))),
                [
                    Assay.factory("assay-1", "[GAT]T[ACA]", {666}),
                    Assay.factory("assay-2", "[GAT](T)[ACA]", {666, 662}),
                    Assay.factory("assay-3", "[GAT](T)[ACA]T", {666, 662}),
                    Assay.factory("assay-4", "[GATTACA]T", {666, 662}),
                    Assay.factory("assay-5", "[GA{T]T}A{[CA}]T", {666, 662}),
                ],
            )
            file.seek(0)
            self.assertListEqual(
                list(parse_assays(file, context=(1, 0))),
                [
                    Assay.factory("assay-1", "[GAT]T[ACA]", {666}),
                    Assay.factory("assay-2", "[GAT](T)[ACA]", {666, 662}),
                    Assay.factory("assay-3", "T[GAT](T)[ACA]", {666, 662}),
                    Assay.factory("assay-4", "T[GATTACA]", {666, 662}),
                    Assay.factory("assay-5", "T[GA{T]T}A{[CA}]", {666, 662}),
                ],
            )
            file.seek(0)
            answer = [
                Assay.factory("assay-1", "[GAT]T[ACA]", {666}),
                Assay.factory("assay-2", "[GAT](T)[ACA]", {666, 662}),
                Assay.factory("assay-3", "CAT[GAT](T)[ACA]TAG", {666, 662}),
                Assay.factory("assay-4", "CAT[GATTACA]TAG", {666, 662}),
                Assay.factory("assay-5", "CAT[GA{T]T}A{[CA}]TAG", {666, 662}),
            ]
            self.assertListEqual(list(parse_assays(file, context=(3, 3))), answer)
            file.seek(0)
            self.assertListEqual(list(parse_assays(file, context=(4, 4))), answer)
