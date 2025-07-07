#!/usr/bin/env python3

import os
import unittest
from collections import defaultdict
from itertools import chain
from operator import attrgetter
from pathlib import Path
from subprocess import PIPE, Popen, check_call
from tempfile import NamedTemporaryFile

from Bio import SearchIO, Seq, SeqIO

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

    return next(assay.hits(list(chain.from_iterable(hsps.values())), dFR=(1, 1000), dF3F2=(1, 100), dF2F1c=(1, 100), dF1cB1c=(1, 100)), None)


class Test(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        os.chdir(root)
        check_call(root / "setup.sh")

    def test_factory(self):
        print("\n")
        print(*Assay.factory("[A]", {1}, "assay").records(), "assay", sep="\n")
        print(Assay.factory("[A]", {1}).components, sep="\n")
        print(Assay.factory("[A]", id="id").components, sep="\n")
        print(Assay.factory("N[A]").components, sep="\n")
        print(Assay.factory("[A]N").components, sep="\n")
        print(Assay.factory("N[A]N").components, sep="\n")
        print(Assay.factory("[A][T]").components, sep="\n")
        print(Assay.factory("N[A][T]").components, sep="\n")
        print(Assay.factory("[A]N[T]").components, sep="\n")
        print(Assay.factory("[A][T]N").components, sep="\n")
        print(Assay.factory("N[A](N)[T]N").components, sep="\n")
        print(Assay.factory("[A][C]G[T][AA]CC[GG][TT]").components, sep="\n")
        print(Assay.factory("[A][C](G)[T][AA]CC[GG][TT]").components, sep="\n")
        print(Assay.factory("[A][C]G[T][AA](CC)[GG][TT]").components, sep="\n")
        print(Assay.factory("[A][C](G)[T][AA](CC)[GG][TT]").components, sep="\n")
        print(Assay.factory("N[A][C](G)[T][AA](CC)[GG][TT]N").components, sep="\n")
        print(Assay.factory("N[A][C](G)[T][AA](CC)[GG][TT]N"))

    def test_lamp_internals(self):
        lamp = Assay.factory("[A][C]G[T][AA]CC[GG][TT]")
        self.assertEqual(lamp.subassays[0].definition, "[A]CGTAACCGG[TT]")
        self.assertEqual(lamp.subassays[1].definition, "[C]G[T]")
        self.assertEqual(lamp.subassays[2].definition, "[AA]CC[GG]")
        lamp = Assay.factory("[A][C](G)[T][AA](CC)[GG][TT]")
        self.assertEqual(lamp.subassays[0].definition, "[A]CGTAACCGG[TT]")
        self.assertEqual(lamp.subassays[1].definition, "[C](G)[T]")
        self.assertEqual(lamp.subassays[2].definition, "[AA](CC)[GG]")
        lamp = Assay.factory("[A][C]G[T][AA](CC)[GG][TT]")
        self.assertEqual(lamp.subassays[0].definition, "[A]CGTAACCGG[TT]")
        self.assertEqual(lamp.subassays[1].definition, "[C]G[T]")
        self.assertEqual(lamp.subassays[2].definition, "[AA](CC)[GG]")
        lamp = Assay.factory("[A][C](G)[T][AA]CC[GG][TT]")
        self.assertEqual(lamp.subassays[0].definition, "[A]CGTAACCGG[TT]")
        self.assertEqual(lamp.subassays[1].definition, "[C](G)[T]")
        self.assertEqual(lamp.subassays[2].definition, "[AA]CC[GG]")

    def test_amplicon(self):
        self.assertEqual(Assay.factory("[AAA]CG[GGGT]").amplicon(), "AAACGGGGT")
        self.assertEqual(Assay.factory("N[AAA]CG[GGGT]N").amplicon(), "AAACGGGGT")

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

    def test_invalid_definitions(self):
        val = [
            "[][]",
            "[]()[]",
            "[GATTACA",
            "GATTACA]",
            "[GATTACA[]",
            "[GATTACA]]",
            "[GATTACA]CAT[0123456789]",
            "[GATTACA[CAT]GATTACA]",
            "(GATTACA)C[ATG]A[TTACA]",
            "[GATTACA]C[AT(G]ATTACA)",
            "[A(T]CG[T)A]",
            "N[A(T]CG[T)A]",
            "N[A(T]CG[T)A]N",
            "[A][C]GT[AA]CC[GG][TT]",
            "[A][C](G)T[AA]CC[GG][TT]",
            "[A][C]G[T][AA](CC)[GGTT]",
            "[A][C](G)[T]AA(CC)[GG][TT]",
            "N[A][C]G[T][AA](C)(C)[GG][TT]N",
        ]
        for idx, ele in enumerate(val, start=1):
            with self.assertRaises(Exception):
                print(ele, Assay.factory(f"assay-{idx}", ele, {1}))

    def test_lamp_fip_bip(self):
        assay = Assay.factory(
            "GCTCCAGAAA[TTCAAGTGGAGCACTTGGGCTA]TCAGCATCTATAAAATACTTCG[ATGTTGAGGCGAAGTTTAGGT]AACCTT(TTACCGCCGCAATACGAGC)CCCAAGGTTAGT[GGCATTTGGGCTTGTCGGATCA]ATATCCATAATG[CAAAACGGGCGACGT]TTAGGCCCC(AACACCTCCTGCATAACTCTTGC)CAT[TGCTCGAGTCTGTCGACGA]GCGGTTACCTTATCTAACAATATGGTACTGG[GTACGCCACTGGTAGCAGAAGA]TTGA",
            {666},
            "Chakraborty-ctxA",
        )
        self.assertEqual(assay.fip(), "TGATCCGACAAGCCCAAATGCCATGTTGAGGCGAAGTTTAGGT")
        self.assertEqual(assay.bip(), "CAAAACGGGCGACGTTCGTCGACAGACTCGAGCA")

    def test_lamp_hits(self):
        assay = Assay.factory(
            "GCTCCAGAAA[TTCAAGTGGAGCACTTGGGCTA]TCAGCATCTATAAAATACTTCG[ATGTTGAGGCGAAGTTTAGGT]AACCTT(TTACCGCCGCAATACGAGC)CCCAAGGTTAGT[GGCATTTGGGCTTGTCGGATCA]ATATCCATAATG[CAAAACGGGCGACGT]TTAGGCCCC(AACACCTCCTGCATAACTCTTGC)CAT[TGCTCGAGTCTGTCGACGA]GCGGTTACCTTATCTAACAATATGGTACTGG[GTACGCCACTGGTAGCAGAAGA]TTGA",
            {666},
            "Chakraborty-ctxA",
        )
        subject = f">sbj\n{''.join(filter(str.isalpha, assay.definition))}\n"

        for subassay in assay.subassays:
            self.assertIsNotNone(perfect_hits(subassay, subject))

        self.assertIsNotNone(perfect_hits(assay, subject))

        # revcomp LoopF only
        subject = "GCTCCAGAAA[TTCAAGTGGAGCACTTGGGCTA]TCAGCATCTATAAAATACTTCG[ATGTTGAGGCGAAGTTTAGGT]AACCTT(GCTCGTATTGCGGCGGTAA)CCCAAGGTTAGT[GGCATTTGGGCTTGTCGGATCA]ATATCCATAATG[CAAAACGGGCGACGT]TTAGGCCCC(AACACCTCCTGCATAACTCTTGC)CAT[TGCTCGAGTCTGTCGACGA]GCGGTTACCTTATCTAACAATATGGTACTGG[GTACGCCACTGGTAGCAGAAGA]TTGA"
        subject = f">sbj\n{''.join(filter(str.isalpha, subject))}\n"
        print(subject)
        self.assertIsNone(perfect_hits(assay, subject))

        # revcomp LoopB only
        subject = "GCTCCAGAAA[TTCAAGTGGAGCACTTGGGCTA]TCAGCATCTATAAAATACTTCG[ATGTTGAGGCGAAGTTTAGGT]AACCTT(TTACCGCCGCAATACGAGC)CCCAAGGTTAGT[GGCATTTGGGCTTGTCGGATCA]ATATCCATAATG[CAAAACGGGCGACGT]TTAGGCCCC(GCAAGAGTTATGCAGGAGGTGTT)CAT[TGCTCGAGTCTGTCGACGA]GCGGTTACCTTATCTAACAATATGGTACTGG[GTACGCCACTGGTAGCAGAAGA]TTGA"
        subject = f">sbj\n{''.join(filter(str.isalpha, subject))}\n"
        self.assertIsNone(perfect_hits(assay, subject))

        # revcomp LoopF & LoopB
        subject = "GCTCCAGAAA[TTCAAGTGGAGCACTTGGGCTA]TCAGCATCTATAAAATACTTCG[ATGTTGAGGCGAAGTTTAGGT]AACCTT(GCTCGTATTGCGGCGGTAA)CCCAAGGTTAGT[GGCATTTGGGCTTGTCGGATCA]ATATCCATAATG[CAAAACGGGCGACGT]TTAGGCCCC(GCAAGAGTTATGCAGGAGGTGTT)CAT[TGCTCGAGTCTGTCGACGA]GCGGTTACCTTATCTAACAATATGGTACTGG[GTACGCCACTGGTAGCAGAAGA]TTGA"
        subject = f">sbj\n{''.join(filter(str.isalpha, subject))}\n"
        self.assertIsNone(perfect_hits(assay, subject))

        # revcomp subject
        subject = f">sbj\n{Seq.reverse_complement(''.join(filter(str.isalpha, assay.definition)))}\n"
        self.assertIsNotNone(perfect_hits(assay, subject))

        # test loop/probe right up-against primer
        assay = Assay.factory(
            "TCGCTG[CCAATAACGGCAAAGATATG]CTGATCATGCTGCACCAAATGGGCAATCACGGGCCTGCGTATTTTAAGCGATATGA[TGAAAAGTTTGCCAAATTCA]CGCCAGTGTGTGAAGGTAATGAGCTTGCC[AAGTGCGAACATCAGTCCTT]GATCAATGCTTATGACAATGCCTTG[CTTGCCACCGATGATTTCAT]CGCTCAAA(GTATCCAGTGGCTGCAGACG)[CACAGCAATGCCTATGAT]GTCTCAATGCTGTATGTCAGCGATCATGGCGAAAGTCTGGGTGAGAACGGT[GTCTATCTACATGGTATGCC]AAATGC",
            {2},
            "NG_056412.1_1129-1448_00001",
        )
        subject = """
            >NG_056412.1 Escherichia coli EC15-101 mcr-1 gene for phosphoethanolamine--lipid A transferase MCR-1.12, complete CDS
            ATGATGCACCATACTTCTGTGTGGTACCGACGCTCGGTCAGTCCGTTTGTTCTTGTGGCGAGTGTTGCCG
            TTTTCTTGACCGCGACCGCCAATCTTACCTTTTTTGATAAAATCAGCCAAACCTATCCCATCGCGGACAA
            TCTCGGCTTTGTGCTGACGATCGCTGTCGTGCTCTTTGGCGCGATGCTACTGATCACCACGCTGTTATCA
            TCGTATCGCTATGTGCTAAAGCCTGTGTTGATTTTGCTATTAATCATGGGCGCGGTGACCAGTTATTTTA
            CTGACACTTATGGCACGGTCTATGATACGACCATGCTCCAAAATGCCCTACAGACCGACCAAGCCGAGAC
            CAAGGATCTATTAAACGCAGCGTTTATCATGCGTATCATTGGTTTGGGTGTGCTACCAAGTTTGCTTGTG
            GCTTTTGTTAAGGTGGATTATCCGACTTGGGGCAAGGGTTTGATGCGCCGATTGGGCTTGATCGTGGCAA
            GTCTTGCGCTGATTTTACTGCCTGTGGTGGCGTTCAGCAGTCATTATGCCAGTTTCTTTCGCGTGCATAA
            GCCGCTGCGTAGCTATGTCAATCCGATCATGCCAATCTACTCGGTGGGTAAGCTTGCCAGTATTGAGTAT
            AAAAAAGCCAGTGCGCCAAAAGATACCATTTATCACGCCAAAGACGCGGTACAAGCAACCAAGCCTGATA
            TGCGTAAGCCACGCCTAGTGGTGTTCGTCGTCGGTGAGACGGCACGCGCCGATCATGTCAGCTTCAATGG
            CTATGAGCGCGATACTTTCCCACAGCTTGCCAAGATCGATGGCGTGACCAATTTTAGCAATGTCACATCG
            TGCGGCACATCGACGGCGTATTCTGTGCCGTGTATGTTCAGCTATCTGGGCGCGGATGAGTATGATGTCG
            ATACCGCCAAATACCAAGAAAATGTGCTGGATACGCTGGATCGCTTGGGCGTAAGTATCTTGTGGCGTGA
            TAATAATTCGGACTCAAAAGGCGTGATGGATAAGCTGCCAAAAGCGCAATTTGCCGATTATAAATCCGCG
            ACCAACAACGCCATCTGCAACACCAATCCTTATAACGAATGCCGCGATGTCGGTATGCTCGTTGGCTTAG
            ATGACTTTGTCGCTGCCAATAACGGCAAAGATATGCTGATCATGCTGCACCAAATGGGCAATCACGGGCC
            TGCGTATTTTAAGCGATATGATGAAAAGTTTGCCAAATTCACGCCAGTGTGTGAAGGTAATGAGCTTGCC
            AAGTGCGAACATCAGTCCTTGATCAATGCTTATGACAATGCCTTGCTTGCCACCGATGATTTCATCGCTC
            AAAGTATCCAGTGGCTGCAGACGCACAGCAATGCCTATGATGTCTCAATGCTGTATGTCAGCGATCATGG
            CGAAAGTCTGGGTGAGAACGGTGTCTATCTACATGGTATGCCAAATGCCTTTGCACCAAAAGAACAGCGC
            AGTGTGCCTGCATTTTTCTGGACGGATAAGCAAACTGGCATCACGCCAATGGCAACCGATACCGTCCTGA
            CCCATGACGCGATCACGCCGACATTATTAAAGCTGTTTGATGTCACCGCGGACAAAGTCAAAGACCGCAC
            CGCATTCATCCGCTGA
        """
        subject = "\n".join((ele.strip() for ele in subject[1:].split("\n")))
        self.assertIsNotNone(perfect_hits(assay, subject))

        assay = Assay.factory(
            "[TTCAAGTGGAGCACTTGGGCTA][ATGTTGAGGCGAAGTTTAGGT](TTACCGCCGCAATACGAGC)[GGCATTTGGGCTTGTCGGATCA]T[CAAAACGGGCGACGT](AACACCTCCTGCATAACTCTTGC)[TGCTCGAGTCTGTCGACGA][GTACGCCACTGGTAGCAGAAGA]",
            {666},
            "Chakraborty-ctxA-modified",
        )
        self.assertIsNotNone(perfect_hits(assay, f">sbj\n{''.join(filter(str.isalpha, assay.definition))}\n"))

    def test_pcr_hits(self):
        #                    5'<-> 3'
        # CCCCGCCGCGGGCCGGCGCG <-> CGCGCCGGCCCGCGGCGGGG
        # TTTGGGTGTGGTGTGGTGGT <-> ACCACCACACCACACCCAAA
        assay = Assay.factory("[CCCCGCCGCGGGCCGGCGCG][TTTGGGTGTGGTGTGGTGGT]")
        self.assertIsNotNone(perfect_hits(assay, ">sub\nCCCCGCCGCGGGCCGGCGCGTTTGGGTGTGGTGTGGTGGT\n"))
        assay = Assay.factory("[CCCCGCCGCGGGCCGGCGCG]AAAAAAAAAATTATATTTATATATTATATTAAAAAAAAAA[TTTGGGTGTGGTGTGGTGGT]")
        self.assertIsNotNone(perfect_hits(assay, ">sub\nCCCCGCCGCGGGCCGGCGCGAAAAAAAAAATTATATTTATATATTATATTAAAAAAAAAATTTGGGTGTGGTGTGGTGGT\n"))
        self.assertIsNotNone(perfect_hits(assay, ">sub\nACCACCACACCACACCCAAATTTTTTTTTTAATATAATATATAAATATAATTTTTTTTTTCGCGCCGGCCCGCGGCGGGG\n"))
        self.assertIsNone(perfect_hits(assay, ">sub\nACCACCACACCACACCCAAAAAAAAAAAAATTATATTTATATATTATATTAAAAAAAAAACCCCGCCGCGGGCCGGCGCG\n"))
        self.assertIsNone(perfect_hits(assay, ">sub\nCGCGCCGGCCCGCGGCGGGGTTTTTTTTTTAATATAATATATAAATATAATTTTTTTTTTTTTGGGTGTGGTGTGGTGGT\n"))
        assay = Assay.factory("[CCCCGCCGCGGGCCGGCGCG](TTATATTTATATATTATATT)AAAAAAAAAA[TTTGGGTGTGGTGTGGTGGT]")
        self.assertIsNotNone(perfect_hits(assay, ">sub\nCCCCGCCGCGGGCCGGCGCGTTATATTTATATATTATATTAAAAAAAAAATTTGGGTGTGGTGTGGTGGT\n"))
        assay = Assay.factory("[CCCCGCCGCGGGCCGGCGCG]AAAAAAAAAA(TTATATTTATATATTATATT)[TTTGGGTGTGGTGTGGTGGT]")
        self.assertIsNotNone(perfect_hits(assay, ">sub\nCCCCGCCGCGGGCCGGCGCGTTATATTTATATATTATATTAAAAAAAAAATTTGGGTGTGGTGTGGTGGT\n"))
        assay = Assay.factory("[CCCCGCCGCGGGCCGGCGCG](TTATATTTATATATTATATT)[TTTGGGTGTGGTGTGGTGGT]")
        self.assertIsNotNone(perfect_hits(assay, ">sub\nCCCCGCCGCGGGCCGGCGCGTTATATTTATATATTATATTTTTGGGTGTGGTGTGGTGGT\n"))
        assay = Assay.factory("[CCCCGCCGCGGGCCGGCGCG]AAAAAAAAAA(TTATATTTATATATTATATT)AAAAAAAAAA[TTTGGGTGTGGTGTGGTGGT]")
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
            "CAGCTG[TGCGTCGGCAAACCAATGCT]ATTGAATCACTAGAAGGTCGAGTAACAACTCTTGA(GGCCAGCTTAAAACCCGTTC)AAGACATGGCAAAGACCATATCATCCCTGAATCGCAGCTGT[GCCGAAATGGTTGCAAAATACG]ACCTAC"
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
                        qry = str(queries[hsp.query_id].seq)[hsp.query_start: hsp.query_end]
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
                            qry = queries[hsp1.query_id].seq[hsp1.query_start: hsp1.query_end]
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

    def test_embrace(self):
        print()

        # 001233345567789AABC
        # C[ATG](AT)TA[CAC]AT
        assay = Assay.factory("C[ATG](AT)TA[CAC]AT")

        # 0123456789ABC    0 123  45 67 89A BC
        # CATGATTACACAT    C[ATG](AT)TA[CAC]AT
        # CAT-ATTACA-AT -> C[AT-](AT)TA[CA-]AT
        qaln, saln = "CATGATTACACAT", "CAT-ATTACA-AT"
        self.assertEqual("".join(assay.embrace(qaln, qaln)), "C[ATG](AT)TA[CAC]AT")
        self.assertEqual("".join(assay.embrace(qaln, saln)), "C[AT-](AT)TA[CA-]AT")

        # 01234566789ABC    0 123  45 667 89A BC
        # CATGATT-ACACAT    C[ATG](AT)T-A[CAC]AT
        # CAT-ATTTACA-AT -> C[AT-](AT)TTA[CA-]AT
        qaln, saln = "CATGATT-ACACAT", "CAT-ATTTACA-AT"
        self.assertEqual("".join(assay.embrace(qaln, qaln)), "C[ATG](AT)T-A[CAC]AT")
        self.assertEqual("".join(assay.embrace(qaln, saln)), "C[AT-](AT)TTA[CA-]AT")

        # 0001234567798ABC    000 123  45 677798A BC
        # --CATGATTA-CACAT    --C[ATG](AT)TA-[CAC]AT
        # CCCATGA--ATC-C-T -> CCC[ATG](A-)-AT[C-C]-T
        qaln, saln = "--CATGATTA-CACAT", "CCCATGA--ATC-C-T"
        self.assertEqual("".join(assay.embrace(qaln, qaln)), "--C[ATG](AT)TA-[CAC]AT")
        self.assertEqual("".join(assay.embrace(qaln, saln)), "CCC[ATG](A-)-AT[C-C]-T")

        # 000122222333456789ABC    000 1111123  3345 67 89A BC
        # C--A----TG--ATTACACAT    C--[A----TG]--(AT)TA[CAC]AT
        # CAAATTTTTGGGATTACACAT -> CAA[ATTTTTG]GG(AT)TA[CAC]AT
        qaln, saln = "C--A----TG--ATTACACAT", "CAAATTTTTGGGATTACACAT"
        self.assertEqual("".join(assay.embrace(qaln, qaln)), "C--[A----TG]--(AT)TA[CAC]AT")
        self.assertEqual("".join(assay.embrace(qaln, saln)), "CAA[ATTTTTG]GG(AT)TA[CAC]AT")

        # 0122222333456789ABC    0 1222223 33 45 67 89A BC
        # CAT----G--ATTACACAT    C[AT----G]--(AT)TA[CAC]AT
        # CATTTTTGGGATTACACAT -> C[ATTTTTG]GG(AT)TA[CAC]AT
        qaln, saln = "CAT----G--ATTACACAT", "CATTTTTGGGATTACACAT"
        self.assertEqual("".join(assay.embrace(qaln, qaln)), "C[AT----G]--(AT)TA[CAC]AT")
        self.assertEqual("".join(assay.embrace(qaln, saln)), "C[ATTTTTG]GG(AT)TA[CAC]AT")

        assay = Assay.factory("[ATG](AT)TA[CAC]")
        qaln, saln = "AT----G--ATTACAC", "ATTTTTGGGATTACAC"
        self.assertEqual("".join(assay.embrace(qaln, qaln)), "[AT----G]--(AT)TA[CAC]")
        self.assertEqual("".join(assay.embrace(qaln, saln)), "[ATTTTTG]GG(AT)TA[CAC]")

    def test_parse_assays(self):
        # id	targets	definition
        # assay-1	666	[GAT]T[ACA]
        # assay-2	666;662	[GAT](T)[ACA]
        # assay-3	666;662	CAT[GAT](T)[ACA]TAG
        # assay-4	666;662	CAT[GATTACA]TAG
        path = root / "assay.tsv"
        with path.open() as file:
            self.assertListEqual(
                [ele for ele, _ in parse_assays(file)],
                [
                    Assay.factory("[GAT]T[ACA]", {"666"}, "assay-1"),
                    Assay.factory("[GAT](T)[ACA]", {"666", "662"}, "assay-2"),
                    Assay.factory("[GAT](T)[ACA]", {"666", "662"}, "assay-3"),
                    Assay.factory("[GATTACA]", {"666", "662"}, "assay-4"),
                ],
            )
            file.seek(0)
            self.assertListEqual(
                [ele for ele, _ in parse_assays(file, context=(0, 1))],
                [
                    Assay.factory("[GAT]T[ACA]", {"666"}, "assay-1"),
                    Assay.factory("[GAT](T)[ACA]", {"666", "662"}, "assay-2"),
                    Assay.factory("[GAT](T)[ACA]T", {"666", "662"}, "assay-3"),
                    Assay.factory("[GATTACA]T", {"666", "662"}, "assay-4"),
                ],
            )
            file.seek(0)
            self.assertListEqual(
                [ele for ele, _ in parse_assays(file, context=(1, 0))],
                [
                    Assay.factory("[GAT]T[ACA]", {"666"}, "assay-1"),
                    Assay.factory("[GAT](T)[ACA]", {"666", "662"}, "assay-2"),
                    Assay.factory("T[GAT](T)[ACA]", {"666", "662"}, "assay-3"),
                    Assay.factory("T[GATTACA]", {"666", "662"}, "assay-4"),
                ],
            )
            file.seek(0)
            answer = [
                Assay.factory("[GAT]T[ACA]", {"666"}, "assay-1"),
                Assay.factory("[GAT](T)[ACA]", {"666", "662"}, "assay-2"),
                Assay.factory("CAT[GAT](T)[ACA]TAG", {"666", "662"}, "assay-3"),
                Assay.factory("CAT[GATTACA]TAG", {"666", "662"}, "assay-4"),
            ]
            self.assertListEqual([ele for ele, _ in parse_assays(file, context=(3, 3))], answer)
            file.seek(0)
            self.assertListEqual([ele for ele, _ in parse_assays(file, context=(4, 4))], answer)
