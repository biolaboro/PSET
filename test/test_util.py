#!/usr/bin/env python3

import unittest
from io import StringIO
from subprocess import PIPE, Popen
from tempfile import NamedTemporaryFile

from Bio import SearchIO

from pset.util import (
    batchify,
    contextify,
    iter_limit,
    iter_lines,
    slice_aln,
    sniff_lines,
)


class Test(unittest.TestCase):
    def test_batchify(self):
        self.assertListEqual(list(batchify(range(10), size=2)), [[0, 1], [2, 3], [4, 5], [6, 7], [8, 9]])
        self.assertListEqual(list(batchify(range(11), size=2)), [[0, 1], [2, 3], [4, 5], [6, 7], [8, 9], [10]])

    def test_contextify(self):
        print()

        def test_with_blastn(qry, sbj, ans):
            cmd1 = ("blastn", "-query", "-", "-subject")
            cmd2 = ("blastn", "-outfmt", "6", "-query", "-", "-subject")
            with NamedTemporaryFile(mode="w") as temp:
                temp.write(sbj)
                temp.flush()
                with Popen((*cmd1, temp.name), universal_newlines=True, stdin=PIPE, stdout=PIPE) as proc:
                    with proc.stdin as file:
                        file.write(qry)
                    with proc.stdout as file:
                        print(*file)
                with Popen((*cmd2, temp.name), universal_newlines=True, stdin=PIPE, stdout=PIPE) as proc:
                    with proc.stdin as file:
                        file.write(qry)
                    for res in SearchIO.parse(proc.stdout, "blast-tab"):
                        for hit in res.hits:
                            for hsp in hit.hsps:
                                self.assertEqual(contextify(hsp, qln), ans)

        #              v                                   vv
        #            0123456789012345678901234567890123456789012
        #            123456789A123456789A123456789A123456789A123
        #              ^                                   ^^
        #            GG|||||||||||||||||||||||||||||||||||||GGGG
        sbj = ">sbj\nAAGGGCATGATTATGATTAAATGATTACACATGATGGGGAAAA\n"

        seq = "GGGGGCATGATTATGATTAAATGATTACACATGATGGGGGGGG"
        # Query  3   GGGCATGATTATGATTAAATGATTACACATGATgggg  39
        #            |||||||||||||||||||||||||||||||||||||
        # Sbjct  3   GGGCATGATTATGATTAAATGATTACACATGATGGGG  39
        qln = len(seq)
        qry = f">qry\n{seq}\n"
        test_with_blastn(qry, sbj, ("sbj", 1, 43, "plus"))

        seq = "CCCCCCCCATCATGTGTAATCATTTAATCATAATCATGCCCCC"
        #  Query  5   ccccATCATGTGTAATCATTTAATCATAATCATGCCC  41
        #             |||||||||||||||||||||||||||||||||||||
        #  Sbjct  39  CCCCATCATGTGTAATCATTTAATCATAATCATGCCC  3
        qln = len(seq)
        qry = f">qry\n{seq}\n"
        test_with_blastn(qry, sbj, ("sbj", 1, 43, "minus"))

    def test_iter_limit(self):
        self.assertListEqual(list(iter_limit(range(10))), list(range(10)))
        self.assertListEqual(list(iter_limit(range(10), 5)), list(range(5)))

    def test_iter_lines(self):
        self.assertListEqual(["foo", "bar", "baz"], list(iter_lines(["foo", "bar", "baz"])))
        self.assertListEqual(["foo", "bar", "baz"], list(iter_lines(["foo", "", "#foo", "bar", "   \t ", "baz"])))

    def test_sniff_lines(self):
        """Test sniff..."""
        delimiters = ",;|\t"
        template = "\n# comment\nfoo,bar,baz\n#,;|\t\n    foo,bar,baz    \n"
        for ele in delimiters:
            data = StringIO(template.replace(",", ele))
            with data as file:
                dialect = sniff_lines(file)
                self.assertEqual(dialect.delimiter, ele)

    def test_slice_aln(self):
        print()

        x = "GATTACA"
        for i in range(len(x)):
            self.assertEqual(list(slice_aln("GATTACA", "GATTACA", 0, 0, i)), list(zip(x[:i], x[:i])))

        # -1 0 1 2 3 4 5
        #  - A T T A C A
        #  G A T T A C A
        self.assertEqual(list(slice_aln("-ATTACA", "GATTACA", 0, 2, 5)), list(zip("TAC", "TAC")))
        self.assertEqual(list(slice_aln("-ATTACA", "GATTACA", 0, 0, 5)), list(zip("ATTAC", "ATTAC")))
        self.assertEqual(list(slice_aln("-ATTACA", "GATTACA", 0, 0, 6)), list(zip("ATTACA", "ATTACA")))
        self.assertEqual(list(slice_aln("-ATTACA", "GATTACA", 0, 0, 7)), list(zip("ATTACA", "ATTACA")))

        #  3 4 5 6 7 8 9
        #  - A T T A C A
        #  G A T T A C A
        self.assertEqual(list(slice_aln("-ATTACA", "GATTACA", 4, 1, 6)), list(zip("-AT", "GAT")))
        self.assertEqual(list(slice_aln("-ATTACA", "GATTACA", 4, 3, 7)), list(zip("-ATT", "GATT")))
        self.assertEqual(list(slice_aln("-ATTACA", "GATTACA", 4, 4, 8)), list(zip("ATTA", "ATTA")))
        self.assertEqual(list(slice_aln("-ATTACA", "GATTACA", 4, 5, 9)), list(zip("TTAC", "TTAC")))
        self.assertEqual(list(slice_aln("-ATTACA", "GATTACA", 4, 6, 10)), list(zip("TACA", "TACA")))

        #  0 1 2 3 4 5 5
        #  G A T T A C -
        #  G A T T A C A
        self.assertEqual(list(slice_aln("GATTAC-", "GATTACA", 0, 0, 5)), list(zip("GATTA", "GATTA")))
        self.assertEqual(list(slice_aln("GATTAC-", "GATTACA", 0, 0, 6)), list(zip("GATTAC", "GATTAC")))
        self.assertEqual(list(slice_aln("GATTAC-", "GATTACA", 0, 0, 7)), list(zip("GATTAC-", "GATTACA")))

        #  0 0 1 2 3 3 4
        #  G - T T A - A
        #  G A T T A C A
        self.assertEqual(list(slice_aln("G-TTA-A", "GATTACA", 0, 0, 4)), list(zip("G-TTA", "GATTA")))
        self.assertEqual(list(slice_aln("G-TTA-A", "GATTACA", 0, 0, 5)), list(zip("G-TTA-A", "GATTACA")))

        #  0 0 1 2 2 2 3
        #  G - T T - - A
        #  G A T T A C A
        self.assertEqual(list(slice_aln("G-TT--A", "GATTACA", 0, 0, 3)), list(zip("G-TT", "GATT")))
        self.assertEqual(list(slice_aln("G-TT--A", "GATTACA", 0, 0, 4)), list(zip("G-TT--A", "GATTACA")))

        #  0 1 2 3 4 5 6
        #  G A T T A C A
        #  . . . . . . .
        self.assertEqual(list(slice_aln("GATTACA", "-ATTACA", 0, 0, 7)), list(zip("GATTACA", "-ATTACA")))
        self.assertEqual(list(slice_aln("GATTACA", "GATTAC-", 0, 0, 7)), list(zip("GATTACA", "GATTAC-")))
        self.assertEqual(list(slice_aln("GATTACA", "GA--ACA", 0, 0, 7)), list(zip("GATTACA", "GA--ACA")))

        #  0 1 1 2 3 4 4
        #  G A - T A C -
        #  G A T - A C A
        self.assertEqual(list(slice_aln("GA-TAC-", "GAT-ACA", 0, 1, 7)), list(zip("A-TAC-", "AT-ACA")))

        #  1 2 2 3 4 5 5
        #  G A - T A C -
        #  G A T - A C A
        self.assertEqual(list(slice_aln("GA-TAC-", "GAT-ACA", 1, 2, 5)), list(zip("A-TA", "AT-A")))

        # -1 0 0 1 2 3 3
        #  - A - T A C -
        #  G A T - A C A
        self.assertEqual(list(slice_aln("-A-TAC-", "GAT-ACA", 0, 1, 7)), list(zip("TAC-", "-ACA")))
