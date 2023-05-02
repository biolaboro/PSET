#!/usr/bin/env python3

import json
import sys
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, FileType


def parse_args(argv):
    parser = ArgumentParser(description="FASTA file -> primers", formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument("json", help="the JSON file of Primer3 output", type=FileType())
    return parser.parse_args(argv)


def main(argv):
    args = parse_args(argv[1:])

    with args.json as file:
        results = json.load(file)

    for result in results:
        template = result["SEQUENCE_TEMPLATE"]
        for num in range(result["PRIMER_PAIR_NUM_RETURNED"]):
            # The selected left primer (the primer to the left in the input sequence).
            # i is the 0-based index of the start base of the primer, and n is t its length.
            idx1, len1 = result[f"PRIMER_LEFT_{num}"]
            # The selected internal oligo.
            # Primer3 outputs this tag if PRIMER_PICK_INTERNAL_OLIGO was non-0.
            # If Primer3 fails to pick a middle oligo upon request, this tag will not be output.
            # i is the 0-based index of start base of the internal oligo, and n is its length.
            idx2, len2 = result.get(f"PRIMER_INTERNAL_{num}", (0, 0))
            # The selected right primer (the primer to the right in the input sequence).
            # i is the 0-based index of the last base of the primer, and n is its length.
            idx3, len3 = result[f"PRIMER_RIGHT_{num}"]
            # adjust from 3'->5' to 5'->3'
            idx3 = idx3 - len3 + 1
            # 5'-context + left
            def1 = template[:idx1] + "[" + template[idx1 : idx1 + len1] + "]"
            if idx2:
                # interprimer + internal + interprimer + 3'-context
                def2 = template[idx1 + len1 : idx2] + "(" + template[idx2 : idx2 + len2] + ")" + template[idx2 + len2 : idx3]
            else:
                # interprimer
                def2 = template[idx1 + len1 : idx3]
            # right + 3'-context
            def3 = "[" + template[idx3 : idx3 + len3] + "]" + template[idx3 + len3 :]
            definition = def1 + def2 + def3
            assert "".join(filter(str.isalpha, definition)) in template
            result.update({f"PSET_DEFINITION_{num}": definition})

    json.dump(results, sys.stdout, indent=4)

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
