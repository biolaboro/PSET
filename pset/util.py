#!/usr/bin/env python3

from csv import DictWriter, Sniffer
from subprocess import PIPE, Popen

fields_std = ("qaccver", "saccver", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")

fields_8CB = (*fields_std, "btop")


def argify(*conf, sep="=", pfx="-"):
    """Generate command-line configuration string entries as tuples.

    -key     -> (key, )
    -key=val -> (key, val)

    Args:
        *conf: the configuration strings
        sep (str, optional): the key-value separator. Defaults to "=".
        pfx (str, optional): the key prefix string. Defaults to "-".

    Yields:
        tuple[str, ...]: the next configuration tuple
    """
    for ele in conf:
        tok = ele.split(sep, maxsplit=1)
        if len(tok) == 2:
            yield (pfx + tok[0], tok[1])
        else:
            yield (pfx + tok[0],)


def batchify(entries, size=10):
    """Generate batches from entries.

    Args:
        entries (iterable): the entries
        size (int, optional): the batch size. Defaults to 10.

    Yields:
        list: the next batch
    """
    batch = []
    for idx, ele in enumerate(entries, start=1):
        batch.append(ele)
        if idx % size == 0:
            yield batch
            batch = []
    if batch:
        yield batch


def blastdbcmd_info(db):
    """Get BLAST+ database metadata.

    Args:
        db (str): the BLAST+ database path

    Yields:
        tuple: the next key-value pair
    """
    cmd = ("blastdbcmd", "-info", "-db", db)
    with Popen(cmd, stdout=PIPE, universal_newlines=True) as pipe:
        with pipe.stdout as file:
            for line in map(str.strip, file):
                # Database: ebola
                if line.startswith("Database:"):
                    yield line.split(": ", maxsplit=1)
                # 3,181 sequences; 52,028,153 total bases
                elif line.endswith("total bases"):
                    tokens = line.split("; ")
                    yield tokens[0].split(" ", maxsplit=1)[::-1]
                    yield tokens[1].split(" ", maxsplit=1)[::-1]
                # Date: Oct 7, 2020  4:14 PM	Longest sequence: 19,897 bases
                elif line.startswith("Date:"):
                    tokens = line.split("\t")
                    yield tokens[0].split(": ", maxsplit=1)
                    yield tokens[1].split(": ", maxsplit=1)
                # BLASTDB Version: 5
                elif line.startswith("BLASTDB Version:"):
                    yield line.split(": ", maxsplit=1)
                # Volumes:
                elif line.startswith("Volumes:"):
                    break
            # /Users/.../Documents/projects/bkp.pset/data/ebola/ebola
            yield "Volumes", list(map(str.strip, file))


def contextify(hsp, qlen):
    sstrand = ("minus", "plus")[hsp.hit_strand >= 0]
    flanks = hsp.query_start, qlen - hsp.query_end
    flanks = flanks if sstrand == "plus" else flanks[::-1]
    sstart, send = hsp.hit_start - flanks[0] + 1, hsp.hit_end + flanks[1]
    sstart = 1 if sstart < 1 else sstart
    return (hsp.hit_id, sstart, send, sstrand)


def iter_hsps(parse):
    yield from (hsp for res in parse for hit in res.hits for hsp in hit.hsps)


def iter_limit(entries, limit=None):
    yield from (ele for idx, ele in enumerate(entries) if limit is None or limit <= 0 or idx < limit)


def iter_lines(lines, comment="#"):
    yield from (line for line in lines if not line.isspace() and line and line.lstrip()[:1] != comment)


def unique(iterable, func=None):
    seen = set()
    for ele in iterable:
        val = func(ele) if func else ele
        if val not in seen:
            yield ele
            seen.add(val)


def slice_aln(qry, sbj, pos, start, end):
    pos -= 1
    for i in range(len(qry)):
        pos += qry[i] != "-"
        if pos == end - 1 and qry[i] == "-" or end <= pos:
            break
        elif start <= pos < end:
            yield qry[i], sbj[i]


def sniff_lines(file, size=1, comment="#"):
    pos = file.tell()
    sample = "".join(iter_limit(iter_lines(file, comment=comment), limit=size))
    dialect = Sniffer().sniff(sample)
    file.seek(pos)
    return dialect


def writerows(rows, file, **kwargs):
    row = next(iter(rows), {})
    writer = DictWriter(file, row.keys(), **kwargs)
    writer.writeheader()
    writer.writerow(row)
    writer.writerows(rows)
