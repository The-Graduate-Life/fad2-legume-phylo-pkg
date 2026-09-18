#!/usr/bin/env python3
"""
Compute pairwise percent identity from an aligned FASTA file.
Identity = (identical, non-gap-in-either-sequence columns) / (columns where
neither sequence has a gap), reported as a percentage.

Usage: python3 pairwise_identity.py aligned.fasta > identity_matrix.tsv
"""
import sys
from collections import OrderedDict


def read_fasta(path):
    seqs = OrderedDict()
    name = None
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                name = line[1:].split()[0]
                seqs[name] = []
            else:
                seqs[name].append(line)
    return {k: "".join(v) for k, v in seqs.items()}


def pct_identity(a, b):
    assert len(a) == len(b)
    compared = 0
    identical = 0
    for x, y in zip(a, b):
        if x == "-" or y == "-":
            continue
        compared += 1
        if x == y:
            identical += 1
    return 100.0 * identical / compared if compared else float("nan")


def main():
    if len(sys.argv) != 2:
        sys.exit("usage: pairwise_identity.py aligned.fasta")
    seqs = read_fasta(sys.argv[1])
    names = list(seqs.keys())
    print("\t" + "\t".join(names))
    for n1 in names:
        row = [n1]
        for n2 in names:
            row.append(f"{pct_identity(seqs[n1], seqs[n2]):.1f}")
        print("\t".join(row))


if __name__ == "__main__":
    main()
