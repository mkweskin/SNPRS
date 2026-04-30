#!/usr/bin/env python3

import sys
import os
from itertools import combinations
from Bio import SeqIO

VALID = set("ACTG-")

def compare_sequences(seq1, seq2):
    cocalled = 0
    identical = 0
    different = 0

    for a, b in zip(seq1, seq2):
        if a in VALID and b in VALID:
            cocalled += 1
            if a == b:
                identical += 1
            else:
                different += 1

    return cocalled, identical, different

def main():
    if len(sys.argv) < 2:
        print("Usage: getPairwise_FASTA.py alignment.fasta", file=sys.stderr)
        sys.exit(1)

    fasta = sys.argv[1]
    seqs = {}
    for record in SeqIO.parse(fasta, "fasta"):
        seqs[record.id] = str(record.seq).upper()
    names = list(seqs.keys())

    lengths = {len(s) for s in seqs.values()}
    if len(lengths) != 1:
        raise ValueError("All sequences must be the same length (aligned FASTA required)")

    base = os.path.splitext(fasta)[0]
    outfile = f"{base}_pw_dist.tsv"

    with open(outfile, "w") as out:
        out.write("Sample1\tSample2\tCocalled\tIdentical\tDifferent\n")

        for n1, n2 in combinations(names, 2):
            cocalled, identical, different = compare_sequences(seqs[n1], seqs[n2])
            out.write(f"{n1}\t{n2}\t{cocalled}\t{identical}\t{different}\n")

if __name__ == "__main__":
    main()