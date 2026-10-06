#!/usr/bin/env python3

from pathlib import Path
from collections import OrderedDict

ALIGN_DIR = Path("/home/tensenl/BATS/FASTAS/OUTGROUPS/loci_aligned")
OUTFILE = Path("/home/tensenl/BATS/FASTAS/OUTGROUPS/fixed_P_vs_D.tsv")

DARK = ["D1", "D2", "D3"]
PALE = ["P1", "P2", "P3", "P4"]

VALID = {"A", "C", "G", "T"}


def read_fasta(path):
    seqs = OrderedDict()
    name = None
    chunks = []

    with open(path) as fh:
        for line in fh:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks).upper()

                name = line[1:].split()[0]
                chunks = []

            else:
                chunks.append(line)

        if name is not None:
            seqs[name] = "".join(chunks).upper()

    return seqs


with open(OUTFILE, "w") as out:

    out.write(
        "gene\talignment_pos\tD_bases\tP_bases\t"
        "D1\tD2\tD3\tP1\tP2\tP3\tP4\n"
    )

    total = 0

    for fasta in sorted(ALIGN_DIR.glob("*_aligned.fasta")):

        gene = fasta.name.replace("_aligned.fasta", "")
        seqs = read_fasta(fasta)

        required = DARK + PALE

        missing = [s for s in required if s not in seqs]

        if missing:
            print(
                f"SKIP {gene}: missing "
                + ", ".join(missing)
            )
            continue

        lengths = {len(seqs[s]) for s in required}

        if len(lengths) != 1:
            print(f"SKIP {gene}: unequal alignment lengths")
            continue

        aln_len = next(iter(lengths))
        n_diff = 0

        for i in range(aln_len):

            d_bases = [seqs[s][i] for s in DARK]
            p_bases = [seqs[s][i] for s in PALE]

            # Keep only unambiguous A/C/G/T bases
            d_valid = [b for b in d_bases if b in VALID]
            p_valid = [b for b in p_bases if b in VALID]

            # Need at least one valid base in each group
            if not d_valid or not p_valid:
                continue

            # Sets of bases present in each group
            d_set = set(d_valid)
            p_set = set(p_valid)

            # Report only if the groups have no nucleotide in common
            if d_set.isdisjoint(p_set):

                d_alleles = "/".join(sorted(d_set))
                p_alleles = "/".join(sorted(p_set))

                out.write(
                    f"{gene}\t{i+1}\t"
                    f"{d_alleles}\t{p_alleles}\t"
                    + "\t".join(d_bases + p_bases)
                    + "\n"
                )

                n_diff += 1
                total += 1

        print(f"{gene}: {n_diff} P-vs-D differences")

print(f"\nTotal P-vs-D differences: {total}")
print(f"Results written to: {OUTFILE}")
