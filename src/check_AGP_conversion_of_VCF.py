#!/usr/bin/env python3
"""
verify_liftover.py

Sanity-checks a lifted-over VCF against the TARGET (superscaffold) assembly
FASTA: for every variant, confirms that the VCF's REF allele matches the
actual base(s) present in the assembly at CHROM:POS.

Why this matters: if your AGP-based liftover got an offset wrong, or missed
a reverse-complement on a minus-oriented contig, REF will almost always
disagree with the assembly sequence at that position. A high match rate
(should be ~100% for SNPs called against this same assembly) is strong
evidence the liftover is correct; mismatches pinpoint exactly which variants
to investigate.

Requires: pysam (pip install pysam --break-system-packages)
Requires: the FASTA must be indexed (samtools faidx assembly.fasta) -- pysam
will auto-create the .fai if missing and writable.

Usage:
    python verify_liftover.py lifted.vcf superscaffold_assembly.fasta
"""

import sys

try:
    import pysam
except ImportError:
    sys.exit(
        "This script requires pysam.\n"
        "Install with: pip install pysam --break-system-packages"
    )


def main(vcf_path, fasta_path):
    fasta = pysam.FastaFile(fasta_path)

    n_total = 0
    n_match = 0
    n_chrom_missing = 0
    mismatches = []

    with open(vcf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            chrom, pos, ref = fields[0], int(fields[1]), fields[3]
            n_total += 1

            if chrom not in fasta.references:
                n_chrom_missing += 1
                mismatches.append((chrom, pos, ref, "CHROM_NOT_IN_FASTA"))
                continue

            # VCF POS is 1-based inclusive; pysam.fetch is 0-based half-open
            actual = fasta.fetch(chrom, pos - 1, pos - 1 + len(ref)).upper()

            if actual == ref.upper():
                n_match += 1
            else:
                mismatches.append((chrom, pos, ref, actual))

    print(f"Total variants checked:     {n_total}")
    if n_total > 0:
        print(f"REF matches assembly:       {n_match} ({100 * n_match / n_total:.2f}%)")
    print(f"Mismatches (incl. missing): {len(mismatches)}")
    print(f"  - of which CHROM missing: {n_chrom_missing}")

    if mismatches:
        print("\nFirst 20 mismatches (CHROM, POS, VCF_REF, ASSEMBLY_BASE):")
        for m in mismatches[:20]:
            print("  ", m)
        print(
            "\nTip: if mismatches cluster on contigs with '-' orientation in "
            "your AGP, the reverse-complement step likely has a bug. If they "
            "cluster near contig boundaries, suspect an off-by-one in the "
            "offset calculation."
        )
    else:
        print("\nAll REF alleles match the target assembly. Liftover looks correct.")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit("Usage: python verify_liftover.py lifted.vcf assembly.fasta")
    main(sys.argv[1], sys.argv[2])