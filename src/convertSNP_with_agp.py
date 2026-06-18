#!/usr/bin/env python3
"""
!!!!! this script was written by claude code !!!!!!

Lift VCF coordinates from contigs to superscaffolds using an AGP file.

Usage:
    python agp_liftover.py contigs.vcf assembly.agp > superscaffold.vcf

Handles:
  - + and - orientation of component contigs
  - Updates CHROM, POS, and the ##contig header lines
  - Skips variants on contigs absent from the AGP, or in gap rows
"""

import sys
from tqdm import tqdm

def parse_agp(agp_path):
    """
    Return:
      mapping[contig] = (super_id, obj_start, obj_end, comp_start, comp_end, orient)
        - obj_start/end : 1-based coords on the superscaffold (object)
        - comp_start/end: 1-based coords on the contig (component)
        - orient        : '+' or '-'
      super_len[super_id] = max object coordinate seen (for header)
    """
    mapping = {}
    super_len = {}
    with open(agp_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            obj_id = f[0]
            obj_start = int(f[1])
            obj_end = int(f[2])
            comp_type = f[4]
            super_len[obj_id] = max(super_len.get(obj_id, 0), obj_end)
            # Gap rows: component_type in {N, U}
            if comp_type in ("N", "U"):
                continue
            contig = f[5]
            comp_start = int(f[6])
            comp_end = int(f[7])
            orient = f[8]
            mapping[contig] = (obj_id, obj_start, obj_end,
                               comp_start, comp_end, orient)
    return mapping, super_len


def lift_pos(pos, rec):
    """Map a 1-based contig position to a 1-based superscaffold position."""
    obj_id, obj_start, obj_end, comp_start, comp_end, orient = rec
    if pos < comp_start or pos > comp_end:
        return None  # outside the placed component span
    if orient == "+":
        new_pos = obj_start + (pos - comp_start)
    else:  # '-' : contig is reverse-complemented onto the superscaffold
        new_pos = obj_end - (pos - comp_start)
    return obj_id, new_pos


def main():
    if len(sys.argv) != 3:
        sys.exit("usage: python PhD_chapter3/src/convertSNP_with_agp.py in.vcf assembly.agp > out.vcf")
    vcf_path, agp_path = sys.argv[1], sys.argv[2]
    mapping, super_len = parse_agp(agp_path)

    # Emit new ##contig headers
    contig_header_written = False
    n_skipped = 0
    n_lifted = 0

    with open(vcf_path) as fh:
        for line in tqdm(fh):
            if line.startswith("##"):
                if line.startswith("##contig"):
                    if not contig_header_written:
                        for sid in sorted(super_len):
                            sys.stdout.write(
                                f"##contig=<ID={sid},length={super_len[sid]}>\n")
                        contig_header_written = True
                    continue  # drop old contig headers
                sys.stdout.write(line)
                continue
            if line.startswith("#CHROM"):
                if not contig_header_written:
                    for sid in sorted(super_len):
                        sys.stdout.write(
                            f"##contig=<ID={sid},length={super_len[sid]}>\n")
                    contig_header_written = True
                sys.stdout.write(line)
                continue

            f = line.rstrip("\n").split("\t")
            chrom, pos = f[0], int(f[1])
            rec = mapping.get(chrom)
            if rec is None:
                n_skipped += 1
                continue
            lifted = lift_pos(pos, rec)
            if lifted is None:
                n_skipped += 1
                continue
            f[0], f[1] = lifted[0], str(lifted[1])
            sys.stdout.write("\t".join(f) + "\n")
            n_lifted += 1

    sys.stderr.write(f"lifted {n_lifted} variants, skipped {n_skipped}\n")
    sys.stderr.write("NOTE: output is unsorted and '-' strand REF/ALT alleles are NOT reverse-complemented.\n")


if __name__ == "__main__":
    main()
