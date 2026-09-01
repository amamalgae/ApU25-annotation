#!/usr/bin/env python3
"""Draw a reproducible random sample from step0_out/unmatched.faa.

Used to time a local InterProScan 6 run on a subset and extrapolate to the whole
unmatched set.  The sample must be reproducible, so the seed is explicit and the
draw is done with random.Random(seed).sample over the record indices in file
order — no dependence on dict or set iteration order.

``unmatched.faa`` carries all 235 sequences with an internal stop, and
InterProScan 6 rejects the whole input with "Invalid character(s) found in the
input FASTA file" when one is present.  ``--mask-internal-stop`` replaces ``*``
with the given residue (default off).  Masking with ``X`` keeps the residue count
unchanged, which matters when the sample is used to extrapolate run time by
residues; it is a choice about the data, not a fix, and the count of masked
sequences is reported so it stays visible.

    python3 scripts/09_sample_unmatched.py --n 100 --seed 20260901 \\
        --mask-internal-stop X --out /path/to/sample_100.faa
"""

from __future__ import annotations

import argparse
import random
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DEFAULT_FAA = ROOT / "step0_out" / "unmatched.faa"


def read_fasta(path: Path):
    header, buf = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(buf)
                header, buf = line[1:], []
            else:
                buf.append(line.strip())
    if header is not None:
        yield header, "".join(buf)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--faa", type=Path, default=DEFAULT_FAA)
    ap.add_argument("--n", type=int, default=100)
    ap.add_argument("--seed", type=int, default=20260901)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--mask-internal-stop", metavar="RESIDUE", default=None,
                    help="replace '*' with this residue (e.g. X); off by default")
    args = ap.parse_args()

    records = list(read_fasta(args.faa))
    if args.n > len(records):
        raise SystemExit(f"{args.faa} has {len(records)} records, fewer than "
                         f"--n {args.n}")
    picked = sorted(random.Random(args.seed).sample(range(len(records)), args.n))

    total = 0
    masked_seqs = masked_residues = 0
    with args.out.open("w") as fh:
        for i in picked:
            header, seq = records[i]
            if "*" in seq:
                masked_seqs += 1
                masked_residues += seq.count("*")
                if args.mask_internal_stop:
                    seq = seq.replace("*", args.mask_internal_stop)
            total += len(seq)
            fh.write(f">{header}\n")
            for j in range(0, len(seq), 60):
                fh.write(seq[j:j + 60] + "\n")

    print(f"source     : {args.faa} ({len(records)} sequences, "
          f"{sum(len(s) for _, s in records)} aa)")
    print(f"seed       : {args.seed}")
    print(f"sampled    : {args.n} sequences, {total} aa -> {args.out}")
    if args.mask_internal_stop:
        print(f"masked     : {masked_seqs} sequence(s), {masked_residues} "
              f"'*' -> {args.mask_internal_stop!r}")
    elif masked_seqs:
        print(f"WARNING    : {masked_seqs} sampled sequence(s) contain '*'; "
              "InterProScan 6 rejects the whole file. Use --mask-internal-stop.")


if __name__ == "__main__":
    main()
