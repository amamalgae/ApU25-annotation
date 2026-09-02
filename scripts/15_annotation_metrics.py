#!/usr/bin/env python3
"""Measure the four annotation coverage metrics, each independently.

A, B, C and D count different things and are each derived from their own
source column.  None of them is computed from any of the others, and none of
them is "the annotation rate" — the definitions are in docs/ANNOTATION.md.

  A  sequences carrying at least one InterPro entry
  B  sequences carrying at least one entry whose type is Family or Domain
  C  sequences whose product is still "hypothetical protein"
  D  sequences whose symbol is still empty

    python3 scripts/15_annotation_metrics.py
"""

from __future__ import annotations

import argparse
import collections
import gzip
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DEFAULT_TABLE = ROOT / "annotation" / "UTEX25_gene_table.tsv"
DEFAULT_LOOKUP = ROOT / "step0_out" / "lookup_status.tsv"
DEFAULT_MATCHES = ROOT / "step0_out" / "matches_raw.json.gz"
DEFAULT_IPS6 = ROOT / "ips6_out" / "utex25_unmatched.json.gz"
HYPOTHETICAL = "hypothetical protein"
# the libraries the local InterProScan run covered, i.e. the ones both
# annotation routes have in common
SHARED_LIBRARIES = {"Pfam", "NCBIFAM", "PROSITE patterns", "PROSITE profiles",
                    "SMART", "CDD", "SUPERFAMILY"}


def read_tsv(path: Path):
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(header, line.rstrip("\n").split("\t")))
                for line in fh if line.strip()]


def entry_types(lookup, matches: Path, ips6: Path):
    with gzip.open(matches, "rt", encoding="utf-8") as fh:
        raw = json.load(fh)
    by_md5 = {}
    for batch in raw["batches"]:
        for result in batch["body"]["results"]:
            by_md5[result["md5"].upper()] = result
    with gzip.open(ips6, "rt", encoding="utf-8") as fh:
        localraw = json.load(fh)
    local = {x["id"].split()[0]: r
             for r in localraw["results"] for x in r.get("xref", [])}
    out = {}
    for row in lookup:
        pid = row["protein_id"]
        matchlist = (by_md5[row["md5"].upper()]["matches"]
                     if row["in_uniparc"] == "yes"
                     else local.get(pid, {"matches": []})["matches"])
        types = collections.Counter()
        shared = 0
        for match in matchlist:
            entry = match["signature"].get("entry")
            if entry:
                types[entry["type"]] += 1
                if (match["signature"]["signatureLibraryRelease"]["library"]
                        in SHARED_LIBRARIES):
                    shared += 1
        types["_shared_entry_matches"] = shared
        out[pid] = types
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gene-table", type=Path, default=DEFAULT_TABLE)
    ap.add_argument("--lookup", type=Path, default=DEFAULT_LOOKUP)
    ap.add_argument("--matches", type=Path, default=DEFAULT_MATCHES)
    ap.add_argument("--ips6", type=Path, default=DEFAULT_IPS6)
    args = ap.parse_args()

    genes = read_tsv(args.gene_table)
    lookup = read_tsv(args.lookup)
    types = entry_types(lookup, args.matches, args.ips6)
    total = len(genes)

    def n_entries(pid):
        t = types[pid]
        return sum(v for k, v in t.items() if not k.startswith("_"))

    a = sum(1 for g in genes if n_entries(g["locus_tag"]) > 0)
    a_shared = sum(1 for g in genes
                   if types[g["locus_tag"]]["_shared_entry_matches"] > 0)
    b = sum(1 for g in genes
            if types[g["locus_tag"]]["Family"] or types[g["locus_tag"]]["Domain"])
    c = sum(1 for g in genes if g["product"] == HYPOTHETICAL)
    d = sum(1 for g in genes if not g["symbol"])

    print(f"total sequences: {total}\n")
    for key, label, value in (
            ("A", "InterPro エントリを1つ以上持つ（全ライブラリ）", a),
            ("A'", "同（共通7ライブラリ由来のみ）", a_shared),
            ("B", "Family または Domain を1つ以上持つ", b),
            ("C", 'product が "hypothetical protein" のまま', c),
            ("D", "symbol が空欄のまま", d)):
        print(f"  指標{key}  {label:<38} {value:>5} / {total}  "
              f"{100 * value / total:5.2f} %")

    print("\nproduct_source:")
    for k, v in sorted(collections.Counter(g["product_source"] for g in genes).items()):
        print(f"  {k:<24} {v:>5}  {100 * v / total:5.2f} %")
    print("symbol_source:")
    for k, v in sorted(collections.Counter(g["symbol_source"] for g in genes).items()):
        print(f"  {k:<24} {v:>5}  {100 * v / total:5.2f} %")

    flags = [("interpro_family_unresolved", "TRUE"),
             ("interpro_structural_only", "TRUE")]
    print("v1.1 への申し送り対象:")
    for column, want in flags:
        n = sum(1 for g in genes if g.get(column) == want)
        print(f"  {column:<28} {n:>5}  {100 * n / total:5.2f} %")
    none_at_all = sum(1 for g in genes if n_entries(g["locus_tag"]) == 0)
    print(f"  {'InterPro エントリなし':<28} {none_at_all:>5}  "
          f"{100 * none_at_all / total:5.2f} %")

    dup = [g for g in genes if g["dup_pair_id"]]
    print(f"\ndup_pair_id: {len(set(g['dup_pair_id'] for g in dup))} groups over "
          f"{len(dup)} rows")


if __name__ == "__main__":
    main()
