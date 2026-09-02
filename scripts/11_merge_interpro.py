#!/usr/bin/env python3
"""Merge the two InterPro sources into one table over all 7,413 proteins.

Two routes produced the domain calls:

  matches_api  the sequence was in UniParc, so InterProScan 6 results came back
               precomputed from the Matches API (step0_out/matches_raw.json.gz)
  ips6_local   the sequence was not in UniParc and was scanned locally

Both are annotated against InterPro 109.0.  That is not assumed: the script
compares the signature-library versions of the two sides and aborts on any
disagreement, and records the ``interpro-version`` the local run stamped into
its own output.

Comparability: the local run covers seven libraries, the Matches API covers
eighteen.  The named per-library columns and ``interpro_accessions`` are
therefore restricted to the seven libraries both sides have, so a column means
the same thing whichever route filled it.  The InterPro entries that only the
Matches API could supply (PANTHER, CATH-Gene3D, and the rest) are not discarded
— they go in ``interpro_accessions_other_libs``, which is empty for every
ips6_local row by construction.

    python3 scripts/11_merge_interpro.py --ips6-json <utex25_unmatched.json>
"""

from __future__ import annotations

import argparse
import collections
import gzip
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DEFAULT_MATCHES = ROOT / "step0_out" / "matches_raw.json.gz"
DEFAULT_LOOKUP = ROOT / "step0_out" / "lookup_status.tsv"
DEFAULT_TABLE = ROOT / "annotation" / "UTEX25_gene_table.tsv"
DEFAULT_OUT = ROOT / "annotation" / "UTEX25_interpro.tsv"

IPR = re.compile(r"^IPR\d{6}$")

# Column name -> the signature libraries that feed it.  PROSITE is one column
# because InterProScan reports patterns and profiles as separate libraries.
LIBRARY_COLUMNS = [
    ("pfam_accessions", ("Pfam",)),
    ("ncbifam", ("NCBIFAM",)),
    ("prosite", ("PROSITE patterns", "PROSITE profiles")),
    ("smart", ("SMART",)),
    ("cdd", ("CDD",)),
    ("superfamily", ("SUPERFAMILY",)),
]
SHARED_LIBRARIES = {lib for _, libs in LIBRARY_COLUMNS for lib in libs}


def open_maybe_gzip(path: Path, mode="rt"):
    return (gzip.open if path.suffix == ".gz" else open)(path, mode,
                                                         encoding="utf-8")


def read_tsv(path: Path):
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        return header, [dict(zip(header, line.rstrip("\n").split("\t")))
                        for line in fh if line.strip()]


# --------------------------------------------------------------------------- #
# The two sources, reduced to the same shape
# --------------------------------------------------------------------------- #

def summarise(matches):
    """One protein's matches -> per-column accession sets + library versions."""
    per_column = {name: set() for name, _ in LIBRARY_COLUMNS}
    interpro_shared, interpro_other = set(), set()
    signatures_shared = set()
    versions = {}
    for match in matches:
        signature = match["signature"]
        release = signature["signatureLibraryRelease"]
        library = release["library"]
        versions[library] = release["version"]
        entry = signature.get("entry")
        accession = signature["accession"]
        shared = library in SHARED_LIBRARIES
        if shared:
            signatures_shared.add(accession)
            for name, libraries in LIBRARY_COLUMNS:
                if library in libraries:
                    per_column[name].add(accession)
        if entry and IPR.match(entry["accession"]):
            (interpro_shared if shared else interpro_other).add(
                entry["accession"])
    # An entry reachable from a shared library is not "other", even if some
    # non-shared library also points at it.
    interpro_other -= interpro_shared
    return per_column, interpro_shared, interpro_other, signatures_shared, versions


def load_matches_api(matches: Path, lookup_rows):
    by_md5 = {}
    with open_maybe_gzip(matches) as fh:
        raw = json.load(fh)
    for batch in raw["batches"]:
        if batch.get("status") != 200:
            raise SystemExit(f"{matches}: batch {batch.get('batch')} is not 200")
        for result in batch["body"]["results"]:
            by_md5[result["md5"].upper()] = result

    out, versions = {}, {}
    for row in lookup_rows:
        if row["in_uniparc"] != "yes":
            continue
        result = by_md5[row["md5"].upper()]
        summary = summarise(result["matches"])
        out[row["protein_id"]] = summary[:4]
        versions.update(summary[4])
    return out, versions


def load_ips6(path: Path):
    with open_maybe_gzip(path) as fh:
        raw = json.load(fh)
    out, versions = {}, {}
    for result in raw["results"]:
        ids = [x["id"] for x in result.get("xref", [])]
        if not ids:
            raise SystemExit(f"{path}: a result carries no xref id")
        summary = summarise(result["matches"])
        for identifier in ids:
            # InterProScan keeps the whole FASTA id; ours have no spaces.
            key = identifier.split()[0]
            if key in out:
                raise SystemExit(f"{path}: duplicate query id {key!r}")
            out[key] = summary[:4]
        versions.update(summary[4])
    return out, versions, raw.get("interpro-version"), raw.get(
        "interproscan-version")


# --------------------------------------------------------------------------- #

def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--matches", type=Path, default=DEFAULT_MATCHES)
    ap.add_argument("--lookup", type=Path, default=DEFAULT_LOOKUP)
    ap.add_argument("--gene-table", type=Path, default=DEFAULT_TABLE)
    ap.add_argument("--ips6-json", type=Path, required=True)
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = ap.parse_args()

    _, lookup_rows = read_tsv(args.lookup)
    _, genes = read_tsv(args.gene_table)
    order = [g["locus_tag"] for g in genes]

    api, api_versions = load_matches_api(args.matches, lookup_rows)
    local, local_versions, interpro_version, ips6_version = load_ips6(
        args.ips6_json)

    # -- the alignment check, not an assumption ----------------------------- #
    shared_libs = sorted(set(api_versions) & set(local_versions))
    mismatches = [(lib, api_versions[lib], local_versions[lib])
                  for lib in shared_libs
                  if api_versions[lib] != local_versions[lib]]
    if mismatches:
        raise SystemExit(
            "signature library versions differ between the two sources; "
            "refusing to merge results annotated on different bases:\n  "
            + "\n  ".join(f"{lib}: matches_api {a} vs ips6_local {b}"
                          for lib, a, b in mismatches))

    expected = {row["protein_id"]: row["in_uniparc"] for row in lookup_rows}
    missing_local = [p for p, v in expected.items()
                     if v == "no" and p not in local]
    if missing_local:
        raise SystemExit(
            f"{len(missing_local)} sequence(s) with in_uniparc=no have no row "
            f"in {args.ips6_json} (first: {missing_local[:5]}); the local run "
            "did not cover the whole input")
    extra_local = sorted(set(local) - set(order))
    if extra_local:
        raise SystemExit(f"{args.ips6_json} has ids absent from "
                         f"{args.gene_table}: {extra_local[:5]}")

    header = ["locus_tag", "source", "interpro_accessions", "pfam_accessions",
              "ncbifam", "prosite", "smart", "cdd", "superfamily", "n_domains",
              "interpro_accessions_other_libs"]

    counts = collections.Counter()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w") as fh:
        fh.write("\t".join(header) + "\n")
        for locus in order:
            if expected[locus] == "yes":
                source, summary = "matches_api", api[locus]
            else:
                source, summary = "ips6_local", local[locus]
            per_column, interpro, other, signatures = summary
            counts[source] += 1
            row = [locus, source, ",".join(sorted(interpro))]
            for name, _ in LIBRARY_COLUMNS:
                row.append(",".join(sorted(per_column[name])))
            row.append(str(len(signatures)))
            row.append(",".join(sorted(other)))
            fh.write("\t".join(row) + "\n")

    print(f"wrote {args.out} ({len(order)} rows + header)")
    print(f"  matches_api : {counts['matches_api']}")
    print(f"  ips6_local  : {counts['ips6_local']}")
    print(f"InterProScan version (local run) : {ips6_version}")
    print(f"InterPro version (local run)     : {interpro_version}")
    print(f"library versions compared        : {len(shared_libs)}, "
          f"{len(mismatches)} mismatch")
    for lib in shared_libs:
        print(f"  {lib:<20} {api_versions[lib]}")
    only_api = sorted(set(api_versions) - set(local_versions))
    only_local = sorted(set(local_versions) - set(api_versions))
    print(f"libraries only in matches_api    : {', '.join(only_api) or '(none)'}")
    print(f"libraries only in ips6_local     : {', '.join(only_local) or '(none)'}")


if __name__ == "__main__":
    main()
