#!/usr/bin/env python3
"""Decide symbol and product for every gene and write them into the gene table.

The rules, and the reasoning they come from, are in docs/ANNOTATION.md.  This
module implements them; it does not invent any.  Where a rule cannot decide, the
value stays empty rather than being guessed.

symbol
  1  an existing UTEX 250-A value is kept, never overwritten   -> utex250a
  2  an empty one may take eggNOG's Preferred_name             -> eggnog
  3  otherwise it stays empty                                  -> none

product
  1    InterPro type=Family: its name.  With several, the InterPro 109.0
       parent/child tree picks the most subordinate one        -> interpro_family
  1x   several Family entries with no ancestor relation between them cannot be
       ranked.  Left as "hypothetical protein" on purpose: those entries come
       overwhelmingly from PANTHER, which only the Matches API route ran, so
       choosing one would make the wording depend on which route annotated the
       gene.  Candidates are kept in interpro_family_candidates -> none
  2    no Family, Domain only, deduplicated by accession and ordered N->C:
       one     "<A> domain-containing protein"                 -> interpro_domain
       two     "<A> and <B> domain-containing protein"         -> interpro_domain_multi
       three+  the first two, same wording                     -> interpro_domain_multi
  3    entries that are only Homologous_superfamily / Repeat / Conserved_site
       are not used: a shared fold is not a shared function    -> none
  4    otherwise "hypothetical protein"                        -> none

There is no eggNOG route into product.  emapper 3.0.0-beta6 writes 22 columns
and none of them holds free text, so eggNOG cannot supply a product at all --
the rule that once tried to has been removed rather than left as dead code.
eggNOG still supplies symbols, from Preferred_name, which is a different thing.

Symbols taken from Preferred_name are dropped again when the value is not a
gene symbol: an NCBI LOC placeholder, another organism's locus id, a bare
number, or a single character.  The dropped value is kept in symbol_rejected
with its reason, because removing it from `symbol` and discarding it are
different things.

Nothing is written unless every stop condition passes: no product over 200
characters, every gene matched by some rule, and not one of the existing
symbols altered.

    python3 scripts/13_finalize_annotation.py
"""

from __future__ import annotations

import argparse
import collections
import gzip
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DEFAULT_TABLE = ROOT / "annotation" / "UTEX25_gene_table.tsv"
DEFAULT_LOOKUP = ROOT / "step0_out" / "lookup_status.tsv"
DEFAULT_MATCHES = ROOT / "step0_out" / "matches_raw.json.gz"
DEFAULT_IPS6 = ROOT / "ips6_out" / "utex25_unmatched.json.gz"
DEFAULT_EGGNOG = ROOT / "raw" / "eggnog_7413.emapper.annotations"
DEFAULT_TREE = ROOT / "raw" / "ParentChildTreeFile.txt"

HYPOTHETICAL = "hypothetical protein"
MAX_PRODUCT = 200
STRUCTURAL_ONLY_TYPES = {"Homologous_superfamily", "Repeat", "Conserved_site",
                         "Active_site", "Binding_site", "PTM"}

NEW_COLUMNS = ["symbol_source", "product_source", "dup_pair_id",
               "interpro_family_candidates", "interpro_family_unresolved",
               "interpro_domain_all", "interpro_structural_only",
               "symbol_rejected", "symbol_rejected_reason", "eggnog_data_row"]

# Longest first within each pair, so "domain-containing" is stripped before
# "domain" and "superfamily" before "family".
DOMAIN_SUFFIXES = (" domain-containing", " domain", " superfamily", " family",
                   " repeat")


def domain_stem(text: str) -> str:
    """Drop a trailing type word so " domain-containing protein" does not stutter."""
    stem = text.strip()
    for _ in range(2):
        lowered = stem.lower()
        for suffix in DOMAIN_SUFFIXES:
            if lowered.endswith(suffix):
                stem = stem[: -len(suffix)].rstrip(" ,")
                break
        else:
            break
    return stem or text.strip()


def symbol_rejection(symbol: str):
    """Why an eggNOG Preferred_name is not usable as a gene symbol, or None.

    These are not rejected for matching a pattern but for not being gene
    symbols: a placeholder, another organism's locus id, a bare number, or a
    single letter.  Values that merely look machine-generated but are real
    approved symbols (C1orf74, hapE_2, U2A') are kept.
    """
    if re.fullmatch(r"LOC\d+", symbol):
        return "loc_placeholder"
    if "\\" in symbol or ":" in symbol:
        return "foreign_locus_id"
    if re.fullmatch(r"\d+", symbol):
        return "numeric_only"
    if len(symbol) == 1:
        return "single_char"
    return None


# --------------------------------------------------------------------------- #
# inputs
# --------------------------------------------------------------------------- #

def read_tsv(path: Path):
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        rows = [line.rstrip("\n").split("\t") for line in fh if line.strip()]
    for i, row in enumerate(rows, 2):
        if len(row) != len(header):
            raise SystemExit(f"{path}:{i}: {len(row)} fields, header has "
                             f"{len(header)}")
    return header, rows


def load_hierarchy(path: Path):
    """child -> set(parents).  An accession recurs in the file, so parents accumulate."""
    parents = collections.defaultdict(set)
    stack: dict[int, str] = {}
    for line in open(path):
        m = re.match(r"^(-*)(IPR\d{6})::", line.rstrip("\n"))
        if not m:
            continue
        depth = len(m.group(1)) // 2
        stack[depth] = m.group(2)
        if depth > 0:
            parents[m.group(2)].add(stack[depth - 1])
    if not parents:
        raise SystemExit(f"{path}: no parent/child relations parsed")
    cache: dict[str, set[str]] = {}

    def ancestors(acc: str) -> set[str]:
        if acc in cache:
            return cache[acc]
        out, todo = set(), list(parents.get(acc, ()))
        while todo:
            p = todo.pop()
            if p in out:
                continue
            out.add(p)
            todo.extend(parents.get(p, ()))
        cache[acc] = out
        return out

    return ancestors


def load_eggnog(path: Path):
    columns, out = None, {}
    for line in open(path):
        if line.startswith("##") or not line.strip():
            continue
        if line.startswith("#"):
            columns = [c.lstrip("#").strip() for c in line.rstrip("\n").split("\t")]
            continue
        out[dict(zip(columns, line.rstrip("\n").split("\t")))[columns[0]]] = dict(
            zip(columns, line.rstrip("\n").split("\t")))
    return columns, out


def load_interpro_per_protein(lookup_rows, matches: Path, ips6: Path):
    """protein_id -> {accession: (entry, earliest start)} plus the source label."""
    with gzip.open(matches, "rt", encoding="utf-8") as fh:
        raw = json.load(fh)
    by_md5 = {}
    for batch in raw["batches"]:
        if batch.get("status") != 200:
            raise SystemExit(f"{matches}: batch {batch.get('batch')} is not 200")
        for result in batch["body"]["results"]:
            by_md5[result["md5"].upper()] = result
    with gzip.open(ips6, "rt", encoding="utf-8") as fh:
        localraw = json.load(fh)
    local = {x["id"].split()[0]: r
             for r in localraw["results"] for x in r.get("xref", [])}

    out = {}
    for row in lookup_rows:
        pid = row["protein_id"]
        if row["in_uniparc"] == "yes":
            source, matchlist = "matches_api", by_md5[row["md5"].upper()]["matches"]
        else:
            source = "ips6_local"
            if pid not in local:
                raise SystemExit(f"{ips6}: no result for {pid}")
            matchlist = local[pid]["matches"]
        entries: dict[str, dict] = {}
        starts: dict[str, int] = {}
        for match in matchlist:
            entry = match["signature"].get("entry")
            if not entry:
                continue
            acc = entry["accession"]
            entries[acc] = entry
            begin = min((loc["start"] for loc in match.get("locations", [])),
                        default=None)
            if begin is not None and (acc not in starts or begin < starts[acc]):
                starts[acc] = begin
        out[pid] = (source, entries, starts)
    return out


# --------------------------------------------------------------------------- #
# the rules
# --------------------------------------------------------------------------- #

def label(entry, field):
    value = (entry.get(field) or "").strip()
    return value or entry["accession"]


def decide_product(entries, starts, ancestors, field):
    """-> (product, product_source, family_candidates, domain_all, structural_only)"""
    families = {a for a, e in entries.items() if e["type"] == "Family"}
    domains = {a for a, e in entries.items() if e["type"] == "Domain"}

    if families:
        if len(families) > 1:
            deepest = {a for a in families
                       if not any(a in ancestors(b) for b in families if b != a)}
        else:
            deepest = set(families)
        if len(deepest) == 1:
            acc = next(iter(deepest))
            return label(entries[acc], field), "interpro_family", "", "", ""
        # rule 1x: not rankable, so nothing is chosen
        candidates = "; ".join(f"{a} {label(entries[a], field)}"
                               for a in sorted(deepest))
        return HYPOTHETICAL, "none", candidates, "", ""

    if domains:
        order = sorted(domains, key=lambda a: (starts.get(a, 10 ** 9), a))
        domain_all = "; ".join(f"{a} {label(entries[a], field)}" for a in order)
        if len(order) == 1:
            product = (f"{domain_stem(label(entries[order[0]], field))} "
                       "domain-containing protein")
            return product, "interpro_domain", "", domain_all, ""
        product = (f"{domain_stem(label(entries[order[0]], field))} and "
                   f"{domain_stem(label(entries[order[1]], field))} "
                   "domain-containing protein")
        return product, "interpro_domain_multi", "", domain_all, ""

    if entries:
        # rule 3 -- only fold/repeat/site level evidence
        if all(e["type"] in STRUCTURAL_ONLY_TYPES for e in entries.values()):
            return HYPOTHETICAL, "none", "", "", "TRUE"
        raise SystemExit(
            "a protein has InterPro entries that no rule covers: types "
            f"{sorted({e['type'] for e in entries.values()})}")

    # rule 4
    return HYPOTHETICAL, "none", "", "", ""


def assign_dup_ids(loci, seed_by_locus):
    """One id per maximal run of adjacent loci sharing a seed ortholog."""
    ids = {}
    pairs = runs = 0
    i = 0
    while i < len(loci):
        j = i
        while (j + 1 < len(loci)
               and seed_by_locus.get(loci[j])
               and seed_by_locus.get(loci[j]) == seed_by_locus.get(loci[j + 1])):
            j += 1
        if j > i:
            runs += 1
            pairs += j - i
            for k in range(i, j + 1):
                ids[loci[k]] = f"DUP{runs:04d}"
            i = j + 1
        else:
            i += 1
    return ids, pairs, runs


# --------------------------------------------------------------------------- #

def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gene-table", type=Path, default=DEFAULT_TABLE)
    ap.add_argument("--lookup", type=Path, default=DEFAULT_LOOKUP)
    ap.add_argument("--matches", type=Path, default=DEFAULT_MATCHES)
    ap.add_argument("--ips6", type=Path, default=DEFAULT_IPS6)
    ap.add_argument("--eggnog", type=Path, default=DEFAULT_EGGNOG)
    ap.add_argument("--tree", type=Path, default=DEFAULT_TREE)
    ap.add_argument("--entry-label", choices=("name", "description"),
                    default="description",
                    help="which InterPro entry field the wording uses")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    header, rows = read_tsv(args.gene_table)
    for column in NEW_COLUMNS:
        if column in header:
            raise SystemExit(f"{args.gene_table} already has a {column!r} column; "
                             "this script is not meant to run twice in place")
    i_lt, i_sym, i_prod = (header.index("locus_tag"), header.index("symbol"),
                           header.index("product"))

    lookup_header, lookup_rows = read_tsv(args.lookup)
    lookup = [dict(zip(lookup_header, r)) for r in lookup_rows]
    ancestors = load_hierarchy(args.tree)
    _, eggnog = load_eggnog(args.eggnog)
    interpro = load_interpro_per_protein(lookup, args.matches, args.ips6)

    loci = [r[i_lt] for r in rows]
    seed = {q: v["seed_ortholog"] for q, v in eggnog.items()
            if v["seed_ortholog"] not in ("", "-")}
    dup_ids, dup_pairs, dup_runs = assign_dup_ids(loci, seed)

    original_symbols = {r[i_lt]: r[i_sym] for r in rows}
    out_rows = []
    sym_src = collections.Counter()
    prod_src = collections.Counter()
    rejected = collections.Counter()
    too_long = []

    for row in rows:
        locus = row[i_lt]
        # -- symbol ---------------------------------------------------------- #
        symbol = row[i_sym]
        symbol_rejected = symbol_reason = ""
        if symbol:
            # rule 1: a hand-curated value is never touched
            symbol_source = "utex250a"
        else:
            preferred = (eggnog.get(locus, {}).get("Preferred_name") or "").strip()
            if preferred and preferred != "-":
                reason = symbol_rejection(preferred)
                if reason:
                    symbol_rejected, symbol_reason = preferred, reason
                    symbol, symbol_source = "", "none"
                    rejected[reason] += 1
                else:
                    symbol, symbol_source = preferred, "eggnog"
            else:
                symbol_source = "none"
        sym_src[symbol_source] += 1

        # -- product --------------------------------------------------------- #
        _source, entries, starts = interpro[locus]
        product, product_source, candidates, domain_all, structural = \
            decide_product(entries, starts, ancestors, args.entry_label)
        prod_src[product_source] += 1
        if len(product) > MAX_PRODUCT:
            too_long.append((locus, len(product), product))

        new = list(row)
        new[i_sym] = symbol
        new[i_prod] = product
        new += [symbol_source, product_source, dup_ids.get(locus, ""),
                candidates, "TRUE" if candidates else "",
                domain_all, structural,
                symbol_rejected, symbol_reason,
                "TRUE" if locus in eggnog else "FALSE"]
        out_rows.append(new)

    # -- stop conditions ----------------------------------------------------- #
    if too_long:
        raise SystemExit(
            f"{len(too_long)} product string(s) exceed {MAX_PRODUCT} characters; "
            "stopping instead of truncating. Examples:\n  "
            + "\n  ".join(f"{l} ({n}) {p}" for l, n, p in too_long[:3]))
    changed = [l for l, s in original_symbols.items()
               if s and out_rows[loci.index(l)][i_sym] != s]
    if changed:
        raise SystemExit(f"{len(changed)} existing symbol(s) were altered "
                         f"(first: {changed[:5]}); overwriting is forbidden")
    if sum(prod_src.values()) != len(rows):
        raise SystemExit("some genes matched no rule")

    if args.dry_run:
        print("dry run: nothing written")
    else:
        with args.gene_table.open("w") as fh:
            fh.write("\t".join(header + NEW_COLUMNS) + "\n")
            for row in out_rows:
                fh.write("\t".join(row) + "\n")
        print(f"wrote {args.gene_table} ({len(out_rows)} rows + header)")

    total = len(rows)
    print(f"entry label field           : {args.entry_label}")
    print(f"symbol_source               : "
          + ", ".join(f"{k} {v}" for k, v in sorted(sym_src.items())))
    if rejected:
        print("  eggNOG Preferred_name rejected as not a gene symbol: "
              f"{sum(rejected.values())}")
        for reason, n in sorted(rejected.items()):
            print(f"    {reason:<20} {n}")
    print(f"eggnog_data_row TRUE        : "
          f"{sum(1 for r in out_rows if r[-1] == 'TRUE')}")
    print(f"product_source              : "
          + ", ".join(f"{k} {v}" for k, v in sorted(prod_src.items())))
    kept = sum(1 for r in out_rows if r[i_prod] == HYPOTHETICAL)
    print(f"product still hypothetical  : {kept} / {total} "
          f"({100 * kept / total:.2f} %)")
    filled = sum(1 for r in out_rows if r[i_sym])
    print(f"symbol filled               : {filled} / {total} "
          f"({100 * filled / total:.2f} %)")
    print(f"dup runs / adjacent pairs   : {dup_runs} runs, {dup_pairs} pairs, "
          f"{sum(1 for r in out_rows if r[header.index('locus_tag')] in dup_ids)} rows")
    print(f"existing symbols altered    : 0 (verified over "
          f"{sum(1 for v in original_symbols.values() if v)})")


if __name__ == "__main__":
    main()
