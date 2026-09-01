#!/usr/bin/env python3
"""Fill in the GO terms the InterPro Matches API does not return.

The Matches API carries no GO field at all — the OpenAPI schema defines none and
the 75 raw batches contain zero ``GO:`` strings — so GO has to come from
InterPro's own ``interpro2go`` mapping, applied to the IPR accessions the API
did return.

What this writes:

  step0_out/lookup_status.tsv   gains ``n_interpro_gos`` and ``interpro_gos``
  step0_out/SUMMARY.md          gains (or has replaced) one marked section with
                                the coverage and the overlap against eggNOG

The eggNOG ``GOs`` column is read only for the comparison.  It lives in
``annotation/UTEX25_gene_table_eggnog.tsv`` and is never touched, and the two
GO sources stay in separate columns in separate files — neither overwrites the
other.

    curl -fSL -o raw/interpro2go \\
      "https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/interpro2go"
    python3 scripts/06_interpro2go.py
"""

from __future__ import annotations

import argparse
import collections
import gzip
import json
import re
import statistics
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DEFAULT_MAP = ROOT / "raw" / "interpro2go"
DEFAULT_MATCHES = ROOT / "step0_out" / "matches_raw.json.gz"
DEFAULT_LOOKUP = ROOT / "step0_out" / "lookup_status.tsv"
DEFAULT_EGGNOG_TABLE = ROOT / "annotation" / "UTEX25_gene_table_eggnog.tsv"
DEFAULT_SUMMARY = ROOT / "step0_out" / "SUMMARY.md"

BEGIN = "<!-- BEGIN interpro2go -->"
END = "<!-- END interpro2go -->"

# InterPro:IPR000003 Retinoid X receptor/HNF4 > GO:DNA binding ; GO:0003677
IP2GO_LINE = re.compile(r"^InterPro:(?P<ipr>IPR\d{6})\s.*;\s*(?P<go>GO:\d{7})\s*$")
IPR = re.compile(r"^IPR\d{6}$")


def load_interpro2go(path: Path) -> dict[str, set[str]]:
    mapping: dict[str, set[str]] = collections.defaultdict(set)
    parsed = skipped = 0
    for line in open(path):
        line = line.rstrip("\n")
        if not line or line.startswith("!"):
            continue
        m = IP2GO_LINE.match(line)
        if not m:
            skipped += 1
            continue
        mapping[m.group("ipr")].add(m.group("go"))
        parsed += 1
    if skipped:
        raise SystemExit(
            f"{path}: {skipped} non-comment line(s) did not parse as an "
            "interpro2go mapping; refusing to run on a file whose format is "
            "not the one this script understands")
    if not parsed:
        raise SystemExit(f"{path}: no mappings parsed")
    return dict(mapping)


def load_ipr_per_protein(matches: Path, lookup: Path) -> dict[str, set[str]]:
    """protein_id -> set of IPR accessions from the Matches API response."""
    with gzip.open(matches, "rt", encoding="utf-8") as fh:
        raw = json.load(fh)
    by_md5 = {}
    for batch in raw["batches"]:
        if batch.get("status") != 200:
            raise SystemExit(f"{matches}: batch {batch.get('batch')} is not 200")
        for result in batch["body"]["results"]:
            by_md5[result["md5"].upper()] = result

    out = {}
    with open(lookup) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        i_pid, i_md5 = header.index("protein_id"), header.index("md5")
        i_up = header.index("in_uniparc")
        for line in fh:
            row = line.rstrip("\n").split("\t")
            if row[i_up] != "yes":
                continue
            accessions = set()
            for match in by_md5[row[i_md5].upper()]["matches"]:
                entry = match["signature"].get("entry")
                if entry and IPR.match(entry["accession"]):
                    accessions.add(entry["accession"])
            out[row[i_pid]] = accessions
    return out


def load_eggnog_gos(path: Path, key="locus_tag", column="eggnog_GOs"):
    """protein_id -> set of GO ids from the merged eggNOG table (read-only)."""
    if not path.exists():
        return None
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        if key not in header or column not in header:
            raise SystemExit(f"{path}: expected {key!r} and {column!r} columns")
        i_key, i_go = header.index(key), header.index(column)
        out = {}
        for line in fh:
            row = line.rstrip("\n").split("\t")
            value = row[i_go].strip() if i_go < len(row) else ""
            out[row[i_key]] = ({g for g in value.split(",") if g.strip()}
                               if value else set())
    return out


def rewrite_lookup(path: Path, gos: dict[str, set[str]]) -> tuple[int, int]:
    """Add n_interpro_gos / interpro_gos, replacing them if already present."""
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        rows = [line.rstrip("\n").split("\t") for line in fh if line.strip()]

    added = [c for c in ("n_interpro_gos", "interpro_gos") if c not in header]
    keep = [i for i, c in enumerate(header)
            if c not in ("n_interpro_gos", "interpro_gos")]
    new_header = [header[i] for i in keep] + ["n_interpro_gos", "interpro_gos"]
    i_pid = header.index("protein_id")

    n_with = 0
    with open(path, "w") as fh:
        fh.write("\t".join(new_header) + "\n")
        for row in rows:
            terms = sorted(gos.get(row[i_pid], ()))
            if terms:
                n_with += 1
            fh.write("\t".join([row[i] for i in keep]
                               + [str(len(terms)) if terms else "",
                                  ",".join(terms)]) + "\n")
    return len(rows), n_with, added


def replace_section(path: Path, body: str) -> None:
    text = path.read_text() if path.exists() else ""
    block = f"{BEGIN}\n{body}{END}\n"
    if BEGIN in text and END in text:
        head, _, rest = text.partition(BEGIN)
        _, _, tail = rest.partition(END)
        path.write_text(head + block + tail.lstrip("\n"))
    else:
        if text and not text.endswith("\n"):
            text += "\n"
        path.write_text(text + "\n" + block)


def q(values):
    ordered = sorted(values)
    if not ordered:
        return None
    if len(ordered) == 1:
        v = float(ordered[0])
        return v, v, v
    return statistics.quantiles(ordered, n=4, method="inclusive")


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--interpro2go", type=Path, default=DEFAULT_MAP)
    ap.add_argument("--matches", type=Path, default=DEFAULT_MATCHES)
    ap.add_argument("--lookup", type=Path, default=DEFAULT_LOOKUP)
    ap.add_argument("--eggnog-table", type=Path, default=DEFAULT_EGGNOG_TABLE)
    ap.add_argument("--summary", type=Path, default=DEFAULT_SUMMARY)
    args = ap.parse_args()

    if not args.interpro2go.exists():
        raise SystemExit(
            f"{args.interpro2go} not found. Fetch it first:\n"
            "  curl -fSL -o raw/interpro2go "
            "\"https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/"
            "interpro2go\"")

    ip2go = load_interpro2go(args.interpro2go)
    per_protein_ipr = load_ipr_per_protein(args.matches, args.lookup)

    gos = {}
    ipr_seen, ipr_mapped = set(), set()
    for pid, accessions in per_protein_ipr.items():
        ipr_seen |= accessions
        terms = set()
        for accession in accessions:
            hit = ip2go.get(accession)
            if hit:
                ipr_mapped.add(accession)
                terms |= hit
        gos[pid] = terms

    n_rows, n_with_go, added = rewrite_lookup(args.lookup, gos)

    eggnog = load_eggnog_gos(args.eggnog_table)
    body = build_section(args, ip2go, per_protein_ipr, gos, eggnog,
                         ipr_seen, ipr_mapped, n_rows, n_with_go)
    replace_section(args.summary, body)

    print(f"interpro2go mappings   : {sum(len(v) for v in ip2go.values())} "
          f"over {len(ip2go)} IPR")
    print(f"proteins with IPR      : {sum(1 for v in per_protein_ipr.values() if v)}")
    print(f"proteins with GO       : {n_with_go} / {n_rows}")
    print(f"{args.lookup}: columns {added or '(already present, refreshed)'}")
    print(f"{args.summary}: interpro2go section written")


def build_section(args, ip2go, per_protein_ipr, gos, eggnog, ipr_seen,
                  ipr_mapped, n_rows, n_with_go) -> str:
    out = []
    w = out.append
    total = n_rows

    w("## 8. GO の補完（interpro2go）")
    w("")
    w("Matches API は GO を返さないため、返ってきた IPR アクセッションに")
    w(f"InterPro 公式の `interpro2go` を適用して補った（`{args.interpro2go.relative_to(ROOT)}`）。")
    w("")
    w(f"- `interpro2go` の対応: {sum(len(v) for v in ip2go.values())} 対応 / "
      f"{len(ip2go)} IPR エントリ")
    w(f"- 本データに出現した IPR: {len(ipr_seen)} 種類、うち GO が引けたもの "
      f"{len(ipr_mapped)} 種類 "
      f"({100 * len(ipr_mapped) / len(ipr_seen):.2f} %)")
    w("")
    w(f"`lookup_status.tsv` に `n_interpro_gos` / `interpro_gos` 列を追加した。")
    w("")

    n_ipr = sum(1 for v in per_protein_ipr.values() if v)
    counts = [len(v) for v in gos.values() if v]
    w("| 指標 | 件数 | 全 " + str(total) + " に対して |")
    w("|---|---|---|")
    w(f"| IPR が 1 つ以上あるタンパク質 | {n_ipr} | {100 * n_ipr / total:.2f} % |")
    w(f"| **GO が 1 つ以上引けたタンパク質** | **{n_with_go}** | "
      f"**{100 * n_with_go / total:.2f} %** |")
    w("")
    if counts:
        q1, med, q3 = q(counts)
        w(f"GO が引けた {len(counts)} タンパク質あたりの GO 数: "
          f"Q1 {q1:g} / 中央値 {med:g} / Q3 {q3:g} / 最大 {max(counts)}")
        w("")

    if eggnog is None:
        w(f"eggNOG の GO 列（`{args.eggnog_table.relative_to(ROOT)}`）が見つからないため、")
        w("重なりの集計は行っていない。")
        w("")
        return "\n".join(out) + "\n"

    ids = sorted(eggnog)
    ip_set = {p: gos.get(p, set()) for p in ids}
    eg_set = eggnog

    n_ip = sum(1 for p in ids if ip_set[p])
    n_eg = sum(1 for p in ids if eg_set[p])
    both = [p for p in ids if ip_set[p] and eg_set[p]]
    only_ip = sum(1 for p in ids if ip_set[p] and not eg_set[p])
    only_eg = sum(1 for p in ids if eg_set[p] and not ip_set[p])
    neither = sum(1 for p in ids if not ip_set[p] and not eg_set[p])

    w("### 8.1 eggNOG の GO 列との比較（別列として併存、上書きしていない）")
    w("")
    w(f"- interpro2go 由来: `step0_out/lookup_status.tsv` の `interpro_gos`")
    w(f"- eggNOG 由来: `{args.eggnog_table.relative_to(ROOT)}` の `eggnog_GOs`")
    w("")
    w("| 区分 | タンパク質数 | 全 " + str(len(ids)) + " に対して |")
    w("|---|---|---|")
    w(f"| interpro2go に GO あり | {n_ip} | {100 * n_ip / len(ids):.2f} % |")
    w(f"| eggNOG に GO あり | {n_eg} | {100 * n_eg / len(ids):.2f} % |")
    w(f"| 両方にあり | {len(both)} | {100 * len(both) / len(ids):.2f} % |")
    w(f"| interpro2go のみ | {only_ip} | {100 * only_ip / len(ids):.2f} % |")
    w(f"| eggNOG のみ | {only_eg} | {100 * only_eg / len(ids):.2f} % |")
    w(f"| どちらにも無し | {neither} | {100 * neither / len(ids):.2f} % |")
    w(f"| **少なくとも一方にあり** | **{len(ids) - neither}** | "
      f"**{100 * (len(ids) - neither) / len(ids):.2f} %** |")
    w("")

    shared = sum(len(ip_set[p] & eg_set[p]) for p in both)
    ip_only = sum(len(ip_set[p] - eg_set[p]) for p in both)
    eg_only = sum(len(eg_set[p] - ip_set[p]) for p in both)
    w(f"両方に GO がある {len(both)} タンパク質での GO 項目単位の重なり:")
    w("")
    w("| 区分 | GO 項目の延べ数 |")
    w("|---|---|")
    w(f"| 両方にある | {shared} |")
    w(f"| interpro2go のみ | {ip_only} |")
    w(f"| eggNOG のみ | {eg_only} |")
    w("")
    jaccard = [len(ip_set[p] & eg_set[p]) / len(ip_set[p] | eg_set[p])
               for p in both]
    if jaccard:
        j1, jm, j3 = q(jaccard)
        w(f"タンパク質ごとの Jaccard 係数: Q1 {j1:.3f} / 中央値 {jm:.3f} / Q3 {j3:.3f}")
        w(f"（完全一致 {sum(1 for p in both if ip_set[p] == eg_set[p])} タンパク質）")
        w("")
    w("2 つの GO 集合は由来が異なる。`interpro2go` は InterPro エントリに対して")
    w("キュレーションされた直接の対応であり、eggNOG の `GOs` は seed ortholog の")
    w("GO を転写したもので、`go_namespace_split=True` の設定で出力されている。")
    w("重なりの数値はこの違いを込みで見る必要がある。")
    w("")
    return "\n".join(out) + "\n"


if __name__ == "__main__":
    main()
