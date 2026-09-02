#!/usr/bin/env python3
"""Coverage report for annotation/UTEX25_interpro.tsv.

Writes docs/INTERPRO_COVERAGE.md.  Numbers only — the report draws no
conclusions.

Two InterPro counts are reported side by side because the two sources do not
cover the same libraries:

  共通7ライブラリ  IPR entries reachable from Pfam, NCBIFAM, PROSITE, SMART,
                   CDD, SUPERFAMILY — the libraries both routes ran, so the
                   figure means the same thing for both
  全ライブラリ      the above plus the entries only the Matches API could
                   supply (PANTHER, CATH-Gene3D, ...); identical to the first
                   for every ips6_local row

    python3 scripts/12_interpro_coverage.py
"""

from __future__ import annotations

import argparse
import collections
import statistics
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DEFAULT_INTERPRO = ROOT / "annotation" / "UTEX25_interpro.tsv"
DEFAULT_LOOKUP = ROOT / "step0_out" / "lookup_status.tsv"
DEFAULT_TABLE = ROOT / "annotation" / "UTEX25_gene_table.tsv"
DEFAULT_EGGNOG = ROOT / "annotation" / "UTEX25_gene_table_eggnog.tsv"
DEFAULT_OUT = ROOT / "docs" / "INTERPRO_COVERAGE.md"

LIBRARY_COLUMNS = ["pfam_accessions", "ncbifam", "prosite", "smart", "cdd",
                   "superfamily"]


def read_tsv(path: Path):
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(header, line.rstrip("\n").split("\t")))
                for line in fh if line.strip()]


def items(value: str) -> list[str]:
    value = (value or "").strip()
    return [v for v in value.split(",") if v] if value else []


def quartiles(values):
    ordered = sorted(values)
    if not ordered:
        return None
    if len(ordered) == 1:
        v = ordered[0]
        return dict(n=1, min=v, q1=v, med=v, q3=v, max=v, mean=float(v))
    q1, med, q3 = statistics.quantiles(ordered, n=4, method="inclusive")
    return dict(n=len(ordered), min=ordered[0], q1=q1, med=med, q3=q3,
                max=ordered[-1], mean=statistics.mean(ordered))


def rate(n, d):
    return f"{100 * n / d:.2f} %" if d else "n/a"


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--interpro", type=Path, default=DEFAULT_INTERPRO)
    ap.add_argument("--lookup", type=Path, default=DEFAULT_LOOKUP)
    ap.add_argument("--gene-table", type=Path, default=DEFAULT_TABLE)
    ap.add_argument("--eggnog-table", type=Path, default=DEFAULT_EGGNOG)
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = ap.parse_args()

    rows = read_tsv(args.interpro)
    lookup = {r["protein_id"]: r for r in read_tsv(args.lookup)}
    genes = {g["locus_tag"]: g for g in read_tsv(args.gene_table)}
    eggnog = {g["locus_tag"]: g for g in read_tsv(args.eggnog_table)} \
        if args.eggnog_table.exists() else {}

    total = len(rows)
    for r in rows:
        r["_ipr"] = set(items(r["interpro_accessions"]))
        r["_ipr_all"] = r["_ipr"] | set(items(
            r.get("interpro_accessions_other_libs", "")))
        r["_pfam"] = set(items(r["pfam_accessions"]))
        r["_n"] = int(r["n_domains"])
        r["_qc"] = genes[r["locus_tag"]]["QC"]
        r["_stop"] = lookup[r["locus_tag"]]["has_internal_stop"]
        r["_len"] = int(lookup[r["locus_tag"]]["length"])

    by_source = collections.defaultdict(list)
    for r in rows:
        by_source[r["source"]].append(r)

    fh = args.out.open("w")
    w = lambda s="": fh.write(s + "\n")

    def cover(subset, key):
        return sum(1 for r in subset if r[key])

    w("# InterPro カバレッジ")
    w()
    w(f"- 入力: `{args.interpro.relative_to(ROOT)}`（{total} 行）")
    w(f"- QC / 長さ: `{args.gene_table.relative_to(ROOT)}`、"
      f"`{args.lookup.relative_to(ROOT)}`")
    w(f"- 生成: `scripts/12_interpro_coverage.py`")
    w()
    w("実測値のみ。解釈は書いていない。")
    w()
    w("**2 通りの InterPro 数を併記している。** ローカル実行は 7 ライブラリ、")
    w("Matches API は 18 ライブラリを走らせているため、両者で意味を揃えた")
    w("「共通 7 ライブラリ」と、API 側だけが持つ分を含めた「全ライブラリ」を分けてある。")
    w("`ips6_local` 行では両者は必ず一致する。")
    w()

    # -- 1. overall --------------------------------------------------------- #
    w(f"## 1. 全体（{total} 配列）")
    w()
    w("| 指標 | 件数 | / " + str(total) + " |")
    w("|---|---|---|")
    n_ipr = cover(rows, "_ipr")
    n_ipr_all = cover(rows, "_ipr_all")
    n_pfam = cover(rows, "_pfam")
    n_any = sum(1 for r in rows if r["_n"] > 0)
    w(f"| InterPro エントリが 1 つ以上（共通 7 ライブラリ） | {n_ipr} | "
      f"{rate(n_ipr, total)} |")
    w(f"| InterPro エントリが 1 つ以上（全ライブラリ） | {n_ipr_all} | "
      f"{rate(n_ipr_all, total)} |")
    w(f"| Pfam ドメインが 1 つ以上 | {n_pfam} | {rate(n_pfam, total)} |")
    w(f"| 共通 7 ライブラリのシグネチャが 1 つ以上 | {n_any} | "
      f"{rate(n_any, total)} |")
    w()

    # -- 2. by source ------------------------------------------------------- #
    w("## 2. 取得経路別")
    w()
    w("| 指標 | matches_api | ips6_local | 合計 |")
    w("|---|---|---|---|")
    api, local = by_source.get("matches_api", []), by_source.get("ips6_local", [])
    w(f"| 配列数 | {len(api)} | {len(local)} | {total} |")
    for label, key in (("InterPro ≥1（共通 7 ライブラリ）", "_ipr"),
                       ("InterPro ≥1（全ライブラリ）", "_ipr_all"),
                       ("Pfam ≥1", "_pfam")):
        a, b = cover(api, key), cover(local, key)
        w(f"| {label} | {a} ({rate(a, len(api))}) | {b} ({rate(b, len(local))}) "
          f"| {a + b} ({rate(a + b, total)}) |")
    a = sum(1 for r in api if r["_n"] > 0)
    b = sum(1 for r in local if r["_n"] > 0)
    w(f"| シグネチャ ≥1（共通 7 ライブラリ） | {a} ({rate(a, len(api))}) | "
      f"{b} ({rate(b, len(local))}) | {a + b} ({rate(a + b, total)}) |")
    w()
    for label, subset in (("matches_api", api), ("ips6_local", local)):
        counts = [r["_n"] for r in subset if r["_n"] > 0]
        q = quartiles(counts)
        if q:
            w(f"- {label}: シグネチャが付いた {q['n']} 配列あたりの数は "
              f"中央値 {q['med']:g}（Q1 {q['q1']:g} / Q3 {q['q3']:g} / "
              f"最大 {q['max']}）")
    w()

    # -- 3. per library ----------------------------------------------------- #
    w("## 3. ライブラリ別（共通 7 ライブラリ）")
    w()
    w("| 列 | matches_api | ips6_local | 合計 | / " + str(total) + " |")
    w("|---|---|---|---|---|")
    for column in LIBRARY_COLUMNS:
        a = sum(1 for r in api if items(r[column]))
        b = sum(1 for r in local if items(r[column]))
        w(f"| `{column}` | {a} | {b} | {a + b} | {rate(a + b, total)} |")
    w()

    # -- 4. QC -------------------------------------------------------------- #
    w("## 4. QC フラグ別")
    w()
    groups = collections.defaultdict(list)
    for r in rows:
        groups[r["_qc"]].append(r)
    w("| QC | 配列数 | InterPro ≥1（共通7） | 率 | Pfam ≥1 | 率 |")
    w("|---|---|---|---|---|---|")
    for name, subset in sorted(groups.items(), key=lambda kv: -len(kv[1])):
        i, p = cover(subset, "_ipr"), cover(subset, "_pfam")
        w(f"| `{name}` | {len(subset)} | {i} | {rate(i, len(subset))} | {p} | "
          f"{rate(p, len(subset))} |")
    w(f"| **計** | **{total}** | {n_ipr} | {rate(n_ipr, total)} | {n_pfam} | "
      f"{rate(n_pfam, total)} |")
    w()

    # -- 5. internal stop --------------------------------------------------- #
    w("## 5. 内部終止コドンの有無別")
    w()
    w("`has_internal_stop = yes` の 235 配列は `*` を `X` に置換して投入した"
      "（`docs/IPS6_FEASIBILITY.md` §5.2）。")
    w()
    w("| has_internal_stop | 配列数 | InterPro ≥1（共通7） | 率 | Pfam ≥1 | 率 |")
    w("|---|---|---|---|---|---|")
    for flag in ("no", "yes"):
        subset = [r for r in rows if r["_stop"] == flag]
        i, p = cover(subset, "_ipr"), cover(subset, "_pfam")
        w(f"| {flag} | {len(subset)} | {i} | {rate(i, len(subset))} | {p} | "
          f"{rate(p, len(subset))} |")
    w()

    # -- 6. vs eggNOG ------------------------------------------------------- #
    if eggnog:
        w("## 6. eggNOG の PFAMs カバレッジとの並置")
        w()
        w("同じ 7,413 配列に対する 2 つの独立した Pfam 付与の件数を並べたもの。"
          "数値の並置のみで、比較の解釈は `docs/PFAM_CONCORDANCE.md` を参照。")
        w()
        eg = sum(1 for g in eggnog.values() if items(g.get("eggnog_PFAMs", "")))
        both = sum(1 for r in rows
                   if r["_pfam"] and items(
                       eggnog.get(r["locus_tag"], {}).get("eggnog_PFAMs", "")))
        only_ip = sum(1 for r in rows
                      if r["_pfam"] and not items(
                          eggnog.get(r["locus_tag"], {}).get("eggnog_PFAMs", "")))
        only_eg = sum(1 for r in rows
                      if not r["_pfam"] and items(
                          eggnog.get(r["locus_tag"], {}).get("eggnog_PFAMs", "")))
        neither = total - both - only_ip - only_eg
        w("| 由来 | Pfam が 1 つ以上付いた配列 | / " + str(total) + " |")
        w("|---|---|---|")
        w(f"| InterPro（本表、`pfam_accessions`） | {n_pfam} | "
          f"{rate(n_pfam, total)} |")
        w(f"| eggNOG（`annotation/UTEX25_gene_table_eggnog.tsv` の "
          f"`eggnog_PFAMs`） | {eg} | {rate(eg, total)} |")
        w()
        w("| 区分 | 配列数 | / " + str(total) + " |")
        w("|---|---|---|")
        w(f"| 両方に Pfam あり | {both} | {rate(both, total)} |")
        w(f"| InterPro のみ | {only_ip} | {rate(only_ip, total)} |")
        w(f"| eggNOG のみ | {only_eg} | {rate(only_eg, total)} |")
        w(f"| どちらにも無し | {neither} | {rate(neither, total)} |")
        w(f"| **少なくとも一方** | **{total - neither}** | "
          f"**{rate(total - neither, total)}** |")
        w()

    # -- 7. still nothing --------------------------------------------------- #
    none = [r for r in rows if r["_n"] == 0 and not r["_ipr_all"]]
    w("## 7. ドメインが 1 つも付かなかった配列")
    w()
    w(f"共通 7 ライブラリのシグネチャも、どのライブラリ由来の InterPro エントリも"
      f"付かなかった配列: **{len(none)} / {total} ({rate(len(none), total)})**")
    w()
    if none:
        w("### QC フラグ別")
        w()
        w("| QC | 件数 | 同 QC の全配列 | 割合 |")
        w("|---|---|---|---|")
        sub = collections.Counter(r["_qc"] for r in none)
        for name, count in sorted(sub.items(), key=lambda kv: (-kv[1], kv[0])):
            w(f"| `{name}` | {count} | {len(groups[name])} | "
              f"{rate(count, len(groups[name]))} |")
        w(f"| **計** | **{len(none)}** | {total} | {rate(len(none), total)} |")
        w()
        w("### 取得経路別")
        w()
        w("| source | 件数 | 同経路の全配列 | 割合 |")
        w("|---|---|---|---|")
        for name in ("matches_api", "ips6_local"):
            count = sum(1 for r in none if r["source"] == name)
            w(f"| `{name}` | {count} | {len(by_source[name])} | "
              f"{rate(count, len(by_source[name]))} |")
        w()
        w("### 長さの分布（アミノ酸残基数）")
        w()
        w("| 集合 | n | 最小 | Q1 | 中央値 | Q3 | 最大 | 平均 |")
        w("|---|---|---|---|---|---|---|---|")
        for label, subset in (("ドメインなし", none), ("全 7,413 配列", rows)):
            q = quartiles([r["_len"] for r in subset])
            w(f"| {label} | {q['n']} | {q['min']} | {q['q1']:.1f} | "
              f"{q['med']:.1f} | {q['q3']:.1f} | {q['max']} | {q['mean']:.1f} |")
        w()
        w("### 内部終止コドンの有無")
        w()
        w("| has_internal_stop | 件数 |")
        w("|---|---|")
        stops = collections.Counter(r["_stop"] for r in none)
        for flag in ("no", "yes"):
            w(f"| {flag} | {stops.get(flag, 0)} |")
        w()
    fh.close()

    print(f"wrote {args.out}")
    print(f"InterPro >=1 (shared 7 libs): {n_ipr}/{total} ({rate(n_ipr, total)})")
    print(f"InterPro >=1 (all libs)    : {n_ipr_all}/{total} "
          f"({rate(n_ipr_all, total)})")
    print(f"Pfam >=1                   : {n_pfam}/{total} ({rate(n_pfam, total)})")
    print(f"nothing at all             : {len(none)}/{total} "
          f"({rate(len(none), total)})")


if __name__ == "__main__":
    main()
