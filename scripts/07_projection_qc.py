#!/usr/bin/env python3
"""Cross-tabulate projection quality (identity, QC flags) against UniParc membership.

Joins ``annotation/UTEX25_gene_table.tsv`` (structural annotation) with
``step0_out/lookup_status.tsv`` (InterPro Matches API result) on the locus tag
and writes ``docs/PROJECTION_QC.md``.

Numbers and tables only — the report states no interpretation.

    python3 scripts/07_projection_qc.py
"""

from __future__ import annotations

import argparse
import collections
import statistics
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DEFAULT_TABLE = ROOT / "annotation" / "UTEX25_gene_table.tsv"
DEFAULT_LOOKUP = ROOT / "step0_out" / "lookup_status.tsv"
DEFAULT_OUT = ROOT / "docs" / "PROJECTION_QC.md"

# Bin edges for the identity distribution, as [low, high) except the last.
BINS = [(1.00, 1.01, "= 1.0000"),
        (0.99, 1.00, "0.9900 – 0.9999"),
        (0.95, 0.99, "0.9500 – 0.9899"),
        (0.90, 0.95, "0.9000 – 0.9499"),
        (0.80, 0.90, "0.8000 – 0.8999"),
        (0.00, 0.80, "< 0.8000")]


def read_tsv(path: Path):
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        rows = [dict(zip(header, line.rstrip("\n").split("\t")))
                for line in fh if line.strip()]
    return header, rows


def stats(values):
    if not values:
        return None
    ordered = sorted(values)
    if len(ordered) == 1:
        v = ordered[0]
        return dict(n=1, min=v, q1=v, med=v, q3=v, max=v, mean=v)
    q1, med, q3 = statistics.quantiles(ordered, n=4, method="inclusive")
    return dict(n=len(ordered), min=ordered[0], q1=q1, med=med, q3=q3,
                max=ordered[-1], mean=statistics.mean(ordered))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gene-table", type=Path, default=DEFAULT_TABLE)
    ap.add_argument("--lookup", type=Path, default=DEFAULT_LOOKUP)
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--key", default="locus_tag")
    args = ap.parse_args()

    _, genes = read_tsv(args.gene_table)
    _, lookup = read_tsv(args.lookup)

    status = {r["protein_id"]: r["in_uniparc"] for r in lookup}
    missing = [g[args.key] for g in genes if g[args.key] not in status]
    if missing:
        raise SystemExit(
            f"{len(missing)} locus_tag(s) in {args.gene_table.name} have no row "
            f"in {args.lookup.name} (first: {missing[:5]})")

    rows = []
    for g in genes:
        rows.append({
            "id": g[args.key],
            "identity": float(g["identity"]),
            "qc": g["QC"],
            "hit": status[g[args.key]] == "yes",
        })
    total = len(rows)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    fh = args.out.open("w")
    w = lambda s="": fh.write(s + "\n")

    w("# 投影品質 × UniParc 収録のクロス集計")
    w()
    w(f"- 構造アノテーション: `{args.gene_table.relative_to(ROOT)}`"
      f"（`identity` 列, `QC` 列）")
    w(f"- UniParc 収録: `{args.lookup.relative_to(ROOT)}`（`in_uniparc` 列）")
    w(f"- 生成: `scripts/07_projection_qc.py`")
    w(f"- 対象: {total} 遺伝子（全件が両表に存在）")
    w()
    w("実測値のみ。解釈は書いていない。")
    w()

    # -- 1. 2x2 --------------------------------------------------------------
    exact = [r for r in rows if r["identity"] >= 1.0]
    below = [r for r in rows if r["identity"] < 1.0]
    cells = {
        ("=100%", "yes"): sum(1 for r in exact if r["hit"]),
        ("=100%", "no"): sum(1 for r in exact if not r["hit"]),
        ("<100%", "yes"): sum(1 for r in below if r["hit"]),
        ("<100%", "no"): sum(1 for r in below if not r["hit"]),
    }
    w("## 1. identity = 100% / < 100% × in_uniparc")
    w()
    w("| identity | in_uniparc = yes | in_uniparc = no | 計 | ヒット率 |")
    w("|---|---|---|---|---|")
    for label, group in (("= 100% (1.0000)", exact), ("< 100%", below)):
        key = "=100%" if group is exact else "<100%"
        y, n = cells[(key, "yes")], cells[(key, "no")]
        rate = f"{100 * y / (y + n):.2f} %" if y + n else "n/a"
        w(f"| {label} | {y} | {n} | {y + n} | {rate} |")
    ty = cells[("=100%", "yes")] + cells[("<100%", "yes")]
    tn = cells[("=100%", "no")] + cells[("<100%", "no")]
    w(f"| **計** | **{ty}** | **{tn}** | **{ty + tn}** | "
      f"**{100 * ty / (ty + tn):.2f} %** |")
    w()
    w("依頼された 2 つの数値:")
    w()
    w(f"- **identity = 100% かつ in_uniparc = no: {cells[('=100%', 'no')]} 件** "
      f"（identity = 100% の {len(exact)} 件の "
      f"{100 * cells[('=100%', 'no')] / len(exact):.2f} %、"
      f"全 {total} 件の {100 * cells[('=100%', 'no')] / total:.2f} %）")
    w(f"- **identity < 100% かつ in_uniparc = yes: {cells[('<100%', 'yes')]} 件** "
      f"（identity < 100% の {len(below)} 件の "
      f"{100 * cells[('<100%', 'yes')] / len(below):.2f} %、"
      f"全 {total} 件の {100 * cells[('<100%', 'yes')] / total:.2f} %）")
    w()

    # -- 2. identity distribution by in_uniparc ------------------------------
    w("## 2. identity の分布（in_uniparc 別）")
    w()
    w("| in_uniparc | n | 最小 | Q1 | 中央値 | Q3 | 最大 | 平均 |")
    w("|---|---|---|---|---|---|---|---|")
    for label, subset in (("yes", [r for r in rows if r["hit"]]),
                          ("no", [r for r in rows if not r["hit"]]),
                          ("全体", rows)):
        s = stats([r["identity"] for r in subset])
        w(f"| {label} | {s['n']} | {s['min']:.4f} | {s['q1']:.4f} | "
          f"{s['med']:.4f} | {s['q3']:.4f} | {s['max']:.4f} | {s['mean']:.4f} |")
    w()

    # -- 3. identity bins ----------------------------------------------------
    w("## 3. identity 区間別の in_uniparc ヒット率")
    w()
    w("| identity | 遺伝子数 | in_uniparc = yes | ヒット率 |")
    w("|---|---|---|---|")
    for low, high, label in BINS:
        subset = [r for r in rows if low <= r["identity"] < high]
        y = sum(1 for r in subset if r["hit"])
        rate = f"{100 * y / len(subset):.2f} %" if subset else "n/a"
        w(f"| {label} | {len(subset)} | {y} | {rate} |")
    w()

    # -- 4. QC combination ---------------------------------------------------
    w("## 4. QC フラグ（記載どおりの組み合わせ）別の in_uniparc ヒット率")
    w()
    groups = collections.defaultdict(list)
    for r in rows:
        groups[r["qc"]].append(r)
    w("| QC | 遺伝子数 | in_uniparc = yes | ヒット率 | identity 中央値 |")
    w("|---|---|---|---|---|")
    for qc, subset in sorted(groups.items(), key=lambda kv: -len(kv[1])):
        y = sum(1 for r in subset if r["hit"])
        med = stats([r["identity"] for r in subset])["med"]
        w(f"| `{qc}` | {len(subset)} | {y} | {100 * y / len(subset):.2f} % | "
          f"{med:.4f} |")
    w(f"| **計** | **{total}** | **{ty}** | **{100 * ty / total:.2f} %** | |")
    w()

    # -- 5. individual QC flags ----------------------------------------------
    w("## 5. QC フラグ（個別、組み合わせを分解）別の in_uniparc ヒット率")
    w()
    w("1 遺伝子が複数フラグを持つため、行の合計は遺伝子数と一致しない。")
    w()
    flags = sorted({f for r in rows for f in r["qc"].split(",") if f})
    w("| フラグ | 遺伝子数 | in_uniparc = yes | ヒット率 |")
    w("|---|---|---|---|")
    for flag in flags:
        subset = [r for r in rows if flag in r["qc"].split(",")]
        y = sum(1 for r in subset if r["hit"])
        w(f"| `{flag}` | {len(subset)} | {y} | {100 * y / len(subset):.2f} % |")
    w()

    # -- 6. QC x identity=100% ------------------------------------------------
    w("## 6. QC フラグ別の「identity = 100% かつ in_uniparc = no」")
    w()
    w("| QC | identity = 100% | うち in_uniparc = no | 率 |")
    w("|---|---|---|---|")
    for qc, subset in sorted(groups.items(), key=lambda kv: -len(kv[1])):
        e = [r for r in subset if r["identity"] >= 1.0]
        n = sum(1 for r in e if not r["hit"])
        rate = f"{100 * n / len(e):.2f} %" if e else "n/a"
        w(f"| `{qc}` | {len(e)} | {n} | {rate} |")
    w(f"| **計** | **{len(exact)}** | **{cells[('=100%', 'no')]}** | "
      f"**{100 * cells[('=100%', 'no')] / len(exact):.2f} %** |")
    w()
    fh.close()

    print(f"wrote {args.out}")
    print(f"identity=100% & in_uniparc=no : {cells[('=100%', 'no')]}")
    print(f"identity<100% & in_uniparc=yes: {cells[('<100%', 'yes')]}")


if __name__ == "__main__":
    main()
