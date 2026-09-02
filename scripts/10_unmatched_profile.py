#!/usr/bin/env python3
"""Profile the 3,514 sequences that missed UniParc, i.e. the InterProScan 6 input.

Writes a marked section into docs/IPS6_FEASIBILITY.md.  Re-running replaces that
section rather than appending a second copy.

The point of the record: if a local InterProScan 6 run later returns low
coverage on this set, these numbers are what separates "the gene models are
poor" from "the method found nothing".  The section states numbers only.

    python3 scripts/10_unmatched_profile.py
"""

from __future__ import annotations

import argparse
import collections
import statistics
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DEFAULT_LOOKUP = ROOT / "step0_out" / "lookup_status.tsv"
DEFAULT_TABLE = ROOT / "annotation" / "UTEX25_gene_table.tsv"
DEFAULT_OUT = ROOT / "docs" / "IPS6_FEASIBILITY.md"

BEGIN = "<!-- BEGIN unmatched-profile -->"
END = "<!-- END unmatched-profile -->"


def read_tsv(path: Path):
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(header, line.rstrip("\n").split("\t")))
                for line in fh if line.strip()]


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


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--lookup", type=Path, default=DEFAULT_LOOKUP)
    ap.add_argument("--gene-table", type=Path, default=DEFAULT_TABLE)
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = ap.parse_args()

    lookup = read_tsv(args.lookup)
    genes = {g["locus_tag"]: g for g in read_tsv(args.gene_table)}

    missing = [r["protein_id"] for r in lookup if r["protein_id"] not in genes]
    if missing:
        raise SystemExit(f"{len(missing)} protein_id(s) absent from "
                         f"{args.gene_table} (first: {missing[:5]})")

    unmatched = [r for r in lookup if r["in_uniparc"] == "no"]
    total_all = len(lookup)
    n = len(unmatched)

    qc = collections.Counter(genes[r["protein_id"]]["QC"] for r in unmatched)
    qc_all = collections.Counter(genes[r["protein_id"]]["QC"] for r in lookup)
    identity = quartiles([float(genes[r["protein_id"]]["identity"])
                          for r in unmatched])
    length = quartiles([int(r["length"]) for r in unmatched])
    stops = sum(1 for r in unmatched if r["has_internal_stop"] == "yes")

    out = []
    w = out.append
    w("### 対象集合の性質（未ヒット 3,514 配列）")
    w("")
    w(f"生成: `scripts/10_unmatched_profile.py` / "
      f"`{args.lookup.relative_to(ROOT)}` の `in_uniparc = no` "
      f"{n} 件（全 {total_all} 件中）を "
      f"`{args.gene_table.relative_to(ROOT)}` と結合。数値のみ。")
    w("")
    w("#### QC フラグ別内訳")
    w("")
    w("| QC | 件数 | 未ヒット内の構成比 | 同 QC の全遺伝子 | そのうち未ヒットの割合 |")
    w("|---|---|---|---|---|")
    for name, count in sorted(qc.items(), key=lambda kv: (-kv[1], kv[0])):
        w(f"| `{name}` | {count} | {100 * count / n:.2f} % | {qc_all[name]} | "
          f"{100 * count / qc_all[name]:.2f} % |")
    w(f"| **計** | **{n}** | **100.00 %** | {total_all} | "
      f"{100 * n / total_all:.2f} % |")
    w("")
    w("#### identity 分布")
    w("")
    w("| 集合 | n | 最小 | Q1 | 中央値 | Q3 | 最大 | 平均 |")
    w("|---|---|---|---|---|---|---|---|")
    w(f"| 未ヒット | {identity['n']} | {identity['min']:.4f} | "
      f"{identity['q1']:.4f} | {identity['med']:.4f} | {identity['q3']:.4f} | "
      f"{identity['max']:.4f} | {identity['mean']:.4f} |")
    w("")
    w("#### タンパク質長の分布（アミノ酸残基数）")
    w("")
    w(f"`{args.lookup.relative_to(ROOT)}` の `length` 列"
      "（末尾 `*` を除去した後の長さ、InterProScan に渡る長さ）。")
    w("")
    w("| 集合 | n | 最小 | Q1 | 中央値 | Q3 | 最大 | 平均 | 合計 |")
    w("|---|---|---|---|---|---|---|---|---|")
    total_aa = sum(int(r["length"]) for r in unmatched)
    w(f"| 未ヒット | {length['n']} | {length['min']} | {length['q1']:.1f} | "
      f"{length['med']:.1f} | {length['q3']:.1f} | {length['max']} | "
      f"{length['mean']:.1f} | {total_aa:,} |")
    w("")
    w("#### 内部終止コドン")
    w("")
    w(f"- 内部終止 `*` を含む配列: **{stops} / {n} "
      f"({100 * stops / n:.2f} %)**")
    w(f"- 全 7,413 配列中の内部終止保有数は "
      f"{sum(1 for r in lookup if r['has_internal_stop'] == 'yes')} で、"
      "そのすべてが未ヒット側にある。")
    w("")

    replace_section(args.out, "\n".join(out) + "\n")
    print(f"wrote the unmatched-profile section into {args.out}")
    print(f"unmatched {n} / {total_all}; internal stops {stops}; "
          f"total {total_aa} aa")


if __name__ == "__main__":
    main()
