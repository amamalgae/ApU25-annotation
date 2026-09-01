#!/usr/bin/env python3
"""Ingest eggNOG-mapper annotations into the UTEX 25 gene table.

The structural table (``annotation/UTEX25_gene_table.tsv``) carries no function:
``product`` is "hypothetical protein" for all 7,413 genes.  This script joins the
eggNOG-mapper output onto it on ``locus_tag`` and writes a new table; the input
table is never modified in place.

Two rules the ingest follows, both deliberate:

1. **No hardcoded column layout.**  The header is located as the one line that
   starts with a single ``#`` and whose first field is a query column
   (``query`` / ``query_name`` / ``#query``); every other column is taken from
   that line by name.  ``--require`` names the columns the caller depends on and
   the run aborts if any is missing, rather than silently emitting empty cells.

2. **No guessing.**  Expected row / PFAM counts are checked before anything is
   written (``--expect-rows`` / ``--expect-pfam``); a mismatch aborts.  Genes
   with no eggNOG row get empty cells, never an invented value, and ``product``
   and ``symbol`` from the structural annotation are left exactly as they are.

Provenance note carried into the report: for this dataset the emapper header's
recorded command line does NOT agree with the settings shown on the web server's
job page (https://eggnog-mapper.cgmlab.org/, v3-beta6).  The header is what the
run actually used, so the header is what gets reported.

Usage:

    python3 scripts/ingest_eggnog.py \\
        --annotations raw/eggnog_7413.emapper.annotations \\
        --reference   raw/eggnog_noPfam1916.emapper.annotations \\
        --expect-rows 6588 --expect-pfam 5497
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DEFAULT_TABLE = ROOT / "annotation" / "UTEX25_gene_table.tsv"
OUT_DIR = ROOT / "eggnog_out"

# emapper writes this single character for "no value".
EMPTY = {"", "-"}

# Recognised spellings of the query column, which is what identifies the header
# line among the comment lines.
QUERY_FIELDS = {"query", "#query", "query_name", "#query_name"}


def is_filled(value: str) -> bool:
    return value.strip() not in EMPTY


# --------------------------------------------------------------------------- #
# Parsing
# --------------------------------------------------------------------------- #

class Emapper:
    """One parsed .emapper.annotations file."""

    def __init__(self, path: Path):
        self.path = path
        self.preamble: list[str] = []   # '##' lines above the header
        self.trailer: list[str] = []    # '##' lines below the last row
        self.columns: list[str] = []
        self.rows: list[dict[str, str]] = []
        self._parse()

    def _parse(self):
        header_line = None
        with open(self.path) as fh:
            for lineno, line in enumerate(fh, 1):
                line = line.rstrip("\n")
                if not line.strip():
                    continue
                if line.startswith("##"):
                    (self.trailer if self.columns else self.preamble).append(
                        line.lstrip("#").strip())
                    continue
                if line.startswith("#"):
                    fields = line.split("\t")
                    if fields[0].strip().lower() not in QUERY_FIELDS:
                        raise ValueError(
                            f"{self.path}:{lineno}: single-'#' line whose first "
                            f"field is {fields[0]!r}, not a query column; "
                            "refusing to guess which line is the header")
                    if self.columns:
                        raise ValueError(
                            f"{self.path}:{lineno}: a second header line")
                    self.columns = [f.lstrip("#").strip() for f in fields]
                    header_line = lineno
                    continue
                if not self.columns:
                    raise ValueError(
                        f"{self.path}:{lineno}: data before any header line")
                fields = line.split("\t")
                if len(fields) != len(self.columns):
                    raise ValueError(
                        f"{self.path}:{lineno}: {len(fields)} fields, header "
                        f"(line {header_line}) has {len(self.columns)}")
                self.rows.append(dict(zip(self.columns, fields)))
        if not self.columns:
            raise ValueError(f"{self.path}: no header line found")

    # -- provenance ------------------------------------------------------- #

    @property
    def query_column(self) -> str:
        return self.columns[0]

    @property
    def version(self) -> str | None:
        return next((l for l in self.preamble if l.startswith("emapper-")), None)

    @property
    def command(self) -> str | None:
        return next((l for l in self.preamble if "emapper.py" in l), None)

    @property
    def filters(self) -> dict[str, str]:
        out = {}
        for line in self.preamble:
            if "=" in line and " " not in line.split("=")[0]:
                key, _, value = line.partition("=")
                out[key.strip()] = value.strip()
        return out

    # -- content ---------------------------------------------------------- #

    def require(self, names) -> None:
        missing = [n for n in names if n not in self.columns]
        if missing:
            raise SystemExit(
                f"{self.path}: required column(s) {missing} not in header "
                f"{self.columns}")

    def filled(self, column: str) -> int:
        return sum(1 for r in self.rows if is_filled(r[column]))

    def by_query(self) -> dict[str, dict[str, str]]:
        out = {}
        for row in self.rows:
            key = row[self.query_column]
            if key in out:
                raise ValueError(f"{self.path}: duplicate query {key!r}")
            out[key] = row
        return out


def check_expectations(ann: Emapper, label: str, rows, pfam, pfam_col) -> None:
    """Abort before writing anything if the file is not what was expected."""
    problems = []
    if rows is not None and len(ann.rows) != rows:
        problems.append(f"data rows: expected {rows}, got {len(ann.rows)}")
    if pfam is not None:
        if pfam_col not in ann.columns:
            problems.append(f"no {pfam_col!r} column to check against")
        else:
            got = ann.filled(pfam_col)
            if got != pfam:
                problems.append(
                    f"{pfam_col} filled: expected {pfam}, got {got}")
    if problems:
        raise SystemExit(
            f"{label} does not match the expected values; stopping before the "
            "merge so nothing is written from a file that is not the one "
            "described:\n  " + "\n  ".join(problems))


# --------------------------------------------------------------------------- #
# Gene table
# --------------------------------------------------------------------------- #

def read_gene_table(path: Path):
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        rows = [line.rstrip("\n").split("\t") for line in fh if line.strip()]
    for i, row in enumerate(rows, 2):
        if len(row) != len(header):
            raise ValueError(f"{path}:{i}: {len(row)} fields, header has "
                             f"{len(header)}")
    return header, rows


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #

def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--annotations", type=Path, required=True,
                    help="the .emapper.annotations file to merge")
    ap.add_argument("--reference", type=Path, action="append", default=[],
                    help="parsed and reported on, never merged (repeatable)")
    ap.add_argument("--gene-table", type=Path, default=DEFAULT_TABLE)
    ap.add_argument("--key", default="locus_tag",
                    help="gene-table column joined against the query column")
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    ap.add_argument("--expect-rows", type=int, default=None)
    ap.add_argument("--expect-pfam", type=int, default=None)
    ap.add_argument("--pfam-column", default="PFAMs")
    ap.add_argument("--require", default="PFAMs,Preferred_name,GOs,EC,KEGG_ko,"
                                         "COG_category,eggNOG_OGs",
                    help="comma-separated columns that must exist in the header")
    args = ap.parse_args()

    ann = Emapper(args.annotations)
    ann.require([c for c in args.require.split(",") if c])
    check_expectations(ann, args.annotations.name,
                       args.expect_rows, args.expect_pfam, args.pfam_column)

    references = []
    for path in args.reference:
        ref = Emapper(path)
        references.append(ref)

    header, rows = read_gene_table(args.gene_table)
    if args.key not in header:
        raise SystemExit(f"{args.gene_table}: no {args.key!r} column")
    key_idx = header.index(args.key)

    by_query = ann.by_query()
    unknown = sorted(set(by_query) - {r[key_idx] for r in rows})
    if unknown:
        raise SystemExit(
            f"{args.annotations.name}: {len(unknown)} query id(s) are not in "
            f"{args.gene_table.name} (first: {unknown[:5]}); refusing to merge "
            "a file that does not belong to this gene set")

    # Everything except the query column is carried over, prefixed so the
    # provenance of each column stays obvious in the merged table.
    carried = [c for c in ann.columns if c != ann.query_column]
    out_header = header + [f"eggnog_{c}" for c in carried]

    args.out_dir.mkdir(exist_ok=True)
    merged_path = args.out_dir / "UTEX25_gene_table_eggnog.tsv"
    n_joined = 0
    per_column = {c: 0 for c in carried}
    with open(merged_path, "w") as fh:
        fh.write("\t".join(out_header) + "\n")
        for row in rows:
            hit = by_query.get(row[key_idx])
            if hit is None:
                extra = [""] * len(carried)
            else:
                n_joined += 1
                extra = []
                for c in carried:
                    value = hit[c]
                    extra.append("" if not is_filled(value) else value)
                    if is_filled(value):
                        per_column[c] += 1
            fh.write("\t".join(row + extra) + "\n")

    write_report(args, ann, references, rows, n_joined, per_column, carried,
                 merged_path)
    print(f"wrote {merged_path} ({len(rows)} rows + header)")
    print(f"joined {n_joined}/{len(rows)} genes "
          f"({100 * n_joined / len(rows):.2f}%)")


def write_report(args, ann, references, rows, n_joined, per_column, carried,
                 merged_path):
    total = len(rows)
    out = args.out_dir / "SUMMARY.md"
    L = out.open("w")
    w = lambda s="": L.write(s + "\n")

    w("# eggNOG-mapper 機能アノテーションの取り込み")
    w()
    def rel(path):
        path = Path(path).resolve()
        try:
            return path.relative_to(ROOT)
        except ValueError:
            return path

    w(f"生成: `scripts/ingest_eggnog.py` / 入力 `{rel(args.annotations)}`")
    w(f"構造アノテーション: `{rel(args.gene_table)}`（{total} 遺伝子）")
    w(f"出力: `{rel(merged_path)}`")
    w()
    w("## 1. 取り込んだファイルのヘッダ（実測、決め打ちなし）")
    w()
    w(f"- 検出したヘッダ行のクエリ列名: `{ann.query_column}`")
    w(f"- 列数: {len(ann.columns)}")
    w(f"- 列名: {', '.join('`' + c + '`' for c in ann.columns)}")
    w(f"- データ行数: **{len(ann.rows)}**")
    w()
    w("## 2. 実行時の記録（`##` 行のまま）")
    w()
    w(f"- version: `{ann.version}`")
    w("- command:")
    w()
    w("```")
    w(ann.command or "(なし)")
    w("```")
    w()
    if ann.filters:
        w("- applied filters:")
        w()
        w("| key | value |")
        w("|---|---|")
        for k, v in ann.filters.items():
            w(f"| `{k}` | `{v}` |")
        w()
    w("> ウェブサーバのジョブページ表示と、この `##` 行に記録された実行コマンドは")
    w("> 一致しないことが確認されている。ここに書いてあるのは **実際に走った設定** の方。")
    w()
    w("## 3. 充填率（データ行 " + str(len(ann.rows)) + " 行に対して）")
    w()
    w("| 列 | 値の入った行 | 率 |")
    w("|---|---|---|")
    for c in carried:
        n = ann.filled(c)
        w(f"| `{c}` | {n} | {100 * n / len(ann.rows):.2f} % |")
    w()
    w(f"## 4. 遺伝子表への結合（{total} 遺伝子に対して）")
    w()
    w(f"- eggNOG 行が付いた遺伝子: **{n_joined} / {total} "
      f"({100 * n_joined / total:.2f} %)**")
    w(f"- eggNOG 行が無い遺伝子: {total - n_joined}")
    w()
    w("| 列 | 値の入った遺伝子 | 全遺伝子に対する率 |")
    w("|---|---|---|")
    for c in carried:
        n = per_column[c]
        w(f"| `eggnog_{c}` | {n} | {100 * n / total:.2f} % |")
    w()
    w("結合は左外部結合。eggNOG 行が無い遺伝子の追加列は**空欄**で、推定値は入れていない。")
    w("既存の `product` と `symbol` は一切書き換えていない。")
    w()
    if references:
        w("## 5. 参照のみ（統合していないファイル）")
        w()
        for ref in references:
            w(f"### `{rel(ref.path)}`")
            w()
            w(f"- データ行数: **{len(ref.rows)}**")
            pf = args.pfam_column
            if pf in ref.columns:
                w(f"- `{pf}` が埋まった行: **{ref.filled(pf)}**")
            w(f"- 記録された command: `{ref.command}`")
            same = (ref.command == ann.command)
            w(f"- 取り込み対象ファイルの command と{'一致' if same else '相違'}")
            w()
    L.close()
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
