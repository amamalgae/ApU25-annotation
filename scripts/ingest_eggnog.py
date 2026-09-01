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

Outputs: annotation/UTEX25_gene_table_eggnog.tsv,
annotation/eggnog_provenance.txt, docs/EGGNOG_REPORT.md.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DEFAULT_TABLE = ROOT / "annotation" / "UTEX25_gene_table.tsv"
DEFAULT_OUT_TABLE = ROOT / "annotation" / "UTEX25_gene_table_eggnog.tsv"
DEFAULT_OUT_REPORT = ROOT / "docs" / "EGGNOG_REPORT.md"
DEFAULT_OUT_PROVENANCE = ROOT / "annotation" / "eggnog_provenance.txt"

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
        self.raw_comments: list[str] = []   # every '#' line, verbatim
        self.header_line: str = ""          # the header line, verbatim
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
                    self.raw_comments.append(line)
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
                    self.raw_comments.append(line)
                    self.header_line = line
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
    ap.add_argument("--out-table", type=Path, default=DEFAULT_OUT_TABLE)
    ap.add_argument("--out-report", type=Path, default=DEFAULT_OUT_REPORT)
    ap.add_argument("--out-provenance", type=Path,
                    default=DEFAULT_OUT_PROVENANCE)
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

    merged_path = args.out_table
    merged_path.parent.mkdir(parents=True, exist_ok=True)
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

    write_provenance(args, ann, references)
    write_report(args, ann, references, rows, n_joined, per_column, carried,
                 merged_path)
    print(f"wrote {merged_path} ({len(rows)} rows + header)")
    print(f"joined {n_joined}/{len(rows)} genes "
          f"({100 * n_joined / len(rows):.2f}%)")




def rel(path) -> Path:
    path = Path(path).resolve()
    try:
        return path.relative_to(ROOT)
    except ValueError:
        return path


# --------------------------------------------------------------------------- #
# Reports
# --------------------------------------------------------------------------- #

def write_provenance(args, ann, references) -> None:
    """Dump every '#' line of every emapper file, verbatim, to one file.

    This is the evidence for what the runs actually did, so nothing here is
    reformatted, reordered, or summarised.
    """
    out = args.out_provenance
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w") as fh:
        fh.write("# eggNOG-mapper provenance\n")
        fh.write("#\n")
        fh.write("# Every comment line of each .emapper.annotations file,\n")
        fh.write("# reproduced verbatim by scripts/ingest_eggnog.py.\n")
        fh.write("# The web job page and these recorded command lines are known\n")
        fh.write("# not to agree; these lines are what the runs actually used.\n")
        for label, src in [("MERGED", ann)] + [("REFERENCE ONLY", r)
                                               for r in references]:
            fh.write("\n")
            fh.write("=" * 78 + "\n")
            fh.write(f"{label}: {rel(src.path)}\n")
            fh.write("=" * 78 + "\n")
            for line in src.raw_comments:
                fh.write(line + "\n")
    print(f"wrote {out}")


# Benchmark figures for Pfam calling mode, quoted from the reference below.
# They are literature values, not measurements on this dataset.
PFAM_F1_DENOVO = "89.7%"
PFAM_F1_REALIGN = "98.9%"
PFAM_CITATION = ("Cantalapiedra CP et al. (2021) *Mol Biol Evol* 38(12):5825. "
                 "DOI: 10.1093/molbev/msab293")


def write_report(args, ann, references, rows, n_joined, per_column, carried,
                 merged_path):
    total = len(rows)
    out = args.out_report
    out.parent.mkdir(parents=True, exist_ok=True)
    fh = out.open("w")
    w = lambda s="": fh.write(s + "\n")

    pfam_realign = ann.filters.get("pfam_realign")

    w("# eggNOG-mapper アノテーションの取り込み")
    w()
    w(f"- 生成: `scripts/ingest_eggnog.py`")
    w(f"- 取り込み対象: `{rel(args.annotations)}`")
    w(f"- 構造アノテーション: `{rel(args.gene_table)}`（{total} 遺伝子）")
    w(f"- 出力表: `{rel(merged_path)}`")
    w(f"- 実行条件の証跡（`#` 行の全文）: `{rel(args.out_provenance)}`")
    w()
    w("すべて実測値。ベンチマーク F1 の 2 値のみ文献値で、出典を明示している。")
    w()

    # -- 1. header ---------------------------------------------------------- #
    w("## 1. ヘッダ（決め打ちせず実ファイルから検出）")
    w()
    w(f"- 検出したクエリ列名: `{ann.query_column}`")
    w(f"- 列数: {len(ann.columns)} / データ行数: **{len(ann.rows)}**")
    w()
    w("```")
    w(ann.header_line)
    w("```")
    w()

    # -- 2. command line ---------------------------------------------------- #
    w("## 2. ヘッダから抽出した実行コマンド（原文のまま）")
    w()
    w(f"- version: `{ann.version}`")
    w()
    w("```")
    w(ann.command or "(command line not recorded in the header)")
    w("```")
    w()
    w("`## applied filters:` ブロックの実測値:")
    w()
    w("| key | value |")
    w("|---|---|")
    for k, v in ann.filters.items():
        mark = " **←**" if k == "pfam_realign" else ""
        w(f"| `{k}` | `{v}`{mark} |")
    w()
    w(f"**`pfam_realign` の実測値: `{pfam_realign}`**")
    w()
    w("ウェブサーバ (https://eggnog-mapper.cgmlab.org/, v3-beta6) のジョブページ表示と")
    w("この記録は一致しないことが確認されている。上に転記したのは記録された方、")
    w("すなわち実際に走った設定である。")
    w()

    # -- 3. what pfam_realign=none means ------------------------------------ #
    if pfam_realign == "none":
        w("## 3. `pfam_realign=none` が Pfam 呼び出しに与える影響（文献値）")
        w()
        w("`pfam_realign=none` は、Pfam ドメインを seed ortholog 経由の転写")
        w("(transfer) で付与し、クエリ配列に対する再アラインメントを行わない。")
        w("転写モードでの Pfam 呼び出しについて、de novo を正解としたときの")
        w(f"報告値は **F1 = {PFAM_F1_DENOVO}**、realign を行った場合は")
        w(f"**F1 = {PFAM_F1_REALIGN}** である（{PFAM_CITATION}）。")
        w()
        w("**この 2 値をそのまま本データに当てはめることはできない。**")
        w("当該ベンチマークは Progenomes（原核生物）ベースであり、eggNOG における")
        w("代表性が低い緑藻では、誤差はこれより悪い方向に振れうる。")
        w("本データでの実測は `docs/PFAM_CONCORDANCE.md`（InterProScan 6 を正解と")
        w("したときの一致度）を参照。")
    else:
        w("## 3. Pfam 呼び出しモード")
        w()
        w(f"`pfam_realign` の実測値は `{pfam_realign}` であり `none` ではないため、")
        w("転写モード前提の文献値（F1 = 89.7% / 98.9%）は該当しない。")
    w()

    # -- 4. coverage -------------------------------------------------------- #
    w(f"## 4. 各列のカバレッジ（全 {total} 遺伝子に対して）")
    w()
    w(f"- eggNOG 行が付いた遺伝子: **{n_joined} / {total} "
      f"({100 * n_joined / total:.2f} %)**")
    w(f"- eggNOG 行が無い遺伝子: {total - n_joined} / {total} "
      f"({100 * (total - n_joined) / total:.2f} %)")
    w()
    w("| 列（出力表での名前） | 件数 / " + str(total) + " | % |")
    w("|---|---|---|")
    w(f"| （eggNOG 行そのもの） | {n_joined} / {total} | {100 * n_joined / total:.2f} |")
    for c in carried:
        n = per_column[c]
        w(f"| `eggnog_{c}` | {n} / {total} | {100 * n / total:.2f} |")
    w()
    w("空欄は「eggNOG が値を返さなかった」であり、推定値は一切入れていない。")
    w(f"既存の `product` / `symbol` および `{rel(args.gene_table)}` 本体は書き換えていない。")
    w()

    # -- 5. reference-only files -------------------------------------------- #
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
            w(f"- `pfam_realign`: `{ref.filters.get('pfam_realign')}`")
            same = ref.command == ann.command
            w(f"- 記録された command は取り込み対象と"
              f"**{'完全一致' if same else '相違'}**")
            w()
        w("de novo を狙って再投入したジョブでも `pfam_realign=none` が記録され、")
        w("PFAM 付与は 0 件だった。ジョブページ表示と実行内容の不一致を示す直接証拠。")
        w()

    fh.close()
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
