#!/usr/bin/env python3
"""Compare Pfam domain calls: InterProScan 6 (via the Matches API) vs eggNOG-mapper.

Scope: the proteins for which BOTH tools produced a result, i.e. the sequence is
in UniParc (so the Matches API returned InterProScan 6 output) AND eggNOG-mapper
emitted an annotation row.  A protein where one side called Pfam domains and the
other called none stays in scope — that is exactly a recall/precision event.

Comparison unit: the PFxxxxx accession, version suffix stripped, compared as a
SET per protein (a protein carries several domains).

One thing about the inputs has to be handled explicitly.  eggNOG-mapper 3.0.0-beta6
does NOT write Pfam accessions in its ``PFAMs`` column; it writes Pfam *names*
with the hit coordinates appended, e.g. ``AAA_lid_3_330_372``.  So:

  * the trailing ``_<start>_<end>`` is stripped (the script aborts if a token
    does not have that shape, rather than guessing);
  * the remaining name is mapped to an accession through Pfam's own
    ``Pfam-A.clans.tsv.gz`` (plus ``Pfam-A.dead.gz`` for killed families);
  * names that map to nothing are NOT silently dropped: they are counted,
    listed in the report, and the direction of the bias they cause is stated.

Because that mapping is lossy, the report also carries a second, complete
comparison keyed on the Pfam NAME (InterPro supplies ``signature.name`` for its
Pfam matches), which needs no mapping at all.  The gap between the two is the
size of the mapping loss.

Usage:

    python3 scripts/compare_pfam.py
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
DEFAULT_EGGNOG = ROOT / "raw" / "eggnog_7413.emapper.annotations"
DEFAULT_CLANS = ROOT / "raw" / "Pfam-A.clans.tsv.gz"
DEFAULT_DEAD = ROOT / "raw" / "Pfam-A.dead.gz"
DEFAULT_OUT = ROOT / "docs" / "PFAM_CONCORDANCE.md"

# eggNOG token: <pfam name>_<start>_<end>.  Names themselves contain digits and
# underscores ("AAA_lid_3"), so only the last two numeric fields are coordinates.
EGGNOG_TOKEN = re.compile(r"^(?P<name>.+?)_(?P<start>\d+)_(?P<end>\d+)$")
ACCESSION = re.compile(r"^PF\d{5}$")


def strip_version(accession: str) -> str:
    return accession.split(".", 1)[0]


# --------------------------------------------------------------------------- #
# Pfam name -> accession
# --------------------------------------------------------------------------- #

def load_pfam_map(clans: Path, dead: Path):
    """name -> accession, from Pfam's own release files."""
    name2acc, source = {}, {}
    with gzip.open(clans, "rt") as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                continue
            acc, name = strip_version(fields[0]), fields[3]
            if ACCESSION.match(acc):
                name2acc.setdefault(name, acc)
                source.setdefault(name, "Pfam-A.clans.tsv")

    # Killed families keep their accession; eggNOG's database may still use the
    # old name.  These are added only where the current release has no such name.
    if dead.exists():
        entry = {}
        with gzip.open(dead, "rt") as fh:
            for line in fh:
                line = line.rstrip("\n")
                if line == "//":
                    name, acc = entry.get("ID"), strip_version(entry.get("AC", ""))
                    if name and ACCESSION.match(acc) and name not in name2acc:
                        name2acc[name] = acc
                        source[name] = "Pfam-A.dead"
                    entry = {}
                elif line.startswith("#=GF "):
                    parts = line[5:].split(None, 1)
                    if len(parts) == 2:
                        entry.setdefault(parts[0], parts[1].strip())
    return name2acc, source


# --------------------------------------------------------------------------- #
# InterPro side
# --------------------------------------------------------------------------- #

def load_interpro(matches: Path, lookup: Path):
    """protein_id -> (accession set, name set), for in-UniParc proteins only."""
    with gzip.open(matches, "rt", encoding="utf-8") as fh:
        raw = json.load(fh)

    by_md5 = {}
    for batch in raw["batches"]:
        if batch.get("status") != 200:
            raise SystemExit(
                f"{matches}: batch {batch.get('batch')} has status "
                f"{batch.get('status')}; re-run `05_interpro_matches.py fetch` "
                "before comparing against a partial result")
        for result in batch["body"]["results"]:
            by_md5[result["md5"].upper()] = result

    accs, names, in_uniparc = {}, {}, set()
    with open(lookup) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        i_pid, i_md5 = header.index("protein_id"), header.index("md5")
        i_up = header.index("in_uniparc")
        for line in fh:
            row = line.rstrip("\n").split("\t")
            if row[i_up] != "yes":
                continue
            pid = row[i_pid]
            in_uniparc.add(pid)
            result = by_md5[row[i_md5].upper()]
            a, n = set(), set()
            for match in result["matches"]:
                signature = match["signature"]
                library = signature["signatureLibraryRelease"]["library"]
                if library != "Pfam":
                    continue
                accession = strip_version(signature["accession"])
                if not ACCESSION.match(accession):
                    raise SystemExit(
                        f"{matches}: Pfam signature accession {accession!r} is "
                        "not PFxxxxx; refusing to normalise it by guesswork")
                a.add(accession)
                if signature.get("name"):
                    n.add(signature["name"])
            accs[pid], names[pid] = a, n
    return accs, names, in_uniparc


# --------------------------------------------------------------------------- #
# eggNOG side
# --------------------------------------------------------------------------- #

def load_eggnog(path: Path, name2acc):
    """protein_id -> (accession set, name set, unmapped name set) for rows present."""
    columns = None
    accs, names, unmapped = {}, {}, {}
    unmapped_counts = collections.Counter()
    stats = {"tokens": 0, "accession_shaped": 0}
    for lineno, line in enumerate(open(path), 1):
        line = line.rstrip("\n")
        if line.startswith("##") or not line.strip():
            continue
        if line.startswith("#"):
            columns = [f.lstrip("#").strip() for f in line.split("\t")]
            if columns[0].lower() not in {"query", "query_name"}:
                raise SystemExit(f"{path}:{lineno}: unexpected query column "
                                 f"{columns[0]!r}")
            continue
        if columns is None:
            raise SystemExit(f"{path}:{lineno}: data before the header line")
        row = dict(zip(columns, line.split("\t")))
        pid = row[columns[0]]
        value = row.get("PFAMs", "-")
        a, n, u = set(), set(), set()
        if value not in ("-", ""):
            for token in value.split(","):
                token = token.strip()
                if not token:
                    continue
                m = EGGNOG_TOKEN.match(token)
                if not m:
                    raise SystemExit(
                        f"{path}:{lineno}: PFAMs token {token!r} is not "
                        "<name>_<start>_<end>; refusing to guess how to parse it")
                stats["tokens"] += 1
                if ACCESSION.match(token.split(".", 1)[0]):
                    stats["accession_shaped"] += 1
                name = m.group("name")
                n.add(name)
                accession = name2acc.get(name)
                if accession is None:
                    u.add(name)
                    unmapped_counts[name] += 1
                else:
                    a.add(accession)
        accs[pid], names[pid], unmapped[pid] = a, n, u
    return accs, names, unmapped, unmapped_counts, stats


# --------------------------------------------------------------------------- #
# Comparison
# --------------------------------------------------------------------------- #

class Concordance:
    """Domain-call level confusion over a fixed set of proteins."""

    def __init__(self, proteins, truth, test):
        self.both = self.truth_only = self.test_only = 0
        self.exact = 0
        self.truth_only_counts = collections.Counter()
        self.test_only_counts = collections.Counter()
        self.both_counts = collections.Counter()
        self.n_proteins = len(proteins)
        self.n_both_empty = 0
        for pid in proteins:
            t, e = truth.get(pid, set()), test.get(pid, set())
            self.both += len(t & e)
            self.truth_only += len(t - e)
            self.test_only += len(e - t)
            self.both_counts.update(t & e)
            self.truth_only_counts.update(t - e)
            self.test_only_counts.update(e - t)
            if t == e:
                self.exact += 1
                if not t:
                    self.n_both_empty += 1

    @property
    def precision(self):
        d = self.both + self.test_only
        return self.both / d if d else None

    @property
    def recall(self):
        d = self.both + self.truth_only
        return self.both / d if d else None

    @property
    def f1(self):
        p, r = self.precision, self.recall
        return 2 * p * r / (p + r) if p and r else None

    @property
    def exact_rate(self):
        return self.exact / self.n_proteins if self.n_proteins else None


def pct(x):
    return "n/a" if x is None else f"{100 * x:.2f} %"


# --------------------------------------------------------------------------- #
# Report
# --------------------------------------------------------------------------- #

def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--matches", type=Path, default=DEFAULT_MATCHES)
    ap.add_argument("--lookup", type=Path, default=DEFAULT_LOOKUP)
    ap.add_argument("--eggnog", type=Path, default=DEFAULT_EGGNOG)
    ap.add_argument("--clans", type=Path, default=DEFAULT_CLANS)
    ap.add_argument("--dead", type=Path, default=DEFAULT_DEAD)
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--top", type=int, default=20)
    args = ap.parse_args()

    name2acc, map_source = load_pfam_map(args.clans, args.dead)
    ip_acc, ip_name, in_uniparc = load_interpro(args.matches, args.lookup)
    eg_acc, eg_name, eg_unmapped, unmapped_counts, eg_stats = load_eggnog(
        args.eggnog, name2acc)

    shared = sorted(in_uniparc & set(eg_acc))
    acc_cmp = Concordance(shared, ip_acc, eg_acc)
    name_cmp = Concordance(shared, ip_name, eg_name)

    # accession -> a readable name, for the top-N tables
    acc2name = {}
    for name, acc in name2acc.items():
        acc2name.setdefault(acc, name)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    fh = args.out.open("w")
    w = lambda s="": fh.write(s + "\n")

    w("# Pfam 呼び出しの突き合わせ — InterProScan 6 vs eggNOG-mapper")
    w()
    w("- 生成: `scripts/compare_pfam.py`")
    w(f"- InterPro 側: `{args.matches.relative_to(ROOT)}`"
      f"（InterPro Matches API の生レスポンス、`signatureLibraryRelease.library == \"Pfam\"` のみ）")
    w(f"- eggNOG 側: `{args.eggnog.relative_to(ROOT)}` の `PFAMs` 列")
    w(f"- Pfam 名 → アクセッション対応: `{args.clans.relative_to(ROOT)}`"
      f" + `{args.dead.relative_to(ROOT)}`")
    w()
    w("数値はすべて実測。文献値は最後の 1 段落のみで、出典を明示している。")
    w()

    # -- 0. the input quirk ------------------------------------------------- #
    w("## 0. 比較キーの正規化（実測にもとづく前提）")
    w()
    w("eggNOG-mapper 3.0.0-beta6 の `PFAMs` 列は **PF アクセッションを出力しない**。")
    w("出力されるのは Pfam 名 + ヒット座標（例 `AAA_lid_3_330_372`）で、")
    w(f"本データの Pfam トークン {eg_stats['tokens']:,} 個のうち PFxxxxx 形の"
      f"ものは **{eg_stats['accession_shaped']} 個**だった。")
    w("そのため以下の正規化を行っている。")
    w()
    w("1. 末尾の `_<start>_<end>` を除去して Pfam 名を得る（形が違えば異常終了）")
    w("2. Pfam 公式の名前→アクセッション表で PFxxxxx に変換")
    w("3. InterPro 側は `signature.accession` からバージョン接尾辞を除去")
    w("4. タンパク質ごとに**集合**として比較")
    w()
    total_tokens = sum(unmapped_counts.values())
    n_unmapped_names = len(unmapped_counts)
    w(f"変換できなかった Pfam 名: **{n_unmapped_names} 種類 / 延べ {total_tokens} 出現**。")
    w("これは eggNOG の参照 Pfam 版と現行 Pfam 版で名前が変わった家系で、")
    w("アクセッション単位の比較からは除外している。除外は eggNOG 側の集合のみを")
    w("小さくするので、**「両方」を減らし「InterPro のみ」を増やす向き**に働く")
    w("（eggNOG の recall を過小評価しうる）。影響の大きさを見るため、")
    w("変換を必要としない **Pfam 名単位**の比較を §5 に併記した。")
    w()
    if unmapped_counts:
        w("変換できなかった名前（上位 20、延べ出現数）:")
        w()
        w("| Pfam 名 | 出現 |")
        w("|---|---|")
        for name, n in unmapped_counts.most_common(20):
            w(f"| `{name}` | {n} |")
        w()

    # -- 1. scope ----------------------------------------------------------- #
    w("## 1. 対象タンパク質")
    w()
    w("| 集合 | 件数 |")
    w("|---|---|")
    w(f"| InterPro 側に結果がある（UniParc ヒット） | {len(in_uniparc)} |")
    w(f"| eggNOG 側に行がある | {len(eg_acc)} |")
    w(f"| **両方に結果がある = 比較対象** | **{len(shared)}** |")
    w(f"| InterPro のみ結果あり | {len(in_uniparc - set(eg_acc))} |")
    w(f"| eggNOG のみ結果あり | {len(set(eg_acc) - in_uniparc)} |")
    w()
    n_ip_with = sum(1 for p in shared if ip_acc[p])
    n_eg_with = sum(1 for p in shared if eg_acc[p])
    w(f"比較対象 {len(shared)} のうち、Pfam を 1 つ以上持つのは "
      f"InterPro 側 {n_ip_with}、eggNOG 側 {n_eg_with}。")
    w(f"両側とも 0 個だったタンパク質は {acc_cmp.n_both_empty} で、")
    w("これは「完全一致」に含めて数えている。")
    w()

    # -- 2. confusion ------------------------------------------------------- #
    w("## 2. ドメイン呼び出し単位の混同行列（PF アクセッション）")
    w()
    w("（タンパク質 × PF アクセッションの組を 1 呼び出しと数える）")
    w()
    w("| 区分 | 呼び出し数 |")
    w("|---|---|")
    w(f"| 両方にある | {acc_cmp.both} |")
    w(f"| InterPro のみ（eggNOG の取りこぼし） | {acc_cmp.truth_only} |")
    w(f"| eggNOG のみ | {acc_cmp.test_only} |")
    w(f"| **合計（和集合）** | **{acc_cmp.both + acc_cmp.truth_only + acc_cmp.test_only}** |")
    w()
    w(f"- InterPro 側の呼び出し総数: {acc_cmp.both + acc_cmp.truth_only}")
    w(f"- eggNOG 側の呼び出し総数: {acc_cmp.both + acc_cmp.test_only}")
    w()

    # -- 3. P/R/F1 ---------------------------------------------------------- #
    w("## 3. InterProScan 6 を正解としたときの eggNOG の性能")
    w()
    w("| 指標 | 値 |")
    w("|---|---|")
    w(f"| precision | **{pct(acc_cmp.precision)}** |")
    w(f"| recall | **{pct(acc_cmp.recall)}** |")
    w(f"| F1 | **{pct(acc_cmp.f1)}** |")
    w()
    w("precision = 両方 / (両方 + eggNOG のみ)、")
    w("recall = 両方 / (両方 + InterPro のみ)。")
    w()

    # -- 4. protein level --------------------------------------------------- #
    w("## 4. タンパク質単位の一致率")
    w()
    w(f"- ドメイン集合が完全一致: **{acc_cmp.exact} / {acc_cmp.n_proteins} "
      f"({pct(acc_cmp.exact_rate)})**")
    w(f"  - うち両側とも 0 ドメイン: {acc_cmp.n_both_empty}")
    nonempty = [p for p in shared if ip_acc[p] or eg_acc[p]]
    exact_nonempty = sum(1 for p in nonempty if ip_acc[p] == eg_acc[p])
    w(f"- 少なくとも一方が 1 ドメイン以上ある {len(nonempty)} に限った完全一致: "
      f"**{exact_nonempty} ({pct(exact_nonempty / len(nonempty)) if nonempty else 'n/a'})**")
    w()

    # -- 5. name-level sensitivity ------------------------------------------ #
    w("## 5. 感度確認 — Pfam 名単位（対応表を経由しない比較）")
    w()
    w("§0 の変換損失が結論を動かしていないかの確認。両側とも Pfam 名を持つので")
    w("対応表が不要で、eggNOG のトークンを 1 つも落とさずに比較できる。")
    w()
    w("| 指標 | アクセッション単位 (§3) | 名前単位 |")
    w("|---|---|---|")
    w(f"| 両方にある | {acc_cmp.both} | {name_cmp.both} |")
    w(f"| InterPro のみ | {acc_cmp.truth_only} | {name_cmp.truth_only} |")
    w(f"| eggNOG のみ | {acc_cmp.test_only} | {name_cmp.test_only} |")
    w(f"| precision | {pct(acc_cmp.precision)} | {pct(name_cmp.precision)} |")
    w(f"| recall | {pct(acc_cmp.recall)} | {pct(name_cmp.recall)} |")
    w(f"| F1 | {pct(acc_cmp.f1)} | {pct(name_cmp.f1)} |")
    w(f"| タンパク質単位の完全一致 | {pct(acc_cmp.exact_rate)} | {pct(name_cmp.exact_rate)} |")
    w()
    if (acc_cmp.both == name_cmp.both
            and acc_cmp.truth_only == name_cmp.truth_only):
        w("「両方にある」と「InterPro のみ」は 2 つの単位で**同数**だった。")
        w("すなわち §0 で変換できなかった eggNOG 名は 1 つも InterPro 側と")
        w("一致していない（InterPro も現行 Pfam 名を使うため当然ではある）。")
        w("したがって変換損失は recall を動かしておらず、影響は eggNOG のみ側の")
        w(f"{name_cmp.test_only - acc_cmp.test_only} 件、すなわち precision "
          f"{pct(acc_cmp.precision)} → {pct(name_cmp.precision)} の差に限られる。")
    else:
        w("2 つの単位で「両方にある」または「InterPro のみ」の件数が異なる。")
        w("変換損失が一致判定そのものに影響しているため、両方の数値を併記する。")
    w()

    # -- 6. top lists ------------------------------------------------------- #
    def top_table(counter, title):
        w(f"### {title}")
        w()
        w("| # | PF アクセッション | Pfam 名 | 件数 |")
        w("|---|---|---|---|")
        for i, (acc, n) in enumerate(counter.most_common(args.top), 1):
            w(f"| {i} | `{acc}` | `{acc2name.get(acc, '?')}` | {n} |")
        w()
        w(f"（異なるアクセッション {len(counter)} 種類、延べ {sum(counter.values())} 件）")
        w()

    w(f"## 6. 上位 {args.top} ドメイン")
    w()
    top_table(acc_cmp.truth_only_counts,
              f"6.1 eggNOG が取りこぼした（InterPro のみ）上位 {args.top}")
    top_table(acc_cmp.test_only_counts,
              f"6.2 eggNOG のみにある上位 {args.top}")

    # -- 7. vs literature --------------------------------------------------- #
    w("## 7. 文献値 F1 = 89.7% との比較")
    w()
    w(f"本データでの実測は F1 = **{pct(acc_cmp.f1)}**（アクセッション単位、"
      f"対象 {len(shared)} タンパク質、precision {pct(acc_cmp.precision)} / "
      f"recall {pct(acc_cmp.recall)}）で、名前単位でも "
      f"**{pct(name_cmp.f1)}** だった。")
    w("Cantalapiedra et al. (2021) *Mol Biol Evol* 38(12):5825 "
      "(DOI: 10.1093/molbev/msab293) が転写モードの Pfam 呼び出しについて報告した "
      "F1 = 89.7%（realign 時 98.9%）と、直接は比較できない。"
      "同ベンチマークの正解は de novo の Pfam 呼び出しであるのに対し、"
      "ここでの正解は InterProScan 6 が UniParc 収録配列に対して算出した結果であり、"
      "対象生物も Progenomes（原核）ではなく緑藻 "
      "*Auxenochlorella protothecoides* である。"
      "eggNOG における緑藻の代表性は低く、転写元となる seed ortholog が"
      "系統的に遠くなるため、誤差は文献値より悪い方向に振れうる。"
      "実測値がその方向に出ているか否かは上の表の数値そのものを参照のこと。"
      "本節は 2 つの数値の由来が違うことの注記であって、"
      "文献値を本データの期待値として採用するものではない。")
    w()
    fh.close()

    print(f"wrote {args.out}")
    print(f"proteins compared : {len(shared)}")
    print(f"both/IPonly/EGonly: {acc_cmp.both} / {acc_cmp.truth_only} / "
          f"{acc_cmp.test_only}")
    print(f"P / R / F1        : {pct(acc_cmp.precision)} / "
          f"{pct(acc_cmp.recall)} / {pct(acc_cmp.f1)}")
    print(f"exact set match   : {acc_cmp.exact}/{acc_cmp.n_proteins} "
          f"({pct(acc_cmp.exact_rate)})")
    print(f"unmapped eggNOG names: {n_unmapped_names} names / "
          f"{total_tokens} occurrences")


if __name__ == "__main__":
    main()
