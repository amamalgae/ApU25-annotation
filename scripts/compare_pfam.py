#!/usr/bin/env python3
"""Compare Pfam domain calls: InterProScan 6 (via the Matches API) vs eggNOG-mapper.

Scope: the proteins for which BOTH tools produced a result, i.e. the sequence is
in UniParc (so the Matches API returned InterProScan 6 output) AND eggNOG-mapper
emitted an annotation row.  A protein where one side called Pfam domains and the
other called none stays in scope — that is exactly a recall/precision event.

What this measures is the AGREEMENT between eggNOG-mapper's transferred Pfam
calls and InterProScan 6's calls.  It is not the transfer-vs-de-novo comparison
that eggNOG's own benchmark reports: the two sides here differ in Pfam release,
in clan-conflict resolution, and in pipeline, and those differences are folded
into the disagreement.  The report says so in its first section.

It reports the comparison two ways: over every eggNOG call (keyed on the Pfam
name) and over only those eggNOG calls whose name resolves to an accession in
the Pfam release table (keyed on the accession).  The second is NOT a
release-matched comparison -- eggNOG's own Pfam version is not recorded anywhere
in its output, so no common release can be defined.  It only removes calls from
the eggNOG side.

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
DEFAULT_PFAM_VERSION = ROOT / "raw" / "Pfam.version.gz"
DEFAULT_GENE_TABLE = ROOT / "annotation" / "UTEX25_gene_table.tsv"
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
    libraries = collections.Counter()
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
                release = signature["signatureLibraryRelease"]
                library = release["library"]
                libraries[(library, release["version"])] += 1
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
    return accs, names, in_uniparc, libraries


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




def read_pfam_release(path: Path):
    """Parse Pfam.version(.gz).  Returns the raw text and a {key: value} map."""
    if not path.exists():
        return None, {}
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as fh:
        text = fh.read().strip()
    fields = {}
    for line in text.splitlines():
        if ":" in line:
            key, _, value = line.partition(":")
            fields[key.strip()] = value.strip()
    return text, fields


def load_qc(path: Path, key="locus_tag", column="QC"):
    if not path.exists():
        return None
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        if key not in header or column not in header:
            raise SystemExit(f"{path}: expected {key!r} and {column!r} columns")
        i_key, i_qc = header.index(key), header.index(column)
        return {r[i_key]: r[i_qc]
                for r in (line.rstrip("\n").split("\t") for line in fh
                          if line.strip())}


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
    ap.add_argument("--pfam-version", type=Path, default=DEFAULT_PFAM_VERSION)
    ap.add_argument("--gene-table", type=Path, default=DEFAULT_GENE_TABLE)
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--top", type=int, default=20)
    args = ap.parse_args()

    name2acc, _map_source = load_pfam_map(args.clans, args.dead)
    ip_acc, ip_name, in_uniparc, libraries = load_interpro(args.matches,
                                                           args.lookup)
    eg_acc, eg_name, _eg_unmapped, unmapped_counts, eg_stats = load_eggnog(
        args.eggnog, name2acc)
    version_text, version_fields = read_pfam_release(args.pfam_version)
    qc = load_qc(args.gene_table)

    # Inventory of the Pfam release the mapping table describes.
    clans_accessions = set()
    with gzip.open(args.clans, "rt") as fh:
        for line in fh:
            acc = strip_version(line.split("\t")[0])
            if ACCESSION.match(acc):
                clans_accessions.add(acc)

    shared = sorted(in_uniparc & set(eg_acc))

    # Two views of the same comparison:
    #   resolved -- keyed on PF accession, so eggNOG calls whose name is not in
    #               the Pfam release table drop out of the eggNOG side only
    #   all_calls -- keyed on the Pfam name, so every eggNOG call is kept
    # `resolved` is NOT release-matched: it shrinks one side, nothing more.
    restricted = Concordance(shared, ip_acc, eg_acc)
    unrestricted = Concordance(shared, ip_name, eg_name)

    # What each side loses when the unresolvable names drop out (measured).
    ip_called = {a for pid in shared for a in ip_acc[pid]}
    ip_dropped = ip_called - clans_accessions
    eg_dropped_calls = unrestricted.test_only - restricted.test_only

    acc2name = {}
    for name, acc in name2acc.items():
        acc2name.setdefault(acc, name)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    fh = args.out.open("w")
    w = lambda s="": fh.write(s + "\n")

    interpro_pfam_versions = sorted({v for (lib, v) in libraries if lib == "Pfam"})
    pfam_release = version_fields.get("Pfam release")

    w("# Pfam 呼び出しの突き合わせ — InterProScan 6 vs eggNOG-mapper")
    w()
    w("- 生成: `scripts/compare_pfam.py`")
    w(f"- InterPro 側: `{args.matches.relative_to(ROOT)}`"
      f"（Matches API の生レスポンス、`signatureLibraryRelease.library == \"Pfam\"` のみ）")
    w(f"- eggNOG 側: `{args.eggnog.relative_to(ROOT)}` の `PFAMs` 列")
    w(f"- Pfam 名 → アクセッション対応: `{args.clans.relative_to(ROOT)}`"
      f" + `{args.dead.relative_to(ROOT)}`（版は `{args.pfam_version.relative_to(ROOT)}`）")
    w()
    w("数値はすべて実測。文献値は §9 の 1 段落のみで、出典を明示している。")
    w()

    # -- 0. positioning ------------------------------------------------------ #
    w("## 0. 本解析の位置づけ")
    w()
    w("**本解析が測っているのは、eggNOG-mapper が転写 (transfer) した Pfam 呼び出しと、")
    w("InterProScan 6 の呼び出しとの一致率である。**")
    w("転写 (transfer) と de novo の直接比較ではない。")
    w()
    w("Cantalapiedra et al. (2021) の F1 = 89.7% は eggNOG **内部**で、")
    w("同一の Pfam 版・同一パイプラインのもとで transfer と de novo を比べた値であり、")
    w("本解析とは比較対象がそもそも異なる。本解析の不一致には、少なくとも次の 3 つが")
    w("交絡として含まれている。")
    w()
    w("1. **Pfam リリース差** — 両側が同じ Pfam 版を使っている保証がない（§2 に実測）")
    w("2. **clan 競合の解決方法の実装差** — 同一 clan 内の重複ドメインをどう1つに")
    w("   絞るかは実装依存で、eggNOG の転写と InterProScan では処理が異なる")
    w("3. **パイプライン差** — 閾値、領域の切り方、転写元 seed ortholog の選び方など")
    w()
    w("**したがって Cantalapiedra 2021 の 89.7% と本値を直接比較することはできず、")
    w("本値は「一致率の下限」とみなすべきである。**")
    w()
    w("§4–§6 では集計を 2 通り並記しているが、**どちらも上記 3 つの交絡を")
    w("除去していない**。とくに 1（Pfam リリース差）については、eggNOG 側の")
    w("Pfam 版が特定不能である以上、**両側の共通集合を定義すること自体ができない**。")
    w("2 通りの違いは、eggNOG 側の呼び出しを一部除外するかどうかだけである（§4）。")
    w()

    # -- 1. key normalisation ------------------------------------------------ #
    w("## 1. 比較キーの正規化（実測にもとづく前提）")
    w()
    w("eggNOG-mapper 3.0.0-beta6 の `PFAMs` 列は **PF アクセッションを出力しない**。")
    w(f"本データの Pfam トークン {eg_stats['tokens']:,} 個のうち PFxxxxx 形のものは "
      f"**{eg_stats['accession_shaped']} 個**で、実際は Pfam 名に開始・終了座標を")
    w("アンダースコアで連結した形（例 `AAA_lid_3_330_372`）だった。そのため:")
    w()
    w("1. 末尾の数値 2 フィールド（座標）を除去して Pfam 名を得る（形が違えば異常終了）")
    w("2. Pfam 公式の名前→アクセッション表で PFxxxxx に変換")
    w("3. InterPro 側は `signature.accession` からバージョン接尾辞を除去")
    w("4. タンパク質ごとに**集合**として比較")
    w()
    n_unmapped_names = len(unmapped_counts)
    total_unmapped = sum(unmapped_counts.values())
    w(f"変換できなかった Pfam 名: **{n_unmapped_names} 種類 / 延べ {total_unmapped} 出現**"
      f"（全 {eg_stats['tokens']:,} トークンの "
      f"{100 * total_unmapped / eg_stats['tokens']:.2f} %）。")
    w("これらは現行 Pfam の名前表に存在しない名前で、eggNOG の参照 Pfam 版が")
    w("現行版と異なることの直接の証拠でもある（§2）。")
    w("なお「表に無い」は「その家系が現行版に存在しない」と同義ではない。")
    w(f"家系が残ったまま改名された場合も表引きは失敗するため、この {n_unmapped_names} 種類には")
    w("**廃止された家系と改名された家系の両方が混在している**。")
    w()
    if unmapped_counts:
        w("変換できなかった名前（上位 20、延べ出現数）:")
        w()
        w("| Pfam 名 | 出現 |")
        w("|---|---|")
        for name, n in sorted(unmapped_counts.items(),
                              key=lambda kv: (-kv[1], kv[0]))[:20]:
            w(f"| `{name}` | {n} |")
        w()

    # -- 2. releases --------------------------------------------------------- #
    w("## 2. 両側のリリース版（実測）")
    w()
    w("### 2.1 InterPro 側 — Matches API のレスポンスから集計")
    w()
    w(f"`{args.matches.relative_to(ROOT)}` の全マッチの "
      "`signature.signatureLibraryRelease` を集計した。")
    w()
    w("| ライブラリ | version | マッチ数 |")
    w("|---|---|---|")
    for (lib, ver), n in sorted(libraries.items()):
        mark = " **←**" if lib == "Pfam" else ""
        w(f"| {lib}{mark} | `{ver}` | {n} |")
    w()
    if len(interpro_pfam_versions) == 1:
        w(f"**Pfam の release 値は `{interpro_pfam_versions[0]}` の 1 種類のみ**"
          "（レスポンス中に他の版は現れない）。")
    else:
        w(f"**Pfam の release 値が複数ある: {interpro_pfam_versions}**。")
    w()
    w("### 2.2 eggNOG 側")
    w()
    w("emapper の出力ヘッダ（`annotation/eggnog_provenance.txt`）には Pfam の版が")
    w("記録されていない。実行コマンドにも `## applied filters:` にも該当項目が無く、")
    w("**eggNOG-mapper が参照した Pfam 版はこの成果物からは特定不能**である。")
    w()
    w("特定できるのは次の 2 点のみ:")
    w()
    w(f"- 名前→アクセッション対応表として使った Pfam のリリース版")
    if version_text:
        w()
        w("```")
        w(version_text)
        w("```")
        w()
    else:
        w(f"  — `{args.pfam_version.relative_to(ROOT)}` が無いため特定不能")
    w(f"- eggNOG が出力した Pfam 名のうち {n_unmapped_names} 種類がこの版の名前表に")
    w("  無い → eggNOG の参照 Pfam 版は上記の版と**同一ではない**")
    w()
    if pfam_release and interpro_pfam_versions == [pfam_release]:
        w(f"対応表の Pfam 版 `{pfam_release}` は、InterProScan 6 が報告した Pfam 版 "
          f"`{interpro_pfam_versions[0]}` と**一致する**。")
        w(f"ただし一致するのは InterPro 側と対応表の間だけである。"
          "eggNOG 側の版は不明のままなので、**「両側に共通する Pfam 版」は"
          "定義できていない**。§4 以降で行っているのは、"
          f"Pfam {pfam_release} の名前表で解決できない eggNOG 側の呼び出しを"
          "除外することだけで、バージョン差の除去ではない。")
    elif pfam_release:
        w(f"対応表の Pfam 版 `{pfam_release}` と InterProScan の Pfam 版 "
          f"{interpro_pfam_versions} は**一致しない**。")
    w()

    # -- 3. scope ------------------------------------------------------------ #
    w("## 3. 対象タンパク質")
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
    w(f"両側とも 0 個だったタンパク質は {restricted.n_both_empty} で、"
      "これは「完全一致」に含めて数えている。")
    w()

    # -- 4. confusion, both variants ----------------------------------------- #
    w("## 4. ドメイン呼び出し単位の混同行列")
    w()
    w("（タンパク質 × ドメインの組を 1 呼び出しと数える）")
    w()
    w("集計は 2 通り。**どちらもバージョン差を除去していない**（§0、§2.2）。")
    w()
    w("- **全呼び出し**: キーは Pfam 名。eggNOG の呼び出しを 1 つも落とさない。")
    w(f"- **対応表で解決できた呼び出しのみ**: キーは PF アクセッション。"
      f"Pfam {pfam_release or '(版不明)'} の名前表で解決できない eggNOG 側の"
      "呼び出しを除外する。**除外されるのは eggNOG 側だけ**である。")
    w()
    w("| 区分 | 全呼び出し | 対応表で解決できた呼び出しのみ |")
    w("|---|---|---|")
    w(f"| 両方にある | {unrestricted.both} | {restricted.both} |")
    w(f"| InterPro のみ（eggNOG の取りこぼし） | {unrestricted.truth_only} | "
      f"{restricted.truth_only} |")
    w(f"| eggNOG のみ | {unrestricted.test_only} | {restricted.test_only} |")
    w(f"| 合計（和集合） | "
      f"{unrestricted.both + unrestricted.truth_only + unrestricted.test_only} | "
      f"{restricted.both + restricted.truth_only + restricted.test_only} |")
    w()
    w("除外によって各側が失う呼び出し（実測）:")
    w()
    w(f"- InterPro 側: **{len(ip_dropped)} 呼び出し**。"
      f"比較対象で呼ばれた {len(ip_called)} 種類の PF アクセッションは"
      f"すべて対応表に存在した。")
    w(f"- eggNOG 側: **{eg_dropped_calls} 呼び出し**"
      f"（§1 の変換不能な {n_unmapped_names} 種類に由来）。")
    w()
    w("**この操作は「両側を同じ Pfam 版に揃える」ものではない。**")
    w("片側（eggNOG）の分母を減らしているだけであり、"
      "InterPro 側は 1 呼び出しも減っていない。")
    w()

    # -- 5. P/R/F1 ------------------------------------------------------------ #
    w("## 5. InterProScan 6 を正解としたときの eggNOG の性能")
    w()
    w("| 指標 | 全呼び出し | 対応表で解決できた呼び出しのみ |")
    w("|---|---|---|")
    w(f"| precision | {pct(unrestricted.precision)} | **{pct(restricted.precision)}** |")
    w(f"| recall | {pct(unrestricted.recall)} | **{pct(restricted.recall)}** |")
    w(f"| F1 | {pct(unrestricted.f1)} | **{pct(restricted.f1)}** |")
    w()
    w("precision = 両方 / (両方 + eggNOG のみ)、"
      "recall = 両方 / (両方 + InterPro のみ)。")
    w()
    if (restricted.both == unrestricted.both
            and restricted.truth_only == unrestricted.truth_only):
        w("「両方にある」と「InterPro のみ」は 2 通りで**同数**である。")
        w("すなわち §1 で変換できなかった eggNOG 名は 1 つも InterPro 側と一致していない。")
        w(f"よって recall は {pct(unrestricted.recall)} のまま動かず、")
        w(f"**F1 が {pct(unrestricted.f1)} → {pct(restricted.f1)} に上がったのは、")
        w(f"precision の分母（eggNOG 側の呼び出し総数）が "
          f"{unrestricted.both + unrestricted.test_only} → "
          f"{restricted.both + restricted.test_only} に減ったためである**。")
        w("バージョン差が除かれた結果ではない。")
        w()

    # -- 6. protein level ----------------------------------------------------- #
    w("## 6. タンパク質単位の一致率")
    w()
    w("| 指標 | 全呼び出し | 対応表で解決できた呼び出しのみ |")
    w("|---|---|---|")
    w(f"| ドメイン集合が完全一致 | {unrestricted.exact} / {unrestricted.n_proteins}"
      f" ({pct(unrestricted.exact_rate)}) | {restricted.exact} / "
      f"{restricted.n_proteins} ({pct(restricted.exact_rate)}) |")
    nonempty_r = [p for p in shared if ip_acc[p] or eg_acc[p]]
    nonempty_u = [p for p in shared if ip_name[p] or eg_name[p]]
    er = sum(1 for p in nonempty_r if ip_acc[p] == eg_acc[p])
    eu = sum(1 for p in nonempty_u if ip_name[p] == eg_name[p])
    w(f"| うち両側とも 0 ドメインを除く | {eu} / {len(nonempty_u)} "
      f"({pct(eu / len(nonempty_u)) if nonempty_u else 'n/a'}) | "
      f"{er} / {len(nonempty_r)} "
      f"({pct(er / len(nonempty_r)) if nonempty_r else 'n/a'}) |")
    w()

    # -- 7. QC strata ---------------------------------------------------------- #
    w("## 7. QC フラグ別の層別集計（対応表で解決できた呼び出しのみ）")
    w()
    if qc is None:
        w(f"`{args.gene_table.relative_to(ROOT)}` が無いため層別できない。")
        w()
    else:
        missing = [p for p in shared if p not in qc]
        if missing:
            raise SystemExit(
                f"{len(missing)} of the compared proteins have no QC value in "
                f"{args.gene_table} (first: {missing[:5]})")
        groups = collections.defaultdict(list)
        for pid in shared:
            groups[qc[pid]].append(pid)
        w(f"比較対象 {len(shared)} タンパク質を `{args.gene_table.relative_to(ROOT)}` の")
        w("`QC` 列で層別した。数値のみ。")
        w()
        all_qc = sorted({v for v in qc.values()})
        empty = [v for v in all_qc if v not in groups]
        if empty:
            w(f"遺伝子表に現れる QC 値は {len(all_qc)} 種類あるが、比較対象に 1 件も"
              f"含まれない値が {len(empty)} 種類ある: "
              + "、".join(f"`{v}`" for v in empty) + "。")
            w("これらは in_uniparc = no のため InterPro 側の結果が無く、"
              "比較対象に入らない（`docs/PROJECTION_QC.md` の QC 別ヒット率を参照）。")
            w()
        w("| QC | タンパク質数 | 両方 | InterProのみ | eggNOGのみ | precision | "
          "recall | F1 | 完全一致 |")
        w("|---|---|---|---|---|---|---|---|---|")
        for name, pids in sorted(groups.items(), key=lambda kv: -len(kv[1])):
            c = Concordance(pids, ip_acc, eg_acc)
            w(f"| `{name}` | {len(pids)} | {c.both} | {c.truth_only} | "
              f"{c.test_only} | {pct(c.precision)} | {pct(c.recall)} | "
              f"{pct(c.f1)} | {c.exact} ({pct(c.exact_rate)}) |")
        w(f"| **全体** | **{len(shared)}** | {restricted.both} | "
          f"{restricted.truth_only} | {restricted.test_only} | "
          f"{pct(restricted.precision)} | {pct(restricted.recall)} | "
          f"{pct(restricted.f1)} | {restricted.exact} "
          f"({pct(restricted.exact_rate)}) |")
        w()
        flags = sorted({f for pid in shared for f in qc[pid].split(",") if f})
        w("個別フラグに分解（1 タンパク質が複数フラグを持つため合計は一致しない）:")
        w()
        w("| フラグ | タンパク質数 | precision | recall | F1 | 完全一致 |")
        w("|---|---|---|---|---|---|")
        for flag in flags:
            pids = [p for p in shared if flag in qc[p].split(",")]
            c = Concordance(pids, ip_acc, eg_acc)
            w(f"| `{flag}` | {len(pids)} | {pct(c.precision)} | {pct(c.recall)} | "
              f"{pct(c.f1)} | {c.exact} ({pct(c.exact_rate)}) |")
        w()

    # -- 8. top lists ---------------------------------------------------------- #
    def top_table(counter, title):
        w(f"### {title}")
        w()
        w("| # | PF アクセッション | Pfam 名 | 件数 |")
        w("|---|---|---|---|")
        ranked = sorted(counter.items(), key=lambda kv: (-kv[1], kv[0]))
        for i, (acc, n) in enumerate(ranked[:args.top], 1):
            w(f"| {i} | `{acc}` | `{acc2name.get(acc, '?')}` | {n} |")
        w()
        w(f"（異なるアクセッション {len(counter)} 種類、延べ {sum(counter.values())} 件）")
        w()

    w(f"## 8. 上位 {args.top} ドメイン（対応表で解決できた呼び出しのみ）")
    w()
    top_table(restricted.truth_only_counts,
              f"8.1 eggNOG が取りこぼした（InterPro のみ）上位 {args.top}")
    top_table(restricted.test_only_counts,
              f"8.2 eggNOG のみにある上位 {args.top}")

    # -- 9. vs literature ------------------------------------------------------ #
    w("## 9. 文献値 F1 = 89.7% との比較")
    w()
    w(f"本データでの実測は、対応表で解決できた呼び出しのみで F1 = **{pct(restricted.f1)}**"
      f"（precision {pct(restricted.precision)} / recall {pct(restricted.recall)}、"
      f"対象 {len(shared)} タンパク質）、"
      f"全呼び出しで F1 = **{pct(unrestricted.f1)}** である。")
    w(f"§4 のとおりこの 2 値の差はバージョン差の除去によるものではないので、"
      f"**{pct(restricted.f1)} も一致率の下限**であり、"
      "以下の理由により文献値との直接比較はできない。")
    w("Cantalapiedra et al. (2021) *Mol Biol Evol* 38(12):5825 "
      "(DOI: 10.1093/molbev/msab293) が転写モードの Pfam 呼び出しについて報告した "
      "F1 = 89.7%（realign 時 98.9%）とは、§0 に述べたとおり比較対象が異なるため、"
      "直接比較できない。同ベンチマークの正解は同一 Pfam 版・同一パイプラインでの "
      "de novo 呼び出しであり、本解析の正解は InterProScan 6 が UniParc 収録配列に対して"
      "算出した結果である。対象生物も Progenomes（原核）ではなく緑藻 "
      "*Auxenochlorella protothecoides* で、eggNOG における代表性は低い。"
      "本節は 2 つの数値の由来が違うことの注記であって、"
      "文献値を本データの期待値として採用するものではない。")
    w()
    fh.close()

    print(f"wrote {args.out}")
    print(f"proteins compared     : {len(shared)}")
    print(f"restricted   P/R/F1   : {pct(restricted.precision)} / "
          f"{pct(restricted.recall)} / {pct(restricted.f1)}")
    print(f"unrestricted P/R/F1   : {pct(unrestricted.precision)} / "
          f"{pct(unrestricted.recall)} / {pct(unrestricted.f1)}")
    print(f"InterPro Pfam release : {interpro_pfam_versions}")
    print(f"mapping table release : {pfam_release}")
    print(f"dropped by restriction: InterPro {len(ip_dropped)} accessions / "
          f"eggNOG {eg_dropped_calls} calls")


if __name__ == "__main__":
    main()
