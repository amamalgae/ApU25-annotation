# UTEX 25 遺伝子アノテーション（相同性投影版）

Auxenochlorella protothecoides UTEX 25（PRJNA1328465, 12染色体, 21,976,416 bp）に対し、
UTEX 250-A 両ハプロタイプのタンパク質配列を miniprot でスプライスアラインし、
遺伝子モデルを投影したもの。作成日 2026-08-31。

**これは実験的に検証されたアノテーションではない。** RNA-seq / IsoSeq による裏づけはなく、
すべて相同性からの投影である。UTEX 250-A に存在しない遺伝子は原理的に検出できない。

---

## 作成手順

| 段階 | 内容 |
|---|---|
| 1 | UTEX 250-A hap1 (PRJNA1195245) / hap2 (PRJNA1195244) の GenBank から CDS 20,011件を抽出 |
| 2 | locus_tag ごとに最長アイソフォームを代表として選択 → hap1 7,508、hap2 7,514 タンパク質 |
| 3 | miniprot 0.18-r281 で UTEX 25 ゲノムにマップ（`--outn=1`、最良アラインメントのみ） |
| 4 | locus_tag 番号で A/B 座位を対応づけ、スコアの高い側を採用 |
| 5 | ゲノム上で 50% 超重複するモデルは高スコア側のみ残す |
| 6 | CDS をゲノムから抽出し翻訳、GenBank / GFF3 / FASTA / TSV に出力 |

## 統計

| 項目 | 値 |
|---|---|
| 遺伝子モデル | **7,413** |
| マップ率（hap1 / hap2） | 7,494/7,508 (99.8%) / 7,500/7,514 (99.8%) |
| 遺伝子シンボル付き | 2,418 (32.6%) |
| QC 完全通過（開始Met あり・内部終止なし・3の倍数） | 6,786 (91.5%) |
| 内部終止コドンあり | 235 |
| 開始 Met なし | 422 |
| 長さが3の倍数でない | 182 |
| 翻訳ラウンドトリップ検証 | 7,413/7,413 (100%) |

QC 6,786 は「開始Met あり・内部終止なし・3の倍数」の3条件を満たす遺伝子モデル数。UTEX25_gene_table.tsv の QC 列が pass の 6,762 とは定義が異なり、差の 24 件は frameshift フラグのみを持つモデルである。

### 由来サブゲノム

| 由来 | 座位数 |
|---|---|
| A/B 同スコア（tie） | 2,846 |
| B 側が高スコア | 2,317 |
| A 側が高スコア | 2,076 |
| B のみ存在 | 93 |
| A のみ存在 | 81 |

A と B のどちらが UTEX 25 に近いかは**系統的な偏りがない**（同一性差の中央値 0.00000、
A優位 1,976 / B優位 2,310 / 同値 3,052、n=7,338ペア）。対応ペアの 99.9% が同一染色体に落ちる。

## ファイル

| ファイル | 内容 |
|---|---|
| `*.gb` ×12 | 染色体ごとの GenBank（UTEX 250-A と同形式。塩基配列 ORIGIN 付き） |
| `UTEX25_annotation.gff3` | GFF3（gene / mRNA / CDS） |
| `UTEX25_proteins.faa` | 投影タンパク質配列 7,413 本 |
| `UTEX25_cds.fna` | CDS 塩基配列 |
| `UTEX25_gene_table.tsv` | 全遺伝子の一覧表（座標・由来座位・同一性・QCフラグ） |

## 各 CDS に付与した情報

- `/locus_tag` : `APU25_NNNNN`（ゲノム座標順、5刻み）
- `/gene` : UTEX 250-A のシンボル（存在する場合のみ）
- `/product` : UTEX 250-A の product
- `/inference` : 由来した UTEX 250-A の locus_tag
- `/note` : 由来座位、由来サブゲノム（A / B / tie / A-only / B-only）、miniprot 同一性、QCフラグ

## 既知の限界

1. UTEX 250-A に存在しない、または 250-A で未モデル化の遺伝子は検出されない
2. product の 67.4% は "hypothetical protein"。**機能記述の欠落は本ファイルでは解消していない**
   （Swiss-Prot / Pfam / eggNOG による機能アノテーションは外部DBが必要で別途実施が要る）
3. UTEX 25 固有の構造変異・遺伝子欠失は、モデルが投影されないことでのみ間接的に示される
4. 選択的スプライシングは扱っていない（1遺伝子1モデル）
5. tRNA / rRNA / ncRNA は含まない

## 妥当性の傍証

- chr3:1,171,807–1,173,904 に `APU25_05450 = CHL27`（同一性 1.0000）。
  pJLM0021/pJLM0022 の 5'HA 末端は 1,171,806、3'HA 始点は 1,173,915 であり、
  投影モデルは**両アームに正確に挟まれる**。
- chr12:300,918–303,165 に `APU25_32315 = PGK1B`（同一性 1.0000）。
  pJLM0022 の PGK1B プロモーターブロックは 300,123–300,920、
  ターミネーターブロックは 303,163–303,415 であり、
  投影モデルは**この2ブロックに正確に挟まれる**。

---

# 機能アノテーション（v1.0.0）

構造アノテーション（上記）とは別工程。symbol と product をここで確定させた。

## 設計思想

**規則より先にこれを読むこと。** 規則は下の3つの用途から導かれている。
規則の文面で判断できないケースが出たら、文面ではなくここに立ち返って考える。

### 用途

1. **DNA修復系・脂質代謝系の遺伝子を検索して見つけられること。**
   このプロジェクトの実際の問いは「形質転換が動かない理由」と
   「FAD2 などの代謝改変ターゲットの同定」である。symbol 検索で遺伝子が
   引けないと調査が止まる。過去に BRCA2 / MSH2 / XRCC2 を
   「シンボル検索で出ないから不在」と誤判定した実例がある。
2. **後から検証・修正できること。** 本アノテーションは投影と自動注釈の産物で、
   実験的裏づけはない。誤りは必ず含まれる。重要なのは誤りをゼロにすることではなく、
   誤りを後から特定して直せる状態にしておくことである。
3. **外部に出せること。** リポジトリは公開されており、将来的に論文または
   規制当局向け資料の根拠になりうる。根拠のない断定を書かないことが、
   詳しく書くことより優先される。

### 原則

| # | 原則 | 由来 |
|---|---|---|
| 1 | **人手の判断 > 自動注釈。上書きしない。** 既存 symbol 2,418件は Craig ら (DOI: 10.1093/plcell/koaf259) が手作業で付けた値であり、自動値で上書きすればより確からしい情報をより不確かな情報で置き換えることになる | 用途2 |
| 2 | **空欄は失敗ではない。誤った断定のほうが高くつく。** `hypothetical protein` は「調べたが分からなかった」を正しく表現している。誤った product はそれを信じた人の数日〜数週間を奪う。埋まっている割合は品質指標ではない | 用途3 |
| 3 | **由来の異なる情報を混ぜない。** `symbol_source` / `product_source` は装飾ではなく、後から検証可能であるための本体である。新しい情報源を足すときは必ず新しい source 値を定義する | 用途2 |
| 4 | **情報源が変わっても粒度が揃っていること。** 同じ品質のはずの遺伝子が、たまたま通った経路の違いで片方は具体名・片方は一般名になる状態は外部に説明できない。一貫性は平均的な詳しさより優先される | 用途3 |
| 5 | **特異的な名前 > 一般的な名前。ただし裏づけの範囲内で。** `SDR family protein` は検索でヒットしても何も分からない。`chlorophyll b reductase` なら調査が前に進む。ただし原則2に劣後する | 用途1 |

## 出典と版

| 項目 | 版 |
|---|---|
| InterProScan | 6.0.1 |
| InterPro | **109.0**（`--interpro 109.0` で固定。`latest` は使っていない） |
| Pfam | 38.2 |
| NCBIFAM | 19.0 |
| PROSITE patterns / profiles | 2026_01 |
| SMART | 9.0 |
| CDD | 3.21 |
| SUPERFAMILY | 1.75 |
| 階層 | InterPro 109.0 `ParentChildTreeFile.txt` |
| eggNOG-mapper | 3.0.0-beta6（symbol の Preferred_name のみ使用） |

**PANTHER と CATH-Gene3D は未導入。** ローカル実行の帯域制約による
（`docs/IPS6_FEASIBILITY.md` §3.3）。この2つが無いことが後述の200件の判断に直結する。

取得経路は2つあり、`annotation/UTEX25_interpro.tsv` の `source` 列で区別できる。
UniParc 収録済みの 3,899件は Matches API（18ライブラリ）、残り 3,514件はローカル実行
（上表の7ライブラリ）。両者が InterPro 109.0 で揃っていることは
`scripts/08_ips6_versions.py` で検証済み（18/18 一致、不一致0）。

## symbol の決定規則

上から順に適用し、最初に該当したもので確定する。

1. 既存値（UTEX 250-A 由来）があれば**保持。上書き禁止**（原則1）→ `utex250a`
2. 空欄の行に限り eggNOG の `Preferred_name` を入れる → `eggnog`
3. なお空欄なら空欄のまま。推測で埋めない（原則2）→ `none`

## product の決定規則

上から順に適用し、最初に該当したもので確定する。

1. InterPro の `type=Family` があればその `name`。複数ある場合は
   `ParentChildTreeFile` で最も下位（子）を採る（原則5）→ `interpro_family`
   - **1x**: Family 同士に祖先・子孫関係がなく順位が付かない場合は
     採用せず `hypothetical protein` のまま（後述）→ `none`
2. Family が無く Domain のみの場合。**まず IPR アクセッションで重複排除する**
   （同一アクセッションの複数マッチは繰り返しドメインであり1種類と数える）。
   配列上の N末端側→C末端側の順（各ドメインの最も N 側のマッチ開始残基で判定、
   同一開始位置なら IPR アクセッション昇順）に並べ:
   - ユニーク1種: `<A> domain-containing protein` → `interpro_domain`
   - ユニーク2種: `<A> and <B> domain-containing protein` → `interpro_domain_multi`
   - ユニーク3種以上: N→C 順の先頭2つで同じ書式。3つ目以降は product に入れない
     → `interpro_domain_multi`
3. InterPro が無く eggNOG の Description があればそれ。
   **本データでは適用不能**（後述）。
4. **3.5**: InterPro エントリが `Homologous_superfamily` / `Repeat` /
   `Conserved_site` のみで構成される場合は product に使用せず
   `hypothetical protein`（後述）→ `none`
5. いずれも無ければ `hypothetical protein`（原則2）→ `none`

`<name>` は InterPro エントリの **`name` フィールド**（短縮名）を使う。
`description`（長い記述文）を使うと `... domain domain-containing protein` の
重複が生じ、最長174文字となって GenBank `/product` として不適切なため採用しない
（`name` 使用時の最長は90文字）。

## 規則で決めなかったもの、およびその理由

### Family 同士に階層関係がない200件（2.70%）— 一貫性による判断

Family エントリが複数あり、かつ相互に祖先・子孫関係がないため規則1で順位が
付かないケースが200件ある。例: `APU25_00260` で `IPR002347 SDR_fam` と
`IPR052625 Chl_b_Red` が並立する。

**全件 `hypothetical protein` / `product_source = none` とした。
これは正確性による判断ではなく、一貫性（原則4）による判断である。**

根拠は実測である。200件の内訳は `matches_api` 188件 / `ips6_local` 12件。
競合の82.00%（164件）が IPR05xxxx 番台を含み、**その164件のうち159件（96.95%）が
`matches_api` 由来**だった。IPR05xxxx 番台は PANTHER サブファミリー由来の
高特異性エントリ群で、本環境では PANTHER をローカル導入していない。
全7,413配列で見ても、IPR05xxxx を持つ割合は `matches_api` 20.98% に対し
`ips6_local` 6.35% と約3.3倍の開きがある。

つまりここで特異的なほうを採ると、**遺伝子の性質ではなく、たまたま通った
取得経路によって product の粒度が変わる**。原則5（特異的な名前を優先）と
原則4（経路によらず粒度が揃っていること）が衝突し、原則4を優先した。

IPR05xxxx を含まない36件も分離せず同一に扱った。**これはコストによる判断である。**
(a) うち29件（80.6%）も `matches_api` 由来で偏りの方向が同じ、
(b) 36件（0.49%）のために別規則・別 source 値・InterPro API 依存を増やすのは
複雑さに見合わない。正確性の判断ではない。

競合した候補は捨てていない。`interpro_family_candidates` 列にアクセッション昇順で
全件保持し、`interpro_family_unresolved = TRUE` を立ててある。GenBank では
`/product="hypothetical protein"` とし、CDS の `/note` に
`unresolved InterPro family assignment: IPRxxxxxx <name>; ...` を記録した。

### 構造的手がかりのみの331件（4.47%）— 規則3.5

InterPro エントリはあるが `Homologous_superfamily` / `Repeat` /
`Conserved_site` のみで、Family も Domain も無い。product には使わない。

- `Homologous_superfamily`（SUPERFAMILY・CATH-Gene3D 由来）は**構造フォールドの
  類似のみを主張し、機能の共有を含意しない**。P-loop containing nucleoside
  triphosphate hydrolase superfamily にはヘリカーゼ・キナーゼ・ABCトランスポーター・
  Gタンパク質・RecA/RAD51 がすべて含まれる。product に書いても調査は前に進まない（用途1）。
- これを `<name> domain-containing protein` と書くのは**型の詐称**である。
  superfamily は domain ではなく、Pfam 相当のドメイン同定があったと読み手に
  誤解させる（原則2・原則3）。
- `Repeat` / `Conserved_site` も同様に機能を決定しない。

該当エントリのアクセッションと name は `annotation/UTEX25_interpro.tsv` の
`interpro_accessions` 列に従来どおり保持している。**product に使わないことと
情報を捨てることは別である。** `interpro_structural_only = TRUE` で再抽出できる。

### 規則3（eggNOG Description）は適用不能

emapper 3.0.0-beta6 の出力には **`Description` 列が存在しない**（22列: query,
seed_ortholog, evalue, score, eggNOG_OGs, tax_ceiling, farthest_donor_lineage,
COG_category, Preferred_name, GOs, EC, KEGG_ko, KEGG_Pathway, KEGG_Module,
KEGG_Reaction, KEGG_rclass, BRITE, KEGG_TC, CAZy, BiGG_Reaction, PFAMs,
annotation_confidence）。したがって規則3は一度も発火せず、
`product_source = eggnog` は **0件**である。

規則3が対象にしえた行（InterPro エントリが無く eggNOG 行がある）は788件あり、
それらは規則4に落ちて `hypothetical protein` のままになっている。

## 指標の定義

**A〜D は別々の指標である。相互に足し引きして導出してはならない。
どれか一つを「アノテーション率」と呼んではならない。** 全て 7,413 配列に対する実測値。

| 指標 | 定義 | 件数 | 割合 |
|---|---|---|---|
| **A** | InterPro エントリを1つ以上持つ配列数（**全ライブラリ**。`matches_api` 側の PANTHER・CATH-Gene3D 等を含む） | **5,926** | 79.94 % |
| **A'** | 同（**共通7ライブラリ由来のみ**。両経路で意味を揃えた値） | **5,833** | 78.69 % |
| **B** | `Family` または `Domain` を1つ以上持つ配列数 | **5,595** | 75.48 % |
| **C** | product が `hypothetical protein` のまま残った配列数 | **2,018** | 27.22 % |
| **D** | symbol が空欄のまま残った配列数 | **2,090** | 28.19 % |

A と A' は母集団のライブラリ範囲が違うだけで、どちらも「InterPro エントリの有無」を
数えている。B はさらに型を Family / Domain に限る。C と D は product / symbol 列を
直接数えたもので、A・B とは母数も対象も異なる。

`docs/INTERPRO_COVERAGE.md` の「ドメインが1つも付かなかった配列 1,411件」は
また別の定義（共通7ライブラリのシグネチャも、どのライブラリ由来の InterPro
エントリも付かない）であり、上表のどれとも一致しない。

### source 別内訳

| `product_source` | 件数 | 割合 |
|---|---|---|
| `interpro_family` | 3,194 | 43.09 % |
| `interpro_domain` | 1,203 | 16.23 % |
| `interpro_domain_multi` | 998 | 13.46 % |
| `eggnog` | **0** | 0.00 % |
| `none` | 2,018 | 27.22 % |

| `symbol_source` | 件数 | 割合 |
|---|---|---|
| `utex250a` | 2,418 | 32.62 % |
| `eggnog` | 2,905 | 39.19 % |
| `none` | 2,090 | 28.19 % |

`eggnog` 由来 2,905件のうち **523件は `LOC` + 数字の形**（例 `LOC107807905`）、
3件はバックスラッシュを含む他生物の識別子（例 `Dgri\GH23632`）である。
規則2どおり `Preferred_name` をそのまま入れた結果であり、遺伝子シンボルとしての
情報量は低い。`symbol_source = eggnog` で一括して除外・再検討できる。

## eggNOG の Pfam の扱い

**eggNOG の Pfam は参考値である。** eggNOG-mapper が参照した Pfam の版は
出力のどこにも記録されておらず特定できない（出力された Pfam 名のうち135種類が
Pfam 38.2 の名前表に存在しないことから、38.2 でないことだけは分かる）。
したがって **InterProScan の Pfam と直接比較できない。**

`docs/PFAM_CONCORDANCE.md` の一致率 **F1 = 83.87% は下限値**である。不一致には
Pfam リリース差・clan 競合の実装差・パイプライン差が交絡として含まれている。

**Cantalapiedra et al. (2021) *Mol Biol Evol* 38(12):5825
(DOI: 10.1093/molbev/msab293) の F1 = 89.7% とは比較対象が異なるため、
並べて論じてはならない。** 同ベンチマークは eggNOG 内部で同一 Pfam 版・
同一パイプラインのもとで transfer と de novo を比べた値であり、本値は
eggNOG の transfer と InterProScan の一致率である。

`eggnog_PFAMs` 列は InterProScan との差分検証のために保持してある。
eggNOG 由来の `KEGG_ko` / `COG_category` / `EC` / `CAZy` / `BiGG_Reaction` は
InterProScan が提供できないため、`annotation/UTEX25_gene_table_eggnog.tsv` に
全て保持している。

## 重複遺伝子モデル 418組は未解消

隣接する locus_tag が同一の eggNOG `seed_ortholog` を持つ組が **418組**ある。
連続する3個・4個の並びを1群にまとめると **397群 / 815行**になる。
`dup_pair_id` 列（`DUP0001`〜`DUP0397`）で識別できる。

**フラグを立てただけで、削除も統合もしていない。** したがって
**実遺伝子数は 7,413 ではなく約 7,000 と見るべきである。**
7,413 は「投影された遺伝子モデル数」であって「遺伝子数」ではない。

## QC 件数の2つの値について

`README.md` の「QC 完全通過 6,786（91.5%）」と `UTEX25_gene_table.tsv` の
`QC == pass` 6,762（91.22%）は**定義が異なるだけで、どちらも誤りではない**。
差の24件は frameshift フラグのみを持つモデルの扱いの違いによる。
**両方の数字とも変更していない。**

## v1.1 への申し送り

以下の3群は Foldseek + ProstT5（Heinzinger 2024,
DOI: 10.1093/nargab/lqae150）で最も改善が見込める対象である。
いずれも列で再抽出できるようにしてある。

| 群 | 列 | 件数 |
|---|---|---|
| Family 競合により未確定 | `interpro_family_unresolved = TRUE` | 200 |
| 構造的手がかりのみ | `interpro_structural_only = TRUE` | 331 |
| InterPro エントリなし | `annotation/UTEX25_interpro.tsv` の `n_domains = 0` かつ InterPro 列が空 | 1,487 |

PANTHER を導入した時点で、`interpro_family_unresolved` の200件は
一括して見直すこと（経路依存が解消されるため原則4の制約が外れる）。

## 再現手順

```bash
python3 scripts/13_finalize_annotation.py     # symbol / product を確定
python3 scripts/14_regenerate_derived.py      # GenBank / GFF3 / FASTA に反映
python3 scripts/15_annotation_metrics.py      # 指標 A〜D
```

`13` は既に新列がある表に対しては実行を拒否する（二重適用の防止）。
`14` は配列本体（ORIGIN・`/translation`・FASTA 本体）が1文字も変わらないことを
実行のたびに検証し、変化があれば異常終了する。
