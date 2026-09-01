# Pfam 呼び出しの突き合わせ — InterProScan 6 vs eggNOG-mapper

- 生成: `scripts/compare_pfam.py`
- InterPro 側: `step0_out/matches_raw.json.gz`（Matches API の生レスポンス、`signatureLibraryRelease.library == "Pfam"` のみ）
- eggNOG 側: `raw/eggnog_7413.emapper.annotations` の `PFAMs` 列
- Pfam 名 → アクセッション対応: `raw/Pfam-A.clans.tsv.gz` + `raw/Pfam-A.dead.gz`（版は `raw/Pfam.version.gz`）

数値はすべて実測。文献値は §9 の 1 段落のみで、出典を明示している。

## 0. 本解析の位置づけ

**本解析が測っているのは、eggNOG-mapper が転写 (transfer) した Pfam 呼び出しと、
InterProScan 6 の呼び出しとの一致率である。**
転写 (transfer) と de novo の直接比較ではない。

Cantalapiedra et al. (2021) の F1 = 89.7% は eggNOG **内部**で、
同一の Pfam 版・同一パイプラインのもとで transfer と de novo を比べた値であり、
本解析とは比較対象がそもそも異なる。本解析の不一致には、少なくとも次の 3 つが
交絡として含まれている。

1. **Pfam リリース差** — 両側が同じ Pfam 版を使っている保証がない（§2 に実測）
2. **clan 競合の解決方法の実装差** — 同一 clan 内の重複ドメインをどう1つに
   絞るかは実装依存で、eggNOG の転写と InterProScan では処理が異なる
3. **パイプライン差** — 閾値、領域の切り方、転写元 seed ortholog の選び方など

**したがって Cantalapiedra 2021 の 89.7% と本値を直接比較することはできず、
本値は「一致率の下限」とみなすべきである。**

§4–§6 では 2 つの版を並記している。**共通リリース限定版**は上記 1 の交絡を
測れる範囲で落としたものだが、**2 と 3 は依然として残る**。
限定版もまた transfer vs de novo の値ではない。

## 1. 比較キーの正規化（実測にもとづく前提）

eggNOG-mapper 3.0.0-beta6 の `PFAMs` 列は **PF アクセッションを出力しない**。
本データの Pfam トークン 9,610 個のうち PFxxxxx 形のものは **0 個**で、実際は Pfam 名に開始・終了座標を
アンダースコアで連結した形（例 `AAA_lid_3_330_372`）だった。そのため:

1. 末尾の数値 2 フィールド（座標）を除去して Pfam 名を得る（形が違えば異常終了）
2. Pfam 公式の名前→アクセッション表で PFxxxxx に変換
3. InterPro 側は `signature.accession` からバージョン接尾辞を除去
4. タンパク質ごとに**集合**として比較

変換できなかった Pfam 名: **135 種類 / 延べ 247 出現**（全 9,610 トークンの 2.57 %）。
これらは現行 Pfam の名前表に存在しない名前で、eggNOG の参照 Pfam 版が
現行版と異なることの直接の証拠でもある（§2）。
なお「表に無い」は「その家系が現行版に存在しない」と同義ではない。
家系が残ったまま改名された場合も表引きは失敗するため、この 135 種類には
**廃止された家系と改名された家系の両方が混在している**。

変換できなかった名前（上位 20、延べ出現数）:

| Pfam 名 | 出現 |
|---|---|
| `HA2` | 12 |
| `SNARE_assoc` | 11 |
| `Exostosin` | 8 |
| `Chorein_N` | 6 |
| `UPF0016` | 6 |
| `Prot_ATP_ID_OB` | 5 |
| `zf-RanBP` | 5 |
| `DUF4217` | 4 |
| `E2F_TDP` | 4 |
| `Glyco_hydro_31` | 4 |
| `Kelch_4` | 4 |
| `Pyrid_oxidase_2` | 4 |
| `TFIIS_C` | 4 |
| `rRNA_proc-arch` | 4 |
| `IKI3` | 3 |
| `MS_channel` | 3 |
| `NOT2_3_5` | 3 |
| `Pectate_lyase_3` | 3 |
| `Rtt106` | 3 |
| `STI1` | 3 |

## 2. 両側のリリース版（実測）

### 2.1 InterPro 側 — Matches API のレスポンスから集計

`step0_out/matches_raw.json.gz` の全マッチの `signature.signatureLibraryRelease` を集計した。

| ライブラリ | version | マッチ数 |
|---|---|---|
| CATH-FunFam | `4.3.0` | 1360 |
| CATH-Gene3D | `4.3.0` | 4151 |
| CDD | `3.21` | 1505 |
| COILS | `2.2.1` | 484 |
| HAMAP | `2026_01` | 297 |
| MobiDB-lite | `4.0` | 1997 |
| NCBIFAM | `19.0` | 793 |
| PANTHER | `19.0` | 3133 |
| PIRSF | `3.10` | 254 |
| PIRSR | `2025_05` | 4668 |
| PRINTS | `42.0` | 544 |
| PROSITE patterns | `2026_01` | 885 |
| PROSITE profiles | `2026_01` | 1553 |
| Pfam **←** | `38.2` | 4279 |
| Phobius | `1.01` | 4914 |
| SFLD | `4` | 95 |
| SMART | `9.0` | 1261 |
| SUPERFAMILY | `1.75` | 3232 |

**Pfam の release 値は `38.2` の 1 種類のみ**（レスポンス中に他の版は現れない）。

### 2.2 eggNOG 側

emapper の出力ヘッダ（`annotation/eggnog_provenance.txt`）には Pfam の版が
記録されていない。実行コマンドにも `## applied filters:` にも該当項目が無く、
**eggNOG-mapper が参照した Pfam 版はこの成果物からは特定不能**である。

特定できるのは次の 2 点のみ:

- 名前→アクセッション対応表として使った Pfam のリリース版

```
Pfam release       : 38.2
Pfam-A families    : 30134
Date               : 2026-01
Based on UniProtKB : 2025_03
```

- eggNOG が出力した Pfam 名のうち 135 種類がこの版の名前表に
  無い → eggNOG の参照 Pfam 版は上記の版と**同一ではない**

対応表の Pfam 版 `38.2` は、InterProScan 6 が報告した Pfam 版 `38.2` と**一致する**。
よって以下の「共通リリース限定」は、「Pfam 38.2 の名前表で解決できる呼び出しに限定する」ことを意味する。

## 3. 対象タンパク質

| 集合 | 件数 |
|---|---|
| InterPro 側に結果がある（UniParc ヒット） | 3899 |
| eggNOG 側に行がある | 6588 |
| **両方に結果がある = 比較対象** | **3580** |
| InterPro のみ結果あり | 319 |
| eggNOG のみ結果あり | 3008 |

比較対象 3580 のうち、Pfam を 1 つ以上持つのは InterPro 側 3045、eggNOG 側 2965。
両側とも 0 個だったタンパク質は 427 で、これは「完全一致」に含めて数えている。

## 4. ドメイン呼び出し単位の混同行列

（タンパク質 × ドメインの組を 1 呼び出しと数える）

- **全体版**: キーは Pfam 名。eggNOG の呼び出しを 1 つも落とさない。
- **共通リリース限定版**: キーは PF アクセッション。Pfam 38.2 の名前表で解決できる呼び出しに限定。

| 区分 | 全体版 | 共通リリース限定版 |
|---|---|---|
| 両方にある | 3544 | 3544 |
| InterPro のみ（eggNOG の取りこぼし） | 700 | 700 |
| eggNOG のみ | 773 | 663 |
| 合計（和集合） | 5017 | 4907 |

限定によって各側が失う呼び出し（実測）:

- InterPro 側: **0 呼び出し**。比較対象で呼ばれた 2418 種類の PF アクセッションはすべて対応表に存在した。
- eggNOG 側: **110 呼び出し**（§1 の変換不能な 135 種類に由来）。

## 5. InterProScan 6 を正解としたときの eggNOG の性能

| 指標 | 全体版 | 共通リリース限定版 |
|---|---|---|
| precision | 82.09 % | **84.24 %** |
| recall | 83.51 % | **83.51 %** |
| F1 | 82.79 % | **83.87 %** |

precision = 両方 / (両方 + eggNOG のみ)、recall = 両方 / (両方 + InterPro のみ)。

「両方にある」と「InterPro のみ」は 2 つの版で**同数**である。
すなわち §1 で変換できなかった eggNOG 名は 1 つも InterPro 側と一致しておらず、
リリース差の除去は recall を動かさない。差は precision 82.09 % → 84.24 % のみに現れる。

## 6. タンパク質単位の一致率

| 指標 | 全体版 | 共通リリース限定版 |
|---|---|---|
| ドメイン集合が完全一致 | 2667 / 3580 (74.50 %) | 2672 / 3580 (74.64 %) |
| うち両側とも 0 ドメインを除く | 2242 / 3155 (71.06 %) | 2245 / 3153 (71.20 %) |

## 7. QC フラグ別の層別集計（共通リリース限定版）

比較対象 3580 タンパク質を `annotation/UTEX25_gene_table.tsv` の
`QC` 列で層別した。数値のみ。

遺伝子表に現れる QC 値は 8 種類あるが、比較対象に 1 件も含まれない値が 6 種類ある: `frameshift`、`internal_stop`、`internal_stop,frameshift`、`internal_stop,no_start_Met`、`internal_stop,no_start_Met,frameshift`、`no_start_Met,frameshift`。
これらは in_uniparc = no のため InterPro 側の結果が無く、比較対象に入らない（`docs/PROJECTION_QC.md` の QC 別ヒット率を参照）。

| QC | タンパク質数 | 両方 | InterProのみ | eggNOGのみ | precision | recall | F1 | 完全一致 |
|---|---|---|---|---|---|---|---|---|
| `pass` | 3579 | 3542 | 700 | 663 | 84.23 % | 83.50 % | 83.86 % | 2671 (74.63 %) |
| `no_start_Met` | 1 | 2 | 0 | 0 | 100.00 % | 100.00 % | 100.00 % | 1 (100.00 %) |
| **全体** | **3580** | 3544 | 700 | 663 | 84.24 % | 83.51 % | 83.87 % | 2672 (74.64 %) |

個別フラグに分解（1 タンパク質が複数フラグを持つため合計は一致しない）:

| フラグ | タンパク質数 | precision | recall | F1 | 完全一致 |
|---|---|---|---|---|---|
| `no_start_Met` | 1 | 100.00 % | 100.00 % | 100.00 % | 1 (100.00 %) |
| `pass` | 3579 | 84.23 % | 83.50 % | 83.86 % | 2671 (74.63 %) |

## 8. 上位 20 ドメイン（共通リリース限定版）

### 8.1 eggNOG が取りこぼした（InterPro のみ）上位 20

| # | PF アクセッション | Pfam 名 | 件数 |
|---|---|---|---|
| 1 | `PF00400` | `WD40` | 9 |
| 2 | `PF12796` | `Ank_2` | 9 |
| 3 | `PF00702` | `Hydrolase` | 7 |
| 4 | `PF13499` | `EF-hand_7` | 7 |
| 5 | `PF25390` | `WD40_RLD` | 6 |
| 6 | `PF04408` | `WHD_HA2` | 5 |
| 7 | `PF16450` | `Prot_ATP_ID_OB_C` | 5 |
| 8 | `PF21010` | `HA2_C` | 5 |
| 9 | `PF01096` | `Zn_ribbon_TFIIS` | 4 |
| 10 | `PF12937` | `F-box-like` | 4 |
| 11 | `PF17177` | `PPR_long` | 4 |
| 12 | `PF23231` | `HAT_Syf1_CNRKL1_C` | 4 |
| 13 | `PF00097` | `zf-C3HC4` | 3 |
| 14 | `PF00098` | `zf-CCHC` | 3 |
| 15 | `PF00149` | `Metallophos` | 3 |
| 16 | `PF00270` | `DEAD` | 3 |
| 17 | `PF00561` | `Abhydrolase_1` | 3 |
| 18 | `PF00924` | `MS_channel_2nd` | 3 |
| 19 | `PF02037` | `SAP` | 3 |
| 20 | `PF05347` | `Complex1_LYR` | 3 |

（異なるアクセッション 549 種類、延べ 700 件）

### 8.2 eggNOG のみにある上位 20

| # | PF アクセッション | Pfam 名 | 件数 |
|---|---|---|---|
| 1 | `PF00076` | `RRM_1` | 13 |
| 2 | `PF00400` | `WD40` | 12 |
| 3 | `PF00271` | `Helicase_C` | 9 |
| 4 | `PF00270` | `DEAD` | 6 |
| 5 | `PF00415` | `RCC1` | 6 |
| 6 | `PF13419` | `HAD_2` | 6 |
| 7 | `PF00004` | `AAA` | 5 |
| 8 | `PF00005` | `ABC_tran` | 4 |
| 9 | `PF00036` | `EF-hand_1` | 4 |
| 10 | `PF00176` | `SNF2-rel_dom` | 4 |
| 11 | `PF00249` | `Myb_DNA-binding` | 4 |
| 12 | `PF02518` | `HATPase_c` | 4 |
| 13 | `PF13202` | `EF-hand_5` | 4 |
| 14 | `PF13637` | `Ank_4` | 4 |
| 15 | `PF13857` | `Ank_5` | 4 |
| 16 | `PF17862` | `AAA_lid_3` | 4 |
| 17 | `PF00023` | `Ank` | 3 |
| 18 | `PF00153` | `Mito_carr` | 3 |
| 19 | `PF00481` | `PP2C` | 3 |
| 20 | `PF00583` | `Acetyltransf_1` | 3 |

（異なるアクセッション 495 種類、延べ 663 件）

## 9. 文献値 F1 = 89.7% との比較

本データでの実測は、共通リリース限定版で F1 = **83.87 %**（precision 84.24 % / recall 83.51 %、対象 3580 タンパク質）、全体版で F1 = **82.79 %** である。
Cantalapiedra et al. (2021) *Mol Biol Evol* 38(12):5825 (DOI: 10.1093/molbev/msab293) が転写モードの Pfam 呼び出しについて報告した F1 = 89.7%（realign 時 98.9%）とは、§0 に述べたとおり比較対象が異なるため、直接比較できない。同ベンチマークの正解は同一 Pfam 版・同一パイプラインでの de novo 呼び出しであり、本解析の正解は InterProScan 6 が UniParc 収録配列に対して算出した結果である。対象生物も Progenomes（原核）ではなく緑藻 *Auxenochlorella protothecoides* で、eggNOG における代表性は低い。本節は 2 つの数値の由来が違うことの注記であって、文献値を本データの期待値として採用するものではない。

