# Pfam 呼び出しの突き合わせ — InterProScan 6 vs eggNOG-mapper

- 生成: `scripts/compare_pfam.py`
- InterPro 側: `step0_out/matches_raw.json.gz`（InterPro Matches API の生レスポンス、`signatureLibraryRelease.library == "Pfam"` のみ）
- eggNOG 側: `raw/eggnog_7413.emapper.annotations` の `PFAMs` 列
- Pfam 名 → アクセッション対応: `raw/Pfam-A.clans.tsv.gz` + `raw/Pfam-A.dead.gz`

数値はすべて実測。文献値は最後の 1 段落のみで、出典を明示している。

## 0. 比較キーの正規化（実測にもとづく前提）

eggNOG-mapper 3.0.0-beta6 の `PFAMs` 列は **PF アクセッションを出力しない**。
出力されるのは Pfam 名 + ヒット座標（例 `AAA_lid_3_330_372`）で、
本データの Pfam トークン 9,610 個のうち PFxxxxx 形のものは **0 個**だった。
そのため以下の正規化を行っている。

1. 末尾の `_<start>_<end>` を除去して Pfam 名を得る（形が違えば異常終了）
2. Pfam 公式の名前→アクセッション表で PFxxxxx に変換
3. InterPro 側は `signature.accession` からバージョン接尾辞を除去
4. タンパク質ごとに**集合**として比較

変換できなかった Pfam 名: **135 種類 / 延べ 247 出現**。
これは eggNOG の参照 Pfam 版と現行 Pfam 版で名前が変わった家系で、
アクセッション単位の比較からは除外している。除外は eggNOG 側の集合のみを
小さくするので、**「両方」を減らし「InterPro のみ」を増やす向き**に働く
（eggNOG の recall を過小評価しうる）。影響の大きさを見るため、
変換を必要としない **Pfam 名単位**の比較を §5 に併記した。

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

## 1. 対象タンパク質

| 集合 | 件数 |
|---|---|
| InterPro 側に結果がある（UniParc ヒット） | 3899 |
| eggNOG 側に行がある | 6588 |
| **両方に結果がある = 比較対象** | **3580** |
| InterPro のみ結果あり | 319 |
| eggNOG のみ結果あり | 3008 |

比較対象 3580 のうち、Pfam を 1 つ以上持つのは InterPro 側 3045、eggNOG 側 2965。
両側とも 0 個だったタンパク質は 427 で、
これは「完全一致」に含めて数えている。

## 2. ドメイン呼び出し単位の混同行列（PF アクセッション）

（タンパク質 × PF アクセッションの組を 1 呼び出しと数える）

| 区分 | 呼び出し数 |
|---|---|
| 両方にある | 3544 |
| InterPro のみ（eggNOG の取りこぼし） | 700 |
| eggNOG のみ | 663 |
| **合計（和集合）** | **4907** |

- InterPro 側の呼び出し総数: 4244
- eggNOG 側の呼び出し総数: 4207

## 3. InterProScan 6 を正解としたときの eggNOG の性能

| 指標 | 値 |
|---|---|
| precision | **84.24 %** |
| recall | **83.51 %** |
| F1 | **83.87 %** |

precision = 両方 / (両方 + eggNOG のみ)、
recall = 両方 / (両方 + InterPro のみ)。

## 4. タンパク質単位の一致率

- ドメイン集合が完全一致: **2672 / 3580 (74.64 %)**
  - うち両側とも 0 ドメイン: 427
- 少なくとも一方が 1 ドメイン以上ある 3153 に限った完全一致: **2245 (71.20 %)**

## 5. 感度確認 — Pfam 名単位（対応表を経由しない比較）

§0 の変換損失が結論を動かしていないかの確認。両側とも Pfam 名を持つので
対応表が不要で、eggNOG のトークンを 1 つも落とさずに比較できる。

| 指標 | アクセッション単位 (§3) | 名前単位 |
|---|---|---|
| 両方にある | 3544 | 3544 |
| InterPro のみ | 700 | 700 |
| eggNOG のみ | 663 | 773 |
| precision | 84.24 % | 82.09 % |
| recall | 83.51 % | 83.51 % |
| F1 | 83.87 % | 82.79 % |
| タンパク質単位の完全一致 | 74.64 % | 74.50 % |

「両方にある」と「InterPro のみ」は 2 つの単位で**同数**だった。
すなわち §0 で変換できなかった eggNOG 名は 1 つも InterPro 側と
一致していない（InterPro も現行 Pfam 名を使うため当然ではある）。
したがって変換損失は recall を動かしておらず、影響は eggNOG のみ側の
110 件、すなわち precision 84.24 % → 82.09 % の差に限られる。

## 6. 上位 20 ドメイン

### 6.1 eggNOG が取りこぼした（InterPro のみ）上位 20

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

### 6.2 eggNOG のみにある上位 20

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

## 7. 文献値 F1 = 89.7% との比較

本データでの実測は F1 = **83.87 %**（アクセッション単位、対象 3580 タンパク質、precision 84.24 % / recall 83.51 %）で、名前単位でも **82.79 %** だった。
Cantalapiedra et al. (2021) *Mol Biol Evol* 38(12):5825 (DOI: 10.1093/molbev/msab293) が転写モードの Pfam 呼び出しについて報告した F1 = 89.7%（realign 時 98.9%）と、直接は比較できない。同ベンチマークの正解は de novo の Pfam 呼び出しであるのに対し、ここでの正解は InterProScan 6 が UniParc 収録配列に対して算出した結果であり、対象生物も Progenomes（原核）ではなく緑藻 *Auxenochlorella protothecoides* である。eggNOG における緑藻の代表性は低く、転写元となる seed ortholog が系統的に遠くなるため、誤差は文献値より悪い方向に振れうる。実測値がその方向に出ているか否かは上の表の数値そのものを参照のこと。本節は 2 つの数値の由来が違うことの注記であって、文献値を本データの期待値として採用するものではない。

