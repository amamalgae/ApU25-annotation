# eggNOG-mapper アノテーションの取り込み

- 生成: `scripts/ingest_eggnog.py`
- 取り込み対象: `raw/eggnog_7413.emapper.annotations`
- 構造アノテーション: `annotation/UTEX25_gene_table.tsv`（7413 遺伝子）
- 出力表: `annotation/UTEX25_gene_table_eggnog.tsv`
- 実行条件の証跡（`#` 行の全文）: `annotation/eggnog_provenance.txt`

すべて実測値。ベンチマーク F1 の 2 値のみ文献値で、出典を明示している。

## 1. ヘッダ（決め打ちせず実ファイルから検出）

- 検出したクエリ列名: `query`
- 列数: 22 / データ行数: **6588**

```
#query	seed_ortholog	evalue	score	eggNOG_OGs	tax_ceiling	farthest_donor_lineage	COG_category	Preferred_name	GOs	EC	KEGG_ko	KEGG_Pathway	KEGG_Module	KEGG_Reaction	KEGG_rclass	BRITE	KEGG_TC	CAZy	BiGG_Reaction	PFAMs	annotation_confidence
```

## 2. ヘッダから抽出した実行コマンド（原文のまま）

- version: `emapper-3.0.0-beta6`

```
/usr/local/bin/emapper.py -i input.fasta -o query --output_dir output --temp_dir output --data_dir /eggnog-data --cpu 10 --override --itype proteins -m diamond --dmnd_iterate no --dmnd_block_size 2 --dmnd_index_chunks 4
```

`## applied filters:` ブロックの実測値:

| key | value |
|---|---|
| `annot_evalue` | `0.001` |
| `annot_score` | `null` |
| `donor_pool` | `closest` |
| `excluded_taxa` | `null` |
| `go_namespace_split` | `True` |
| `lazy_cascade` | `True` |
| `pfam_realign` | `none` **←** |
| `sort_entries` | `True` |
| `target_orthologs` | `all` |
| `target_taxa` | `null` |
| `tax_scope` | `auto` |

**`pfam_realign` の実測値: `none`**

ウェブサーバ (https://eggnog-mapper.cgmlab.org/, v3-beta6) のジョブページ表示と
この記録は一致しないことが確認されている。上に転記したのは記録された方、
すなわち実際に走った設定である。

## 3. `pfam_realign=none` が Pfam 呼び出しに与える影響（文献値）

`pfam_realign=none` は、Pfam ドメインを seed ortholog 経由の転写
(transfer) で付与し、クエリ配列に対する再アラインメントを行わない。
転写モードでの Pfam 呼び出しについて、de novo を正解としたときの
報告値は **F1 = 89.7%**、realign を行った場合は
**F1 = 98.9%** である（Cantalapiedra CP et al. (2021) *Mol Biol Evol* 38(12):5825. DOI: 10.1093/molbev/msab293）。

**この 2 値をそのまま本データに当てはめることはできない。**
当該ベンチマークは Progenomes（原核生物）ベースであり、eggNOG における
代表性が低い緑藻では、誤差はこれより悪い方向に振れうる。
本データでの実測は `docs/PFAM_CONCORDANCE.md`（InterProScan 6 を正解と
したときの一致度）を参照。

## 4. 各列のカバレッジ（全 7413 遺伝子に対して）

- eggNOG 行が付いた遺伝子: **6588 / 7413 (88.87 %)**
- eggNOG 行が無い遺伝子: 825 / 7413 (11.13 %)

| 列（出力表での名前） | 件数 / 7413 | % |
|---|---|---|
| （eggNOG 行そのもの） | 6588 / 7413 | 88.87 |
| `eggnog_seed_ortholog` | 6588 / 7413 | 88.87 |
| `eggnog_evalue` | 6588 / 7413 | 88.87 |
| `eggnog_score` | 6588 / 7413 | 88.87 |
| `eggnog_eggNOG_OGs` | 6516 / 7413 | 87.90 |
| `eggnog_tax_ceiling` | 6588 / 7413 | 88.87 |
| `eggnog_farthest_donor_lineage` | 6588 / 7413 | 88.87 |
| `eggnog_COG_category` | 6240 / 7413 | 84.18 |
| `eggnog_Preferred_name` | 5040 / 7413 | 67.99 |
| `eggnog_GOs` | 4491 / 7413 | 60.58 |
| `eggnog_EC` | 2521 / 7413 | 34.01 |
| `eggnog_KEGG_ko` | 4627 / 7413 | 62.42 |
| `eggnog_KEGG_Pathway` | 3354 / 7413 | 45.24 |
| `eggnog_KEGG_Module` | 970 / 7413 | 13.09 |
| `eggnog_KEGG_Reaction` | 0 / 7413 | 0.00 |
| `eggnog_KEGG_rclass` | 0 / 7413 | 0.00 |
| `eggnog_BRITE` | 4626 / 7413 | 62.40 |
| `eggnog_KEGG_TC` | 0 / 7413 | 0.00 |
| `eggnog_CAZy` | 133 / 7413 | 1.79 |
| `eggnog_BiGG_Reaction` | 348 / 7413 | 4.69 |
| `eggnog_PFAMs` | 5497 / 7413 | 74.15 |
| `eggnog_annotation_confidence` | 6588 / 7413 | 88.87 |

空欄は「eggNOG が値を返さなかった」であり、推定値は一切入れていない。
既存の `product` / `symbol` および `annotation/UTEX25_gene_table.tsv` 本体は書き換えていない。

## 5. 参照のみ（統合していないファイル）

### `raw/eggnog_noPfam1916.emapper.annotations`

- データ行数: **1091**
- `PFAMs` が埋まった行: **0**
- `pfam_realign`: `none`
- 記録された command は取り込み対象と**完全一致**

de novo を狙って再投入したジョブでも `pfam_realign=none` が記録され、
PFAM 付与は 0 件だった。ジョブページ表示と実行内容の不一致を示す直接証拠。

