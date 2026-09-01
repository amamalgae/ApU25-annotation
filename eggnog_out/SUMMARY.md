# eggNOG-mapper 機能アノテーションの取り込み

生成: `scripts/ingest_eggnog.py` / 入力 `raw/eggnog_7413.emapper.annotations`
構造アノテーション: `annotation/UTEX25_gene_table.tsv`（7413 遺伝子）
出力: `eggnog_out/UTEX25_gene_table_eggnog.tsv`

## 1. 取り込んだファイルのヘッダ（実測、決め打ちなし）

- 検出したヘッダ行のクエリ列名: `query`
- 列数: 22
- 列名: `query`, `seed_ortholog`, `evalue`, `score`, `eggNOG_OGs`, `tax_ceiling`, `farthest_donor_lineage`, `COG_category`, `Preferred_name`, `GOs`, `EC`, `KEGG_ko`, `KEGG_Pathway`, `KEGG_Module`, `KEGG_Reaction`, `KEGG_rclass`, `BRITE`, `KEGG_TC`, `CAZy`, `BiGG_Reaction`, `PFAMs`, `annotation_confidence`
- データ行数: **6588**

## 2. 実行時の記録（`##` 行のまま）

- version: `emapper-3.0.0-beta6`
- command:

```
/usr/local/bin/emapper.py -i input.fasta -o query --output_dir output --temp_dir output --data_dir /eggnog-data --cpu 10 --override --itype proteins -m diamond --dmnd_iterate no --dmnd_block_size 2 --dmnd_index_chunks 4
```

- applied filters:

| key | value |
|---|---|
| `annot_evalue` | `0.001` |
| `annot_score` | `null` |
| `donor_pool` | `closest` |
| `excluded_taxa` | `null` |
| `go_namespace_split` | `True` |
| `lazy_cascade` | `True` |
| `pfam_realign` | `none` |
| `sort_entries` | `True` |
| `target_orthologs` | `all` |
| `target_taxa` | `null` |
| `tax_scope` | `auto` |

> ウェブサーバのジョブページ表示と、この `##` 行に記録された実行コマンドは
> 一致しないことが確認されている。ここに書いてあるのは **実際に走った設定** の方。

## 3. 充填率（データ行 6588 行に対して）

| 列 | 値の入った行 | 率 |
|---|---|---|
| `seed_ortholog` | 6588 | 100.00 % |
| `evalue` | 6588 | 100.00 % |
| `score` | 6588 | 100.00 % |
| `eggNOG_OGs` | 6516 | 98.91 % |
| `tax_ceiling` | 6588 | 100.00 % |
| `farthest_donor_lineage` | 6588 | 100.00 % |
| `COG_category` | 6240 | 94.72 % |
| `Preferred_name` | 5040 | 76.50 % |
| `GOs` | 4491 | 68.17 % |
| `EC` | 2521 | 38.27 % |
| `KEGG_ko` | 4627 | 70.23 % |
| `KEGG_Pathway` | 3354 | 50.91 % |
| `KEGG_Module` | 970 | 14.72 % |
| `KEGG_Reaction` | 0 | 0.00 % |
| `KEGG_rclass` | 0 | 0.00 % |
| `BRITE` | 4626 | 70.22 % |
| `KEGG_TC` | 0 | 0.00 % |
| `CAZy` | 133 | 2.02 % |
| `BiGG_Reaction` | 348 | 5.28 % |
| `PFAMs` | 5497 | 83.44 % |
| `annotation_confidence` | 6588 | 100.00 % |

## 4. 遺伝子表への結合（7413 遺伝子に対して）

- eggNOG 行が付いた遺伝子: **6588 / 7413 (88.87 %)**
- eggNOG 行が無い遺伝子: 825

| 列 | 値の入った遺伝子 | 全遺伝子に対する率 |
|---|---|---|
| `eggnog_seed_ortholog` | 6588 | 88.87 % |
| `eggnog_evalue` | 6588 | 88.87 % |
| `eggnog_score` | 6588 | 88.87 % |
| `eggnog_eggNOG_OGs` | 6516 | 87.90 % |
| `eggnog_tax_ceiling` | 6588 | 88.87 % |
| `eggnog_farthest_donor_lineage` | 6588 | 88.87 % |
| `eggnog_COG_category` | 6240 | 84.18 % |
| `eggnog_Preferred_name` | 5040 | 67.99 % |
| `eggnog_GOs` | 4491 | 60.58 % |
| `eggnog_EC` | 2521 | 34.01 % |
| `eggnog_KEGG_ko` | 4627 | 62.42 % |
| `eggnog_KEGG_Pathway` | 3354 | 45.24 % |
| `eggnog_KEGG_Module` | 970 | 13.09 % |
| `eggnog_KEGG_Reaction` | 0 | 0.00 % |
| `eggnog_KEGG_rclass` | 0 | 0.00 % |
| `eggnog_BRITE` | 4626 | 62.40 % |
| `eggnog_KEGG_TC` | 0 | 0.00 % |
| `eggnog_CAZy` | 133 | 1.79 % |
| `eggnog_BiGG_Reaction` | 348 | 4.69 % |
| `eggnog_PFAMs` | 5497 | 74.15 % |
| `eggnog_annotation_confidence` | 6588 | 88.87 % |

結合は左外部結合。eggNOG 行が無い遺伝子の追加列は**空欄**で、推定値は入れていない。
既存の `product` と `symbol` は一切書き換えていない。

## 5. 参照のみ（統合していないファイル）

### `raw/eggnog_noPfam1916.emapper.annotations`

- データ行数: **1091**
- `PFAMs` が埋まった行: **0**
- 記録された command: `/usr/local/bin/emapper.py -i input.fasta -o query --output_dir output --temp_dir output --data_dir /eggnog-data --cpu 10 --override --itype proteins -m diamond --dmnd_iterate no --dmnd_block_size 2 --dmnd_index_chunks 4`
- 取り込み対象ファイルの command と一致

