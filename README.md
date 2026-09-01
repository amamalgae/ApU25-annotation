# apu25-annotation

*Auxenochlorella protothecoides* UTEX 25 の遺伝子アノテーションを、
UTEX 250-A のタンパク質配列から miniprot で投影するパイプライン。

生成物（GenBank / GFF3 / FASTA / TSV）は git には含めず、GitHub Release に添付する。
入力はすべて公開データなので、このリポジトリだけで再現できる。

## 入力

| データ | 出所 |
|---|---|
| UTEX 25 ゲノム（12染色体, 21,976,416 bp） | BioProject PRJNA1328465 |
| UTEX 250-A haplotype 1 アノテーション | BioProject PRJNA1195245 |
| UTEX 250-A haplotype 2 アノテーション | BioProject PRJNA1195244 |

原著: Craig RJ, Dueñas MA, Camacho DJ, Gallaher SD, Blaby-Haas CE, Moseley JL, Merchant SS (2025)
*Targeted genetic manipulation and yeast-like evolutionary genomics in the green alga Auxenochlorella.*
Plant Cell 37(11):koaf259. DOI: 10.1093/plcell/koaf259

## 実行

```bash
git clone --depth 1 https://github.com/lh3/miniprot.git && (cd miniprot && make -j4)
python3 scripts/01_parse_250A.py     # 250-A GenBank -> hap1.faa / hap2.faa
bash    scripts/02_map_miniprot.sh   # miniprot マッピング（約2分）
python3 scripts/03_build_models.py   # A/B 統合・重複解消・翻訳
python3 scripts/04_write_genbank.py  # GenBank 出力
```

所要 4 分程度、ピークメモリ約 1.4 GB。外部データベースは不要。

## 出力

7,413 遺伝子モデル。QC 完全通過 6,786（91.5%）、翻訳ラウンドトリップ 100%。

## 注意

- 相同性投影であり、RNA-seq / IsoSeq による実験的裏づけはない
- UTEX 250-A に存在しない遺伝子は検出されない
- product の 67.4% は "hypothetical protein" のまま（機能アノテーションは別途）
- 構造アノテーションの著作性は原著者に帰属する。再配布時は必ず上記を引用すること

## ライセンス

スクリプト: MIT。生成物: 原データが INSDC 公開データであるため CC0 相当だが、
再配布時は Craig et al. 2025 の引用を必須とする。

## 機能アノテーション（構造アノテーションとは別工程）

### step 0 — InterPro Matches API

UniParc に既収録の配列について、InterProScan 6 の計算済み結果を MD5 で引く。

```bash
python3 scripts/05_interpro_matches.py probe
python3 scripts/05_interpro_matches.py fetch \
    --endpoint https://www.ebi.ac.uk/interpro/matches/api/matches
python3 scripts/05_interpro_matches.py report
```

7,413 配列中 3,899 (52.60%) が UniParc にヒット。未ヒット 3,514 配列は
`step0_out/unmatched.faa`（ローカル InterProScan 6 の入力）。
数値と実際のレスポンススキーマは `step0_out/SUMMARY.md`。

### eggNOG-mapper 結果の取り込み

```bash
python3 scripts/ingest_eggnog.py \
    --annotations raw/eggnog_7413.emapper.annotations \
    --reference   raw/eggnog_noPfam1916.emapper.annotations \
    --expect-rows 6588 --expect-pfam 5497
```

7,413 遺伝子中 6,588 (88.87%) に eggNOG 行が付く。結合結果は
`eggnog_out/UTEX25_gene_table_eggnog.tsv`（既存 14 列 + `eggnog_*` 21 列）で、
`annotation/UTEX25_gene_table.tsv` 自体は書き換えていない。
充填率と実行時の記録は `eggnog_out/SUMMARY.md`。

## リポジトリ構成

| パス | 内容 |
|---|---|
| `annotation/` | 生成物（GenBank 12染色体・GFF3・faa・fna・TSV） |
| `raw/` | 外部サービスから取得した生の機能アノテーション（eggNOG-mapper） |
| `scripts/` | 再現用パイプライン 4本 + 機能アノテーション 2本 |
| `step0_out/` | InterPro Matches API の結果 |
| `eggnog_out/` | eggNOG-mapper 取り込み結果 |
| `docs/ANNOTATION.md` | 生成物の詳細・統計・限界 |
| `docs/PUSH.md` | GitHub 運用メモ |
