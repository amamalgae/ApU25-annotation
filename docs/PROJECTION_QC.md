# 投影品質 × UniParc 収録のクロス集計

- 構造アノテーション: `annotation/UTEX25_gene_table.tsv`（`identity` 列, `QC` 列）
- UniParc 収録: `step0_out/lookup_status.tsv`（`in_uniparc` 列）
- 生成: `scripts/07_projection_qc.py`
- 対象: 7413 遺伝子（全件が両表に存在）

実測値のみ。解釈は書いていない。

## 1. identity = 100% / < 100% × in_uniparc

| identity | in_uniparc = yes | in_uniparc = no | 計 | ヒット率 |
|---|---|---|---|---|
| = 100% (1.0000) | 2163 | 22 | 2185 | 98.99 % |
| < 100% | 1736 | 3492 | 5228 | 33.21 % |
| **計** | **3899** | **3514** | **7413** | **52.60 %** |

依頼された 2 つの数値:

- **identity = 100% かつ in_uniparc = no: 22 件** （identity = 100% の 2185 件の 1.01 %、全 7413 件の 0.30 %）
- **identity < 100% かつ in_uniparc = yes: 1736 件** （identity < 100% の 5228 件の 33.21 %、全 7413 件の 23.42 %）

## 2. identity の分布（in_uniparc 別）

| in_uniparc | n | 最小 | Q1 | 中央値 | Q3 | 最大 | 平均 |
|---|---|---|---|---|---|---|---|
| yes | 3899 | 0.8081 | 0.9955 | 1.0000 | 1.0000 | 1.0000 | 0.9965 |
| no | 3514 | 0.4623 | 0.9825 | 0.9918 | 0.9959 | 1.0000 | 0.9810 |
| 全体 | 7413 | 0.4623 | 0.9900 | 0.9959 | 1.0000 | 1.0000 | 0.9892 |

## 3. identity 区間別の in_uniparc ヒット率

| identity | 遺伝子数 | in_uniparc = yes | ヒット率 |
|---|---|---|---|
| = 1.0000 | 2185 | 2163 | 98.99 % |
| 0.9900 – 0.9999 | 3376 | 1364 | 40.40 % |
| 0.9500 – 0.9899 | 1614 | 365 | 22.61 % |
| 0.9000 – 0.9499 | 126 | 4 | 3.17 % |
| 0.8000 – 0.8999 | 71 | 3 | 4.23 % |
| < 0.8000 | 41 | 0 | 0.00 % |

## 4. QC フラグ（記載どおりの組み合わせ）別の in_uniparc ヒット率

| QC | 遺伝子数 | in_uniparc = yes | ヒット率 | identity 中央値 |
|---|---|---|---|---|
| `pass` | 6762 | 3898 | 57.65 % | 0.9964 |
| `no_start_Met` | 391 | 1 | 0.26 % | 0.9879 |
| `internal_stop,frameshift` | 143 | 0 | 0.00 % | 0.9798 |
| `internal_stop` | 62 | 0 | 0.00 % | 0.9640 |
| `frameshift` | 24 | 0 | 0.00 % | 0.9858 |
| `internal_stop,no_start_Met` | 16 | 0 | 0.00 % | 0.8598 |
| `internal_stop,no_start_Met,frameshift` | 14 | 0 | 0.00 % | 0.8505 |
| `no_start_Met,frameshift` | 1 | 0 | 0.00 % | 0.7135 |
| **計** | **7413** | **3899** | **52.60 %** | |

## 5. QC フラグ（個別、組み合わせを分解）別の in_uniparc ヒット率

1 遺伝子が複数フラグを持つため、行の合計は遺伝子数と一致しない。

| フラグ | 遺伝子数 | in_uniparc = yes | ヒット率 |
|---|---|---|---|
| `frameshift` | 182 | 0 | 0.00 % |
| `internal_stop` | 235 | 0 | 0.00 % |
| `no_start_Met` | 422 | 1 | 0.24 % |
| `pass` | 6762 | 3898 | 57.65 % |

## 6. QC フラグ別の「identity = 100% かつ in_uniparc = no」

| QC | identity = 100% | うち in_uniparc = no | 率 |
|---|---|---|---|
| `pass` | 2168 | 6 | 0.28 % |
| `no_start_Met` | 17 | 16 | 94.12 % |
| `internal_stop,frameshift` | 0 | 0 | n/a |
| `internal_stop` | 0 | 0 | n/a |
| `frameshift` | 0 | 0 | n/a |
| `internal_stop,no_start_Met` | 0 | 0 | n/a |
| `internal_stop,no_start_Met,frameshift` | 0 | 0 | n/a |
| `no_start_Met,frameshift` | 0 | 0 | n/a |
| **計** | **2185** | **22** | **1.01 %** |

