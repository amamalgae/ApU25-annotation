# step 0 — InterPro Matches API lookup (UTEX 25 projected proteome)

実行日: 2026-09-01 / 入力: `annotation/UTEX25_proteins.faa` (7,413 配列)
スクリプト: `scripts/05_interpro_matches.py` (`probe` → `fetch` → `report`)

すべて実測値。取得できなかった項目は「取得できなかった」と明記している。

## 1. エンドポイント（未確定事項 (a) の確定）

`https://www.ebi.ac.uk/interpro/matches/api/docs` は **HTTP 404**。
`https://www.ebi.ac.uk/interpro/matches/api/openapi.json` が **HTTP 200** で、
そこに実際のパスが書かれていた（`step0_out/openapi.json` に保存）。

| 項目 | 実測 |
|---|---|
| **200 を返した URL** | `https://www.ebi.ac.uk/interpro/matches/api/matches` (POST) |
| OpenAPI `servers[0].url` | `/interpro/matches/api` |
| 定義されているパス | `POST /matches`（最大100件）, `GET /matches/{md5}`（1件） |
| API バージョン (`info.version`) | `0.6.0` — "InterPro Matches API" |
| リクエスト形式 | `Content-Type: application/json`, body `{"md5": [...]}`（大文字16進、1〜100件） |

`probe` が候補の第1 URL で 200 を得たため、`https://www.ebi.ac.uk/interpro/matches/api/`
（末尾スラッシュ）は試行していない。

## 2. レスポンス JSON の実際のキー名（未確定事項 (b) の確定）

1件だけ POST した生レスポンスは `step0_out/probe_response.json` に保存済み。
実際に返ってきたキー名は以下（OpenAPI の `BatchResponse` / `BatchResult` / `Match`
/ `Signature` / `Entry` / `Location` と一致することも確認）。

```
{"results": [ {"md5", "found", "matches"} ]}
```

| 階層 | 実際のキー名 |
|---|---|
| トップレベル | `results`（配列。これ1つのみ） |
| `results[]` | `md5`, `found`, `matches` |
| `matches[]` | `signature`, `model-ac`, `source`, `locations`, `score`, `evalue`（PANTHER のみ `ancestralNode`, `graphscan` が付く） |
| `signature` | `accession`, `name`, `description`, `type`, `signatureLibraryRelease`, `entry` |
| `signatureLibraryRelease` | `library`, `version` |
| `signature.entry` | `accession`(=IPR), `name`, `description`, `type`, `parent` — 未対応シグネチャでは `null` |
| `locations[]` | `start`, `end`, `location-fragments`, `hmmStart`, `hmmEnd`, `hmmLength`, `hmmBounds`, `envelopeStart`, `envelopeEnd`, `evalue`, `score` |

**未ヒットの扱い（実測）**: 送った MD5 は必ず全件返る。UniParc に無い MD5 は
`"found": false` かつ `"matches": []` で返る（32桁のゼロを混ぜて実測確認）。
したがって `found` が UniParc 収録の判定に使える唯一の権威ある値であり、
集計はこれを使っている（`matches` の空/非空ではない）。

**GO について**: このレスポンスには GO 項目が一切含まれない。
OpenAPI スキーマに GO 用のフィールドが定義されておらず（`goXRefs` 等は存在しない）、
全75バッチの生レスポンス 3.1 MB を走査しても `GO:` で始まる文字列は **0 件**だった。
よって下の GO 数は「API が返さないので全て 0」であり、
「この遺伝子群に GO が無い」という意味ではない。GO は別途 InterPro entry → GO の
対応表（`interpro2go`）を引く必要がある。

## 3. 取得

| 項目 | 値 |
|---|---|
| ユニーク MD5 | 7,413（完全重複配列 0 のため配列数と同数） |
| リクエスト数 | 75（100件/回） |
| レート | 2.5 req/s（≦3 req/s） |
| 429 / 5xx による再試行 | 0 回 |
| 失敗バッチ | **0** |
| 生レスポンス | `step0_out/matches_raw.json.gz`（gzip, 3.1 MB） |

## 4. ヒット率

| 母数 | ヒット | 率 |
|---|---|---|
| 全 7,413 配列 | 3,899 | **52.60 %** |
| 内部終止 235 を除く 7,178 配列 | 3,899 | **54.32 %** |

内部終止 `*` を含む 235 配列のヒットは **0 件**（実測）。
そのため 2 つの母数でヒット数は同じ 3,899 になる。
未ヒットは 3,514 配列で、`step0_out/unmatched.faa` に出力した（= `in_uniparc` が `no` の行数と一致）。

## 5. ヒットした 3,899 配列の統計

| 指標 | 最小 | Q1 | 中央値 | Q3 | 最大 | 平均 |
|---|---|---|---|---|---|---|
| InterPro エントリ数 (IPR) | 0 | 1 | **3** | 4 | 20 | 3.05 |
| Pfam 数 (PF) | 0 | 1 | **1** | 1 | 8 | 1.10 |
| GO 数 | 0 | 0 | **0** | 0 | 0 | 0.00 |
| シグネチャマッチ総数 | 0 | 5 | 7 | 11 | 95 | 9.08 |

四分位は `statistics.quantiles(n=4, method="inclusive")`、母集団はヒットした
3,899 配列すべて（0 件のものを含む）。

- InterPro エントリが 1 つ以上付いたもの: 3,294 / 3,899 (84.48 %)
- Pfam が 1 つ以上付いたもの: 3,077 / 3,899 (78.92 %)
- UniParc にはあるがマッチ 0 件: 80 配列

### シグネチャライブラリ別のマッチ数（ヒット配列の合計）

| ライブラリ | version | マッチ数 |
|---|---|---|
| Phobius | 1.01 | 4,914 |
| PIRSR | 2025_05 | 4,668 |
| Pfam | 38.2 | 4,279 |
| CATH-Gene3D | 4.3.0 | 4,151 |
| SUPERFAMILY | 1.75 | 3,232 |
| PANTHER | 19.0 | 3,133 |
| MobiDB-lite | 4.0 | 1,997 |
| PROSITE profiles | 2026_01 | 1,553 |
| CDD | 3.21 | 1,505 |
| CATH-FunFam | 4.3.0 | 1,360 |
| SMART | 9.0 | 1,261 |
| PROSITE patterns | 2026_01 | 885 |
| NCBIFAM | 19.0 | 793 |
| PRINTS | 42.0 | 544 |
| COILS | 2.2.1 | 484 |
| HAMAP | 2026_01 | 297 |
| PIRSF | 3.10 | 254 |
| SFLD | 4 | 95 |

## 6. 成果物

| ファイル | 内容 |
|---|---|
| `lookup_status.tsv` | 7,413 行 + ヘッダ = 7,414 行。`protein_id` / `md5` / `length` / `has_internal_stop` / `in_uniparc` / `n_interpro` / `n_pfam` / `n_go`。未ヒット行の計数3列は空欄 |
| `unmatched.faa` | 3,514 配列（`in_uniparc` = `no` の行数と一致） |
| `matches_raw.json.gz` | 75 バッチの生レスポンス（gzip） |
| `probe_response.json` | 1 件だけ POST した生レスポンス（スキーマ確定の根拠） |
| `openapi.json` | 取得した OpenAPI 定義（エンドポイント確定の根拠） |
| `SUMMARY.md` | 本ファイル |

`md5_table.tsv` は `lookup_status.tsv` に `length` / `has_internal_stop` 列として
統合したため削除した。

## 7. 次段階

未ヒット 3,514 配列（内部終止 235 を含む）はローカル InterProScan 6 にかける必要がある。
GO は API から取れないため、`interpro2go` を別途引くこと。
