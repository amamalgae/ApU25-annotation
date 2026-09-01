# step 0 — InterPro Matches API による UniParc 事前計算結果の取得

ApU25（*Auxenochlorella protothecoides* UTEX 25）投影プロテオームに対し、
InterPro Matches API から InterProScan 6 の事前計算結果を引き当て、
未ヒット分だけをローカル InterProScan 6 に回すための前段。

## 状態: 未完了 — API へ到達できず

この実行環境の egress ポリシーが `ebi.ac.uk` 全体を遮断しているため、
API 問い合わせを 1 回も実行できていない。
**オフラインで確定できる部分（MD5 の算出まで）のみを本 PR に含める。**
ヒット率を含む API 依存の数値は、推測値を書かず未取得として明示する。

### 遮断の実測記録

egress プロキシは TLS の CONNECT 段階で 403 を返すため、パスに関係なく全遮断となる。

| URL | 結果 |
|---|---|
| `https://www.ebi.ac.uk/interpro/matches/api/openapi.json` | `000`（CONNECT に 403 / `connect_rejected`） |
| `https://www.ebi.ac.uk/interpro/matches/api/docs` | `000`（同上） |
| `https://www.ebi.ac.uk/` | `000`（同上） |
| `https://ebi.ac.uk/` | `000`（同上） |
| `https://ftp.ebi.ac.uk/pub/databases/interpro/releases/latest/` | `000`（同上） |
| `https://github.com/`（対照） | `400`（到達可・遮断されていない） |

プロキシの `recentRelayFailures` に記録された理由:
`gateway answered 403 to CONNECT (policy denial or upstream failure)` — `www.ebi.ac.uk:443`。

セルフホスト（README の Docker/Singularity 手順）も、データ書庫の配布元が
同じ `ftp.ebi.ac.uk` であるため、この環境では代替にならない。

## 実測できた項目

| 項目 | 値 |
|---|---|
| 使用 FASTA | `annotation/UTEX25_proteins.faa` |
| 配列数 | **7,413** |
| ユニーク MD5 数 | **7,413** |
| 完全重複配列（配列数 − ユニーク MD5 数） | **0** |
| 重複する protein_id | 0 |
| 内部終止コドン `*` を含む配列 | **235** |
| 内部終止を含まない配列 | 7,178 |
| 配列長（正規化後） | 最短 19 aa / 最長 17,905 aa |
| 必要リクエスト数（100 件/回） | **75** |

完全重複が 0 件のため、MD5 集約によるリクエスト削減は起きない（事前の見積りどおり）。

MD5 は仕様どおり「大文字化し末尾の `*` を除去した配列」に対して計算し、16 進大文字で保持している。
本 FASTA には末尾 `*` を持つ配列は 1 件も無く、`*` は計 816 文字すべてが配列内部に出現する。

## 未取得の項目（API 到達後に確定するもの）

| 項目 | 状態 |
|---|---|
| 実際に 200 を返した URL | **未取得** — 候補 2 件を試行できていない |
| レスポンス JSON の実際のキー名 | **未取得** — 生レスポンスを 1 件も観測していない |
| UniParc ヒット数・百分率（全 7,413 に対する率） | **未取得** |
| UniParc ヒット数・百分率（内部終止を除く 7,178 に対する率） | **未取得** |
| ヒット配列あたり InterPro エントリ数の中央値・四分位 | **未取得** |
| ヒット配列あたり Pfam 数の中央値・四分位 | **未取得** |
| ヒット配列あたり GO 数の中央値・四分位 | **未取得** |
| 未ヒット件数 | **未取得** |

内部終止を含む 235 配列が UniParc に存在せず必ずミスになる、というのは
依頼者から与えられた前提であり、本作業で検証したものではない。
実測したのは「内部終止を含む配列が 235 件ある」ことのみ。

## 生成物

| ファイル | 状態 |
|---|---|
| `md5_table.tsv` | **生成済み** — 7,413 行 + ヘッダ。`protein_id` / `md5` / `length` / `has_internal_stop` |
| `lookup_status.tsv` | 未生成（`in_uniparc` が API 依存のため） |
| `unmatched.faa` | 未生成（未ヒット集合が未確定のため） |
| `matches_raw.json.gz` | 未生成（生レスポンスが無いため） |

`lookup_status.tsv` を `in_uniparc` 未確定のまま置くと実測値と誤認されうるため、
オフラインで確定する 2 列のみを別名 `md5_table.tsv` として出力してある。

## API 到達可能な環境での実行手順

`scripts/05_interpro_matches.py` に全工程を実装済み。`ebi.ac.uk` に出られる環境で:

```sh
# 1. エンドポイント特定 + 生レスポンス 1 件の観測（未確定事項 (a) と (b)）
python3 scripts/05_interpro_matches.py probe

# 2. 全 75 リクエスト（2.5 req/s、429/5xx は 2/4/8/16 秒で最大 4 回再試行）
python3 scripts/05_interpro_matches.py fetch --endpoint <probe が 200 を返した URL>

# 3. 集計（lookup_status.tsv / unmatched.faa と SUMMARY 用の数値）
python3 scripts/05_interpro_matches.py report
```

`probe` は候補 2 件を順に POST し、どちらも 200 でなければ `/openapi.json` と `/docs` を
GET して保存する。200 を返した URL と生レスポンスは `probe_response.json` に残るので、
上表の「実際に 200 を返した URL」「レスポンス JSON の実際のキー名」はそこから確定できる。

### レスポンススキーマの扱い

スキーマ未観測のまま集計を書くと当て推量になるため、`resolve_by_md5()` は
想定しうる数通りの形（`{"results":[…]}` / `{"MD5":{…}}` / 素の配列 など）を受理し、
**どれにも当てはまらなければ観測されたキー名を添えて例外送出して停止する**。
黙って誤集計しない。合成データでの動作は確認済み。

InterPro / Pfam / GO の計数はペイロード内の全文字列を走査してアクセッション表記
（`IPR\d{6}` / `PF\d{5}` / `GO:\d{7}`）の異なり数を数える方式で、入れ子構造に依存しない。

なお `probe` 実行後に確認が要る点が 1 つある: UniParc に存在するがマッチ 0 件の配列に対し
API がエントリを返すのか、そもそも省くのか。現状は「エントリが返れば UniParc にあり」と
解釈している（マッチ 0 件でも `in_uniparc=yes`）。`probe_response.json` で要確認。
