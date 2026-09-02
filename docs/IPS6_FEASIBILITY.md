# InterProScan 6 ローカル実行の可能性判定

本書はこの環境（ディスク 30 GB / 4 vCPU / 16 GB RAM）で実行できるかの判定と、
そのための実測値の記録である。§0〜§7 は本番実行の**前**に行った判定で、
そのまま残してある。**本番実行は §9 に記録した（実施済み、exit 0）。**

数値はすべて実測。計測していない項目は「未計測」と明記し、推定値は書いていない。

## 0. 判定

| 問い | 判定 |
|---|---|
| この環境で InterProScan 6 6.0.1 は動くか | **動く**。Pfam のみ / 7 アプリケーション構成のいずれも exit 0 で完走 |
| 30 GB に収まる applications の組み合わせ | **Pfam, NCBIFAM, PROSITE patterns, PROSITE profiles, SMART, CDD, SUPERFAMILY**（+ InterPro コア）で **8.7 GB**（実測） |
| PANTHER / CATH-Gene3D を加えて収まるか | **未計測**（理由は §3.3） |
| Matches API と同じ InterPro 版に固定できるか | **できる**。`--interpro 109.0`。18 ライブラリの版が完全一致（§4） |
| `unmatched.faa` をそのまま投入できるか | **できない**。内部終止 `*` を含むため入力検証で落ちる（§5） |
| 3,514 配列の所要時間（Pfam のみ） | 実測 3 点からの外挿で **約 18 分**（§6） |
| 同（7 ライブラリ、本番実行の実測） | **4,320 s = 72 分**（§9） |

そのまま実行に移せる状態ではあるが、事前に 2 つの前提（§2）と入力の前処理（§5）が要る。

## 1. 環境（実測）

| 項目 | 実測値 | 取得方法 |
|---|---|---|
| 空きディスク | 30 GB（開始時） | `df -h /` の Avail |
| CPU | 4 | `nproc` |
| メモリ | 15.7 GiB | `docker info` の Total Memory |
| docker | client / daemon とも 29.3.1、storage-driver overlayfs | `docker info` |
| java | openjdk 21.0.10 | `java -version` |
| nextflow | **未インストール**。26.04.6 を `get.nextflow.io` から導入 | `command -v nextflow` |

**docker デーモンは起動していなかった**（`/var/run/docker.sock` が存在しない）。
`dockerd` を手動起動して使用した。cgroup v1 の非推奨警告は出るが動作する。

## 2. この環境で必要になった前提（実測）

### 2.1 パイプラインを GitHub から直接取れない

```
nextflow run ebi-pf-team/interproscan6 -r 6.0.1 -profile docker,test ...
  -> WARN: Cannot read project manifest -- Cause: Forbidden
     Forbidden -- Provide your GitHub user name and password to access this repository
```

nextflow が叩く GitHub API がセッションの egress プロキシに拒否される。
`git clone` は通るので、ローカル clone に対して実行すれば回避できる:

```sh
git clone --depth 1 -b 6.0.1 https://github.com/ebi-pf-team/interproscan6.git
nextflow run ./interproscan6 ...
```

### 2.2 コンテナがプロキシの CA を信頼しない

回避しないと `INIT_DATABASES:DOWNLOAD_INTERPRO` が全リトライで落ちる:

```
curl: (60) SSL certificate OpenSSL verify result: self-signed certificate in certificate chain (19)
```

コンテナに CA を渡し、プロキシに到達できるようにする設定を `-c` で足す:

```groovy
docker {
    enabled    = true
    runOptions = '--network host -v /root/.ccr/ca-bundle.crt:/ccr-ca.crt:ro'
}
env {
    CURL_CA_BUNDLE     = '/ccr-ca.crt'
    SSL_CERT_FILE      = '/ccr-ca.crt'
    REQUESTS_CA_BUNDLE = '/ccr-ca.crt'
    HTTPS_PROXY        = 'http://127.0.0.1:41551'
    https_proxy        = 'http://127.0.0.1:41551'
}
```

`--network host` が要るのは、プロキシが `127.0.0.1:41551` にあり、
bridge ネットワークのコンテナからは届かないため。

## 3. データ量

### 3.1 配布物のサイズ（FTP の一覧、実測）

データ元: `https://ftp.ebi.ac.uk/pub/software/unix/iprscan/6/6.0/`

| ディレクトリ | ファイル | サイズ |
|---|---|---|
| antifam | antifam-8.0.tar.gz | 2.5M |
| cath | cath-4.3.0.tar.gz | **2.5G** |
| cdd | cdd-3.21.tar.gz | 481M |
| hamap | hamap-2026_01.tar.gz | 18M |
| interpro | interpro-109.0.tar.gz | 11M |
| ncbifam | ncbifam-19.0.tar.gz | 402M |
| panther | panther-19.0.tar.gz | **1.0G** |
| pfam | pfam-38.2.tar.gz | 391M |
| pirsf | pirsf-3.10.tar.gz | 95M |
| pirsr | pirsr-2025_05.tar.gz | 63M |
| prints | prints-42.0.tar.gz | 7.4M |
| prosite | prosite-2026_01.tar.gz | 7.5M |
| sfld | sfld-4.tar.gz | 8.4M |
| smart | smart-9.0.tar.gz | 30M |
| superfamily | superfamily-1.75.tar.gz | 947M |

（各ディレクトリには旧版も置かれている。上表は 109.0 が参照する版のみ）

### 3.2 展開後のディスク使用量（`du -sh data`、実測）

| 構成 | `--applications` | data 合計 | 所要 | 終了 |
|---|---|---|---|---|
| Pfam のみ | `pfam` | **2.3 G** | 373 s | exit 0 |
| + 中規模 6 種 | `pfam,ncbifam,prositepatterns,prositeprofiles,smart,cdd,superfamily` | **8.7 G** | 1119 s | exit 0 |

内訳（8.7 G 構成、`du -sh data/*`）:

| ディレクトリ | 展開後 | 配布物 |
|---|---|---|
| interpro | 112M | 11M |
| prosite | 35M | 7.5M |
| smart | 160M | 30M |
| cdd | 1.1G | 481M |
| pfam | 2.1G | 391M |
| superfamily | 2.6G | 947M |
| ncbifam | 2.7G | 402M |
| **合計** | **8.7G** | 2.3G |

この構成の実行後、空きは 19 GB。

### 3.3 PANTHER / CATH-Gene3D — 未計測

**計測していない。理由は帯域制約。**

プロキシ経由の実効ダウンロード速度が約 0.7 MB/s で、
PANTHER (1.0G) と CATH (2.5G) の取得だけで数時間かかる見込みだったため、
得られる情報が制約判定のみであることに対して費用対効果が見合わないと判断し、中止した。

中止時点の実測:

- PANTHER: **266,223,745 バイト（254 MB）/ 1.0 G** をダウンロードした時点で中断
- CATH-Gene3D: **ダウンロードを開始していない**

したがって以下は**いずれも未計測**であり、本書に推定値は書かない:

- PANTHER を加えた場合の `data` 合計
- CATH-Gene3D を加えた場合の `data` 合計
- 両者を加えて 30 GB に収まるか否か

判定するには、この 2 つを実際に取得して `du -sh data` を測る必要がある。

## 4. InterPro リリースの一致（実測）

Matches API はリリース番号を返さず、シグネチャライブラリごとの版だけを返す。
InterProScan 6 はリリースごとに `databases.json` を同梱しており、両者を突き合わせられる。

- `versions.json`: `{"interpro": [105.0, 106.0, 107.0, 108.0, "109.0"]}` → `--interpro latest` は **109.0**
  （FTP には `interpro-110.0.tar.gz` も置かれているが `versions.json` には載っていない）
- **`data/interpro/109.0/databases.json` と Matches API のレスポンスは 18 ライブラリすべてで一致**
  （Pfam 38.2, PANTHER 19.0, CATH-Gene3D 4.3.0, NCBIFAM 19.0, PROSITE 2026_01, CDD 3.21,
  SUPERFAMILY 1.75, SMART 9.0, PIRSF 3.10, PIRSR 2025_05, PRINTS 42.0, SFLD 4,
  HAMAP 2026_01, COILS 2.2.1, MobiDB-lite 4.0, Phobius 1.01, CATH-FunFam 4.3.0,
  PROSITE patterns 2026_01）。不一致 0、未掲載 0。

再検証:

```sh
python3 scripts/08_ips6_versions.py
```

メンバ DB の版が実際に固定されることも確認した。`--interpro 109.0` を指定すると、
FTP に `ncbifam-20.0.tar.gz` があるにもかかわらず **19.0 が取得される**。

**したがって `--interpro 109.0` を指定すれば、ヒット分（Matches API）と未ヒット分
（ローカル実行）で注釈の基準が揃う。** `--interpro latest` は現時点でたまたま 109.0 だが、
110.0 が `versions.json` に載った時点でずれるため、版は明示的に固定すべきである。

## 5. 対象集合

### 5.1 そのままでは投入できない（実測）

`step0_out/unmatched.faa` をそのまま渡すと入力検証で落ちる:

```
[ERROR] ERROR ~ Invalid character(s) found in the input FASTA file.
```

不正文字は **`*` のみ**（配列中に現れる文字は `*` と 20 種の標準残基だけ）。
内部終止を持つ 235 配列は**全件が未ヒット側にある**ため、必ず衝突する。

### 5.2 採用した方式: `*` を `X` に置換

本番入力は `scripts/09_sample_unmatched.py --all --mask-internal-stop X` で
`*` を `X` に置換したものとする（3,514 配列 / 1,993,896 aa、置換 235 配列 / 816 文字、
残基数は不変）。理由は 3 点。

第一に、**内部終止を持つ 235 配列が 1 遺伝子なのか 2 遺伝子なのかは未確定**であり、
`*` の位置で配列を分割することは「2 遺伝子である」という判断を先に下すことになる。
本プロジェクトは相同性投影であって実験的裏づけを持たないので、その判断の根拠がない。
`X` 置換は 1 配列 1 エントリのまま扱い、この判断を保留する。
第二に、`X` 置換は eggNOG-mapper に投入した `UTEX25_proteins_for_eggnog.faa` と
同じ前処理であり、両者の比較可能性が保たれる（ただし当該ファイルは本リポジトリに
含まれていないため、同一処理であることは依頼者の申告によるもので、ここでは検証していない）。
第三に、除外方式と違って 235 配列が結果表から落ちない。
`X` 置換では終止をまたぐドメインが誤って 1 つに連結される可能性が残るが、
その影響は `step0_out/lookup_status.tsv` の `has_internal_stop` 列でいつでも層別できる
（カバレッジの層別は `docs/INTERPRO_COVERAGE.md` §5）。

置換後のファイルは `annotation/` には置かず作業用ディレクトリに置いた
（`step0_out/unmatched.faa` から決定的に再生成できるため）。

<!-- BEGIN unmatched-profile -->
### 対象集合の性質（未ヒット 3,514 配列）

生成: `scripts/10_unmatched_profile.py` / `step0_out/lookup_status.tsv` の `in_uniparc = no` 3514 件（全 7413 件中）を `annotation/UTEX25_gene_table.tsv` と結合。数値のみ。

#### QC フラグ別内訳

| QC | 件数 | 未ヒット内の構成比 | 同 QC の全遺伝子 | そのうち未ヒットの割合 |
|---|---|---|---|---|
| `pass` | 2864 | 81.50 % | 6762 | 42.35 % |
| `no_start_Met` | 390 | 11.10 % | 391 | 99.74 % |
| `internal_stop,frameshift` | 143 | 4.07 % | 143 | 100.00 % |
| `internal_stop` | 62 | 1.76 % | 62 | 100.00 % |
| `frameshift` | 24 | 0.68 % | 24 | 100.00 % |
| `internal_stop,no_start_Met` | 16 | 0.46 % | 16 | 100.00 % |
| `internal_stop,no_start_Met,frameshift` | 14 | 0.40 % | 14 | 100.00 % |
| `no_start_Met,frameshift` | 1 | 0.03 % | 1 | 100.00 % |
| **計** | **3514** | **100.00 %** | 7413 | 47.40 % |

#### identity 分布

| 集合 | n | 最小 | Q1 | 中央値 | Q3 | 最大 | 平均 |
|---|---|---|---|---|---|---|---|
| 未ヒット | 3514 | 0.4623 | 0.9825 | 0.9918 | 0.9959 | 1.0000 | 0.9810 |

#### タンパク質長の分布（アミノ酸残基数）

`step0_out/lookup_status.tsv` の `length` 列（末尾 `*` を除去した後の長さ、InterProScan に渡る長さ）。

| 集合 | n | 最小 | Q1 | 中央値 | Q3 | 最大 | 平均 | 合計 |
|---|---|---|---|---|---|---|---|---|
| 未ヒット | 3514 | 19 | 274.0 | 454.0 | 693.8 | 17905 | 567.4 | 1,993,896 |

#### 内部終止コドン

- 内部終止 `*` を含む配列: **235 / 3514 (6.69 %)**
- 全 7,413 配列中の内部終止保有数は 235 で、そのすべてが未ヒット側にある。

<!-- END unmatched-profile -->
## 6. 所要時間（実測と外挿）

### 6.1 test プロファイルは外挿に使えない

`-profile test` の入力 `assets/test.faa` は 103 配列 / 68,238 aa だが、
**102 配列が Matches API のルックアップで解決され、ローカル走査されたのは 1 配列
（780 aa）だけ**だった。この所要時間から 3,514 配列を外挿することはできない。

### 6.2 計測方法

`step0_out/unmatched.faa` から seed 固定で無作為抽出し、Matches API を無効化して
Pfam のみを走らせた。

```sh
python3 scripts/09_sample_unmatched.py --n 100 --seed 20260901 \
    --mask-internal-stop X --out sample_100.faa

nextflow run ./interproscan6 -profile docker -c ccr.config \
    --input sample_100.faa --datadir data --interpro 109.0 \
    --applications pfam --no-matches-api --cpus 4 --outprefix t100
```

抽出は `random.Random(seed).sample` をファイル順のインデックスに対して行うので、
同じ seed なら同じ集合になる（`--seed 20260901` で 2 回生成して一致を確認済み）。

計測はパイプライン全体の wall time（入力検証・分割・hmmsearch・パース・
出力書き出しを含む）。他のジョブが動いていない状態で実行した。
n=100 は 2 回繰り返して 133 s / 133 s（別途 136 s）と再現した。

### 6.3 実測値

| n | 残基数 | wall time | マスクした配列 |
|---|---|---|---|
| 100 | 50,828 aa | **133 s** | 3 |
| 400 | 197,388 aa | **207 s** | 17 |
| 1,000 | 546,418 aa | **381 s** | 48 |

いずれも exit 0、`[SUCCESS] completed=15 failed=0`。

### 6.4 外挿

3 点の最小二乗直線（残基数に対して）:

```
t(x) = 107.9 + 5.0007e-4 * x        （x = 残基数、t = 秒）
```

| n | 実測 | 直線の値 | 残差 |
|---|---|---|---|
| 100 | 133 s | 133.3 s | -0.3 s |
| 400 | 207 s | 206.6 s | +0.4 s |
| 1,000 | 381 s | 381.1 s | -0.1 s |

区間ごとの傾きは 5.049e-4 / 4.985e-4 / 5.004e-4 s/aa で、この範囲では直線とみなせる。

**3,514 配列 / 1,993,896 aa への外挿: 1,105 s ≒ 18.4 分**（Pfam のみ、4 vCPU）。

比較のため、n=100 の 1 点だけから比例で外挿した場合は次のようになる。
切片 107.9 s（固定費）を無視するため大きく過大になる。

| 外挿方法 | 結果 |
|---|---|
| 3 点の最小二乗（採用） | **1,105 s = 18.4 分** |
| n=100 から残基数比例 | 5,217 s = 87.0 分 |
| n=100 から配列数比例 | 4,674 s = 77.9 分 |

### 6.5 外挿の前提

以下を仮定している。いずれも 3,514 配列での実測ではない。

1. **線形性** — 50,828〜546,418 aa の範囲で直線が成り立つことは実測したが、
   1,993,896 aa まで同じ傾きが続くことは検証していない。
2. **バッチ分割の regime が変わらない** — `batchSize = 5000` / `subBatchSize = 500`
   のため 3,514 配列は 1 バッチに収まる。計測した 3 点も同じ 1 バッチ regime だが、
   サブバッチ数は増える。
3. **4 vCPU、`--cpus 4`** — hmmsearch は `--cpu 4` で起動される。
   他のジョブが同時に動けば伸びる。
4. **データはダウンロード済み** — 上の時間に取得時間は含まない（§3.2 の 373 s / 1119 s が別途要る）。
5. **Pfam のみ** — 他の applications を足した場合の時間は**未計測**。
6. **Matches API を使わない前提** — 未ヒット 3,514 配列は定義上 UniParc に無いので、
   本番でもルックアップは 1 件も当たらない。

## 7. セットアップスクリプト案

セッション中に取得したデータは次のセッションに残らない。同じ待ち時間を繰り返さないために、
環境のセットアップスクリプトで先にデータを用意しておく案を以下に示す。
**これは案であり、設定は未実施。**

`docs` の記述にしたがい、非 critical なコマンドには `|| true` を付けて、
5 分以内に終わらない場合でもセッション起動をブロックしないようにしてある。
データ取得はどうしても 5 分を超えるため（§3.2 の実測で 1119 s）、
**取得はバックグラウンドに逃がし、完了フラグで判定する**構成にしている。

```sh
#!/usr/bin/env bash
# InterProScan 6 のローカル実行準備。
# 5 分以内に終わらない処理はバックグラウンドに逃がし、セッション起動を止めない。
set -u

IPS6_HOME="${HOME}/ips6"
IPS6_REV="6.0.1"
IPS6_INTERPRO="109.0"          # Matches API と一致する版（docs/IPS6_FEASIBILITY.md §4）
IPS6_APPS="pfam,ncbifam,prositepatterns,prositeprofiles,smart,cdd,superfamily"

mkdir -p "$IPS6_HOME" || true

# 1. docker デーモン（この環境では自動起動しない）
if ! docker info >/dev/null 2>&1; then
    rm -f /var/run/docker.pid || true
    setsid nohup dockerd > "$IPS6_HOME/dockerd.log" 2>&1 < /dev/null &
    for _ in $(seq 1 30); do docker info >/dev/null 2>&1 && break; sleep 1; done
fi

# 2. nextflow（プリインストールされていない）
if [ ! -x "$IPS6_HOME/nextflow" ]; then
    curl -fsSL https://get.nextflow.io -o "$IPS6_HOME/nextflow" || true
    chmod +x "$IPS6_HOME/nextflow" || true
    NXF_HOME="$IPS6_HOME/.nextflow" "$IPS6_HOME/nextflow" -version >/dev/null 2>&1 || true
fi

# 3. パイプライン本体
#    `nextflow run ebi-pf-team/interproscan6` は GitHub API がプロキシに拒否されるため使えない。
#    clone してローカルパスに対して実行する。
if [ ! -d "$IPS6_HOME/interproscan6" ]; then
    git clone --depth 1 -b "$IPS6_REV" \
        https://github.com/ebi-pf-team/interproscan6.git \
        "$IPS6_HOME/interproscan6" || true
fi

# 4. コンテナにプロキシ CA を渡す設定（無いと全ダウンロードが curl exit 60 で落ちる）
cat > "$IPS6_HOME/ccr.config" <<'CFG' || true
docker {
    enabled    = true
    runOptions = '--network host -v /root/.ccr/ca-bundle.crt:/ccr-ca.crt:ro'
}
env {
    CURL_CA_BUNDLE     = '/ccr-ca.crt'
    SSL_CERT_FILE      = '/ccr-ca.crt'
    REQUESTS_CA_BUNDLE = '/ccr-ca.crt'
    HTTPS_PROXY        = 'http://127.0.0.1:41551'
    https_proxy        = 'http://127.0.0.1:41551'
}
CFG

# 5. データ取得。実測 1119 s かかるのでバックグラウンドで走らせ、
#    完了したら .ready を置く。セッションはここでブロックしない。
if [ ! -f "$IPS6_HOME/data/.ready" ] && [ ! -f "$IPS6_HOME/data.fetching" ]; then
    touch "$IPS6_HOME/data.fetching" || true
    (
        cd "$IPS6_HOME" &&
        NXF_HOME="$IPS6_HOME/.nextflow" "$IPS6_HOME/nextflow" run \
            "$IPS6_HOME/interproscan6" \
            -profile docker,test -c ccr.config \
            --datadir "$IPS6_HOME/data" \
            --interpro "$IPS6_INTERPRO" \
            --applications "$IPS6_APPS" \
            --outprefix warmup > "$IPS6_HOME/warmup.log" 2>&1 &&
        touch "$IPS6_HOME/data/.ready"
        rm -f "$IPS6_HOME/data.fetching"
    ) > /dev/null 2>&1 < /dev/null &
fi

exit 0
```

意図:

- **`-profile docker,test` で warm-up する** — `test` プロファイルの入力は 103 配列なので、
  データ一式を取得しつつ実行が通ることまで確認できる。所要は §3.2 の 1119 s。
- **`--interpro 109.0` を明示** — `latest` は 110.0 が `versions.json` に載った時点でずれる。
- **`--applications` を絞る** — 指定した 7 種のデータだけを取得する（8.7 GB）。
  PANTHER / CATH-Gene3D を足す場合は §3.3 のとおり容量が未計測なので、
  先に `du -sh data` で確認する必要がある。
- **`data/.ready`** — 本番実行の前にこのフラグを見れば、取得が完了しているか判定できる。

本番実行はこの後、別途（実際に実行した内容は §9）:

```sh
python3 scripts/09_sample_unmatched.py \
    --all --mask-internal-stop X --out ips6_input.faa   # §5.2

NXF_HOME=~/ips6/.nextflow ~/ips6/nextflow run ~/ips6/interproscan6 \
    -profile docker -c ~/ips6/ccr.config \
    --input ips6_input.faa --datadir ~/ips6/data \
    --interpro 109.0 --applications pfam --no-matches-api --cpus 4 \
    --outprefix utex25_unmatched
```

## 8. 未計測・未実施のまとめ

| 項目 | 状態 |
|---|---|
| PANTHER のデータ量 | **未計測**（帯域制約、254 MB / 1.0 G で中断） |
| CATH-Gene3D のデータ量 | **未計測**（未着手） |
| PANTHER + CATH-Gene3D で 30 GB に収まるか | **未計測** |
| Pfam 以外を含めた場合の所要時間 | **未計測** |
| 3,514 配列の本番実行 | **実施済み**（§9）。§6.4 の外挿は Pfam のみの値 |
| 内部終止 `*` の本番での扱い | **決定済み**: `X` 置換（§5.2） |

## 9. 本番実行の記録（実施済み）

`docs/IPS6_FEASIBILITY.md` §0〜§8 の判定にもとづいて実行した。

| 項目 | 実測 |
|---|---|
| 入力 | `step0_out/unmatched.faa` を `X` 置換したもの（3,514 配列 / 1,993,896 aa） |
| applications | `pfam,ncbifam,prositepatterns,prositeprofiles,smart,cdd,superfamily` |
| `--interpro` | `109.0`（`latest` は使っていない） |
| `--cpus` | 4 |
| Matches API ルックアップ | 有効のまま |
| 開始 | 2026-09-02T00:29:29Z |
| **所要** | **4,320 s = 72.0 分** |
| 終了状態 | **exit 0** / `[SUCCESS] completed=45 failed=0 cached=0` |
| 途中で落ちた地点 | なし |

データは前段で取得済みだったため、この 72 分にダウンロード時間は含まれない
（取得は §3.2 の 1119 s）。

出力（`ips6_out/`）:

| ファイル | 内容 |
|---|---|
| `utex25_unmatched.json.gz` | 生の JSON 出力（24 MB → gzip 3.2 MB） |
| `utex25_unmatched.tsv.gz` | 同 TSV |
| `input_prep.log` | 入力生成のログ（置換件数・残基数） |
| `run_provenance.txt` | 実行コマンド・開始時刻・所要・終了状態 |

出力の検証:

- `results` は **3,514 件**、xref id も 3,514 件でユニーク（入力と一致、取りこぼしなし）
- `interproscan-version: 6.0.1` / **`interpro-version: 109.0`**
- ライブラリ版は Pfam 38.2 / NCBIFAM 19.0 / PROSITE patterns 2026_01 /
  PROSITE profiles 2026_01 / SMART 9.0 / CDD 3.21 / SUPERFAMILY 1.75
- マッチが 1 つ以上付いた配列: 2,687 / 3,514

### 9.1 Pfam のみの外挿（§6.4）との関係

§6.4 の 18.4 分は **Pfam のみ**の外挿値で、本番は 7 ライブラリなので直接は比較できない。
両者を並べると次のとおり。7 ライブラリ構成での所要時間は事前に計測していなかったため、
事前の予測値は存在しない。

| 構成 | 所要 |
|---|---|
| Pfam のみ（§6.4、外挿） | 1,105 s = 18.4 分 |
| 7 ライブラリ（本番、実測） | 4,320 s = 72.0 分 |
