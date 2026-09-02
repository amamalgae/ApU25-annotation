# 環境メモ（この実行環境で踏んだ問題と対処）

Claude Code のリモート実行環境で作業したときに判明した、環境固有の事項の記録。
将来同じ問題を踏まないための覚え書きで、いずれも実測にもとづく。

## 1. コンテナにプロキシ CA を渡さないと全ダウンロードが落ちる

**症状** — docker コンテナの中から HTTPS を叩く処理が例外なく失敗する。
InterProScan 6 では `INIT_DATABASES:DOWNLOAD_INTERPRO` がリトライを使い切って落ちた。

```
curl: (60) SSL certificate OpenSSL verify result: self-signed certificate in certificate chain (19)
```

**原因** — このセッションの外向き HTTPS は `127.0.0.1:41551` のプロキシを経由し、
そこで TLS が再終端される。ホスト側には CA バンドル `/root/.ccr/ca-bundle.crt` が
配置され各種の環境変数も設定済みだが、**コンテナの中にはどちらも入らない**。
さらにプロキシは `127.0.0.1` にあるので、bridge ネットワークのコンテナからは届かない。

**対処** — CA をマウントし、ホストネットワークを使い、CA の場所を環境変数で教える。
nextflow なら次の内容を `-c` で足す:

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

`docker run` を直接使う場合も同じ 3 点（`--network host`、CA のマウント、
`CURL_CA_BUNDLE` などの環境変数）が要る。

**やってはいけないこと** — TLS 検証の無効化と `HTTPS_PROXY` の unset。
`/root/.ccr/README.md` に明記されている。

なお **docker デーモンは自動起動していない**（`/var/run/docker.sock` が無い）。
`dockerd` を手動で起動する必要がある。cgroup v1 の非推奨警告は出るが動作する。

```sh
rm -f /var/run/docker.pid
setsid nohup dockerd > dockerd.log 2>&1 < /dev/null &
for _ in $(seq 1 30); do docker info >/dev/null 2>&1 && break; sleep 1; done
```

## 2. `nextflow run <github org>/<repo>` は GitHub API 拒否で失敗する

**症状**

```
$ nextflow run ebi-pf-team/interproscan6 -r 6.0.1 -profile docker,test ...
WARN: Cannot read project manifest -- Cause: Forbidden
Forbidden -- Provide your GitHub user name and password to access this repository
```

**原因** — nextflow はパイプラインを取得する前に GitHub API でマニフェストを読む。
その API 呼び出しがセッションの egress プロキシに拒否される。
リポジトリが公開かどうかとは関係なく、認証情報を足しても解決しない。

**対処** — `git clone` は通るので、clone してローカルパスに対して実行する。

```sh
git clone --depth 1 -b 6.0.1 https://github.com/ebi-pf-team/interproscan6.git
nextflow run ./interproscan6 -profile docker -c ccr.config ...
```

`-r <revision>` の代わりに `git clone -b <tag>` でリビジョンを固定する。

## 3. プロキシ帯域は約 0.7 MB/s（実測）

大きなデータを落とす作業では、この帯域が支配的な律速になる。実測例:

| 対象 | サイズ | 実測 |
|---|---|---|
| InterProScan 6 データ（Pfam のみ） | 391M | 373 s（展開・スキャン込み） |
| 同（7 ライブラリ） | 約 2.3G | 1119 s（同上） |
| PANTHER | 1.0G | 途中で中止（254 MB 時点） |

**作業計画への影響** — GB 級のダウンロードを伴う計測は、開始前に所要時間を見積もること。
数 GB あれば時間単位でかかる。実際、PANTHER と CATH-Gene3D（計 3.5G）の計測は
費用対効果が見合わないと判断して中止した（`docs/IPS6_FEASIBILITY.md` §3.3）。

## 4. セッションをまたいで残るもの・残らないもの

- **セットアップスクリプトが書いたものは環境キャッシュに保持される。**
  繰り返し必要になる重いデータは、セットアップスクリプトで用意しておくのが確実。
- **セッション中に落としたデータは次のセッションに残らない可能性が高い。**
  今回は同一セッション内で `data/` 8.7 GB が保持されたが、これは前提にできない。

大きなデータ取得を伴う作業では、セットアップスクリプトに逃がす案を先に検討する。
具体案は `docs/IPS6_FEASIBILITY.md` §7 に記載した。

セットアップスクリプトは 5 分以内に終わらないとセッション起動を阻害するため、
長時間かかる処理はバックグラウンドに逃がし、完了フラグで判定する。
非 critical なコマンドには `|| true` を付ける。

## 5. InterProScan 6 は入力 FASTA の `*` を不正文字として拒否する

**症状**

```
[ERROR] ERROR ~ Invalid character(s) found in the input FASTA file.
```

**原因** — 内部終止コドンを表す `*` が 1 文字でもあると、ファイル全体が弾かれる。
`step0_out/unmatched.faa` では 235 配列 / 計 816 文字が該当した。

**対処** — 投入前に置換するか、該当配列を除外する。本プロジェクトでは
`scripts/09_sample_unmatched.py --mask-internal-stop X` で `X` に置換した。
残基数が変わらないので所要時間の外挿や座標が壊れない。
どう扱うかはデータの解釈に関わる判断なので、採用理由は
`docs/IPS6_FEASIBILITY.md` §5.2 に記載してある。

## 6. ディスク

書き込み可能な容量はセッションごとの割り当てで、`df` の表示は実態と異なる。
`Avail` が 0 で `Used` が小さい場合は「割り当てを使い切った」であって故障ではない。
今回の割り当ては約 30 GB だった。

`No space left on device` が出たら、不要になった大きいファイル
（ビルド生成物、キャッシュ、古い clone）を消す。書き込みが失敗する状態でも削除は通り、
解放した分はすぐ書けるようになる。
