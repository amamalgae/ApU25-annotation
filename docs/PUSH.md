# リポジトリ作成と push

```bash
cd apu25-annotation
git init -b main && git add . && git commit -m "UTEX 25 annotation pipeline"
gh repo create apu25-annotation --private --source=. --push

# 生成物は git に入れず Release に添付
tar czf apu25-annotation-v1.0.0.tar.gz -C /path/to UTEX25_annotation
gh release create v1.0.0 apu25-annotation-v1.0.0.tar.gz \
  --title "UTEX 25 projected annotation v1.0.0" \
  --notes "7,413 gene models projected from UTEX 250-A onto UTEX 25 with miniprot 0.18-r281."
```

## サイズ

| | |
|---|---|
| GenBank 12ファイル（配列込み） | 40.6 MB |
| GFF3 + faa + fna + TSV | 21.9 MB |
| 全部 tar.gz | 約 12 MB |
| 最大単一ファイル | 5.4 MB (chr12.gb) |

GitHub の単一ファイル上限 100 MB（50 MB で警告）、Release アセット 2 GB/ファイル。
Git LFS は無料枠がストレージ 1 GB・帯域 1 GB/月と狭いので使わない。
