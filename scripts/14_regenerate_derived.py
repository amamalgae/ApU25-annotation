#!/usr/bin/env python3
"""Push the finalised symbol / product from the gene table into the derived files.

Touches only the annotation fields.  Coordinates, feature structure and — the
point of the exercise — every sequence byte are left exactly as they were, and
the script verifies that afterwards rather than trusting it.

  GenBank x12   /gene on gene and CDS, /product on mRNA and CDS, and an extra
                /note on CDS for the unresolved-family and 3+-domain cases
  GFF3          Name= on gene, product= on mRNA
  proteins.faa  header line
  CDS.fna       header line

Qualifier wrapping follows scripts/04_write_genbank.py: textwrap at width 58
with a 21-space indent, so re-written lines are indistinguishable from the
originals.

    python3 scripts/14_regenerate_derived.py
"""

from __future__ import annotations

import argparse
import hashlib
import re
import textwrap
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
ANN = ROOT / "annotation"
DEFAULT_TABLE = ANN / "UTEX25_gene_table.tsv"

QUAL_START = re.compile(r"^ {21}/(\w+)=")
FEATURE_START = re.compile(r"^ {5}(\S+)\s+\S")


def qual(key, value, quote=True):
    text = f"/{key}=" + (f'"{value}"' if quote else value)
    return textwrap.wrap(text, 58, initial_indent=" " * 21,
                         subsequent_indent=" " * 21, break_long_words=True,
                         break_on_hyphens=False)


def read_tsv(path: Path):
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(header, line.rstrip("\n").split("\t")))
                for line in fh if line.strip()]


def extra_notes(gene):
    """The /note lines the rules ask for, or None."""
    if gene["interpro_family_unresolved"] == "TRUE":
        return ("unresolved InterPro family assignment: "
                + gene["interpro_family_candidates"])
    if gene["product_source"] == "interpro_domain_multi":
        parts = [p.strip() for p in gene["interpro_domain_all"].split(";") if p.strip()]
        if len(parts) >= 3:
            return "additional InterPro domains: " + "; ".join(parts[2:])
    return None


# --------------------------------------------------------------------------- #
# GenBank
# --------------------------------------------------------------------------- #

def sequence_fingerprint_gb(lines):
    """MD5 over everything from ORIGIN to the end — the sequence body."""
    md5 = hashlib.md5()
    inside = False
    for line in lines:
        if line.startswith("ORIGIN"):
            inside = True
        if inside:
            md5.update(line.encode())
    return md5.hexdigest()


def rewrite_genbank(path: Path, genes: dict, stats):
    lines = path.read_text().splitlines(keepends=True)
    before = sequence_fingerprint_gb(lines)

    out = []
    i = 0
    feature = None
    locus = None
    while i < len(lines):
        line = lines[i]
        if line.startswith("ORIGIN"):
            out.extend(lines[i:])
            break
        m = FEATURE_START.match(line)
        if m:
            feature = m.group(1)
            if feature not in ("gene", "mRNA", "CDS"):
                locus = None
            out.append(line)
            i += 1
            continue
        q = QUAL_START.match(line)
        if not q:
            out.append(line)
            i += 1
            continue

        key = q.group(1)
        block = [line]
        j = i + 1
        while (j < len(lines) and not QUAL_START.match(lines[j])
               and not FEATURE_START.match(lines[j])
               and not lines[j].startswith("ORIGIN")):
            block.append(lines[j])
            j += 1

        if key == "locus_tag":
            locus = re.search(r'/locus_tag="([^"]+)"', block[0]).group(1)
            out.extend(block)
            # a gene that has newly acquired a symbol needs /gene inserting here
            gene = genes.get(locus)
            if gene and gene["symbol"] and feature in ("gene", "CDS"):
                following = lines[j] if j < len(lines) else ""
                if not following.startswith(" " * 21 + "/gene="):
                    out.extend(l + "\n" for l in qual("gene", gene["symbol"]))
                    stats["gene_inserted"] += 1
            i = j
            continue

        gene = genes.get(locus) if locus else None
        if gene and key == "gene" and feature in ("gene", "CDS"):
            out.extend(l + "\n" for l in qual("gene", gene["symbol"]))
            stats["gene_rewritten"] += 1
            i = j
            continue
        if gene and key == "product" and feature in ("mRNA", "CDS"):
            out.extend(l + "\n" for l in qual("product", gene["product"]))
            stats["product_rewritten"] += 1
            i = j
            continue
        if gene and key == "note" and feature == "CDS":
            out.extend(block)
            note = extra_notes(gene)
            already = "".join(block)
            if note and note.split(":")[0] in already:
                note = None      # already added by an earlier run
            if note:
                out.extend(l + "\n" for l in qual("note", note))
                stats["note_added"] += 1
            i = j
            continue

        out.extend(block)
        i = j

    path.write_text("".join(out))
    after = sequence_fingerprint_gb(path.read_text().splitlines(keepends=True))
    if before != after:
        raise SystemExit(f"{path}: the ORIGIN block changed; aborting")
    return before


# --------------------------------------------------------------------------- #
# GFF3 / FASTA
# --------------------------------------------------------------------------- #

def rewrite_gff3(path: Path, genes: dict, stats):
    out = []
    for line in path.read_text().splitlines(keepends=True):
        if line.startswith("#") or "\t" not in line:
            out.append(line)
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) != 9:
            out.append(line)
            continue
        if f[2] == "gene":
            locus = re.search(r"ID=gene-([^;]+)", f[8])
            gene = genes.get(locus.group(1)) if locus else None
            if gene:
                name = gene["symbol"] or gene["locus_tag"]
                # a lambda replacement: symbol values contain backslashes
                f[8] = re.sub(r"Name=[^;]*", lambda _m: "Name=" + name, f[8])
                stats["gff_name"] += 1
        elif f[2] == "mRNA":
            locus = re.search(r"ID=rna-([^;]+)", f[8])
            gene = genes.get(locus.group(1)) if locus else None
            if gene:
                f[8] = re.sub(r"product=[^;]*",
                              lambda _m, g=gene: "product=" + g["product"], f[8])
                stats["gff_product"] += 1
        out.append("\t".join(f) + "\n")
    path.write_text("".join(out))


def rewrite_fasta(path: Path, genes: dict, stats, key):
    out, sequences_before, current, buf = [], {}, None, []
    for line in path.read_text().splitlines(keepends=True):
        if line.startswith(">"):
            if current is not None:
                sequences_before[current] = "".join(buf)
            fields = line[1:].rstrip("\n").split(None, 1)
            locus = fields[0]
            current, buf = locus, []
            gene = genes.get(locus)
            if gene:
                rest = fields[1] if len(fields) > 1 else ""
                coords = rest.split(" ", 2)[1] if len(rest.split(" ", 2)) > 1 else ""
                symbol = gene["symbol"] or "-"
                out.append(f">{locus} {symbol} {coords} {gene['product']}\n")
                stats[key] += 1
            else:
                out.append(line)
        else:
            buf.append(line)
            out.append(line)
    if current is not None:
        sequences_before[current] = "".join(buf)
    path.write_text("".join(out))
    return sequences_before


def fasta_sequences(path: Path):
    out, current, buf = {}, None, []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if current is not None:
                out[current] = "".join(buf)
            current, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if current is not None:
        out[current] = "".join(buf)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gene-table", type=Path, default=DEFAULT_TABLE)
    ap.add_argument("--annotation-dir", type=Path, default=ANN)
    args = ap.parse_args()

    genes = {g["locus_tag"]: g for g in read_tsv(args.gene_table)}
    stats = {k: 0 for k in ("gene_rewritten", "gene_inserted", "product_rewritten",
                            "note_added", "gff_name", "gff_product",
                            "faa_header", "fna_header")}

    faa = args.annotation_dir / "UTEX25_proteins.faa"
    fna = args.annotation_dir / "UTEX25_cds.fna"
    faa_before, fna_before = fasta_sequences(faa), fasta_sequences(fna)

    gb_files = sorted(args.annotation_dir.glob("*.gb"))
    if len(gb_files) != 12:
        raise SystemExit(f"expected 12 GenBank files, found {len(gb_files)}")
    for path in gb_files:
        rewrite_genbank(path, genes, stats)
    rewrite_gff3(args.annotation_dir / "UTEX25_annotation.gff3", genes, stats)
    rewrite_fasta(faa, genes, stats, "faa_header")
    rewrite_fasta(fna, genes, stats, "fna_header")

    for label, path, before in (("proteins.faa", faa, faa_before),
                                ("CDS.fna", fna, fna_before)):
        after = fasta_sequences(path)
        if before != after:
            differing = [k for k in before if before[k] != after.get(k)]
            raise SystemExit(f"{label}: {len(differing)} sequence(s) changed "
                             f"(first: {differing[:3]})")
        print(f"{label}: {len(after)} sequences, all identical to before")

    for key, value in stats.items():
        print(f"{key:<20} {value}")


if __name__ == "__main__":
    main()
