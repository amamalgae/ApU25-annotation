#!/usr/bin/env bash
# UTEX 250-A proteins -> UTEX 25 genome (miniprot 0.18-r281)
set -euo pipefail
MP=${MP:-./miniprot/miniprot}
GENOME=${GENOME:-U25.fa}
"$MP" -t1 -d U25.mpi "$GENOME"
for h in hap1 hap2; do
  "$MP" -t1 -K 1M --outn=1 --gff U25.mpi "$h.faa" > "$h.gff" 2> "$h.log"
  echo "$h: $(grep -cP '\tmRNA\t' "$h.gff") models"
done
