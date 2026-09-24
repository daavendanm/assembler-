#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<EOF
Usage: $(basename "$0") -1 R1.fastq.gz -2 R2.fastq.gz -g GENOME_SIZE_BP -o OUTPUT_DIR [-x TARGET_COVERAGE] [-t THREADS] [-s SEED]

De novo assembly of paired-end Illumina reads. Reads are first downsampled
(with seqtk) to a target coverage, since raw coverage far above what an
assembler needs mostly adds runtime/memory without improving quality, then
assembled with Unicycler (falls back to SPAdes --isolate if Unicycler is
not installed). Both are short-read bacterial-genome assemblers; no long
reads are used here.

Required:
  -1  Forward reads, FASTQ(.gz) (clean/trimmed reads)
  -2  Reverse reads, FASTQ(.gz)
  -g  Approximate genome size in bp (e.g. 4376788, from the earlier
      reference-based consensus length) -- used only to compute the
      downsampling fraction
  -o  Output directory

Optional:
  -x  Target coverage after downsampling (default: 100)
  -t  Threads (default: 4)
  -s  seqtk random seed, must be identical for R1/R2 to keep pairing
      (default: 100)

Requires 'seqtk' and either 'unicycler' or 'spades.py' on PATH.
EOF
    exit 1
}

TARGET_COV=100
THREADS=4
SEED=100

while getopts "1:2:g:o:x:t:s:h" opt; do
    case "$opt" in
        1) R1="$OPTARG" ;;
        2) R2="$OPTARG" ;;
        g) GENOME_SIZE="$OPTARG" ;;
        o) OUTDIR="$OPTARG" ;;
        x) TARGET_COV="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        s) SEED="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done

if [[ -z "${R1:-}" || -z "${R2:-}" || -z "${GENOME_SIZE:-}" || -z "${OUTDIR:-}" ]]; then
    usage
fi

for f in "$R1" "$R2"; do
    [[ -f "$f" ]] || { echo "Error: file not found: $f" >&2; exit 1; }
done

command -v seqtk >/dev/null 2>&1 || { echo "Error: seqtk not found in PATH" >&2; exit 1; }

ASSEMBLER=""
if command -v unicycler >/dev/null 2>&1; then
    ASSEMBLER="unicycler"
elif command -v spades.py >/dev/null 2>&1; then
    ASSEMBLER="spades"
else
    echo "Error: neither unicycler nor spades.py found in PATH" >&2
    exit 1
fi

mkdir -p "$OUTDIR"

count_bases() {
    local f="$1"
    if [[ "$f" == *.gz ]]; then
        zcat "$f" | awk 'NR%4==2 { total += length($0) } END { print total+0 }'
    else
        awk 'NR%4==2 { total += length($0) } END { print total+0 }' "$f"
    fi
}

echo "[1/3] Counting input bases..." >&2
BASES_R1=$(count_bases "$R1")
BASES_R2=$(count_bases "$R2")
TOTAL_BASES=$((BASES_R1 + BASES_R2))
CURRENT_COV=$(awk -v b="$TOTAL_BASES" -v g="$GENOME_SIZE" 'BEGIN { printf "%.2f", b/g }')
echo "  Total input bases: $TOTAL_BASES (R1=$BASES_R1, R2=$BASES_R2)" >&2
echo "  Estimated current coverage: ${CURRENT_COV}x (genome size ${GENOME_SIZE} bp)" >&2

FRACTION=$(awk -v cur="$CURRENT_COV" -v tgt="$TARGET_COV" 'BEGIN { f = tgt/cur; if (f > 1) f = 1; printf "%.6f", f }')

DS_R1="$OUTDIR/reads_ds_1.fastq.gz"
DS_R2="$OUTDIR/reads_ds_2.fastq.gz"

if awk -v f="$FRACTION" 'BEGIN { exit !(f >= 1.0) }'; then
    echo "[2/3] Current coverage already at or below target (${TARGET_COV}x); skipping downsampling." >&2
    cp "$R1" "$DS_R1"
    cp "$R2" "$DS_R2"
else
    echo "[2/3] Downsampling to ~${TARGET_COV}x (fraction=$FRACTION, seed=$SEED)..." >&2
    seqtk sample -s"$SEED" "$R1" "$FRACTION" | gzip > "$DS_R1"
    seqtk sample -s"$SEED" "$R2" "$FRACTION" | gzip > "$DS_R2"
fi

ASM_OUT="$OUTDIR/assembly"

echo "[3/3] Running $ASSEMBLER..." >&2
if [[ "$ASSEMBLER" == "unicycler" ]]; then
    unicycler -1 "$DS_R1" -2 "$DS_R2" -o "$ASM_OUT" -t "$THREADS"
    FINAL_FASTA="$ASM_OUT/assembly.fasta"
else
    spades.py --isolate -1 "$DS_R1" -2 "$DS_R2" -o "$ASM_OUT" -t "$THREADS"
    FINAL_FASTA="$ASM_OUT/contigs.fasta"
fi

echo "" >&2
if [[ -f "$FINAL_FASTA" ]]; then
    N_CONTIGS=$(grep -c '^>' "$FINAL_FASTA")
    TOTAL_LEN=$(awk '/^>/ { next } { total += length($0) } END { print total+0 }' "$FINAL_FASTA")
    echo "Done. Assembly: $FINAL_FASTA" >&2
    echo "Contigs: $N_CONTIGS, total length: $TOTAL_LEN bp" >&2
else
    echo "Done. Check $ASM_OUT for assembler output (expected FASTA not found at $FINAL_FASTA)." >&2
fi
