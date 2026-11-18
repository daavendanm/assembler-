#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<EOF
Usage: $(basename "$0") -1 R1.fastq.gz -2 R2.fastq.gz -g GENOME_SIZE_BP -o OUTPUT_DIR [-x TARGET_COVERAGE] [-t THREADS] [-s SEED] [-r REFERENCE.fasta] [-k KRAKEN2_DB] [-l BUSCO_LINEAGE] [-T]

De novo assembly of paired-end Illumina reads, with read trimming and
assembly QC. Reads are adapter/quality trimmed with fastp, downsampled
(with seqtk) to a target coverage, since raw coverage far above what an
assembler needs mostly adds runtime/memory without improving quality, then
assembled with Unicycler (falls back to SPAdes --isolate if Unicycler is
not installed). Both are short-read bacterial-genome assemblers; no long
reads are used here.

The assembly is then checked with QUAST (contiguity, optional reference
comparison), BUSCO (completeness and duplication) and, if a Kraken2
database is given, Kraken2 on both reads and contigs (contamination).
A per-contig length/GC/depth table is always written: contigs whose GC or
depth is far from the rest are the first place to look for contamination
or plasmids. Results are collected in OUTPUT_DIR/assembly_qc_summary.tsv.

This assembly, not the reference-based consensus, is the sequence to use
for ANI, dDDH and Mashtree.

Required:
  -1  Forward reads, FASTQ(.gz) (raw; trimmed here unless -T is given)
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
  -r  Reference FASTA for QUAST comparison (e.g. the top ANI hit). Optional;
      without it QUAST reports contiguity only.
  -k  Kraken2 database directory. Without it the contamination screen is
      skipped and the summary says so -- do not report the assembly as
      contamination-checked in that case.
  -l  BUSCO lineage (default: bacillales_odb10)
  -T  Reads are already trimmed; skip fastp

Requires 'seqtk', 'fastp' (unless -T), 'quast.py', 'busco', either
'unicycler' or 'spades.py', and 'kraken2' if -k is given. All are in
environment_denovo.yml.
EOF
    exit 1
}

TARGET_COV=100
THREADS=4
SEED=100
REFERENCE=""
KRAKEN_DB=""
BUSCO_LINEAGE="bacillales_odb10"
SKIP_TRIM=0

while getopts "1:2:g:o:x:t:s:r:k:l:Th" opt; do
    case "$opt" in
        1) R1="$OPTARG" ;;
        2) R2="$OPTARG" ;;
        g) GENOME_SIZE="$OPTARG" ;;
        o) OUTDIR="$OPTARG" ;;
        x) TARGET_COV="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        s) SEED="$OPTARG" ;;
        r) REFERENCE="$OPTARG" ;;
        k) KRAKEN_DB="$OPTARG" ;;
        l) BUSCO_LINEAGE="$OPTARG" ;;
        T) SKIP_TRIM=1 ;;
        h) usage ;;
        *) usage ;;
    esac
done

if [[ -z "${R1:-}" || -z "${R2:-}" || -z "${GENOME_SIZE:-}" || -z "${OUTDIR:-}" ]]; then
    usage
fi

for f in "$R1" "$R2" ${REFERENCE:+"$REFERENCE"}; do
    [[ -f "$f" ]] || { echo "Error: file not found: $f" >&2; exit 1; }
done
if [[ -n "$KRAKEN_DB" && ! -d "$KRAKEN_DB" ]]; then
    echo "Error: Kraken2 database directory not found: $KRAKEN_DB" >&2
    exit 1
fi

REQUIRED_TOOLS=(seqtk quast.py busco)
[[ "$SKIP_TRIM" -eq 1 ]] || REQUIRED_TOOLS+=(fastp)
[[ -z "$KRAKEN_DB" ]] || REQUIRED_TOOLS+=(kraken2)
for tool in "${REQUIRED_TOOLS[@]}"; do
    command -v "$tool" >/dev/null 2>&1 || { echo "Error: $tool not found in PATH" >&2; exit 1; }
done

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

if [[ "$SKIP_TRIM" -eq 1 ]]; then
    echo "[1/5] -T given; using reads as supplied (no trimming)." >&2
    IN_R1="$R1"
    IN_R2="$R2"
else
    echo "[1/5] Trimming adapters and low-quality tails with fastp..." >&2
    IN_R1="$OUTDIR/reads_trimmed_1.fastq.gz"
    IN_R2="$OUTDIR/reads_trimmed_2.fastq.gz"
    fastp -i "$R1" -I "$R2" -o "$IN_R1" -O "$IN_R2" \
        --detect_adapter_for_pe --cut_tail --cut_mean_quality 20 \
        --length_required 50 -w "$THREADS" \
        -j "$OUTDIR/fastp.json" -h "$OUTDIR/fastp.html"
fi

echo "[2/5] Counting input bases..." >&2
BASES_R1=$(count_bases "$IN_R1")
BASES_R2=$(count_bases "$IN_R2")
TOTAL_BASES=$((BASES_R1 + BASES_R2))
CURRENT_COV=$(awk -v b="$TOTAL_BASES" -v g="$GENOME_SIZE" 'BEGIN { printf "%.2f", b/g }')
echo "  Total input bases: $TOTAL_BASES (R1=$BASES_R1, R2=$BASES_R2)" >&2
echo "  Estimated current coverage: ${CURRENT_COV}x (genome size ${GENOME_SIZE} bp)" >&2

FRACTION=$(awk -v cur="$CURRENT_COV" -v tgt="$TARGET_COV" 'BEGIN { f = tgt/cur; if (f > 1) f = 1; printf "%.6f", f }')

DS_R1="$OUTDIR/reads_ds_1.fastq.gz"
DS_R2="$OUTDIR/reads_ds_2.fastq.gz"

if awk -v f="$FRACTION" 'BEGIN { exit !(f >= 1.0) }'; then
    echo "  Current coverage already at or below target (${TARGET_COV}x); skipping downsampling." >&2
    cp "$IN_R1" "$DS_R1"
    cp "$IN_R2" "$DS_R2"
else
    echo "  Downsampling to ~${TARGET_COV}x (fraction=$FRACTION, seed=$SEED)..." >&2
    seqtk sample -s"$SEED" "$IN_R1" "$FRACTION" | gzip > "$DS_R1"
    seqtk sample -s"$SEED" "$IN_R2" "$FRACTION" | gzip > "$DS_R2"
fi

ASM_OUT="$OUTDIR/assembly"

echo "[3/5] Running $ASSEMBLER..." >&2
if [[ "$ASSEMBLER" == "unicycler" ]]; then
    unicycler -1 "$DS_R1" -2 "$DS_R2" -o "$ASM_OUT" -t "$THREADS"
    FINAL_FASTA="$ASM_OUT/assembly.fasta"
else
    spades.py --isolate -1 "$DS_R1" -2 "$DS_R2" -o "$ASM_OUT" -t "$THREADS"
    FINAL_FASTA="$ASM_OUT/contigs.fasta"
fi

if [[ ! -s "$FINAL_FASTA" ]]; then
    echo "Error: expected assembly not found at $FINAL_FASTA; check $ASM_OUT" >&2
    exit 1
fi

echo "[4/5] Per-contig length, GC and depth..." >&2
# Depth comes from the assembler's headers: SPAdes 'cov_N' (k-mer coverage)
# or Unicycler 'depth=Nx' (relative to the median contig, so 1.0 = typical).
CONTIG_TABLE="$OUTDIR/contigs_length_gc_depth.tsv"
awk 'function out() { printf "%s\t%d\t%.2f\t%s\n", name, len, len ? 100 * gc / len : 0, d }
    BEGIN { print "contig\tlength\tGC_pct\tdepth" }
    /^>/ {
        if (name != "") out()
        name = substr($1, 2); d = "NA"; len = 0; gc = 0
        if (match($0, /cov_[0-9.]+/)) d = substr($0, RSTART + 4, RLENGTH - 4)
        else if (match($0, /depth=[0-9.]+/)) d = substr($0, RSTART + 6, RLENGTH - 6)
        next
    }
    { len += length($0); gc += gsub(/[GgCc]/, "") }
    END { if (name != "") out() }' "$FINAL_FASTA" > "$CONTIG_TABLE"

echo "[5/5] Assembly QC (QUAST, BUSCO${KRAKEN_DB:+, Kraken2})..." >&2
quast.py "$FINAL_FASTA" -o "$OUTDIR/quast" -t "$THREADS" ${REFERENCE:+-r "$REFERENCE"}

busco -i "$FINAL_FASTA" -m genome -l "$BUSCO_LINEAGE" -c "$THREADS" \
    -o busco --out_path "$OUTDIR" -f
BUSCO_SUMMARY=$(find "$OUTDIR/busco" -maxdepth 1 -name 'short_summary*.txt' | sort | head -n 1)

top_species() {
    # Top 5 species-level assignments from a Kraken2 report: "pct% name"
    awk -F'\t' '$4 == "S" { gsub(/^ +/, "", $6); printf "%s%% %s\n", $1 + 0, $6 }' "$1" \
        | sort -rn | head -n 5 | paste -sd ';' -
}

if [[ -n "$KRAKEN_DB" ]]; then
    kraken2 --db "$KRAKEN_DB" --threads "$THREADS" --paired --gzip-compressed \
        --report "$OUTDIR/kraken2_reads.report" --output /dev/null "$DS_R1" "$DS_R2"
    kraken2 --db "$KRAKEN_DB" --threads "$THREADS" \
        --report "$OUTDIR/kraken2_contigs.report" \
        --output "$OUTDIR/kraken2_contigs.tsv" "$FINAL_FASTA"
fi

SUMMARY="$OUTDIR/assembly_qc_summary.tsv"
{
    printf 'assembly\t%s\n' "$FINAL_FASTA"
    printf 'assembler\t%s\n' "$ASSEMBLER"
    printf 'reads_trimmed_with_fastp\t%s\n' "$([[ "$SKIP_TRIM" -eq 1 ]] && echo no || echo yes)"
    awk -F'\t' 'NR > 1 { n++; t += $2; gc += $2 * $3 / 100 }
        END {
            printf "contigs\t%d\n", n
            printf "total_length_bp\t%d\n", t
            printf "GC_pct\t%.2f\n", t ? 100 * gc / t : 0
        }' "$CONTIG_TABLE"
    tail -n +2 "$CONTIG_TABLE" | cut -f2 | sort -rn | awk '
        { l[NR] = $1; t += $1 }
        END {
            for (i = 1; i <= NR; i++) { c += l[i]; if (c >= t / 2) { printf "N50_bp\t%d\n", l[i]; break } }
            printf "largest_contig_bp\t%d\n", l[1]
        }'
    printf 'quast_report\t%s\n' "$OUTDIR/quast/report.tsv"
    if [[ -n "$BUSCO_SUMMARY" ]]; then
        printf 'busco_%s\t%s\n' "$BUSCO_LINEAGE" "$(grep -m1 -o 'C:[^ ]*' "$BUSCO_SUMMARY")"
    else
        printf 'busco_%s\tNOT FOUND\n' "$BUSCO_LINEAGE"
    fi
    if [[ -n "$KRAKEN_DB" ]]; then
        printf 'kraken2_reads_top_species\t%s\n' "$(top_species "$OUTDIR/kraken2_reads.report")"
        printf 'kraken2_contigs_top_species\t%s\n' "$(top_species "$OUTDIR/kraken2_contigs.report")"
    else
        printf 'contamination_screen\tNOT RUN (no -k Kraken2 database)\n'
    fi
} > "$SUMMARY"

echo "" >&2
echo "Done. Assembly: $FINAL_FASTA" >&2
echo "  Per-contig table: $CONTIG_TABLE" >&2
echo "  QUAST:            $OUTDIR/quast/report.html" >&2
echo "  BUSCO:            ${BUSCO_SUMMARY:-not found}" >&2
[[ "$SKIP_TRIM" -eq 1 ]] || echo "  fastp report:     $OUTDIR/fastp.html" >&2
echo "  QC summary:       $SUMMARY" >&2
echo >&2
sed 's/^/    /' "$SUMMARY" >&2
