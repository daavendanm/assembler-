#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<EOF
Usage: $(basename "$0") -r REFERENCE.fasta -1 R1.fastq.gz -2 R2.fastq.gz -o OUTPUT_PREFIX [-t THREADS] [-q MIN_QUAL] [-d MIN_DEPTH] [-k MASK_DEPTH]

Reference-based mapping, variant calling, and consensus generation for
paired-end Illumina reads (bwa-mem + bcftools). No de novo assembly.

Variants are called with a haploid model (--ploidy 1). Reference positions
with read depth below MASK_DEPTH are written as N in the consensus instead
of silently inheriting the reference base, so the consensus contains only
sequence the reads actually support. Breadth of coverage and the number of
masked bases are written to the summary file.

Required:
  -r  Reference genome, FASTA (e.g. closest public Bacillus licheniformis genome)
  -1  Forward reads, FASTQ(.gz)
  -2  Reverse reads, FASTQ(.gz)
  -o  Output prefix (path/basename for all generated files)

Optional:
  -t  Threads (default: 4)
  -q  Minimum variant QUAL to keep (default: 20) -- review against your own
      coverage before trusting this default
  -d  Minimum read depth (DP) to keep a variant (default: 10) -- same caveat
  -k  Positions with depth below this are masked with N in the consensus
      (default: 1, i.e. only zero-coverage bases). Deliberately separate
      from -d: masking removes sequence, so it should fire only where there
      is genuinely no read support.
EOF
    exit 1
}

THREADS=4
MIN_QUAL=20
MIN_DEPTH=10
MASK_DEPTH=1

while getopts "r:1:2:o:t:q:d:k:h" opt; do
    case "$opt" in
        r) REF="$OPTARG" ;;
        1) R1="$OPTARG" ;;
        2) R2="$OPTARG" ;;
        o) OUT="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        q) MIN_QUAL="$OPTARG" ;;
        d) MIN_DEPTH="$OPTARG" ;;
        k) MASK_DEPTH="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done

if [[ -z "${REF:-}" || -z "${R1:-}" || -z "${R2:-}" || -z "${OUT:-}" ]]; then
    usage
fi

for tool in bwa samtools bcftools; do
    command -v "$tool" >/dev/null 2>&1 || { echo "Error: $tool not found in PATH" >&2; exit 1; }
done

for f in "$REF" "$R1" "$R2"; do
    [[ -f "$f" ]] || { echo "Error: file not found: $f" >&2; exit 1; }
done

outdir="$(dirname "$OUT")"
[[ "$outdir" == "." ]] || mkdir -p "$outdir"

if [[ ! -f "${REF}.bwt" ]]; then
    echo "[1/6] Indexing reference..." >&2
    bwa index "$REF"
fi
[[ -f "${REF}.fai" ]] || samtools faidx "$REF"

echo "[2/6] Mapping reads with bwa-mem..." >&2
bwa mem -t "$THREADS" "$REF" "$R1" "$R2" \
    | samtools sort -@ "$THREADS" -o "${OUT}.sorted.bam" -
samtools index "${OUT}.sorted.bam"

echo "[3/6] Calling variants with bcftools..." >&2
bcftools mpileup -f "$REF" -Ou "${OUT}.sorted.bam" \
    | bcftools call -mv --ploidy 1 -Oz -o "${OUT}.raw.vcf.gz"
bcftools index -f "${OUT}.raw.vcf.gz"

echo "[4/6] Filtering variants (QUAL>=${MIN_QUAL}, DP>=${MIN_DEPTH})..." >&2
bcftools filter -e "QUAL<${MIN_QUAL} || DP<${MIN_DEPTH}" \
    "${OUT}.raw.vcf.gz" -Oz -o "${OUT}.filtered.vcf.gz"
bcftools index -f "${OUT}.filtered.vcf.gz"

echo "[5/6] Building consensus sequence (depth < ${MASK_DEPTH} masked as N)..." >&2
# -aa: report every reference position, including zero-depth ones and
# contigs with no reads at all. -J: count deletions as covered, so a
# called deletion is not also masked.
samtools depth -aa -J "${OUT}.sorted.bam" > "${OUT}.depth.tsv"
awk -v min="$MASK_DEPTH" 'BEGIN { OFS = "\t" }
    $3 < min {
        if ($1 == chr && $2 == end + 1) { end = $2; next }
        if (chr != "") print chr, start - 1, end
        chr = $1; start = $2; end = $2; next
    }
    { if (chr != "") print chr, start - 1, end; chr = "" }
    END { if (chr != "") print chr, start - 1, end }' \
    "${OUT}.depth.tsv" > "${OUT}.low_coverage.bed"
bcftools consensus -f "$REF" -m "${OUT}.low_coverage.bed" \
    "${OUT}.filtered.vcf.gz" > "${OUT}.consensus.fasta"

echo "[6/6] Writing variant call and coverage statistics..." >&2
bcftools stats "${OUT}.filtered.vcf.gz" > "${OUT}.stats.txt"
awk -v mask="$MASK_DEPTH" '
    { n++; if ($3 >= 1) c1++; if ($3 >= 10) c10++; if ($3 < mask) m++ }
    END {
        if (n == 0) { print "Error: no positions in depth table" > "/dev/stderr"; exit 1 }
        printf "reference_bases\t%d\n", n
        printf "breadth_ge1x_pct\t%.4f\n", 100 * c1 / n
        printf "breadth_ge10x_pct\t%.4f\n", 100 * c10 / n
        printf "masked_bases_depth_lt_%s\t%d\n", mask, m
    }' "${OUT}.depth.tsv" > "${OUT}.coverage_summary.tsv"
grep -v '^>' "${OUT}.consensus.fasta" | tr -d '\n' | awk '
    { len = length($0); n = gsub(/[Nn]/, ""); gc = gsub(/[GgCc]/, "") }
    END {
        printf "consensus_length\t%d\n", len
        printf "consensus_N_bases\t%d\n", n
        printf "consensus_GC_pct_excl_N\t%.4f\n", (len - n) ? 100 * gc / (len - n) : 0
    }' >> "${OUT}.coverage_summary.tsv"

echo "Done." >&2
echo "  Sorted BAM:        ${OUT}.sorted.bam" >&2
echo "  Filtered variants: ${OUT}.filtered.vcf.gz" >&2
echo "  Consensus genome:  ${OUT}.consensus.fasta (N = no read support)" >&2
echo "  Masked regions:    ${OUT}.low_coverage.bed" >&2
echo "  Per-base depth:    ${OUT}.depth.tsv" >&2
echo "  Coverage summary:  ${OUT}.coverage_summary.tsv" >&2
echo "  VCF statistics:    ${OUT}.stats.txt" >&2
echo >&2
sed 's/^/    /' "${OUT}.coverage_summary.tsv" >&2
