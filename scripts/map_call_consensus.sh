#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<EOF
Usage: $(basename "$0") -r REFERENCE.fasta -1 R1.fastq.gz -2 R2.fastq.gz -o OUTPUT_PREFIX [-t THREADS] [-q MIN_QUAL] [-d MIN_DEPTH]

Reference-based mapping, variant calling, and consensus generation for
paired-end Illumina reads (bwa-mem + bcftools). No de novo assembly.

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
EOF
    exit 1
}

THREADS=4
MIN_QUAL=20
MIN_DEPTH=10

while getopts "r:1:2:o:t:q:d:h" opt; do
    case "$opt" in
        r) REF="$OPTARG" ;;
        1) R1="$OPTARG" ;;
        2) R2="$OPTARG" ;;
        o) OUT="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        q) MIN_QUAL="$OPTARG" ;;
        d) MIN_DEPTH="$OPTARG" ;;
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
    | bcftools call -mv -Oz -o "${OUT}.raw.vcf.gz"
bcftools index -f "${OUT}.raw.vcf.gz"

echo "[4/6] Filtering variants (QUAL>=${MIN_QUAL}, DP>=${MIN_DEPTH})..." >&2
bcftools filter -e "QUAL<${MIN_QUAL} || DP<${MIN_DEPTH}" \
    "${OUT}.raw.vcf.gz" -Oz -o "${OUT}.filtered.vcf.gz"
bcftools index -f "${OUT}.filtered.vcf.gz"

echo "[5/6] Building consensus sequence..." >&2
bcftools consensus -f "$REF" "${OUT}.filtered.vcf.gz" > "${OUT}.consensus.fasta"

echo "[6/6] Writing variant call statistics..." >&2
bcftools stats "${OUT}.filtered.vcf.gz" > "${OUT}.stats.txt"

echo "Done." >&2
echo "  Sorted BAM:        ${OUT}.sorted.bam" >&2
echo "  Filtered variants: ${OUT}.filtered.vcf.gz" >&2
echo "  Consensus genome:  ${OUT}.consensus.fasta" >&2
echo "  VCF statistics:    ${OUT}.stats.txt" >&2
