#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<EOF
Usage: $(basename "$0") -r REFERENCE.fasta -i READS.fastq[.gz] -o OUTPUT_PREFIX [options]

Reference-based mapping, variant calling, and consensus generation for a
plasmid sequenced with Oxford Nanopore (ONT) long reads (minimap2 +
bcftools). This is NOT de novo assembly: it reconstructs the sample's
plasmid sequence by mapping reads onto the closest reference you supply
and calling the consensus, so it can only find differences from that
reference (SNPs/small indels/coverage gaps) -- it cannot discover larger
structural rearrangements, insertions of unrelated sequence, or confirm
that the plasmid is genuinely circular in the sample. Cross-check the
result (e.g. with scripts/ncbi_strain_search.py) before treating it as a
finished sequence.

Required:
  -r  Reference plasmid, FASTA, single record (closest known/public backbone)
  -i  ONT reads, FASTQ(.gz), single file
  -o  Output prefix (path/basename for all generated files)

Optional:
  -c  Treat the plasmid as circular (default: off). Enables the
      reference-doubling step so reads spanning the reference's
      start/end junction map contiguously instead of being clipped.
      Leave this off only if you know the input is a linear molecule.
  -p  Pad length in bp for circular handling (default: 2000). Should be
      at least as long as the reads you expect to span the junction --
      check your read length distribution (e.g. with NanoPlot) and raise
      this if many reads are longer than the default. Ignored without -c.
  -t  Threads (default: 4)
  -q  Minimum variant QUAL to keep (default: 20) -- provisional default,
      review against your own coverage and error profile before trusting it
  -d  Minimum read depth (DP) to keep a variant (default: 20) -- ONT
      plasmid runs are often very high depth; raise this if yours is
      -- provisional default, review before trusting it
  -D  Max depth passed to bcftools mpileup (default: 4000) -- raise if
      your per-base depth exceeds this, or pileup will be capped
  -Q  Minimum base quality for mpileup (default: 7) -- ONT base qualities
      are on a different scale than Illumina's; this is a lenient
      provisional default, review before trusting it

Notes:
  - bcftools mpileup/call was designed primarily for short, high base-
    quality reads. It is used here for consistency with
    scripts/map_call_consensus.sh and because it requires no extra
    services, but it is not the most accurate consensus caller for ONT
    data, especially in homopolymers. For higher-accuracy consensus
    polishing consider Oxford Nanopore's own medaka
    (https://github.com/nanoporetech/medaka) as a follow-up step.
  - BAQ (per-base alignment quality recalibration) is disabled (-B) in
    mpileup, since it assumes a low-indel-rate short-read error model
    that does not fit ONT's higher indel rate.
EOF
    exit 1
}

THREADS=4
MIN_QUAL=20
MIN_DEPTH=20
MAX_DEPTH=4000
MIN_BASEQ=7
PAD=2000
CIRCULAR=0

while getopts "r:i:o:t:q:d:D:Q:p:ch" opt; do
    case "$opt" in
        r) REF="$OPTARG" ;;
        i) READS="$OPTARG" ;;
        o) OUT="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        q) MIN_QUAL="$OPTARG" ;;
        d) MIN_DEPTH="$OPTARG" ;;
        D) MAX_DEPTH="$OPTARG" ;;
        Q) MIN_BASEQ="$OPTARG" ;;
        p) PAD="$OPTARG" ;;
        c) CIRCULAR=1 ;;
        h) usage ;;
        *) usage ;;
    esac
done

if [[ -z "${REF:-}" || -z "${READS:-}" || -z "${OUT:-}" ]]; then
    usage
fi

for tool in minimap2 samtools bcftools python3; do
    command -v "$tool" >/dev/null 2>&1 || { echo "Error: $tool not found in PATH" >&2; exit 1; }
done

for f in "$REF" "$READS"; do
    [[ -f "$f" ]] || { echo "Error: file not found: $f" >&2; exit 1; }
done

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

outdir="$(dirname "$OUT")"
[[ "$outdir" == "." ]] || mkdir -p "$outdir"

STEP_TOTAL=7
STEP=1

if [[ "$CIRCULAR" -eq 1 ]]; then
    echo "[${STEP}/${STEP_TOTAL}] Building doubled reference for circular handling (pad=${PAD}bp)..." >&2
    MAPREF="${OUT}.ref_extended.fasta"
    python3 "${SCRIPT_DIR}/rotate_circular_fasta.py" extend \
        -i "$REF" -o "$MAPREF" --pad "$PAD"
else
    echo "[${STEP}/${STEP_TOTAL}] Circular handling disabled (-c not set); using reference as-is." >&2
    MAPREF="$REF"
fi
STEP=$((STEP + 1))

samtools faidx "$MAPREF"

echo "[${STEP}/${STEP_TOTAL}] Mapping ONT reads with minimap2 (map-ont)..." >&2
minimap2 -ax map-ont -t "$THREADS" "$MAPREF" "$READS" \
    | samtools sort -@ "$THREADS" -o "${OUT}.sorted.bam" -
samtools index "${OUT}.sorted.bam"
STEP=$((STEP + 1))

echo "[${STEP}/${STEP_TOTAL}] Calling variants with bcftools (haploid, BAQ off)..." >&2
bcftools mpileup -f "$MAPREF" -B -Q "$MIN_BASEQ" --max-depth "$MAX_DEPTH" -Ou "${OUT}.sorted.bam" \
    | bcftools call -mv --ploidy 1 -Oz -o "${OUT}.raw.vcf.gz"
bcftools index -f "${OUT}.raw.vcf.gz"
STEP=$((STEP + 1))

echo "[${STEP}/${STEP_TOTAL}] Filtering variants (QUAL>=${MIN_QUAL}, DP>=${MIN_DEPTH})..." >&2
bcftools filter -e "QUAL<${MIN_QUAL} || DP<${MIN_DEPTH}" \
    "${OUT}.raw.vcf.gz" -Oz -o "${OUT}.filtered.vcf.gz"
bcftools index -f "${OUT}.filtered.vcf.gz"
STEP=$((STEP + 1))

echo "[${STEP}/${STEP_TOTAL}] Building consensus sequence..." >&2
if [[ "$CIRCULAR" -eq 1 ]]; then
    bcftools consensus -f "$MAPREF" "${OUT}.filtered.vcf.gz" > "${OUT}.consensus_extended.fasta"
    python3 "${SCRIPT_DIR}/rotate_circular_fasta.py" reconcile \
        -i "${OUT}.consensus_extended.fasta" --meta "${MAPREF}.meta" \
        -o "${OUT}.consensus.fasta"
else
    bcftools consensus -f "$MAPREF" "${OUT}.filtered.vcf.gz" > "${OUT}.consensus.fasta"
fi
STEP=$((STEP + 1))

echo "[${STEP}/${STEP_TOTAL}] Writing coverage and variant call statistics..." >&2
samtools depth -a "${OUT}.sorted.bam" > "${OUT}.depth.txt"
awk -v mind="$MIN_DEPTH" '
    { sum += $3; n++; if ($3 >= mind) covered++; if ($3 == 0) zero++ }
    END {
        if (n == 0) { print "No positions in depth file."; exit }
        printf "Reference positions (mapping coordinates): %d\n", n
        printf "Mean depth: %.1f\n", sum / n
        printf "Positions with depth >= %d: %d (%.2f%%)\n", mind, covered, 100.0 * covered / n
        printf "Positions with zero depth: %d (%.2f%%)\n", zero, 100.0 * zero / n
    }
' "${OUT}.depth.txt" > "${OUT}.coverage_summary.txt"
bcftools stats "${OUT}.filtered.vcf.gz" > "${OUT}.stats.txt"
STEP=$((STEP + 1))

echo "Done." >&2
echo "  Sorted BAM:         ${OUT}.sorted.bam" >&2
echo "  Filtered variants:  ${OUT}.filtered.vcf.gz" >&2
echo "  Consensus plasmid:  ${OUT}.consensus.fasta" >&2
echo "  Coverage summary:   ${OUT}.coverage_summary.txt" >&2
echo "  VCF statistics:     ${OUT}.stats.txt" >&2
echo "" >&2
echo "Zero-depth or low-depth regions in the coverage summary mean the" >&2
echo "reference backbone may not match the sample there -- inspect" >&2
echo "${OUT}.depth.txt before trusting the consensus in those stretches." >&2
