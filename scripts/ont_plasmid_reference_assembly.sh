#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<EOF
Usage: $(basename "$0") -r REFERENCE.fasta -i READS.fastq[.gz] -o OUTPUT_PREFIX [options]

Reference-based mapping, variant calling, and consensus generation for a
plasmid sequenced with Oxford Nanopore (ONT) long reads (minimap2 +
bcftools), plus a structural difference report that accounts for large
deletions and insertions relative to the reference.

This is NOT de novo assembly: it reconstructs the sample's plasmid by
mapping reads onto the reference you supply. It can therefore measure how
much of the reference is present in the sample and where the sample
departs from it, but it can never reconstruct sequence that is absent
from the reference -- inserted sequence is reported and extracted from
the reads, not assembled. Cross-check the result (e.g. with
scripts/ncbi_strain_search.py) before treating it as a finished sequence.

Regions with no read support are masked with N in the consensus rather
than silently inheriting reference bases, and a trimmed consensus with
those regions removed is written alongside it. Use this to confirm, for
example, that a cassette was excised and the rest re-cloned intact.

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
  -s  Minimum size in bp for a structural event (deletion/insertion) to be
      reported (default: 50). bcftools only calls small variants, so
      anything at this scale comes from the alignment analysis instead.
  -k  Depth below which a base counts as having NO read support and is
      masked with N (default: 1, i.e. only zero-coverage bases). This is
      deliberately not the same as -d: masking removes sequence from the
      supported consensus, so it should only fire where there is genuinely
      no evidence, not merely low coverage.

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
MIN_SV_SIZE=50
MASK_DEPTH=1

while getopts "r:i:o:t:q:d:D:Q:p:s:k:ch" opt; do
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
        s) MIN_SV_SIZE="$OPTARG" ;;
        k) MASK_DEPTH="$OPTARG" ;;
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

STEP_TOTAL=8
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

echo "[${STEP}/${STEP_TOTAL}] Analyzing structural differences (large indels, clips, coverage)..." >&2
REPORT_ARGS=(report -a "${OUT}.sorted.bam" -r "$MAPREF" -o "$OUT"
             --min-sv-size "$MIN_SV_SIZE" --min-depth "$MIN_DEPTH"
             --mask-depth "$MASK_DEPTH")
if [[ "$CIRCULAR" -eq 1 ]]; then
    REPORT_ARGS+=(--circular-meta "${MAPREF}.meta")
fi
python3 "${SCRIPT_DIR}/indel_report.py" "${REPORT_ARGS[@]}"
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

echo "[${STEP}/${STEP_TOTAL}] Building consensus sequence (regions without read support masked as N)..." >&2
if [[ "$CIRCULAR" -eq 1 ]]; then
    bcftools consensus -f "$MAPREF" -m "${OUT}.low_coverage.bed" \
        "${OUT}.filtered.vcf.gz" > "${OUT}.consensus_extended.fasta"
    python3 "${SCRIPT_DIR}/rotate_circular_fasta.py" reconcile \
        -i "${OUT}.consensus_extended.fasta" --meta "${MAPREF}.meta" \
        -o "${OUT}.consensus.fasta"
else
    bcftools consensus -f "$MAPREF" -m "${OUT}.low_coverage.bed" \
        "${OUT}.filtered.vcf.gz" > "${OUT}.consensus.fasta"
fi
STEP=$((STEP + 1))

echo "[${STEP}/${STEP_TOTAL}] Writing supported-only consensus and variant statistics..." >&2
python3 "${SCRIPT_DIR}/indel_report.py" trim \
    -c "${OUT}.consensus.fasta" -o "${OUT}.consensus_supported.fasta"
bcftools stats "${OUT}.filtered.vcf.gz" > "${OUT}.stats.txt"
STEP=$((STEP + 1))

echo "Done." >&2
echo "  Sorted BAM:            ${OUT}.sorted.bam" >&2
echo "  Filtered variants:     ${OUT}.filtered.vcf.gz" >&2
echo "  Consensus (N-masked):  ${OUT}.consensus.fasta" >&2
echo "  Consensus (supported): ${OUT}.consensus_supported.fasta" >&2
echo "  Structural report:     ${OUT}.indel_report.txt" >&2
echo "  Structural variants:   ${OUT}.structural_variants.tsv" >&2
echo "  Unsupported regions:   ${OUT}.low_coverage.bed" >&2
echo "  Inserted sequences:    ${OUT}.inserted_sequences.fasta (if any)" >&2
echo "  Per-base depth:        ${OUT}.depth.tsv" >&2
echo "  VCF statistics:        ${OUT}.stats.txt" >&2
echo "" >&2
echo "Read ${OUT}.indel_report.txt first: it states how much of the reference" >&2
echo "the reads support, what is missing, and what the sample carries that the" >&2
echo "reference does not." >&2
