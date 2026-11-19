#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<EOF
Usage: $(basename "$0") -r REFERENCE.fasta -1 R1.fastq.gz -2 R2.fastq.gz -o OUTPUT_PREFIX [-t THREADS] [-q MIN_QUAL] [-d MIN_DEPTH] [-k MASK_DEPTH] [-m MIN_MAPQ] [-b MIN_BASEQ] [-D MAX_DEPTH] [-M SORT_MEM] [-K KRAKEN2_DB] [-T]

Reference-based mapping, variant calling, and consensus generation for
paired-end Illumina reads (bwa-mem + bcftools). No de novo assembly.

Variants are called with a haploid model (--ploidy 1). Reference positions
with read depth below MASK_DEPTH are written as N in the consensus instead
of silently inheriting the reference base, so the consensus contains only
sequence the reads actually support. Breadth of coverage and the number of
masked bases are written to the summary file. PCR/optical duplicates are
marked (samtools fixmate + markdup) and are therefore ignored by mpileup;
the duplicate rate and mapping statistics are written next to the BAM.

Reads are adapter/quality trimmed with fastp before mapping. Indels are
left-aligned and normalised (bcftools norm) before filtering, and SNPs
within 3 bp of an indel or indels within 10 bp of each other are thinned
(bcftools filter -g/-G), the classic source of false calls around
homopolymers. Reads that do not map to the reference are written to FASTQ
and, if a Kraken2 database is given, classified to show whether they are
contamination or strain-specific sequence.

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
  -m  Minimum mapping quality for mpileup (default: 20) -- excludes
      multi-mapping reads in rRNA operons, IS elements and other repeats
  -b  Minimum base quality for mpileup (default: 20)
  -D  Max per-file depth for mpileup (default: 2000). bcftools' own default
      is 250, which silently discards most reads at typical bacterial
      Illumina depth; set this to a few times your expected coverage.
  -M  Memory per samtools sort thread (default: 512M). Sort uses roughly
      THREADS x SORT_MEM; lower this or -t on low-memory machines (WSL,
      laptops) if sort fails with "couldn't allocate memory".
  -K  Kraken2 database directory, used to classify unmapped reads. Without
      it the unmapped reads are still extracted, but the summary records
      the contamination screen as NOT RUN.
  -T  Reads are already trimmed; skip fastp
EOF
    exit 1
}

THREADS=4
MIN_QUAL=20
MIN_DEPTH=10
MASK_DEPTH=1
MIN_MAPQ=20
MIN_BASEQ=20
MAX_DEPTH=2000
SORT_MEM=512M
KRAKEN_DB=""
SKIP_TRIM=0

while getopts "r:1:2:o:t:q:d:k:m:b:D:M:K:Th" opt; do
    case "$opt" in
        r) REF="$OPTARG" ;;
        1) R1="$OPTARG" ;;
        2) R2="$OPTARG" ;;
        o) OUT="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        q) MIN_QUAL="$OPTARG" ;;
        d) MIN_DEPTH="$OPTARG" ;;
        k) MASK_DEPTH="$OPTARG" ;;
        m) MIN_MAPQ="$OPTARG" ;;
        b) MIN_BASEQ="$OPTARG" ;;
        D) MAX_DEPTH="$OPTARG" ;;
        M) SORT_MEM="$OPTARG" ;;
        K) KRAKEN_DB="$OPTARG" ;;
        T) SKIP_TRIM=1 ;;
        h) usage ;;
        *) usage ;;
    esac
done

if [[ -z "${REF:-}" || -z "${R1:-}" || -z "${R2:-}" || -z "${OUT:-}" ]]; then
    usage
fi

REQUIRED_TOOLS=(bwa samtools bcftools)
[[ "$SKIP_TRIM" -eq 1 ]] || REQUIRED_TOOLS+=(fastp)
[[ -z "$KRAKEN_DB" ]] || REQUIRED_TOOLS+=(kraken2)
for tool in "${REQUIRED_TOOLS[@]}"; do
    command -v "$tool" >/dev/null 2>&1 || { echo "Error: $tool not found in PATH" >&2; exit 1; }
done

for f in "$REF" "$R1" "$R2"; do
    [[ -f "$f" ]] || { echo "Error: file not found: $f" >&2; exit 1; }
done
if [[ -n "$KRAKEN_DB" && ! -d "$KRAKEN_DB" ]]; then
    echo "Error: Kraken2 database directory not found: $KRAKEN_DB" >&2
    exit 1
fi

outdir="$(dirname "$OUT")"
[[ "$outdir" == "." ]] || mkdir -p "$outdir"

if [[ ! -f "${REF}.bwt" ]]; then
    echo "[1/8] Indexing reference..." >&2
    bwa index "$REF"
fi
[[ -f "${REF}.fai" ]] || samtools faidx "$REF"

if [[ "$SKIP_TRIM" -eq 1 ]]; then
    echo "[2/8] -T given; mapping reads as supplied (no trimming)." >&2
    MAP_R1="$R1"
    MAP_R2="$R2"
else
    echo "[2/8] Trimming adapters and low-quality tails with fastp..." >&2
    MAP_R1="${OUT}.trimmed_1.fastq.gz"
    MAP_R2="${OUT}.trimmed_2.fastq.gz"
    fastp -i "$R1" -I "$R2" -o "$MAP_R1" -O "$MAP_R2" \
        --detect_adapter_for_pe --cut_tail --cut_mean_quality 20 \
        --length_required 50 -w "$THREADS" \
        -j "${OUT}.fastp.json" -h "${OUT}.fastp.html"
fi

echo "[3/8] Mapping reads with bwa-mem and marking duplicates..." >&2
# bwa-mem output is grouped by read name, which is what fixmate needs.
# fixmate -m adds the mate-score tags markdup uses; duplicates are flagged
# (0x400), not removed, and mpileup's default filter then skips them.
bwa mem -t "$THREADS" "$REF" "$MAP_R1" "$MAP_R2" \
    | samtools fixmate -u -m -@ "$THREADS" - - \
    | samtools sort -u -@ "$THREADS" -m "$SORT_MEM" -T "${OUT}.sorttmp" - \
    | samtools markdup -@ "$THREADS" -f "${OUT}.markdup.txt" - "${OUT}.sorted.bam"
samtools index "${OUT}.sorted.bam"
samtools flagstat -@ "$THREADS" "${OUT}.sorted.bam" > "${OUT}.flagstat.txt"

echo "[4/8] Calling and normalising variants (MAPQ>=${MIN_MAPQ}, BQ>=${MIN_BASEQ}, max depth ${MAX_DEPTH})..." >&2
# norm left-aligns indels against the reference so the same event in a
# homopolymer/repeat always gets the same coordinates; -c e aborts if a
# VCF REF allele disagrees with the FASTA (e.g. a stale index).
bcftools mpileup -f "$REF" -q "$MIN_MAPQ" -Q "$MIN_BASEQ" -d "$MAX_DEPTH" \
    -a AD,DP,ADF,ADR -Ou "${OUT}.sorted.bam" \
    | bcftools call -mv --ploidy 1 -Ou \
    | bcftools norm -f "$REF" -c e -Oz -o "${OUT}.raw.vcf.gz"
bcftools index -f "${OUT}.raw.vcf.gz"

echo "[5/8] Filtering variants (QUAL>=${MIN_QUAL}, DP>=${MIN_DEPTH}, SnpGap 3, IndelGap 10)..." >&2
# INFO/DP (site depth) must be named explicitly: mpileup -a DP also adds
# FORMAT/DP, and bcftools >= 1.10 refuses a bare "DP" as ambiguous.
# -e drops failing records outright, but -g/-G only set FILTER=SnpGap /
# IndelGap; consensus would still apply those, so keep PASS records only.
bcftools filter -e "QUAL<${MIN_QUAL} || INFO/DP<${MIN_DEPTH}" -g 3 -G 10 \
    "${OUT}.raw.vcf.gz" -Ou \
    | bcftools view -f 'PASS,.' -Oz -o "${OUT}.filtered.vcf.gz"
bcftools index -f "${OUT}.filtered.vcf.gz"

echo "[6/8] Building consensus sequence (depth < ${MASK_DEPTH} masked as N)..." >&2
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

echo "[7/8] Extracting unmapped reads${KRAKEN_DB:+ and classifying with Kraken2}..." >&2
# -f 4: reads that are themselves unmapped. Pairs where both mates are
# unmapped go to _1/_2; an unmapped read whose mate mapped goes to _single.
samtools collate -u -O -@ "$THREADS" "${OUT}.sorted.bam" "${OUT}.collatetmp" \
    | samtools fastq -f 4 -@ "$THREADS" \
        -1 "${OUT}.unmapped_1.fastq.gz" -2 "${OUT}.unmapped_2.fastq.gz" \
        -s "${OUT}.unmapped_single.fastq.gz" -0 /dev/null -
if [[ -n "$KRAKEN_DB" ]]; then
    # Classified per read (not --paired) so singletons are included.
    kraken2 --db "$KRAKEN_DB" --threads "$THREADS" --gzip-compressed \
        --report "${OUT}.unmapped.kraken2.report" --output /dev/null \
        "${OUT}.unmapped_1.fastq.gz" "${OUT}.unmapped_2.fastq.gz" \
        "${OUT}.unmapped_single.fastq.gz"
fi

echo "[8/8] Writing variant call and coverage statistics..." >&2
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
{
    for v in raw filtered; do
        printf '%s_snps\t%d\n' "$v" "$(bcftools view -H -v snps "${OUT}.${v}.vcf.gz" | wc -l)"
        printf '%s_indels\t%d\n' "$v" "$(bcftools view -H -v indels "${OUT}.${v}.vcf.gz" | wc -l)"
    done
    printf 'reads_unmapped\t%d\n' "$(samtools view -c -f 4 "${OUT}.sorted.bam")"
    if [[ -n "$KRAKEN_DB" ]]; then
        # Top 5 species among unmapped reads, from the Kraken2 report.
        printf 'unmapped_kraken2_unclassified_pct\t%s\n' \
            "$(awk -F'\t' '$4 == "U" { print $1 + 0 }' "${OUT}.unmapped.kraken2.report")"
        printf 'unmapped_kraken2_top_species\t%s\n' \
            "$(awk -F'\t' '$4 == "S" { gsub(/^ +/, "", $6); printf "%s%% %s\n", $1 + 0, $6 }' \
                "${OUT}.unmapped.kraken2.report" | sort -rn | head -n 5 | paste -sd ';' -)"
    else
        printf 'unmapped_contamination_screen\tNOT RUN (no -K Kraken2 database)\n'
    fi
} >> "${OUT}.coverage_summary.tsv"

echo "Done." >&2
echo "  Sorted BAM:        ${OUT}.sorted.bam (duplicates marked)" >&2
echo "  Mapping stats:     ${OUT}.flagstat.txt" >&2
[[ "$SKIP_TRIM" -eq 1 ]] || echo "  fastp report:      ${OUT}.fastp.html" >&2
echo "  Unmapped reads:    ${OUT}.unmapped_{1,2,single}.fastq.gz" >&2
[[ -z "$KRAKEN_DB" ]] || echo "  Unmapped Kraken2:  ${OUT}.unmapped.kraken2.report" >&2
echo "  Duplicate stats:   ${OUT}.markdup.txt" >&2
echo "  Filtered variants: ${OUT}.filtered.vcf.gz" >&2
echo "  Consensus genome:  ${OUT}.consensus.fasta (N = no read support)" >&2
echo "  Masked regions:    ${OUT}.low_coverage.bed" >&2
echo "  Per-base depth:    ${OUT}.depth.tsv" >&2
echo "  Coverage summary:  ${OUT}.coverage_summary.tsv" >&2
echo "  VCF statistics:    ${OUT}.stats.txt" >&2
echo >&2
sed 's/^/    /' "${OUT}.coverage_summary.tsv" >&2
