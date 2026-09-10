#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<EOF
Usage: $(basename "$0") -q CONSENSUS.fasta -o OUTPUT_DIR [-t "Taxon1,Taxon2,..."]

Downloads one reference genome per taxon from NCBI and computes whole-genome
ANI (FastANI) against your consensus sequence -- use this to resolve
species-level identity before trusting a single BLAST hit's percent identity.

Required:
  -q  Your consensus genome, FASTA
  -o  Output directory (downloaded reference genomes + results go here)

Optional:
  -t  Comma-separated taxa to compare against (default:
      "Bacillus licheniformis,Bacillus paralicheniformis,Bacillus sonorensis,
      Bacillus subtilis,Bacillus haynesii")

Requires the 'datasets' CLI (NCBI Datasets) and 'fastANI' on PATH, plus
'unzip'. Install via conda:
  conda install -c conda-forge ncbi-datasets-cli
  conda install -c bioconda fastani
  sudo apt-get install -y unzip
EOF
    exit 1
}

TAXA="Bacillus licheniformis,Bacillus paralicheniformis,Bacillus sonorensis,Bacillus subtilis,Bacillus haynesii"

while getopts "q:o:t:h" opt; do
    case "$opt" in
        q) QUERY="$OPTARG" ;;
        o) OUTDIR="$OPTARG" ;;
        t) TAXA="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done

if [[ -z "${QUERY:-}" || -z "${OUTDIR:-}" ]]; then
    usage
fi

[[ -f "$QUERY" ]] || { echo "Error: file not found: $QUERY" >&2; exit 1; }

for tool in datasets unzip fastANI; do
    command -v "$tool" >/dev/null 2>&1 || { echo "Error: $tool not found in PATH" >&2; exit 1; }
done

mkdir -p "$OUTDIR/genomes"
REF_LIST="$OUTDIR/ref_list.txt"
: > "$REF_LIST"

IFS=',' read -ra TAXA_ARR <<< "$TAXA"
for taxon in "${TAXA_ARR[@]}"; do
    slug=$(echo "$taxon" | tr ' ' '_')
    zipfile="$OUTDIR/${slug}.zip"

    echo "Downloading reference genome for: $taxon" >&2
    if ! datasets download genome taxon "$taxon" --reference --include genome --filename "$zipfile"; then
        echo "Warning: no reference-flagged genome for '$taxon'; retrying without --reference" >&2
        datasets download genome taxon "$taxon" --include genome --filename "$zipfile"
    fi

    unzip -oq "$zipfile" -d "$OUTDIR/genomes/$slug"
    genome_fasta=$(find "$OUTDIR/genomes/$slug" -name "*.fna" | head -n 1)
    if [[ -z "$genome_fasta" ]]; then
        echo "Warning: no genome FASTA found for $taxon, skipping" >&2
        continue
    fi
    cp "$genome_fasta" "$OUTDIR/genomes/${slug}.fasta"
    echo "$OUTDIR/genomes/${slug}.fasta" >> "$REF_LIST"
done

if [[ ! -s "$REF_LIST" ]]; then
    echo "Error: no reference genomes were downloaded for any taxon" >&2
    exit 1
fi

echo "Running FastANI ($(wc -l < "$REF_LIST") reference genome(s))..." >&2
fastANI -q "$QUERY" --rl "$REF_LIST" -o "$OUTDIR/ani_results.tsv"

echo "" >&2
echo "Results, sorted by ANI (highest first):" >&2
sort -t$'\t' -k3,3 -rn "$OUTDIR/ani_results.tsv" \
    | awk -F'\t' 'BEGIN{OFS="\t"; print "query","reference","ANI_percent","fragments_matched","fragments_total"} {print}' \
    | column -t -s$'\t'
