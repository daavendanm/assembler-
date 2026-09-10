#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<EOF
Usage: $(basename "$0") -q CONSENSUS.fasta -o OUTPUT_DIR [-t "Taxon"] [-n MAX_GENOMES] [-l ASSEMBLY_LEVEL]

Downloads multiple genome assemblies for a single species from NCBI and
ranks them by whole-genome ANI (FastANI) against your consensus sequence,
to find the closest known strain within that species.

Required:
  -q  Your consensus genome, FASTA
  -o  Output directory (downloaded genomes + results go here)

Optional:
  -t  Taxon to search within (default: "Bacillus paralicheniformis")
  -n  Maximum number of genomes to download and compare (default: 25)
  -l  Assembly level filter passed to 'datasets': one of
      complete, chromosome, scaffold, contig (default: complete --
      keeps the download small and the genomes high-quality; loosen
      this if too few complete genomes exist). Check your installed
      CLI's accepted values with:
      datasets summary genome taxon --help

Requires 'datasets' (NCBI Datasets CLI), 'jq', 'unzip', and 'fastANI' on
PATH. Install via conda:
  conda install -c conda-forge ncbi-datasets-cli jq
  conda install -c bioconda fastani
  sudo apt-get install -y unzip

Note: this script relies on the installed 'datasets' CLI's flags and JSON
field names (--assembly-level values, the per-assembly "accession" field
in --as-json-lines output). These have changed across CLI versions -- if
-l's default is rejected, pass the value your version's --help lists; if
the accession lookup returns nothing, run:
  datasets summary genome taxon "Bacillus paralicheniformis" --as-json-lines | head -n 1 | jq .
to confirm the field name, and adjust the jq filter below if needed.
EOF
    exit 1
}

TAXON="Bacillus paralicheniformis"
MAX_GENOMES=25
ASSEMBLY_LEVEL="complete"

while getopts "q:o:t:n:l:h" opt; do
    case "$opt" in
        q) QUERY="$OPTARG" ;;
        o) OUTDIR="$OPTARG" ;;
        t) TAXON="$OPTARG" ;;
        n) MAX_GENOMES="$OPTARG" ;;
        l) ASSEMBLY_LEVEL="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done

if [[ -z "${QUERY:-}" || -z "${OUTDIR:-}" ]]; then
    usage
fi

[[ -f "$QUERY" ]] || { echo "Error: file not found: $QUERY" >&2; exit 1; }

for tool in datasets jq unzip fastANI; do
    command -v "$tool" >/dev/null 2>&1 || { echo "Error: $tool not found in PATH" >&2; exit 1; }
done

mkdir -p "$OUTDIR/genomes/raw"

echo "Looking up ${ASSEMBLY_LEVEL} assemblies for: $TAXON" >&2
ACCESSIONS_FILE="$OUTDIR/accessions.txt"
datasets summary genome taxon "$TAXON" --assembly-level "$ASSEMBLY_LEVEL" --as-json-lines \
    | jq -r '.accession' \
    | sort -u \
    | head -n "$MAX_GENOMES" > "$ACCESSIONS_FILE"

N_ACC=$(wc -l < "$ACCESSIONS_FILE")
if [[ "$N_ACC" -eq 0 ]]; then
    echo "Error: no ${ASSEMBLY_LEVEL} assemblies found for '$TAXON'." >&2
    echo "Try a different -l value (chromosome, scaffold, contig) or check the taxon name." >&2
    exit 1
fi
echo "Found $N_ACC assemblies (capped at $MAX_GENOMES); downloading..." >&2

ZIPFILE="$OUTDIR/genomes.zip"
datasets download genome accession --inputfile "$ACCESSIONS_FILE" --include genome --filename "$ZIPFILE"
unzip -oq "$ZIPFILE" -d "$OUTDIR/genomes/raw"

REF_LIST="$OUTDIR/ref_list.txt"
: > "$REF_LIST"
ACC_NAME_MAP="$OUTDIR/accession_names.tsv"
: > "$ACC_NAME_MAP"

while IFS= read -r acc; do
    genome_fasta=$(find "$OUTDIR/genomes/raw" -path "*${acc}*" -name "*.fna" | head -n 1)
    if [[ -z "$genome_fasta" ]]; then
        echo "Warning: no FASTA found for $acc, skipping" >&2
        continue
    fi
    dest="$OUTDIR/genomes/${acc}.fasta"
    cp "$genome_fasta" "$dest"
    echo "$dest" >> "$REF_LIST"
    header=$(head -n 1 "$dest")
    printf '%s\t%s\n' "$acc" "${header#>}" >> "$ACC_NAME_MAP"
done < "$ACCESSIONS_FILE"

if [[ ! -s "$REF_LIST" ]]; then
    echo "Error: no genomes were successfully prepared" >&2
    exit 1
fi

echo "Running FastANI against $(wc -l < "$REF_LIST") genome(s)..." >&2
fastANI -q "$QUERY" --rl "$REF_LIST" -o "$OUTDIR/ani_results.tsv"

echo "" >&2
echo "Closest strains (sorted by ANI, highest first):" >&2
sort -t$'\t' -k3,3 -rn "$OUTDIR/ani_results.tsv" \
    | awk -F'\t' 'BEGIN{OFS="\t"; print "query","reference","ANI_percent","fragments_matched","fragments_total"} {print}' \
    | column -t -s$'\t'

echo "" >&2
echo "Accession -> assembly description map: $ACC_NAME_MAP" >&2
