#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<EOF
Usage: $(basename "$0") -q CONSENSUS.fasta -o OUTPUT_DIR [-t "Taxon"] [-n MAX_GENOMES] [-l ASSEMBLY_LEVEL]

Downloads multiple genome assemblies for a single species from NCBI and
ranks them by whole-genome ANI (FastANI) against your consensus sequence,
to find the closest known strain within that species.

The candidate set is every assembly NCBI returns for the taxon at the
chosen level(s), with atypical assemblies excluded and GenBank/RefSeq
pairs (GCA_x / GCF_x of the same assembly) collapsed to one entry, keeping
the RefSeq copy. Nothing is silently truncated: if -n is set and more
assemblies are found, the script stops instead of picking an arbitrary
subset. The query date, datasets version, query and counts are written to
OUTPUT_DIR/candidate_set_manifest.tsv; cite them in the methods.

Required:
  -q  Your consensus genome, FASTA
  -o  Output directory (downloaded genomes + results go here)

Optional:
  -t  Taxon to search within (default: "Bacillus paralicheniformis")
  -n  Safety limit on the number of genomes (default: 0 = no limit). If
      more assemblies than this are found, the script exits with an error
      rather than comparing an arbitrary subset; narrow -l or raise -n.
  -l  Assembly level filter passed to 'datasets', comma-separated, from
      complete, chromosome, scaffold, contig (default: all four).
      Completeness is unrelated to relatedness, so restricting to complete
      genomes can drop the closest strain. Check your installed CLI's
      accepted values with:
      datasets summary genome taxon --help

Requires 'datasets' (NCBI Datasets CLI), 'jq', 'unzip', and 'fastANI' on
PATH. Install via conda:
  conda install -c conda-forge ncbi-datasets-cli jq
  conda install -c bioconda fastani
  sudo apt-get install -y unzip

Note: this script relies on the installed 'datasets' CLI's flags and JSON
field names (--assembly-level values, --exclude-atypical, and the
per-assembly "accession" and "paired_accession" fields in --as-json-lines
output). These have changed across CLI versions -- if
-l's default is rejected, pass the value your version's --help lists; if
the accession lookup returns nothing, run:
  datasets summary genome taxon "Bacillus paralicheniformis" --as-json-lines | head -n 1 | jq .
to confirm the field name, and adjust the jq filter below if needed.
EOF
    exit 1
}

TAXON="Bacillus paralicheniformis"
MAX_GENOMES=0
ASSEMBLY_LEVEL="complete,chromosome,scaffold,contig"

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
[[ "$MAX_GENOMES" =~ ^[0-9]+$ ]] || { echo "Error: -n must be a non-negative integer" >&2; exit 1; }

for tool in datasets jq unzip fastANI; do
    command -v "$tool" >/dev/null 2>&1 || { echo "Error: $tool not found in PATH" >&2; exit 1; }
done

mkdir -p "$OUTDIR/genomes/raw"

echo "Looking up ${ASSEMBLY_LEVEL} assemblies for: $TAXON" >&2
SUMMARY_JSONL="$OUTDIR/assembly_summary.jsonl"
ACCESSIONS_FILE="$OUTDIR/accessions.txt"
datasets summary genome taxon "$TAXON" --assembly-level "$ASSEMBLY_LEVEL" \
    --exclude-atypical --as-json-lines > "$SUMMARY_JSONL"

N_RAW=$(jq -r '.accession' "$SUMMARY_JSONL" | sort -u | wc -l)
# Collapse GCA_/GCF_ pairs: keep every GCF_, and a GCA_ only if its paired
# RefSeq accession is not itself in the result set.
jq -rs '
    (map(.accession) | map({(.): true}) | add // {}) as $have
    | .[]
    | select((.accession | startswith("GCF_"))
             or ((.paired_accession // "") as $p | ($p == "" or ($have[$p] | not))))
    | .accession' "$SUMMARY_JSONL" | sort -u > "$ACCESSIONS_FILE"

N_ACC=$(wc -l < "$ACCESSIONS_FILE")
if [[ "$N_ACC" -eq 0 ]]; then
    echo "Error: no ${ASSEMBLY_LEVEL} assemblies found for '$TAXON'." >&2
    echo "Try a different -l value or check the taxon name." >&2
    exit 1
fi
if [[ "$MAX_GENOMES" -gt 0 && "$N_ACC" -gt "$MAX_GENOMES" ]]; then
    echo "Error: $N_ACC assemblies found, more than the -n limit of $MAX_GENOMES." >&2
    echo "Refusing to compare an arbitrary subset; narrow -l or raise/unset -n." >&2
    exit 1
fi

MANIFEST="$OUTDIR/candidate_set_manifest.tsv"
{
    printf 'query_date_utc\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    printf 'datasets_version\t%s\n' "$(datasets --version 2>&1 | head -n 1)"
    printf 'taxon\t%s\n' "$TAXON"
    printf 'assembly_level\t%s\n' "$ASSEMBLY_LEVEL"
    printf 'exclude_atypical\ttrue\n'
    printf 'assemblies_returned\t%s\n' "$N_RAW"
    printf 'assemblies_after_GCA_GCF_dereplication\t%s\n' "$N_ACC"
    printf 'accession_list\t%s\n' "$ACCESSIONS_FILE"
} > "$MANIFEST"
echo "Found $N_RAW assemblies, $N_ACC after GCA/GCF dereplication; downloading..." >&2

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
echo "Candidate-set manifest (cite in methods): $MANIFEST" >&2
