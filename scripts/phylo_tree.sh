#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<EOF
Usage: $(basename "$0") -q ISOLATE_CONSENSUS.fasta -g GENOMES_DIR -b BACKBONE.fasta -o OUTPUT_DIR [-t THREADS]

Builds a bootstrap-supported phylogenetic tree placing your isolate among
a set of reference genomes: Parsnp aligns the core genome shared across all
inputs (fast, designed for closely related bacterial genomes), then IQ-TREE
builds a maximum-likelihood tree with ultrafast bootstrap support from that
alignment.

Required:
  -q  Isolate consensus genome, FASTA
  -g  Directory of reference genome FASTA files to include (e.g. the
      'genomes' folder already produced by fastani_strain_search.sh or
      fastani_species_check.sh)
  -b  One reference genome FASTA to use as Parsnp's alignment backbone
      (typically your closest strain, e.g. the BL-09 genome)
  -o  Output directory

Optional:
  -t  Threads (default: 4)

Requires 'parsnp', 'harvesttools', and 'iqtree2' (or 'iqtree') on PATH.
Install via conda:
  conda install -c bioconda parsnp harvesttools iqtree
EOF
    exit 1
}

THREADS=4

while getopts "q:g:b:o:t:h" opt; do
    case "$opt" in
        q) QUERY="$OPTARG" ;;
        g) GENOMES_DIR="$OPTARG" ;;
        b) BACKBONE="$OPTARG" ;;
        o) OUTDIR="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done

if [[ -z "${QUERY:-}" || -z "${GENOMES_DIR:-}" || -z "${BACKBONE:-}" || -z "${OUTDIR:-}" ]]; then
    usage
fi

[[ -f "$QUERY" ]] || { echo "Error: file not found: $QUERY" >&2; exit 1; }
[[ -d "$GENOMES_DIR" ]] || { echo "Error: directory not found: $GENOMES_DIR" >&2; exit 1; }
[[ -f "$BACKBONE" ]] || { echo "Error: file not found: $BACKBONE" >&2; exit 1; }

command -v parsnp >/dev/null 2>&1 || { echo "Error: parsnp not found in PATH" >&2; exit 1; }
command -v harvesttools >/dev/null 2>&1 || { echo "Error: harvesttools not found in PATH" >&2; exit 1; }

IQTREE_BIN=""
for cand in iqtree2 iqtree; do
    if command -v "$cand" >/dev/null 2>&1; then
        IQTREE_BIN="$cand"
        break
    fi
done
[[ -n "$IQTREE_BIN" ]] || { echo "Error: neither iqtree2 nor iqtree found in PATH" >&2; exit 1; }

mkdir -p "$OUTDIR"
INPUT_DIR="$OUTDIR/genomes_for_tree"
mkdir -p "$INPUT_DIR"

BACKBONE_REAL=$(realpath "$BACKBONE")
for f in "$GENOMES_DIR"/*.fasta; do
    [[ -e "$f" ]] || continue
    if [[ "$(realpath "$f")" == "$BACKBONE_REAL" ]]; then
        continue
    fi
    cp "$f" "$INPUT_DIR"/
done
awk '/^>/ { print ">isolate_query"; next } { print }' "$QUERY" > "$INPUT_DIR/isolate_query.fasta"

DUPES=$(grep -h '^>' "$INPUT_DIR"/*.fasta "$BACKBONE" | sort | uniq -d || true)
if [[ -n "$DUPES" ]]; then
    echo "Error: duplicate sequence header(s) found across input genomes -- Parsnp indexes by header and will fail or crash on these:" >&2
    echo "$DUPES" >&2
    exit 1
fi

N_GENOMES=$(find "$INPUT_DIR" -name "*.fasta" | wc -l)
if [[ "$N_GENOMES" -lt 3 ]]; then
    echo "Error: only $N_GENOMES genome(s) found in $INPUT_DIR (need at least 3 for a meaningful tree)" >&2
    exit 1
fi
echo "Building tree from $N_GENOMES genomes (including the isolate and the backbone)..." >&2

PARSNP_OUT="$OUTDIR/parsnp_run"
rm -rf "$PARSNP_OUT"
parsnp -r "$BACKBONE" -d "$INPUT_DIR" -o "$PARSNP_OUT" -c -p "$THREADS"

echo "Extracting core-genome alignment..." >&2
harvesttools -i "$PARSNP_OUT/parsnp.ggr" -M "$OUTDIR/core_alignment.fasta"

echo "Running IQ-TREE (GTR+G, 1000 ultrafast bootstrap replicates)..." >&2
"$IQTREE_BIN" -s "$OUTDIR/core_alignment.fasta" -m GTR+G -bb 1000 -nt AUTO -pre "$OUTDIR/tree" -redo

echo "" >&2
echo "Done. Tree file (Newick, with bootstrap support): $OUTDIR/tree.treefile" >&2
echo "View it by uploading that file to https://itol.embl.de (or open in FigTree)." >&2
echo "Ultrafast bootstrap values are percentages out of 100; conventionally, >=95 is considered strong support." >&2
