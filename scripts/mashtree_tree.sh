#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<EOF
Usage: $(basename "$0") -q ISOLATE_CONSENSUS.fasta -g GENOMES_DIR -o OUTPUT_DIR [-r BOOTSTRAP_REPS] [-t THREADS]

Builds a bootstrap-supported phylogenetic tree placing your isolate among
a set of reference genomes, using Mashtree (Mash-distance-based
neighbor-joining). This is a lighter-weight, more robust alternative to
Parsnp+IQ-TREE (phylo_tree.sh) when that toolchain fails in a given
environment -- it has no compiled core-alignment binary to crash, and it
labels tree tips by filename rather than FASTA header, so header
collisions are a non-issue.

Required:
  -q  Isolate consensus genome, FASTA
  -g  Directory of reference genome FASTA files to include (e.g. the
      'genomes' folder produced by fastani_strain_search.sh or
      fastani_species_check.sh)
  -o  Output directory

Optional:
  -r  Bootstrap replicates (default: 100)
  -t  Threads (default: 4)

Requires 'mashtree_bootstrap.pl' on PATH. Install via conda:
  conda install -c bioconda mashtree
EOF
    exit 1
}

REPS=100
THREADS=4

while getopts "q:g:o:r:t:h" opt; do
    case "$opt" in
        q) QUERY="$OPTARG" ;;
        g) GENOMES_DIR="$OPTARG" ;;
        o) OUTDIR="$OPTARG" ;;
        r) REPS="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done

if [[ -z "${QUERY:-}" || -z "${GENOMES_DIR:-}" || -z "${OUTDIR:-}" ]]; then
    usage
fi

[[ -f "$QUERY" ]] || { echo "Error: file not found: $QUERY" >&2; exit 1; }
[[ -d "$GENOMES_DIR" ]] || { echo "Error: directory not found: $GENOMES_DIR" >&2; exit 1; }

command -v mashtree_bootstrap.pl >/dev/null 2>&1 || { echo "Error: mashtree_bootstrap.pl not found in PATH" >&2; exit 1; }

mkdir -p "$OUTDIR"
INPUT_DIR="$OUTDIR/genomes_for_tree"
rm -rf "$INPUT_DIR"
mkdir -p "$INPUT_DIR"

cp "$GENOMES_DIR"/*.fasta "$INPUT_DIR"/
cp "$QUERY" "$INPUT_DIR/isolate_query.fasta"

GENOME_FILES=("$INPUT_DIR"/*.fasta)
N_GENOMES=${#GENOME_FILES[@]}
if [[ "$N_GENOMES" -lt 3 ]]; then
    echo "Error: only $N_GENOMES genome(s) found (need at least 3 for a meaningful tree)" >&2
    exit 1
fi
echo "Building tree from $N_GENOMES genomes ($REPS bootstrap replicates)..." >&2

mashtree_bootstrap.pl --reps "$REPS" --numcpus "$THREADS" "${GENOME_FILES[@]}" > "$OUTDIR/tree.dnd"

echo "" >&2
echo "Done. Tree file (Newick, with bootstrap support): $OUTDIR/tree.dnd" >&2
echo "View it by uploading that file to https://itol.embl.de (or open in FigTree)." >&2
echo "Tips are labeled by input filename (e.g. isolate_query.fasta, GCA_000876525.1.fasta)." >&2
echo "Bootstrap values are fractions out of 1.0; conventionally, >=0.95 is considered strong support." >&2
