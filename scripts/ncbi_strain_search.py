#!/usr/bin/env python3
"""Search NCBI (BLASTN) for public genomes matching a reference-based
consensus sequence, to identify the closest known strain.

Note: submitting a whole bacterial genome (a few Mb) as a single query to
NCBI's web BLAST is slow and not what the service is meant for. For
genome-scale comparison prefer a local BLAST+ database or FastANI/Mash
against downloaded genomes; this script is best suited to individual
contigs/regions of interest, or used sparingly with --delay respected.
"""

import argparse
import csv
import sys
import time

from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", help="Consensus FASTA file to search")
    parser.add_argument("-o", "--output", default="ncbi_hits.csv",
                         help="Output CSV path (default: ncbi_hits.csv)")
    parser.add_argument("--database", default="nt",
                         help="NCBI database to search (default: nt)")
    parser.add_argument("--program", default="blastn",
                         help="BLAST program (default: blastn)")
    parser.add_argument("--megablast", action="store_true",
                         help="Use megablast for highly similar sequences (faster)")
    parser.add_argument("--entrez-query", default="Bacillus licheniformis[Organism]",
                         help="Entrez query to restrict the search "
                              "(default: 'Bacillus licheniformis[Organism]'; "
                              "pass '' to search unrestricted)")
    parser.add_argument("--hitlist-size", type=int, default=10,
                         help="Number of hits to keep per query sequence (default: 10)")
    parser.add_argument("--delay", type=float, default=10.0,
                         help="Seconds to wait between successive NCBI requests "
                              "when the input has multiple sequences (default: 10)")
    return parser.parse_args()


def run_blast(record, args):
    entrez_query = args.entrez_query or None
    handle = NCBIWWW.qblast(
        program=args.program,
        database=args.database,
        sequence=record.format("fasta"),
        entrez_query=entrez_query,
        hitlist_size=args.hitlist_size,
        megablast=args.megablast,
    )
    blast_record = NCBIXML.read(handle)
    handle.close()
    return blast_record


def collect_hits(query_id, blast_record, query_length):
    rows = []
    for alignment in blast_record.alignments:
        hsp = max(alignment.hsps, key=lambda h: h.identities)
        percent_identity = 100.0 * hsp.identities / hsp.align_length
        query_coverage = 100.0 * hsp.align_length / query_length
        rows.append({
            "query_id": query_id,
            "hit_accession": alignment.accession,
            "hit_description": alignment.hit_def,
            "percent_identity": round(percent_identity, 2),
            "query_coverage_pct": round(query_coverage, 2),
            "evalue": hsp.expect,
        })
    rows.sort(key=lambda r: (-r["percent_identity"], -r["query_coverage_pct"]))
    return rows


def main():
    args = parse_args()
    records = list(SeqIO.parse(args.input, "fasta"))
    if not records:
        sys.exit(f"No sequences found in {args.input}")

    all_rows = []
    for i, record in enumerate(records):
        print(f"Submitting {record.id} to NCBI BLAST "
              f"({args.program} vs {args.database})...", file=sys.stderr)
        blast_record = run_blast(record, args)
        rows = collect_hits(record.id, blast_record, len(record.seq))
        all_rows.extend(rows)
        if rows:
            top = rows[0]
            print(f"  Closest match: {top['hit_accession']} "
                  f"{top['hit_description']} "
                  f"({top['percent_identity']}% identity, "
                  f"{top['query_coverage_pct']}% coverage)", file=sys.stderr)
        else:
            print("  No hits returned.", file=sys.stderr)
        if i < len(records) - 1:
            time.sleep(args.delay)

    if not all_rows:
        sys.exit("No hits found for any query sequence.")

    with open(args.output, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(all_rows[0].keys()))
        writer.writeheader()
        writer.writerows(all_rows)
    print(f"Wrote {len(all_rows)} hit(s) to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
