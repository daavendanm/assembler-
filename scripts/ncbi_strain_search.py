#!/usr/bin/env python3
"""Search NCBI (BLASTN) for public genomes matching a reference-based
consensus sequence, to identify the closest known strain.

NCBI's web BLAST rejects any single query over 1,000,000 bases, so a
whole bacterial genome (a few Mb) is split into chunks below that limit
and each chunk is submitted separately. This makes a whole-genome search
slow (one NCBI queue wait per chunk); for genome-scale comparison a local
BLAST+ database or FastANI/Mash against downloaded genomes is faster and
more appropriate. This script remains a reasonable option when you want
results straight from NCBI's own database with no local setup.
"""

import argparse
import csv
import os
import sys
import time

from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

NCBI_MAX_QUERY_LENGTH = 1_000_000


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
                              "(default: 10)")
    parser.add_argument("--chunk-size", type=int, default=900_000,
                         help="Split any sequence longer than this into "
                              "consecutive chunks before submitting to NCBI "
                              "(default: 900000; NCBI's hard limit is 1000000)")
    return parser.parse_args()


def chunk_record(record, chunk_size):
    seq = str(record.seq)
    if len(seq) <= chunk_size:
        yield record
        return
    for start in range(0, len(seq), chunk_size):
        end = min(start + chunk_size, len(seq))
        sub_seq = Seq(seq[start:end])
        yield SeqRecord(sub_seq, id=f"{record.id}_{start + 1}-{end}", description="")


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

    queries = [chunk for record in records for chunk in chunk_record(record, args.chunk_size)]
    if len(queries) > len(records):
        print(f"Input split into {len(queries)} chunks of up to "
              f"{args.chunk_size} bases (NCBI's limit is "
              f"{NCBI_MAX_QUERY_LENGTH}). This will take a while: each "
              f"chunk is a separate NCBI submission.", file=sys.stderr)

    fieldnames = ["query_id", "hit_accession", "hit_description",
                  "percent_identity", "query_coverage_pct", "evalue"]
    out_dir = os.path.dirname(args.output)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    total_rows = 0
    with open(args.output, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for i, query in enumerate(queries):
            print(f"Submitting {query.id} to NCBI BLAST "
                  f"({args.program} vs {args.database})...", file=sys.stderr)
            blast_record = run_blast(query, args)
            rows = collect_hits(query.id, blast_record, len(query.seq))
            writer.writerows(rows)
            fh.flush()
            total_rows += len(rows)
            if rows:
                top = rows[0]
                print(f"  Closest match: {top['hit_accession']} "
                      f"{top['hit_description']} "
                      f"({top['percent_identity']}% identity, "
                      f"{top['query_coverage_pct']}% coverage)", file=sys.stderr)
            else:
                print("  No hits returned.", file=sys.stderr)
            if i < len(queries) - 1:
                time.sleep(args.delay)

    if total_rows == 0:
        sys.exit("No hits found for any query sequence.")
    print(f"Wrote {total_rows} hit(s) to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
