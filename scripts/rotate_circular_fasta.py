#!/usr/bin/env python3
"""Helper for reference-based consensus building on circular plasmids.

A plasmid reference is a linearized representation of a circular molecule.
When ONT reads are mapped to it with a linear aligner, any read whose
insert spans the artificial start/end junction of the reference gets
soft-clipped (or split) instead of aligning as one contiguous read. That
locally depresses coverage and consensus quality right at the junction,
which is exactly the kind of "invisible" assembly error you don't want in
a plasmid backbone.

The standard workaround is reference doubling: append the first PAD bases
of the reference to its own end before mapping, so a read crossing the
junction still aligns contiguously (into the appended copy). After
consensus calling on that extended reference, the true single-copy
circular sequence is reconstructed by taking, for each region, whichever
copy (original vs appended) actually had contiguous read support across
the junction, and discarding the other.

Two subcommands implement the two ends of this workflow:

  extend     ref.fasta -> ref.fasta + first PAD bp appended, plus a
             sidecar .meta file recording the original length and pad
             (consumed by `reconcile`).
  reconcile  the consensus called on the extended reference -> the
             final single-copy circular consensus, using the sidecar
             .meta file to know where the appended copy starts.

Both subcommands assume the input FASTA contains exactly one record
(a single plasmid reference/consensus), which is the normal case for
this pipeline. If your reference file holds more than one sequence,
split it first.
"""

import argparse
import sys

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def read_single_record(path):
    records = list(SeqIO.parse(path, "fasta"))
    if len(records) != 1:
        sys.exit(
            f"Error: {path} must contain exactly one FASTA record for "
            f"circular handling (found {len(records)})."
        )
    return records[0]


def write_meta(meta_path, original_length, pad):
    with open(meta_path, "w") as fh:
        fh.write(f"original_length={original_length}\n")
        fh.write(f"pad={pad}\n")


def read_meta(meta_path):
    values = {}
    with open(meta_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or "=" not in line:
                continue
            key, val = line.split("=", 1)
            values[key] = int(val)
    for key in ("original_length", "pad"):
        if key not in values:
            sys.exit(f"Error: {meta_path} is missing '{key}'.")
    return values["original_length"], values["pad"]


def cmd_extend(args):
    record = read_single_record(args.input)
    seq = str(record.seq)
    length = len(seq)
    if args.pad <= 0:
        sys.exit("Error: --pad must be a positive integer.")
    if args.pad >= length:
        sys.exit(
            f"Error: --pad ({args.pad}) must be smaller than the reference "
            f"length ({length}) -- it should only be as long as the reads "
            f"you expect to span the junction, not the whole plasmid."
        )
    extended_seq = seq + seq[: args.pad]
    out_record = SeqRecord(
        Seq(extended_seq),
        id=record.id,
        description=f"{record.description} [circular-extended pad={args.pad}]",
    )
    SeqIO.write([out_record], args.output, "fasta")

    meta_path = args.meta or f"{args.output}.meta"
    write_meta(meta_path, length, args.pad)

    print(
        f"Wrote extended reference ({length} + {args.pad} = "
        f"{len(extended_seq)} bp) to {args.output}",
        file=sys.stderr,
    )
    print(f"Wrote metadata to {meta_path}", file=sys.stderr)


def cmd_reconcile(args):
    original_length, pad = read_meta(args.meta)
    record = read_single_record(args.input)
    seq = str(record.seq)
    consensus_length = len(seq)

    expected = original_length + pad
    drift = consensus_length - expected
    if consensus_length <= original_length:
        sys.exit(
            f"Error: consensus length ({consensus_length} bp) is not "
            f"longer than the original reference ({original_length} bp). "
            "This extended consensus does not look like it was built from "
            "the doubled reference -- check the pipeline steps before "
            "trusting this output."
        )
    if abs(drift) > 0.01 * expected:
        print(
            f"Warning: extended consensus length ({consensus_length} bp) "
            f"differs from the expected {expected} bp (original "
            f"{original_length} + pad {pad}) by {drift} bp -- more than "
            "1%. This can happen with real indels near the junction, but "
            "double-check the junction region of the output by eye before "
            "trusting it.",
            file=sys.stderr,
        )
    elif drift != 0:
        print(
            f"Note: extended consensus length differs from the expected "
            f"{expected} bp by {drift} bp (small indel(s) called near the "
            "junction or elsewhere).",
            file=sys.stderr,
        )

    # The appended copy (originally ref[0:pad]) carries the consensus for
    # positions that reads crossing the junction actually supported
    # contiguously; the rest comes from the untouched body of the
    # extended reference. Anchoring on the original coordinates (not a
    # length-scaled split) is only approximate if indels shifted the
    # junction itself, hence the drift check/warning above.
    junction_supported_head = seq[original_length:]
    body = seq[pad:original_length]
    final_seq = junction_supported_head + body

    out_record = SeqRecord(
        Seq(final_seq),
        id=record.id.replace("_extended", ""),
        description="circularized reference-based consensus",
    )
    SeqIO.write([out_record], args.output, "fasta")
    print(
        f"Wrote circularized consensus ({len(final_seq)} bp, reference was "
        f"{original_length} bp) to {args.output}",
        file=sys.stderr,
    )


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)

    p_extend = sub.add_parser("extend", help="Build a doubled reference for mapping.")
    p_extend.add_argument("-i", "--input", required=True, help="Reference FASTA (one record).")
    p_extend.add_argument("-o", "--output", required=True, help="Path for the extended reference FASTA.")
    p_extend.add_argument("--pad", type=int, default=2000,
                           help="Bases to duplicate from the start of the reference onto its end "
                                "(default: 2000). Should be at least as long as the reads you expect "
                                "to span the plasmid's start/end junction -- check your read length "
                                "distribution and raise this if needed.")
    p_extend.add_argument("--meta", default=None,
                           help="Path for the sidecar metadata file (default: <output>.meta).")
    p_extend.set_defaults(func=cmd_extend)

    p_reconcile = sub.add_parser("reconcile", help="Rebuild the single-copy circular consensus.")
    p_reconcile.add_argument("-i", "--input", required=True,
                              help="Consensus FASTA called against the extended reference (one record).")
    p_reconcile.add_argument("-o", "--output", required=True, help="Path for the final circular consensus FASTA.")
    p_reconcile.add_argument("--meta", required=True, help="Metadata file produced by `extend`.")
    p_reconcile.set_defaults(func=cmd_reconcile)

    return parser.parse_args()


def main():
    args = parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
