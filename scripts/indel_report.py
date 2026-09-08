#!/usr/bin/env python3
"""Account for every base of a reference and of the reads mapped to it.

Reference-based consensus callers (bcftools included) only handle small
variants: they will not call a kilobase-scale deletion, and
`bcftools consensus` silently copies reference bases wherever nothing was
called -- so an excised cassette can vanish from the report and reappear
in the consensus FASTA as if it were still there. Sequence present in the
sample but absent from the reference is worse: it is soft-clipped away and
leaves no trace in either the VCF or the coverage percentage.

This script closes both gaps by reading the alignments directly:

  report  Walks every alignment record and reports (a) how much of the
          reference is actually covered by reads and what is missing,
          (b) large deletions relative to the reference, from long CIGAR
          D operations AND from split (supplementary) alignments, since
          minimap2 usually splits reads rather than emitting one long
          gap, (c) large insertions relative to the reference, with the
          inserted sequence pulled out of the reads, and (d) soft-clip
          clusters, which is what an insertion or a rearrangement looks
          like when it cannot be represented in the alignment at all.
          Also writes the low-coverage BED used to mask the consensus, so
          uncovered regions become N instead of silently inheriting
          reference bases.

  trim    Removes long N runs from a masked consensus, giving the
          sequence that the reads actually support (e.g. the product of
          an excision). The junction it creates is only as good as the
          reads spanning it -- verify it against the deletion breakpoints
          in the report before treating it as final.

Coordinates in the report are 1-based inclusive and refer to the
reference that the BAM was mapped against. If that reference was built by
`rotate_circular_fasta.py extend`, pass its .meta file to --circular-meta
so the report says where the duplicated region starts.
"""

import argparse
import os
import re
import shutil
import statistics
import subprocess
import sys
from collections import defaultdict

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")
FLAG_UNMAPPED = 0x4
FLAG_REVERSE = 0x10
FLAG_SECONDARY = 0x100
FLAG_SUPPLEMENTARY = 0x800


def parse_cigar(cigar):
    if cigar == "*":
        return []
    ops = [(int(length), op) for length, op in CIGAR_RE.findall(cigar)]
    consumed = sum(len(str(length)) + 1 for length, _ in ops)
    if consumed != len(cigar):
        raise ValueError(f"Could not fully parse CIGAR: {cigar}")
    return ops


def open_alignments(path, samtools_bin):
    """Yield SAM records as split field lists, from a SAM or BAM file."""
    if path.endswith(".sam"):
        with open(path) as fh:
            for line in fh:
                if not line.startswith("@"):
                    yield line.rstrip("\n").split("\t")
        return

    if not shutil.which(samtools_bin):
        sys.exit(f"Error: {samtools_bin} not found in PATH (needed to read {path}).")
    proc = subprocess.Popen(
        [samtools_bin, "view", path], stdout=subprocess.PIPE, text=True
    )
    try:
        for line in proc.stdout:
            if not line.startswith("@"):
                yield line.rstrip("\n").split("\t")
    finally:
        proc.stdout.close()
        if proc.wait() != 0:
            sys.exit(f"Error: {samtools_bin} view failed on {path}.")


def read_reference_lengths(path):
    lengths = {}
    for record in SeqIO.parse(path, "fasta"):
        lengths[record.id] = len(record.seq)
    if not lengths:
        sys.exit(f"Error: no sequences found in {path}")
    return lengths


def query_span(ops):
    """Full read length implied by a CIGAR, including clipped bases."""
    return sum(length for length, op in ops if op in "MIS=XH")


class AlignmentStats:
    def __init__(self, ref_lengths, args):
        self.args = args
        self.ref_lengths = ref_lengths
        # Difference arrays: O(1) per aligned block instead of per base.
        self.cov_diff = {name: [0] * (length + 1) for name, length in ref_lengths.items()}
        self.deletions = []
        self.insertions = []
        self.clips = []
        self.blocks = defaultdict(list)  # split reads only, keyed by (qname, rname)
        self.total_records = 0
        self.secondary = 0
        self.supplementary = 0
        self.unmapped_reads = 0
        self.low_mapq = 0
        self.primary_mapped = 0
        self.total_read_bases = 0
        self.total_aligned_bases = 0
        self.total_clipped_bases = 0
        self.insert_seq_budget = defaultdict(int)

    def consume(self, fields):
        self.total_records += 1
        flag = int(fields[1])
        if flag & FLAG_SECONDARY:
            self.secondary += 1
            return
        if flag & FLAG_UNMAPPED:
            self.unmapped_reads += 1
            if fields[9] != "*":
                self.total_read_bases += len(fields[9])
            return

        qname, rname, pos, mapq, cigar, seq = (
            fields[0], fields[2], int(fields[3]), int(fields[4]), fields[5], fields[9],
        )
        if mapq < self.args.min_mapq:
            self.low_mapq += 1
            return
        if rname not in self.ref_lengths:
            return

        ops = parse_cigar(cigar)
        if not ops:
            return
        is_supplementary = bool(flag & FLAG_SUPPLEMENTARY)
        if is_supplementary:
            self.supplementary += 1
        else:
            self.primary_mapped += 1
            # Count each read's length once, from its primary alignment.
            self.total_read_bases += query_span(ops)

        ref_pos = pos - 1
        query_pos = 0
        ref_start = ref_pos
        diff = self.cov_diff[rname]
        ref_len = self.ref_lengths[rname]

        for length, op in ops:
            if op in "M=X":
                start = min(ref_pos, ref_len)
                end = min(ref_pos + length, ref_len)
                if end > start:
                    diff[start] += 1
                    diff[end] -= 1
                self.total_aligned_bases += length
                ref_pos += length
                query_pos += length
            elif op in "DN":
                if length >= self.args.min_sv_size:
                    self.deletions.append({
                        "qname": qname, "rname": rname, "start": ref_pos,
                        "end": ref_pos + length, "length": length,
                        "evidence": "cigar_deletion",
                    })
                ref_pos += length
            elif op == "I":
                if length >= self.args.min_sv_size:
                    bucket = (rname, ref_pos // max(self.args.cluster_window, 1))
                    inserted = ""
                    if seq != "*" and self.insert_seq_budget[bucket] < self.args.max_insert_seqs:
                        inserted = seq[query_pos:query_pos + length]
                        self.insert_seq_budget[bucket] += 1
                    self.insertions.append({
                        "qname": qname, "rname": rname, "start": ref_pos,
                        "end": ref_pos, "length": length, "sequence": inserted,
                        "evidence": "cigar_insertion",
                    })
                query_pos += length
            elif op == "S":
                self.total_clipped_bases += length
                query_pos += length
            elif op == "H":
                pass
            elif op == "P":
                pass

        self.record_clips(qname, rname, ops, pos - 1, ref_pos)
        if is_supplementary or self.has_sa_tag(fields):
            self.record_block(qname, rname, flag, ops, pos - 1, ref_pos)

    def record_clips(self, qname, rname, ops, ref_start, ref_end):
        leading = ops[0]
        if leading[1] == "H" and len(ops) > 1:
            leading = ops[1]
        if leading[1] == "S" and leading[0] >= self.args.min_clip:
            self.clips.append({"rname": rname, "pos": ref_start, "side": "left",
                               "length": leading[0], "qname": qname})
        trailing = ops[-1]
        if trailing[1] == "H" and len(ops) > 1:
            trailing = ops[-2]
        if trailing[1] == "S" and trailing[0] >= self.args.min_clip:
            self.clips.append({"rname": rname, "pos": ref_end, "side": "right",
                               "length": trailing[0], "qname": qname})

    @staticmethod
    def has_sa_tag(fields):
        return any(field.startswith("SA:Z:") for field in fields[11:])

    def record_block(self, qname, rname, flag, ops, ref_start, ref_end):
        read_length = query_span(ops)
        lead = ops[0][0] if ops[0][1] in "SH" else 0
        tail = ops[-1][0] if ops[-1][1] in "SH" else 0
        if flag & FLAG_REVERSE:
            q_start, q_end = tail, read_length - lead
        else:
            q_start, q_end = lead, read_length - tail
        self.blocks[(qname, rname)].append({
            "ref_start": ref_start, "ref_end": ref_end,
            "q_start": q_start, "q_end": q_end,
            "reverse": bool(flag & FLAG_REVERSE),
        })

    def split_read_gaps(self):
        """Classify the gap between two aligned blocks of the same read.

        Comparing the reference-side gap with the read-side gap is what
        separates the cases: reference gap only means sequence is missing from
        the sample, read gap only means the sample carries sequence the
        reference lacks, and both at once is a replaced segment, which this
        does not try to size.
        """
        deletions, insertions, junctions = [], [], []
        min_sv = self.args.min_sv_size
        for (qname, rname), blocks in self.blocks.items():
            if len(blocks) < 2:
                continue
            blocks.sort(key=lambda b: b["ref_start"])
            for left, right in zip(blocks, blocks[1:]):
                ref_gap = right["ref_start"] - left["ref_end"]
                query_gap = right["q_start"] - left["q_end"]
                colinear = (
                    left["reverse"] == right["reverse"]
                    and left["q_start"] <= right["q_start"]
                )
                base = {"qname": qname, "rname": rname}
                if not colinear:
                    if ref_gap >= min_sv:
                        junctions.append({**base, "start": left["ref_end"],
                                          "end": right["ref_start"], "length": ref_gap,
                                          "evidence": "split_read_noncolinear"})
                elif ref_gap >= min_sv and query_gap < min_sv:
                    deletions.append({**base, "start": left["ref_end"],
                                      "end": right["ref_start"], "length": ref_gap,
                                      "evidence": "split_read"})
                elif query_gap >= min_sv and ref_gap < min_sv:
                    insertions.append({**base, "start": left["ref_end"],
                                       "end": left["ref_end"], "length": query_gap,
                                       "sequence": "", "evidence": "split_read"})
                elif ref_gap >= min_sv and query_gap >= min_sv:
                    junctions.append({**base, "start": left["ref_end"],
                                      "end": right["ref_start"], "length": ref_gap,
                                      "evidence": "split_read_replacement"})
        return deletions, insertions, junctions

    def coverage(self):
        result = {}
        for name, length in self.ref_lengths.items():
            diff = self.cov_diff[name]
            depths = [0] * length
            running = 0
            for i in range(length):
                running += diff[i]
                depths[i] = running
            result[name] = depths
        return result


def cluster_events(events, window):
    """Group events that agree on breakpoint and size within `window` bp."""
    clusters = []
    for event in sorted(events, key=lambda e: (e["rname"], e["start"], e["length"])):
        placed = False
        for cluster in clusters:
            if (cluster["rname"] == event["rname"]
                    and abs(cluster["start_repr"] - event["start"]) <= window
                    and abs(cluster["length_repr"] - event["length"]) <= window):
                cluster["events"].append(event)
                placed = True
                break
        if not placed:
            clusters.append({
                "rname": event["rname"], "start_repr": event["start"],
                "length_repr": event["length"], "events": [event],
            })

    summaries = []
    for cluster in clusters:
        events_in = cluster["events"]
        evidence = defaultdict(int)
        for event in events_in:
            evidence[event["evidence"]] += 1
        summaries.append({
            "rname": cluster["rname"],
            "start": int(statistics.median(e["start"] for e in events_in)),
            "end": int(statistics.median(e["end"] for e in events_in)),
            "length": int(statistics.median(e["length"] for e in events_in)),
            "reads": len({e["qname"] for e in events_in}),
            "evidence": ", ".join(f"{k}={v}" for k, v in sorted(evidence.items())),
            "events": events_in,
        })
    summaries.sort(key=lambda c: (-c["reads"], -c["length"]))
    return summaries


def cluster_clips(clips, window, min_support):
    grouped = defaultdict(list)
    for clip in clips:
        grouped[(clip["rname"], clip["side"], clip["pos"] // max(window, 1))].append(clip)
    summaries = []
    for (rname, side, _), members in grouped.items():
        if len(members) < min_support:
            continue
        summaries.append({
            "rname": rname, "side": side,
            "pos": int(statistics.median(c["pos"] for c in members)),
            "reads": len({c["qname"] for c in members}),
            "median_length": int(statistics.median(c["length"] for c in members)),
            "max_length": max(c["length"] for c in members),
        })
    summaries.sort(key=lambda c: -c["reads"])
    return summaries


def low_coverage_blocks(depths, min_depth, min_block):
    blocks = []
    start = None
    for i, depth in enumerate(depths):
        if depth < min_depth:
            if start is None:
                start = i
        elif start is not None:
            if i - start >= min_block:
                blocks.append((start, i))
            start = None
    if start is not None and len(depths) - start >= min_block:
        blocks.append((start, len(depths)))
    return blocks


def percent(part, total):
    return 100.0 * part / total if total else 0.0


def interior_depth(coverage, event):
    """Mean depth inside a candidate deletion: near 0 confirms it is real."""
    depths = coverage.get(event["rname"], [])
    window = depths[event["start"]:event["end"]]
    return sum(window) / len(window) if window else 0.0


def cmd_report(args):
    ref_lengths = read_reference_lengths(args.reference)
    stats = AlignmentStats(ref_lengths, args)
    for fields in open_alignments(args.alignments, args.samtools):
        stats.consume(fields)

    coverage = stats.coverage()
    split_del, split_ins, junctions = stats.split_read_gaps()

    deletions = [d for d in cluster_events(stats.deletions + split_del, args.cluster_window)
                 if d["reads"] >= args.min_support]
    rearrangements = [r for r in cluster_events(junctions, args.cluster_window)
                      if r["reads"] >= args.min_support]
    insertions = [i for i in cluster_events(stats.insertions + split_ins, args.cluster_window)
                  if i["reads"] >= args.min_support]
    clips = cluster_clips(stats.clips, args.cluster_window, args.min_support)

    circular = read_circular_meta(args.circular_meta)

    write_depth(args.out_prefix + ".depth.tsv", coverage)
    mask_path = args.out_prefix + ".low_coverage.bed"
    write_mask_bed(mask_path, coverage, args.mask_depth, args.mask_min_block)
    write_variants_tsv(args.out_prefix + ".structural_variants.tsv",
                       deletions, rearrangements, insertions, clips, coverage)
    n_written = write_inserted_sequences(
        args.out_prefix + ".inserted_sequences.fasta", insertions
    )
    report_path = args.out_prefix + ".indel_report.txt"
    write_report(report_path, stats, coverage, deletions, rearrangements, insertions,
                 clips, circular, n_written, args)

    with open(report_path) as fh:
        sys.stdout.write(fh.read())


def read_circular_meta(path):
    """Length/pad of a reference built by rotate_circular_fasta.py extend."""
    if not path or not os.path.exists(path):
        return None
    meta = {}
    with open(path) as fh:
        for line in fh:
            if "=" in line:
                key, value = line.strip().split("=", 1)
                meta[key] = value
    if "original_length" not in meta or "pad" not in meta:
        return None
    return {"original_length": int(meta["original_length"]), "pad": int(meta["pad"])}


def annotate_position(pos_1based, circular):
    """Flag coordinates that fall in the duplicated tail of an extended reference."""
    if circular and pos_1based > circular["original_length"]:
        return f"{pos_1based} [dup of {pos_1based - circular['original_length']}]"
    return str(pos_1based)


def biological_length(depths, circular):
    """Length to report percentages against: the real plasmid, not the padding."""
    if circular:
        return min(len(depths), circular["original_length"])
    return len(depths)


def write_depth(path, coverage):
    with open(path, "w") as fh:
        for name, depths in coverage.items():
            for i, depth in enumerate(depths, start=1):
                fh.write(f"{name}\t{i}\t{depth}\n")


def write_mask_bed(path, coverage, min_depth, min_block):
    with open(path, "w") as fh:
        for name, depths in coverage.items():
            for start, end in low_coverage_blocks(depths, min_depth, min_block):
                fh.write(f"{name}\t{start}\t{end}\n")


def write_variants_tsv(path, deletions, rearrangements, insertions, clips, coverage):
    with open(path, "w") as fh:
        fh.write("type\treference\tstart\tend\tlength\tsupporting_reads\t"
                 "mean_depth_inside\tevidence\n")
        for event in deletions:
            fh.write(f"deletion\t{event['rname']}\t{event['start'] + 1}\t{event['end']}\t"
                     f"{event['length']}\t{event['reads']}\t"
                     f"{interior_depth(coverage, event):.1f}\t{event['evidence']}\n")
        for event in rearrangements:
            fh.write(f"complex_junction\t{event['rname']}\t{event['start'] + 1}\t"
                     f"{event['end']}\t{event['length']}\t{event['reads']}\t"
                     f"{interior_depth(coverage, event):.1f}\t{event['evidence']}\n")
        for event in insertions:
            fh.write(f"insertion\t{event['rname']}\t{event['start'] + 1}\t{event['start'] + 1}\t"
                     f"{event['length']}\t{event['reads']}\tNA\t{event['evidence']}\n")
        for clip in clips:
            fh.write(f"softclip_{clip['side']}\t{clip['rname']}\t{clip['pos'] + 1}\t"
                     f"{clip['pos'] + 1}\t{clip['median_length']}\t{clip['reads']}\tNA\t"
                     f"max_clip={clip['max_length']}\n")


def write_inserted_sequences(path, insertions):
    records = []
    for i, event in enumerate(insertions, start=1):
        candidates = [e for e in event["events"] if e["sequence"]]
        if not candidates:
            continue
        target = event["length"]
        best = min(candidates, key=lambda e: abs(len(e["sequence"]) - target))
        records.append(SeqRecord(
            Seq(best["sequence"]),
            id=f"insertion_{i}_{event['rname']}_{event['start'] + 1}",
            description=(f"length={len(best['sequence'])} supporting_reads={event['reads']} "
                         f"raw_read_sequence_uncorrected"),
        ))
    if records:
        SeqIO.write(records, path, "fasta")
    return len(records)


def write_report(path, stats, coverage, deletions, rearrangements, insertions, clips,
                 circular, n_inserted_seqs, args):
    lines = []
    add = lines.append

    add("=" * 72)
    add("STRUCTURAL DIFFERENCE REPORT (reads vs reference)")
    add("=" * 72)
    add("")
    if circular:
        original = circular["original_length"]
        add(f"NOTE: reads were mapped to a circular-extended reference. The real "
            f"plasmid is {original} bp; positions {original + 1}-"
            f"{original + circular['pad']} are a duplicate of positions "
            f"1-{circular['pad']}, added so reads crossing the start/end junction "
            f"align contiguously. Percentages below are against the real "
            f"{original} bp, and coordinates in the duplicated tail are marked "
            f"[dup of N].")
        add("")

    max_depth = max((max(d) if d else 0) for d in coverage.values()) if coverage else 0
    if max_depth < args.min_support or max_depth < args.min_depth:
        add("!" * 72)
        add("WARNING: this dataset cannot satisfy the thresholds in use.")
        add(f"  Maximum depth anywhere on the reference: {max_depth}")
        if max_depth < args.min_support:
            add(f"  --min-support is {args.min_support}, so any event supported by fewer")
            add("  than that many reads is ABSENT from the tables below, even when the")
            add("  alignments show it. The coverage section still reports it.")
        if max_depth < args.min_depth:
            add(f"  --min-depth is {args.min_depth}, so the 'covered at >= {args.min_depth}x'")
            add("  figure will read 0% regardless of how well the reference is covered.")
        add("  If you are comparing a single assembled sequence against the reference")
        add("  rather than a read set, rerun with --min-support 1 --min-depth 1.")
        add("!" * 72)
        add("")

    add("--- Reference ---")
    total_ref = 0
    for name, depths in coverage.items():
        length = biological_length(depths, circular)
        suffix = f" (+{len(depths) - length} bp duplicated padding)" if len(depths) > length else ""
        add(f"  {name}: {length} bp{suffix}")
        total_ref += length
    add("")

    add("--- Read mapping ---")
    add(f"  Alignment records read:        {stats.total_records}")
    add(f"  Primary alignments (mapped):   {stats.primary_mapped}")
    add(f"  Supplementary alignments:      {stats.supplementary}  "
        "(one read split across several reference positions)")
    add(f"  Secondary alignments (skipped):{stats.secondary}")
    add(f"  Unmapped reads:                {stats.unmapped_reads}")
    if stats.low_mapq:
        add(f"  Excluded by MAPQ < {args.min_mapq}:        {stats.low_mapq}")
    add("")
    add(f"  Total read bases:              {stats.total_read_bases}")
    add(f"  Bases aligned to reference:    {stats.total_aligned_bases} "
        f"({percent(stats.total_aligned_bases, stats.total_read_bases):.2f}% of read bases)")
    add(f"  Soft-clipped bases:            {stats.total_clipped_bases} "
        f"({percent(stats.total_clipped_bases, stats.total_read_bases):.2f}% of read bases)")
    add("  A large clipped fraction means the sample carries sequence that is not")
    add("  in the reference; reference-based consensus cannot recover it.")
    add("")

    add("--- Reference coverage ---")
    covered = missing = 0
    for name, depths in coverage.items():
        n = biological_length(depths, circular)
        window = depths[:n]
        at_threshold = sum(1 for d in window if d >= args.min_depth)
        zero = sum(1 for d in window if d == 0)
        covered += at_threshold
        missing += zero
        mean_depth = sum(window) / n if n else 0.0
        add(f"  {name}")
        add(f"    Mean depth:                  {mean_depth:.1f}")
        add(f"    Covered at >= {args.min_depth}x:{'':<15}{at_threshold} bp "
            f"({percent(at_threshold, n):.2f}%)")
        add(f"    Zero coverage:               {zero} bp ({percent(zero, n):.2f}%)")
    add("")
    add(f"  ALIGNS:  {covered} of {total_ref} reference bp "
        f"({percent(covered, total_ref):.2f}%) are covered at >= {args.min_depth}x")
    add(f"  MISSING: {missing} of {total_ref} reference bp "
        f"({percent(missing, total_ref):.2f}%) have no read support at all")
    add("")

    add(f"--- Reference regions absent from the sample (depth < {args.mask_depth}, "
        f">= {args.mask_min_block} bp) ---")
    add("  These are masked with N in the consensus and removed from the")
    add("  .consensus_supported.fasta file.")
    any_block = False
    for name, depths in coverage.items():
        limit = biological_length(depths, circular)
        for start, end in low_coverage_blocks(depths[:limit], args.mask_depth, args.mask_min_block):
            any_block = True
            window = depths[start:end]
            add(f"  {name}:{annotate_position(start + 1, circular)}-{end}  "
                f"length={end - start} bp  mean_depth={sum(window) / len(window):.1f}")
    if not any_block:
        add("  (none)")
    add("")

    add(f"--- Deletions relative to the reference (>= {args.min_sv_size} bp, "
        f">= {args.min_support} reads) ---")
    if deletions:
        add(f"  {'position':<30}{'length':>10}{'reads':>8}{'depth_in':>10}  evidence")
        for event in deletions:
            position = (f"{event['rname']}:"
                        f"{annotate_position(event['start'] + 1, circular)}-{event['end']}")
            add(f"  {position:<30}{event['length']:>10}{event['reads']:>8}"
                f"{interior_depth(coverage, event):>10.1f}  {event['evidence']}")
        add("")
        add("  'depth_in' is the mean read depth inside the deleted interval. A real deletion")
        add("  has no reads there, so depth_in near 0 confirms it; a high depth_in means the")
        add("  reference sequence IS present in the sample and the event is an alignment")
        add("  artifact, which reads longer than a circular plasmid readily produce.")
        add("  'cigar_deletion' is a gap inside one alignment; 'split_read' is a read that")
        add("  aligned in two pieces either side of the gap, which is how minimap2 usually")
        add("  represents a kilobase-scale deletion. Both are evidence for the same event.")
    else:
        add("  (none)")
    add("")

    add(f"--- Complex or non-colinear junctions (>= {args.min_support} reads) ---")
    if rearrangements:
        add(f"  {'position':<30}{'span':>10}{'reads':>8}  evidence")
        for event in rearrangements:
            position = (f"{event['rname']}:"
                        f"{annotate_position(event['start'] + 1, circular)}-{event['end']}")
            add(f"  {position:<30}{event['length']:>10}{event['reads']:>8}  {event['evidence']}")
        add("")
        add("  'split_read_noncolinear': the two pieces of the read are not in the same order")
        add("  as in the reference. That is the expected signature of a read crossing the")
        add("  junction of a circular molecule, and also of a genuine rearrangement.")
        add("  'split_read_replacement': reference sequence is missing AND the sample carries")
        add("  sequence the reference lacks, at the same junction. Neither is counted as a")
        add("  plain deletion. On a circular plasmid, expect some of these.")
    else:
        add("  (none)")
    add("")

    add(f"--- Insertions relative to the reference (>= {args.min_sv_size} bp, "
        f">= {args.min_support} reads) ---")
    if insertions:
        add(f"  {'position':<30}{'length':>10}{'reads':>8}  evidence")
        for event in insertions:
            position = f"{event['rname']}:{annotate_position(event['start'] + 1, circular)}"
            add(f"  {position:<30}{event['length']:>10}{event['reads']:>8}  {event['evidence']}")
        add("")
        add(f"  {n_inserted_seqs} inserted sequence(s) written to the .inserted_sequences.fasta")
        add("  file. These are raw read sequences and still carry ONT errors -- use them to")
        add("  identify what was inserted (e.g. by BLAST), not as a final sequence.")
    else:
        add("  (none)")
    add("")

    add(f"--- Soft-clip clusters (>= {args.min_clip} bp, >= {args.min_support} reads) ---")
    if clips:
        add(f"  {'position':<30}{'side':>8}{'reads':>8}{'median':>10}{'max':>10}")
        for clip in clips:
            position = f"{clip['rname']}:{annotate_position(clip['pos'] + 1, circular)}"
            add(f"  {position:<30}{clip['side']:>8}{clip['reads']:>8}"
                f"{clip['median_length']:>10}{clip['max_length']:>10}")
        add("")
        add("  Many reads clipped at the same position mark a breakpoint: sequence that")
        add("  continues in the sample but not in the reference. Clip clusters at a deletion")
        add("  breakpoint are expected -- they are the same event seen from reads that were")
        add("  clipped instead of split.")
        if circular:
            add("  Clusters at the very first and very last position of the mapping reference")
            add("  are edge artifacts of the linear representation, not biology.")
    else:
        add("  (none)")
    add("")

    add("--- Caveats ---")
    add(f"  - Only bases with depth < {args.mask_depth} are masked and trimmed. Regions")
    add(f"    covered below the {args.min_depth}x confidence threshold but above that are")
    add("    still called from few reads: see the coverage percentages above before")
    add("    trusting them.")
    add("  - The trimmed consensus cuts where coverage is zero, so a few reference bases")
    add("    can survive at a deletion breakpoint wherever reads misalign into the deleted")
    add("    region. Take the junction from the deletion coordinates above, not from the")
    add("    length of the trimmed sequence.")
    add("  - This compares reads against the reference you supplied. It reports what")
    add("    differs; it cannot confirm that the rest of the construct is correct beyond")
    add("    the small variants in the VCF.")
    add("  - Inserted sequence is reported but never reconstructed: a reference-based")
    add("    consensus can only ever return reference-shaped sequence. Use a de novo")
    add("    assembler for the inserted part.")
    add("  - Breakpoint coordinates come from read alignments and can shift by a few bp")
    add("    in repeats or homopolymers. Confirm exact junctions by eye in a viewer (IGV).")
    add("")

    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")


def cmd_trim(args):
    records = list(SeqIO.parse(args.consensus, "fasta"))
    if not records:
        sys.exit(f"Error: no sequences found in {args.consensus}")

    trimmed_records = []
    total_removed = 0
    for record in records:
        seq = str(record.seq)
        pieces = []
        removed = 0
        index = 0
        for match in re.finditer(r"[Nn]{%d,}" % args.min_block, seq):
            pieces.append(seq[index:match.start()])
            removed += match.end() - match.start()
            index = match.end()
        pieces.append(seq[index:])
        new_seq = "".join(pieces)
        total_removed += removed
        trimmed_records.append(SeqRecord(
            Seq(new_seq), id=record.id,
            description=f"unsupported_blocks_removed={removed}bp length={len(new_seq)}",
        ))
        print(f"{record.id}: {len(seq)} bp -> {len(new_seq)} bp "
              f"({removed} bp of unsupported sequence removed)", file=sys.stderr)

    SeqIO.write(trimmed_records, args.output, "fasta")
    if total_removed:
        print("Note: the junction(s) created by removing unsupported blocks are only as "
              "reliable as the reads spanning them -- check the deletion breakpoints in "
              "the indel report before treating this as a final sequence.", file=sys.stderr)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)

    p_report = sub.add_parser("report", help="Analyze alignments for structural differences.")
    p_report.add_argument("-a", "--alignments", required=True,
                          help="Sorted BAM (or SAM) of reads mapped to the reference.")
    p_report.add_argument("-r", "--reference", required=True,
                          help="Reference FASTA the reads were mapped against.")
    p_report.add_argument("-o", "--out-prefix", required=True, help="Prefix for output files.")
    p_report.add_argument("--min-sv-size", type=int, default=50,
                          help="Smallest deletion/insertion to report, bp (default: 50).")
    p_report.add_argument("--min-clip", type=int, default=100,
                          help="Smallest soft-clip to report, bp (default: 100).")
    p_report.add_argument("--min-depth", type=int, default=10,
                          help="Depth at which a reference base counts as covered "
                               "(default: 10) -- provisional, review against your run.")
    p_report.add_argument("--mask-depth", type=int, default=1,
                          help="Depth below which a base is treated as having no read "
                               "support at all and is masked with N (default: 1, i.e. only "
                               "zero-coverage bases). Kept separate from --min-depth on "
                               "purpose: --min-depth is a confidence threshold for "
                               "statistics, while masking (and the trimming that follows "
                               "it) removes sequence, so it must only fire where there is "
                               "genuinely no evidence.")
    p_report.add_argument("--mask-min-block", type=int, default=20,
                          help="Smallest unsupported block to mask/report, bp (default: 20).")
    p_report.add_argument("--min-support", type=int, default=2,
                          help="Reads required to report an event (default: 2).")
    p_report.add_argument("--cluster-window", type=int, default=100,
                          help="Breakpoint/size tolerance when grouping events, bp "
                               "(default: 100).")
    p_report.add_argument("--min-mapq", type=int, default=0,
                          help="Skip alignments below this MAPQ (default: 0, keep all).")
    p_report.add_argument("--max-insert-seqs", type=int, default=50,
                          help="Inserted read sequences to keep per site (default: 50).")
    p_report.add_argument("--circular-meta", default=None,
                          help="Optional .meta file from rotate_circular_fasta.py extend, "
                               "to annotate coordinates in the duplicated region.")
    p_report.add_argument("--samtools", default="samtools", help="samtools binary (default: samtools).")
    p_report.set_defaults(func=cmd_report)

    p_trim = sub.add_parser("trim", help="Drop long N runs from a masked consensus.")
    p_trim.add_argument("-c", "--consensus", required=True, help="Masked consensus FASTA.")
    p_trim.add_argument("-o", "--output", required=True, help="Path for the trimmed FASTA.")
    p_trim.add_argument("--min-block", type=int, default=20,
                        help="Smallest N run to remove, bp (default: 20).")
    p_trim.set_defaults(func=cmd_trim)

    return parser.parse_args()


def main():
    args = parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
