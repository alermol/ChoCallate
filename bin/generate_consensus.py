#!/usr/bin/env python3

"""Generate the consensus BCF of one sample from its caller BCFs.

Every position of the (bgzipped, tabixed) coverage BED is written to the output,
so the consensus contains all sites, monomorphic ones included.  At each
position the callers vote with their genotype; the most frequent call wins if at
least ``--consensus_threshold`` callers support it.  Positions without a caller
record - and hom-reference calls of any ploidy - vote for the reference allele
with the ALT ``.`` placeholder.  Records with a missing genotype (``./.``,
``0/.``) are not evidence of anything and abstain from the vote.  FORMAT/AD and
FORMAT/DP are recovered from the sample BAM: overlapping mate pairs count once
(the observation with the highest base quality), and DP counts every covering
read except reference skips, matching ``samtools depth -J`` used for the
coverage BED.

Performance notes
-----------------
The output has one record per BED position, so the per-position cost decides the
runtime:

* Genotypes, pileups and records are only touched where they are needed.
  Positions without a caller record (the vast majority) take a fast path that
  writes the placeholder record directly, without building genotype tuples or a
  ``Counter``.
* FORMAT/AD support counts come from pysam's bulk pileup accessors
  (``get_query_names``/``get_query_sequences``) and overlapping mate pairs are
  only collapsed when a column actually repeats a query name.  If the optional
  Cython extension ``consensus_accel`` is built (``bin/build_consensus_accel.py``)
  the pileup is driven and tallied in C and the records are written through
  htslib, which together is another ~5x faster.
* BED intervals are grouped into spans (gaps up to ``SPAN_MAX_GAP`` bases) that
  share one pileup iterator and one caller-record scan.  A BAM region query
  costs ~2 ms because htslib decompresses a whole index bin, which dominates the
  runtime on fragmented coverage BEDs.

Set ``CHOCALLATE_NO_ACCEL=1`` in the environment to force the pure-pysam path
even when the Cython extension is importable (used by the test suite to check
that both implementations agree).
"""

import argparse
import os
from collections import namedtuple
from contextlib import ExitStack

import pysam

if os.environ.get('CHOCALLATE_NO_ACCEL'):
    BcfLineWriter = None
    PileupTally = None
else:
    try:
        # Optional Cython helpers built by bin/build_consensus_accel.py.
        from consensus_accel import BcfLineWriter, PileupTally
    except ImportError:  # pragma: no cover - the pure-pysam fallback is equivalent
        BcfLineWriter = None
        PileupTally = None

PileupSite = namedtuple('PileupSite', ['snp', 'ins', 'dele', 'dp'])


# BED intervals separated by at most this many bases share a pileup iterator.
# Bridging a gap costs roughly 0.2us per covered column (advanced in C when the
# accelerator is built), a new region query ~1.6ms, so gaps of several kb are
# much cheaper to walk through than to re-query.
SPAN_MAX_GAP = 10000

# Sentinel for "no caller record left in this span".
_NO_RECORD = 1 << 60

# Pre-encoded alleles for the monomorphic fast path (REF + ALT='.').
_MONOMORPHIC_ALLELES = {base: (base + ',.').encode('ascii') for base in 'ACGTN'}


def batched_intervals(bed_file, max_gap):
    """Group BED intervals into spans that share one pileup iterator.

    Args:
        bed_file: pysam.TabixFile of the (bgzipped) coverage BED.
        max_gap: Maximum gap between intervals of the same span.

    Yields:
        Tuples of (contig, [(start, end), ...]) where start/end are the 0-based
        half-open BED coordinates.
    """
    contig = None
    span = []
    span_end = None
    for line in bed_file.fetch():
        chrom, start, end = line.split('\t')
        start = int(start)
        end = int(end)
        if span and (chrom != contig or start - span_end > max_gap):
            yield contig, span
            span = []
        span.append((start, end))
        contig = chrom
        span_end = end
    if span:
        yield contig, span


def genotype_key(record, ref_base):
    """Consensus voting key of one caller record.

    Hom-reference calls (any ploidy) are mapped onto the position's reference
    placeholder, exactly like positions without a record, so they vote for the
    reference allele.  Records with a missing genotype (``./.``, ``0/.``) are
    not evidence of anything: they abstain from the vote.

    Args:
        record: pysam.VariantRecord of the caller.
        ref_base: Reference base at the position, upper case.

    Returns:
        (reference allele, alternate alleles, genotype) tuple, or ``None`` when
        the record has no genotyped call and abstains.
    """
    gt = record.samples[0].get('GT')
    if gt is None or any(allele is None for allele in gt):
        return None
    gt = tuple(sorted(gt))
    if all(allele == 0 for allele in gt):
        # any ploidy: (0,), (0, 0), ... all vote on the reference placeholder
        return (ref_base, ('.',), (0, 0))
    return (record.ref.upper(), tuple(sorted(record.alts)), gt)


def consensus_genotype(keys, consensus_threshold):
    """Most frequent caller genotype, or None if it is below the threshold.

    Ties are broken towards the caller that comes first, like
    ``collections.Counter.most_common(1)`` did.  ``keys`` are the voting keys in
    caller order, one entry per caller, placeholders included; ``None`` entries
    are callers that abstained and cast no vote.

    Args:
        keys: List of (ref, alt, gt) tuples and ``None`` for abstentions.
        consensus_threshold: Minimum number of matching calls.

    Returns:
        The winning (ref, alt, gt) tuple, or None.
    """
    counts = {}
    for index, key in enumerate(keys):
        if key is None:
            continue
        entry = counts.get(key)
        if entry is None:
            counts[key] = [1, index]
        else:
            entry[0] += 1
    best = None
    best_count = 0
    for key, entry in counts.items():
        if entry[0] > best_count:
            best = key
            best_count = entry[0]
    if best_count < consensus_threshold:
        return None
    return best


def site_from_entries(names, seqs, quals, dp):
    """Per-allele counters from materialised pileup column entries.

    Overlapping PE mates are deduplicated for allele counts (one observation per
    query name, keeping the highest base quality).

    Args:
        names: Query names at the column.
        seqs: Bases with samtools-style indel annotations ('A+2TT', 'A-2',
            '*', '>', '<').
        quals: Base qualities (255 means missing).
        dp: Depth of the column.

    Returns:
        PileupSite(snp Counter, insertion Counter, deletion-length Counter, depth).
    """
    best = {}
    for i, name in enumerate(names):
        seq = seqs[i]
        if seq == '>' or seq == '<' or seq == '*':
            # reference skip or deleted base: no allele at this position
            continue
        qual = quals[i]
        if qual == 255:
            # htslib's missing-quality sentinel; the original code used 0
            qual = 0
        prev = best.get(name)
        if prev is None or qual > prev[0]:
            best[name] = (qual, seq)
    snp = {}
    ins = {}
    dele = {}
    for _, allele in best.values():
        if '+' in allele:
            sequence_start = allele.index('+') + 1
            length_end = sequence_start
            while length_end < len(allele) and allele[length_end].isdigit():
                length_end += 1
            inserted = allele[length_end:].upper()
            ins[inserted] = ins.get(inserted, 0) + 1
        elif '-' in allele:
            length_start = allele.index('-') + 1
            length_end = length_start
            while length_end < len(allele) and allele[length_end].isdigit():
                length_end += 1
            length = int(allele[length_start:length_end])
            dele[length] = dele.get(length, 0) + 1
        else:
            base = allele.upper()
            snp[base] = snp.get(base, 0) + 1
    return PileupSite(snp=snp, ins=ins, dele=dele, dp=dp)


def build_pileup_site(column):
    """Summarise a single pileup column into per-allele counters.

    DP counts every covering read so it matches samtools depth used for the
    coverage BED.  The per-read data comes from pysam's bulk accessors, which
    are much cheaper than iterating PileupRead/AlignedSegment objects.

    Args:
        column: pysam.PileupColumn.

    Returns:
        PileupSite(snp Counter, insertion Counter, deletion-length Counter, depth).
    """
    names = column.get_query_names()
    seqs = column.get_query_sequences(add_indels=True)
    # reference skips are not part of DP
    dp = len(seqs) - seqs.count('>') - seqs.count('<')
    return site_from_entries(names, seqs, column.get_query_qualities(), dp)


def count_allele_depth(site, ref, allele):
    """
    Count reads supporting a single allele at a pileup site.

    Args:
        site: PileupSite or None.
        ref: Reference allele string.
        allele: Allele string to count ('.' is a non-allele placeholder for hom-ref).

    Returns:
        Integer read count supporting the allele.
    """
    if site is None:
        return 0
    # ALT=. is not a real allele; keep a 0 placeholder so Number=R matches pysam's
    # alleles=[REF, '.'] representation. Merge drops this slot and fill_missing_ad.py
    # writes 0 for the real ALT(s) afterward.
    if allele == '.':
        return 0
    if allele == ref:
        return site.snp.get(ref[0], 0)
    if len(ref) == 1 and len(allele) == 1:
        return site.snp.get(allele, 0)
    if len(allele) > len(ref) and allele.startswith(ref):
        return site.ins.get(allele[len(ref):], 0)
    if len(allele) < len(ref) and ref.startswith(allele):
        return site.dele.get(len(ref) - len(allele), 0)
    if len(allele) == len(ref):
        return site.snp.get(allele[0], 0)
    return 0


def recover_ad_dp(site, alleles):
    """
    Recover FORMAT/AD and FORMAT/DP for consensus alleles from a pileup site.

    Args:
        site: PileupSite for the position, or None for a position without coverage.
        alleles: List/tuple of alleles (REF followed by ALTs).

    Returns:
        Tuple of (AD tuple, DP int).
    """
    ref = alleles[0]
    ad = tuple(count_allele_depth(site, ref, allele) for allele in alleles)
    dp = site.dp if site is not None else 0
    return ad, dp


def reference_support(column, ref_base):
    """Reference-base support and depth of a pileup column.

    Equivalent to ``build_pileup_site`` followed by
    ``count_allele_depth(site, ref_base, ref_base)``, but avoids building the
    insertion/deletion counters and the allele tuples that the monomorphic fast
    path never needs.

    Args:
        column: pysam.PileupColumn.
        ref_base: Reference base at the position (upper case).

    Returns:
        (dp, number of collapsed reads supporting the reference base).
    """
    names = column.get_query_names()
    if not names:
        return 0, 0
    seqs = column.get_query_sequences(add_indels=True)
    # reference skips are not part of DP (see build_pileup_site)
    dp = len(seqs) - seqs.count('>') - seqs.count('<')
    if len(set(names)) == len(names):
        # no overlapping mate pair at this column: nothing to collapse
        return dp, seqs.count(ref_base) + seqs.count(ref_base.lower())
    return dp, _collapsed_reference_support(names, seqs, column.get_query_qualities(), ref_base)


def _collapsed_reference_support(names, seqs, quals, ref_base):
    best = {}
    for i, name in enumerate(names):
        seq = seqs[i]
        if seq == '>' or seq == '<' or seq == '*':
            continue
        qual = quals[i]
        if qual == 255:
            qual = 0
        prev = best.get(name)
        if prev is None or qual > prev[0]:
            best[name] = (qual, seq)
    ref_lower = ref_base.lower()
    support = 0
    for _, seq in best.values():
        if seq == ref_base or seq == ref_lower:
            support += 1
    return support


class ColumnWalker:
    """Walk a pysam pileup column by column (no Cython accelerator)."""

    __slots__ = ('pileup', 'column')

    def __init__(self, pileup):
        self.pileup = pileup
        self.column = next(pileup, None)

    def sync(self, vcf_start):
        """Advance to vcf_start; True if a pileup column is there."""
        column = self.column
        pileup = self.pileup
        while column is not None and column.reference_pos + 1 < vcf_start:
            column = next(pileup, None)
        self.column = column
        return column is not None and column.reference_pos + 1 == vcf_start

    def reference_support(self, ref_base):
        # type: (str) -> tuple
        return reference_support(self.column, ref_base)

    def site(self):
        return build_pileup_site(self.column)


class AccelWalker:
    """Walk a pileup through the Cython tally.

    No PileupColumn object is created per reference position; columns without
    coverage are skipped inside C.
    """

    __slots__ = ('tally',)

    def __init__(self, tally):
        self.tally = tally

    def sync(self, vcf_start):
        return self.tally.sync(vcf_start)

    def reference_support(self, ref_base):
        counts = self.tally.simple(ref_base)
        if counts is None:
            # pathologically deep column: the C scratch table is too small
            names, seqs, quals, dp = self.tally.entries()
            return dp, _collapsed_reference_support(names, seqs, quals, ref_base)
        return counts

    def site(self):
        return site_from_entries(*self.tally.entries())


def write_monomorphic_record(output_file, contig, start0, ref_base, ad_ref, dp):
    """Write the placeholder record of a position without caller records."""
    new_record = output_file.new_record(contig=contig, start=start0, stop=start0 + 1,
                                        alleles=(ref_base, '.'), qual=None)
    sample = new_record.samples[0]
    sample['GT'] = (0, 0)
    sample['AD'] = (ad_ref, 0)
    sample['DP'] = dp
    output_file.write(new_record)


def write_consensus_record(output_file, contig, start0, stop0, alleles, gt, ad, dp):
    """Write a consensus record from winning caller genotypes."""
    new_record = output_file.new_record(contig=contig, start=start0, stop=stop0,
                                        alleles=alleles, qual=None)
    sample = new_record.samples[0]
    sample['GT'] = gt
    sample['AD'] = ad
    sample['DP'] = dp
    output_file.write(new_record)


def generate_mininimal_header(sample_name,
                              reference_file,
                              consensus_threshold,
                              split_multiallelic,
                              remove_invariant,
                              version):
    """
    Generate a minimal VCF header for the consensus file.

    Args:
        sample_name: String sample name.
        reference_file: pysam.FastaFile object.
        consensus_threshold: Integer threshold for the consensus genotype.
        split_multiallelic: Whether multiallelic sites were split.
        remove_invariant: Whether monomorphic sites were removed.
        version: String version of the ChoCallate pipeline.

    Returns:
        A pysam.VariantHeader object.
    """
    header = pysam.VariantHeader()
    header.add_sample(sample_name)
    header.add_line(f'##tool=ChoCallate {version}')
    header.add_line(f'##consensusThreshold={consensus_threshold}')
    header.add_line(f'##multialleleSplit={split_multiallelic}')
    header.add_line(f'##noMonomorphic={remove_invariant}')
    for ref in reference_file.references:
        header.add_line(f'##contig=<ID={ref},length={reference_file.get_reference_length(ref)}>')
    header.add_line('##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of samples with data">')
    header.add_line('##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">')
    header.add_line('##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes">')
    header.add_line('##INFO=<ID=AC_Hom,Number=A,Type=Integer,Description="Allele counts in homozygous genotypes">')
    header.add_line('##INFO=<ID=AC_Het,Number=A,Type=Integer,Description="Allele counts in heterozygous genotypes">')
    header.add_line('##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">')
    header.add_line('##INFO=<ID=MAF,Number=1,Type=Float,Description="Frequency of the second most common allele">')
    header.add_line('##INFO=<ID=TYPE,Number=.,Type=String,Description="Variant type">')
    header.add_line('##INFO=<ID=F_MISSING,Number=1,Type=Float,Description="Fraction of missing genotypes">')
    header.add_line('##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth summed across samples">')
    header.add_line('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    header.add_line('##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">')
    header.add_line('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">')
    return header


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True, nargs='+', help='Input BCF/VCF files')
    parser.add_argument('--bed', required=True, help='Input BED file (gzipped and tabixed)')
    parser.add_argument('--bam', required=True, help='Input BAM file (indexed)')
    parser.add_argument('--output', required=True, help='Output BCF file')
    parser.add_argument('--sample_name', required=True, help='Sample name')
    parser.add_argument('--reference', required=True, help='Reference genome assembly (fasta file and fai index)')
    parser.add_argument('--consensus_threshold', required=True, type=int, help='Number of matching calls for consensus')
    parser.add_argument('--version', required=True, type=str, help='Version of ChoCallate pipeline')
    parser.add_argument('--split_multiallelic', action='store_true')
    parser.add_argument('--remove_invariant', action='store_true')
    args = parser.parse_args()

    threshold = args.consensus_threshold
    use_accel = BcfLineWriter is not None

    with ExitStack() as stack:
        bed_file = stack.enter_context(pysam.TabixFile(args.bed, index=args.bed + '.csi'))
        reference_file = stack.enter_context(pysam.FastaFile(args.reference, filepath_index=args.reference + '.fai'))
        bam_file = stack.enter_context(pysam.AlignmentFile(args.bam, 'rb'))
        header = generate_mininimal_header(args.sample_name,
                                           reference_file,
                                           threshold,
                                           args.split_multiallelic,
                                           args.remove_invariant,
                                           args.version)
        contig_ids = {name: header.contigs[name].id for name in header.contigs}
        variant_files = [stack.enter_context(pysam.VariantFile(bcf, index_filename=bcf + '.csi'))
                         for bcf in args.input]
        n_files = len(variant_files)
        # with at least `threshold` callers every position has a consensus (the
        # placeholder wins where no caller has a record)
        write_all_positions = threshold <= n_files
        output_file = None
        writer = None
        if use_accel:
            writer = BcfLineWriter(header, args.output)
        else:
            output_file = stack.enter_context(pysam.VariantFile(args.output, mode='wb', header=header))

        try:
            for contig, intervals in batched_intervals(bed_file, SPAN_MAX_GAP):
                span_start = intervals[0][0]
                span_end = intervals[-1][1]
                rid = contig_ids[contig]
                # Caller records are read once per span (not per BED interval)
                # and walked with a per-file cursor; records between intervals
                # are skipped when the cursor catches up with the next position.
                records_per_file = [list(variant_file.fetch(contig, span_start, span_end))
                                    for variant_file in variant_files]
                indices = [0] * n_files
                min_next = _NO_RECORD
                for records in records_per_file:
                    if records and records[0].pos < min_next:
                        min_next = records[0].pos
                pileup = bam_file.pileup(contig, span_start, span_end,
                                         truncate=True, stepper='all',
                                         min_base_quality=0, max_depth=1000000)
                if use_accel:
                    walker = AccelWalker(PileupTally(pileup))
                else:
                    walker = ColumnWalker(pileup)
                sync = walker.sync
                reference_support_of = walker.reference_support
                site_of = walker.site

                for start0, end0 in intervals:
                    ref_seq = reference_file.fetch(contig, start0, end0).upper()
                    rel = 0
                    # BED 0-based half-open [start0, end0) covers 1-based positions
                    # start0+1 .. end0, so the walk stops after end0.
                    vcf_start = start0 + 1
                    while vcf_start <= end0:
                        # keep the pileup in step with the walk
                        has_column = sync(vcf_start)

                        if min_next < vcf_start:
                            # a multi-base consensus shift or the start of the
                            # next interval moved past some caller records
                            min_next = _NO_RECORD
                            for fi in range(n_files):
                                records = records_per_file[fi]
                                i = indices[fi]
                                n_records = len(records)
                                while i < n_records and records[i].pos < vcf_start:
                                    i += 1
                                indices[fi] = i
                                if i < n_records and records[i].pos < min_next:
                                    min_next = records[i].pos

                        if min_next != vcf_start:
                            # fast path: no caller record here, every caller
                            # votes for the placeholder reference allele
                            if write_all_positions:
                                ref_base = ref_seq[rel]
                                if has_column:
                                    dp, ad_ref = reference_support_of(ref_base)
                                else:
                                    dp = 0
                                    ad_ref = 0
                                if use_accel:
                                    alleles = _MONOMORPHIC_ALLELES.get(ref_base)
                                    if alleles is None:
                                        alleles = (ref_base + ',.').encode('ascii')
                                    writer.write_monomorphic(rid, vcf_start - 1,
                                                             alleles, ad_ref, dp)
                                else:
                                    write_monomorphic_record(output_file, contig,
                                                             vcf_start - 1, ref_base,
                                                             ad_ref, dp)
                            vcf_start += 1
                            rel += 1
                            continue

                        # some caller has a record here: collect the votes
                        ref_base = ref_seq[rel]
                        keys = [(ref_base, ('.',), (0, 0))] * n_files
                        for fi in range(n_files):
                            records = records_per_file[fi]
                            i = indices[fi]
                            if i < len(records) and records[i].pos == vcf_start:
                                keys[fi] = genotype_key(records[i], ref_base)
                                indices[fi] = i + 1
                        min_next = _NO_RECORD
                        for fi in range(n_files):
                            records = records_per_file[fi]
                            i = indices[fi]
                            if i < len(records) and records[i].pos < min_next:
                                min_next = records[i].pos

                        consensus = consensus_genotype(keys, threshold)
                        if consensus is None:
                            vcf_start += 1
                            rel += 1
                            continue

                        alleles = [consensus[0]] + list(consensus[1])
                        if has_column:
                            site = site_of()
                        else:
                            site = None
                        ad, dp = recover_ad_dp(site, alleles)
                        if use_accel:
                            writer.write_record(rid, vcf_start - 1, alleles,
                                                consensus[2], ad, dp)
                        else:
                            write_consensus_record(output_file, contig, vcf_start - 1,
                                                   vcf_start, alleles, consensus[2],
                                                   ad, dp)
                        shift = len(consensus[0])  # shift on reference genotype length
                        vcf_start += shift
                        rel += shift
        finally:
            if writer is not None:
                writer.close()


if __name__ == '__main__':
    main()
