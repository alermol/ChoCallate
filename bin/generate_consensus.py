#!/usr/bin/env python3

"""Generate the consensus BCF of one sample from its caller BCFs.

Every position of the (bgzipped, tabixed) coverage BED is written to the output,
so the consensus contains all sites, monomorphic ones included.  At each
position the callers vote with their genotype in two stages: an ALT allele must
first be supported by at least ``--consensus_threshold`` callers, then the
majority zygosity among the supporters of the winning ALT decides between
heterozygous and homozygous ALT.  Ties are resolved by ``--tie_policy``:
``nocall`` (default, independent of caller order), ``first`` (the earliest
supporter in configured caller order) or ``qual`` (largest summed QUAL).  A
position without ALT support is hom-reference when at least
``--consensus_threshold`` callers vote reference, otherwise it is a no-call.

A caller without a record votes reference by default (``--absent_vote ref``).
With ``--absent_vote abstain`` it only votes reference when the sample pileup
has at least ``--min_coverage`` reads at the position; a caller record with a
missing genotype (``./.``, ``0/.``) always abstains.  Hom-reference calls of any
ploidy vote reference.  Positions without a consensus variant are written as
gVCF-style records with the symbolic ALT ``<*>``: the reference allele plus
"any other allele", whose FORMAT/AD slot carries the depth of reads that do not
support the reference base.  ``bcftools merge --gvcf`` maps that slot onto the
real ALT allele of another sample's variant, so merged hom-reference samples
report their actual non-reference depth instead of a fabricated 0.

FORMAT/AD and FORMAT/DP are recovered from the sample BAM: overlapping mate
pairs count once (the observation with the highest base quality), and DP counts
every covering read except reference skips, matching ``samtools depth -J`` used
for the coverage BED.

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
import sys
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

# Voting key of a reference call (no ALT allele).
_REFERENCE_ALTS = ('.',)

# Pre-encoded alleles for the monomorphic fast path (REF + ALT='<*>').
_MONOMORPHIC_ALLELES = {base: (base + ',<*>').encode('ascii') for base in 'ACGTN'}


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


def genotype_key(records, ref_base):
    """Consensus voting key of one caller at a position.

    ``records`` are all records of the caller at the position: normally one,
    but several after multi-allelic splitting or allele atomisation.  Their ALT
    alleles are merged into one sorted allele set and the genotype indices are
    remapped to that order.  Hom-reference calls (any ploidy) become the
    reference placeholder; records with a missing genotype abstain.

    Args:
        records: Non-empty list of pysam.VariantRecord objects at one position.
        ref_base: Reference base at the position, upper case.

    Returns:
        (reference allele, alternate alleles, genotype) tuple, or ``None`` when
        the records do not form a usable vote.
    """
    alt_copies = {}
    ploidy = 0
    variant_ref = None
    for record in records:
        gt = record.samples[0].get('GT')
        if gt is None or any(allele is None for allele in gt):
            continue
        alts = tuple(allele.upper() for allele in (record.alts or ()))
        if not alts or alts == _REFERENCE_ALTS:
            ploidy = max(ploidy, len(gt))
            continue
        ref = record.ref.upper()
        if variant_ref is None:
            variant_ref = ref
        elif ref != variant_ref:
            # two records at the same position with different REF alleles
            # cannot be represented by one voting key
            return None
        ploidy = max(ploidy, len(gt))
        for allele_index, alt in enumerate(alts, 1):
            copies = sum(1 for allele in gt if allele == allele_index)
            if copies:
                alt_copies[alt] = max(alt_copies.get(alt, 0), copies)
    if variant_ref is None or not alt_copies:
        if ploidy:
            return (ref_base, _REFERENCE_ALTS, (0, 0))
        return None
    alts = tuple(sorted(alt_copies))
    total_copies = sum(alt_copies[alt] for alt in alts)
    if total_copies > ploidy:
        # inconsistent split genotypes (for example 1/1 plus 0/1 in a diploid
        # caller): abstain instead of inventing a genotype
        return None
    gt = [0] * (ploidy - total_copies)
    for allele_index, alt in enumerate(alts, 1):
        gt.extend([allele_index] * alt_copies[alt])
    return (variant_ref, alts, tuple(gt))


def _zygosity(key, alt):
    """Classify a vote supporting ``alt`` as ``het`` or ``hom`` ALT.

    Haploid ALT calls count as homozygous (one ALT copy out of one).  A
    multi-allelic or split genotype is reduced to the winning ALT allele.
    """
    _, alts, gt = key
    indices = {index + 1 for index, allele in enumerate(alts) if allele == alt}
    copies = sum(1 for allele in gt if allele in indices)
    if copies == 0:
        return None
    if copies == len(gt):
        return 'hom'
    return 'het'


def _genotype_for(zygosity, ploidy):
    """Output GT for a zygosity and ploidy (ALT index 1)."""
    if zygosity == 'hom':
        return (1,) * ploidy
    return (0,) * (ploidy - 1) + (1,)


def consensus_calls(keys, consensus_threshold, qualities=None,
                    tie_policy='nocall'):
    """Consensus records for one position.

    Stage 1 (site): every ALT allele supported by at least
    ``consensus_threshold`` callers becomes a consensus candidate.  Stage 2
    (genotype): for each candidate, the majority zygosity (heterozygous vs
    homozygous ALT) among its supporters wins.  A genotype tie is resolved by
    ``tie_policy``:

    ``nocall``
        Skip that ALT (default; independent of caller order).
    ``first``
        Use the earliest supporter in the configured caller order.
    ``qual``
        Use the largest summed QUAL; a remaining tie skips that ALT.

    Every qualifying ALT produces its own split consensus record; the output
    stage re-joins them unless split output was requested.  When no ALT reaches
    the threshold, the position is hom-reference if at least
    ``consensus_threshold`` reference votes remain, otherwise the list is empty.

    Args:
        keys: List of (ref, alt tuple, gt tuple) keys and ``None`` for
            abstentions, one entry per caller.
        consensus_threshold: Minimum number of callers supporting an ALT (or
            the reference) at the site stage.
        qualities: Per-caller QUAL values used by the ``qual`` tie policy.
        tie_policy: ``nocall``, ``first`` or ``qual``.

    Returns:
        A list of (ref, alt tuple, gt tuple) records; normally zero, one or one
        split record per qualifying ALT.
    """
    if qualities is None:
        qualities = [0.0] * len(keys)
    reference_key = None
    reference_count = 0
    alt_support = {}
    for index, key in enumerate(keys):
        if key is None:
            continue
        ref, alts, gt = key
        if not alts or alts == _REFERENCE_ALTS:
            reference_count += 1
            if reference_key is None:
                reference_key = (ref, _REFERENCE_ALTS, (0, 0))
            continue
        for allele_index in sorted({allele for allele in gt if allele > 0}):
            if allele_index <= len(alts):
                alt_support.setdefault(alts[allele_index - 1], []).append(
                    (index, key))

    candidates = sorted(
        (alt, supporters) for alt, supporters in alt_support.items()
        if len(supporters) >= consensus_threshold)
    if not candidates:
        if reference_count >= consensus_threshold:
            return [reference_key]
        return []

    calls = []
    for alt, supporters in candidates:
        classes = {}
        for index, key in supporters:
            zygosity = _zygosity(key, alt)
            if zygosity is not None:
                classes.setdefault(zygosity, []).append((index, key))
        if not classes:
            continue
        max_class = max(len(members) for members in classes.values())
        winners = [(zygosity, members) for zygosity, members in classes.items()
                   if len(members) == max_class]
        if len(winners) > 1:
            if tie_policy == 'first':
                winners = [min(winners,
                               key=lambda item: min(i for i, _ in item[1]))]
            elif tie_policy == 'qual':
                scored = [(sum(qualities[i] for i, _ in members), zygosity, members)
                          for zygosity, members in winners]
                best = max(score for score, _, _ in scored)
                winners = [(zygosity, members)
                           for score, zygosity, members in scored
                           if score == best]
            if len(winners) > 1:
                # unresolved genotype tie: this ALT gets no consensus record
                continue
        zygosity, members = winners[0]
        ploidy = max(len(key[2]) for _, key in members)
        ref = members[0][1][0]
        calls.append((ref, (alt,), _genotype_for(zygosity, ploidy)))
    return calls


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
    # ALT="." is not a real allele.  Hom-reference consensus is written as a
    # gVCF-style <*> record, so this is only a defensive fallback.
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
    """Reference-base support, observation count and depth of a pileup column.

    Equivalent to ``build_pileup_site`` followed by
    ``count_allele_depth(site, ref_base, ref_base)``, but avoids building the
    insertion/deletion counters and the allele tuples that the monomorphic fast
    path never needs.

    Args:
        column: pysam.PileupColumn.
        ref_base: Reference base at the position (upper case).

    Returns:
        (dp, number of collapsed reads supporting the reference base, number of
        collapsed observations that have an allele at this position).  The last
        value feeds the ``<*>`` AD slot; reference skips are not observations.
    """
    names = column.get_query_names()
    if not names:
        return 0, 0, 0
    seqs = column.get_query_sequences(add_indels=True)
    # reference skips are not part of DP (see build_pileup_site)
    dp = len(seqs) - seqs.count('>') - seqs.count('<')
    if len(set(names)) == len(names):
        # no overlapping mate pair at this column: nothing to collapse; every
        # read that is not a reference skip is one observation
        return dp, seqs.count(ref_base) + seqs.count(ref_base.lower()), dp
    support, observations = _collapsed_support(
        names, seqs, column.get_query_qualities(), ref_base)
    return dp, support, observations


def _collapsed_support(names, seqs, quals, ref_base):
    """Collapse overlapping mates and count reference/other observations.

    Args:
        names: Query names at the column.
        seqs: Bases with samtools-style indel annotations.
        quals: Base qualities (255 means missing).
        ref_base: Reference base at the position (upper case).

    Returns:
        (number of reads supporting the reference base, number of collapsed
        observations).  Reads with a base, an insertion or a deletion at the
        position are one observation each; reference skips are ignored.
    """
    best = {}
    for i, name in enumerate(names):
        seq = seqs[i]
        if seq == '>' or seq == '<':
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
    return support, len(best)


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
            support, observations = _collapsed_support(names, seqs, quals, ref_base)
            return dp, support, observations
        return counts

    def site(self):
        return site_from_entries(*self.tally.entries())


def write_monomorphic_record(output_file, contig, start0, ref_base, ad_ref, ad_other, dp):
    """Write the gVCF-style reference record of a position without a call.

    The record carries the symbolic ALT ``<*>``; the second AD slot holds the
    depth of reads that do not support the reference base, so that
    ``bcftools merge --gvcf`` reports honest ALT depth for merged samples that
    are hom-reference.
    """
    new_record = output_file.new_record(contig=contig, start=start0, stop=start0 + 1,
                                        alleles=(ref_base, '<*>'), qual=None)
    sample = new_record.samples[0]
    sample['GT'] = (0, 0)
    sample['AD'] = (ad_ref, ad_other)
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


def _header_value(value):
    """Escape a value for a quoted VCF meta line."""
    return (str(value).replace('\\', '\\\\')
            .replace('"', '\\"')
            .replace('\t', ' ')
            .replace('\r', ' ')
            .replace('\n', ' '))


def generate_mininimal_header(sample_name,
                              reference_file,
                              consensus_threshold,
                              split_multiallelic,
                              remove_invariant,
                              version,
                              callers,
                              callers_failed,
                              tie_policy,
                              absent_vote,
                              min_coverage,
                              command_line=None,
                              caller_versions=None):
    """
    Generate a minimal VCF header for the consensus file.

    Args:
        sample_name: String sample name.
        reference_file: pysam.FastaFile object.
        consensus_threshold: Integer threshold for the site-level vote.
        split_multiallelic: Whether multiallelic sites were split.
        remove_invariant: Whether monomorphic sites were removed.
        version: String version of the ChoCallate pipeline.
        callers: Callers that voted for this sample, in voting order.
        callers_failed: Configured callers without a usable record set.
        tie_policy: Tie policy applied by the voting model.
        absent_vote: Policy for callers without a record at a position.
        min_coverage: Coverage threshold used by ``absent_vote=abstain``.
        command_line: Full pipeline command line, when known.
        caller_versions: Comma-separated caller ``tool=version`` pairs.

    Returns:
        A pysam.VariantHeader object.
    """
    header = pysam.VariantHeader()
    header.add_sample(sample_name)
    meta = [f'Version={version}',
            f'Threshold={consensus_threshold}',
            f'Callers="{",".join(callers)}"',
            'Rule=two-stage',
            f'Tie={tie_policy}',
            f'AbsentVote={absent_vote}']
    if absent_vote == 'abstain':
        meta.append(f'MinCoverage={min_coverage}')
    meta.append(f'RemoveInvariant={remove_invariant}')
    meta.append(f'SplitMultiallelic={split_multiallelic}')
    if callers_failed:
        meta.append(f'CallersFailed="{",".join(callers_failed)}"')
    header.add_line('##ChoCallate=<' + ','.join(meta) + '>')
    if command_line:
        header.add_line(
            '##ChoCallateCommand="{}"'.format(
                _header_value(command_line)))
    if caller_versions:
        header.add_line(
            '##ChoCallateCallerVersions="{}"'.format(
                _header_value(caller_versions)))
    header.add_line('##ALT=<ID=*,Description="Represents allele(s) other than observed.">')
    header.add_line('##INFO=<ID=END,Number=1,Type=Integer,Description="Stop position of the interval">')
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


def resolve_inputs(input_specs, configured_callers):
    """Map the configured callers onto the caller BCF files that exist.

    Args:
        input_specs: ``--input`` values, either ``caller=path`` or a plain
            path (the caller name is then taken from the file stem).
        configured_callers: Comma-separated configured caller list, or None to
            derive it from ``input_specs``.

    Returns:
        (received callers in configured order, failed callers, paths).  A
        configured caller whose file is missing counts as failed; this is how
        a caller that was skipped in opt-in degraded mode is detected.
    """
    specs = []
    for spec in input_specs:
        caller, separator, path = spec.partition('=')
        if not separator:
            caller = os.path.basename(spec)
            if caller.endswith('.bcf'):
                caller = caller[:-4]
            path = spec
        specs.append((caller, path))
    if configured_callers:
        configured = [name.strip() for name in configured_callers.split(',')
                      if name.strip()]
    else:
        configured = [caller for caller, _ in specs]
    present = {caller: path for caller, path in specs if os.path.exists(path)}
    received = [caller for caller in configured if caller in present]
    failed = [caller for caller in configured if caller not in present]
    return received, failed, [present[caller] for caller in received]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True, action='append',
                        metavar='[CALLER=]PATH',
                        help='Caller BCF; repeat once per caller')
    parser.add_argument('--callers', default=None,
                        help='Configured callers in voting order (comma-separated)')
    parser.add_argument('--bed', required=True, help='Input BED file (gzipped and tabixed)')
    parser.add_argument('--bam', required=True, help='Input BAM file (indexed)')
    parser.add_argument('--output', required=True, help='Output BCF file')
    parser.add_argument('--sample_name', required=True, help='Sample name')
    parser.add_argument('--reference', required=True, help='Reference genome assembly (fasta file and fai index)')
    parser.add_argument('--consensus_threshold', required=True, type=int,
                        help='Minimum caller support at the site level')
    parser.add_argument('--version', required=True, type=str, help='Version of ChoCallate pipeline')
    parser.add_argument('--absent_vote', choices=('ref', 'abstain'), default='ref',
                        help='Vote of a caller without a record at a position')
    parser.add_argument('--tie_policy', choices=('nocall', 'first', 'qual'),
                        default='nocall', help='How to break voting ties')
    parser.add_argument('--min_coverage', type=int, default=0,
                        help='Coverage required for an absent caller to vote '
                             'reference when --absent_vote is abstain')
    parser.add_argument('--command_line_file', default=None,
                        help='File containing the pipeline command line')
    parser.add_argument('--caller_versions', default=None,
                        help='Caller provenance as tool=version pairs')
    parser.add_argument('--caller_versions_file', default=None,
                        help='File with tool=version lines')
    parser.add_argument('--split_multiallelic', action='store_true')
    parser.add_argument('--remove_invariant', action='store_true')
    args = parser.parse_args()

    threshold = args.consensus_threshold
    absent_vote = args.absent_vote
    tie_policy = args.tie_policy
    min_coverage = args.min_coverage
    caller_versions = args.caller_versions
    if args.caller_versions_file:
        with open(args.caller_versions_file, encoding='utf-8') as handle:
            caller_versions = ','.join(
                line.strip() for line in handle if line.strip())
    command_line = None
    if args.command_line_file:
        with open(args.command_line_file, encoding='utf-8') as handle:
            command_line = handle.read().strip()
    use_accel = BcfLineWriter is not None
    received, failed, input_paths = resolve_inputs(args.input, args.callers)
    if not received:
        sys.exit('ERROR: no caller BCF available for the configured callers')
    if failed:
        print('[WARN] sample {}: caller(s) without output: {}'.format(
            args.sample_name, ','.join(failed)), file=sys.stderr)
        if threshold > len(received):
            sys.exit(
                'ERROR: consensus threshold {} exceeds the {} surviving '
                'caller(s) ({}) after failure(s) of {}'.format(
                    threshold, len(received), ','.join(received),
                    ','.join(failed)))

    with ExitStack() as stack:
        bed_file = stack.enter_context(pysam.TabixFile(args.bed, index=args.bed + '.csi'))
        reference_file = stack.enter_context(pysam.FastaFile(args.reference, filepath_index=args.reference + '.fai'))
        bam_file = stack.enter_context(pysam.AlignmentFile(args.bam, 'rb'))
        header = generate_mininimal_header(args.sample_name,
                                           reference_file,
                                           threshold,
                                           args.split_multiallelic,
                                           args.remove_invariant,
                                           args.version,
                                           received,
                                           failed,
                                           tie_policy,
                                           absent_vote,
                                           min_coverage,
                                           command_line,
                                           caller_versions)
        contig_ids = {name: header.contigs[name].id for name in header.contigs}
        variant_files = [stack.enter_context(pysam.VariantFile(bcf, index_filename=bcf + '.csi'))
                         for bcf in input_paths]
        n_files = len(variant_files)
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
                            # fast path: no caller record here
                            ref_base = ref_seq[rel]
                            if has_column:
                                dp, ad_ref, observations = reference_support_of(ref_base)
                            else:
                                dp = 0
                                ad_ref = 0
                                observations = 0
                            # Absent callers vote reference by default.  In
                            # abstain mode they need at least min_coverage reads
                            # at the position; a position without a pileup
                            # column has no evidence.
                            if absent_vote == 'ref' or (
                                    has_column and dp >= min_coverage):
                                reference_votes = n_files
                            else:
                                reference_votes = 0
                            if reference_votes >= threshold:
                                ad_other = observations - ad_ref
                                if use_accel:
                                    alleles = _MONOMORPHIC_ALLELES.get(ref_base)
                                    if alleles is None:
                                        alleles = (ref_base + ',<*>').encode('ascii')
                                    writer.write_monomorphic(rid, vcf_start - 1,
                                                             alleles, ad_ref, ad_other, dp)
                                else:
                                    write_monomorphic_record(output_file, contig,
                                                             vcf_start - 1, ref_base,
                                                             ad_ref, ad_other, dp)
                            vcf_start += 1
                            rel += 1
                            continue

                        # some caller has a record here: collect all of its
                        # records at this position (one per split ALT after
                        # normalisation) and build one vote per caller
                        ref_base = ref_seq[rel]
                        keys = [None] * n_files
                        qualities = [0.0] * n_files
                        present = [False] * n_files
                        for fi in range(n_files):
                            records = records_per_file[fi]
                            i = indices[fi]
                            caller_records = []
                            while (i < len(records)
                                   and records[i].pos == vcf_start):
                                caller_records.append(records[i])
                                i += 1
                            if caller_records:
                                key = genotype_key(caller_records, ref_base)
                                keys[fi] = key
                                if key is not None:
                                    qualities[fi] = max(
                                        (float(record.qual)
                                         for record in caller_records
                                         if record.qual is not None),
                                        default=0.0)
                                present[fi] = True
                                indices[fi] = i
                        min_next = _NO_RECORD
                        for fi in range(n_files):
                            records = records_per_file[fi]
                            i = indices[fi]
                            if i < len(records) and records[i].pos < min_next:
                                min_next = records[i].pos

                        reference_detail = None
                        if absent_vote == 'ref':
                            absent_key = (ref_base, _REFERENCE_ALTS, (0, 0))
                        elif any(not seen for seen in present):
                            if has_column:
                                reference_detail = reference_support_of(ref_base)
                                absent_key = (
                                    (ref_base, _REFERENCE_ALTS, (0, 0))
                                    if reference_detail[0] >= min_coverage
                                    else None)
                            else:
                                absent_key = None
                        else:
                            absent_key = None
                        for fi in range(n_files):
                            if not present[fi]:
                                keys[fi] = absent_key

                        calls = consensus_calls(keys, threshold, qualities,
                                                tie_policy)
                        if not calls:
                            vcf_start += 1
                            rel += 1
                            continue

                        if len(calls) == 1 and calls[0][1] == _REFERENCE_ALTS \
                                and all(allele == 0 for allele in calls[0][2]):
                            # hom-reference consensus: write the same gVCF-style
                            # record as the fast path, with honest <*> depth
                            if reference_detail is None:
                                if has_column:
                                    reference_detail = reference_support_of(ref_base)
                                else:
                                    reference_detail = (0, 0, 0)
                            dp, ad_ref, observations = reference_detail
                            ad_other = observations - ad_ref
                            if use_accel:
                                alleles = _MONOMORPHIC_ALLELES.get(ref_base)
                                if alleles is None:
                                    alleles = (ref_base + ',<*>').encode('ascii')
                                writer.write_monomorphic(rid, vcf_start - 1,
                                                         alleles, ad_ref, ad_other, dp)
                            else:
                                write_monomorphic_record(output_file, contig,
                                                         vcf_start - 1, ref_base,
                                                         ad_ref, ad_other, dp)
                            vcf_start += 1
                            rel += 1
                            continue

                        # one split record per qualifying ALT
                        site = site_of() if has_column else None
                        shift = max(len(call[0]) for call in calls)
                        for ref, alts, gt in calls:
                            alleles = [ref] + list(alts)
                            ad, dp = recover_ad_dp(site, alleles)
                            if use_accel:
                                writer.write_record(rid, vcf_start - 1, alleles,
                                                    gt, ad, dp)
                            else:
                                write_consensus_record(output_file, contig,
                                                       vcf_start - 1,
                                                       vcf_start, alleles, gt,
                                                       ad, dp)
                        vcf_start += shift
                        rel += shift
        finally:
            if writer is not None:
                writer.close()


if __name__ == '__main__':
    main()
