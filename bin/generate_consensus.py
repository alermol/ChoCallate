#!/usr/bin/env python3
"""Generate the consensus BCF of one sample from its caller BCFs.

By default every position of the (bgzipped, tabixed) coverage BED that
resolves to a reference or variant consensus is represented in the output, so
the consensus contains all sites, monomorphic ones included.  Positions whose
reference base is an IUPAC ambiguity code (N included) and positions that
resolve to a no-call are omitted: the former cannot be written without an
ambiguity code, and the latter has no consensus genotype to report.  Caller
alleles carrying an ambiguity code are omitted as well, so output sequence
alleles contain only A, C, G and T.  A called deletion's REF span represents
the positions it covers, and ``--variants_only`` omits reference records
instead.  At each position the callers vote with their
genotype in two stages: an ALT allele must first be supported by at least
``--consensus_threshold`` callers, then the majority zygosity among the
supporters of the winning ALT decides between heterozygous and homozygous
ALT.  Ties are resolved by ``--tie_policy``: ``nocall`` (default, independent
of caller order) or ``first`` (the earliest supporter in configured caller
order).  When several ALT alleles qualify, each becomes a split consensus
record, but their combined ALT copies must fit the site ploidy; otherwise the
position is a no-call.  A position without ALT support is hom-reference when
at least ``--consensus_threshold`` callers vote reference, otherwise it is a
no-call.

A caller without a record votes reference by default (``--absent_vote ref``).
With ``--absent_vote abstain`` it only votes reference when the sample pileup
has at least ``--min_coverage`` reads at the position; a caller record with a
missing genotype (``./.``, ``0/.``) always abstains.  A caller record whose
FILTER is not ``PASS``/``.`` abstains as well by default
(``--filtered_vote abstain``); with ``--filtered_vote ref`` it votes reference
instead.  Hom-reference calls of any ploidy vote reference.  Positions
resolved as hom-reference are written as gVCF-style records with the
symbolic ALT ``<*>``: the reference allele plus "any other allele", whose
FORMAT/AD slot carries the depth of reads that do not support the reference
base.  ``bcftools merge --gvcf`` maps that slot onto the real ALT allele of
another sample's variant, so merged hom-reference samples report their actual
non-reference depth instead of a fabricated 0.

With ``--compress_reference_blocks``, runs of consecutive hom-reference
positions without caller records are written as one gVCF block carrying
``INFO/END``.  FORMAT/DP is the minimum depth of the block; FORMAT/AD is omitted
because one value cannot honestly describe every position of the block and
``bcftools merge --gvcf`` then reports a missing depth at variant positions of
other samples instead of a fabricated zero.  Single-base reference records keep
their per-position AD.  Blocks are split at caller records, at positions that
cross ``--min_coverage``, at IUPAC-ambiguous reference positions and wherever
the BED coverage is not contiguous.

With ``--variants_only``, reference records are not written at all: only
positions with a consensus ALT allele are emitted.  The pipeline sets this flag
when per-sample output would discard reference records during export anyway;
the merged output still needs them for ``bcftools merge --gvcf``.

FORMAT/AD and FORMAT/DP are recovered from the sample BAM: overlapping mate
pairs count once (the observation with the highest base quality), and DP counts
the collapsed observations, matching the coverage BED produced by mosdepth,
which also corrects for overlapping mates.

Performance notes
-----------------
The default output has one record per BED position (one per reference run with
``--compress_reference_blocks``), so the per-position cost decides the runtime:

* Genotypes, pileups and records are only touched where they are needed.
  Positions without a caller record (the vast majority) take a fast path that
  writes the reference record directly, without building genotype tuples or
  the per-allele dictionaries.
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


def _warn_accel_unavailable(error):
    print('[WARN] Cython accelerator unavailable, using pure pysam: '
          '{}'.format(error), file=sys.stderr)


if os.environ.get('CHOCALLATE_NO_ACCEL'):
    BcfLineWriter = None
    PileupTally = None
else:
    try:
        # Optional Cython helpers built by bin/build_consensus_accel.py.
        from consensus_accel import BcfLineWriter, PileupTally
    except ModuleNotFoundError as error:
        # A checkout without the optional extension is the documented
        # normal state, so stay silent: this module is imported once per
        # consensus chunk and a warning would repeat for every chunk.  A
        # missing module other than the extension itself means the
        # extension exists but its dependency does not, which is reported.
        if error.name != 'consensus_accel':
            _warn_accel_unavailable(error)
        BcfLineWriter = None
        PileupTally = None
    except Exception as error:
        # A stale extension built against another pysam raises ValueError
        # ("size changed, may indicate binary incompatibility"), not
        # ImportError; an extension that exists but fails to load is
        # reported.
        _warn_accel_unavailable(error)
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

# ALT placeholder of a reference call (no ALT allele).
_REFERENCE_ALTS = ('.',)

# VCF allows only A, C, G, T and N in REF; IUPAC ambiguity codes from the
# reference FASTA (R, Y, ...) or caller records become N.  N is itself an
# ambiguity code, so records carrying one are not written at all.
_VCF_REF_BASES = frozenset('ACGTN')

# Output sequence alleles may not contain IUPAC ambiguity codes, N included.
_UNAMBIGUOUS_BASES = frozenset('ACGT')

# Pre-encoded alleles for the monomorphic fast path (REF + ALT='<*>').
_MONOMORPHIC_ALLELES = {
    base: (base + ',<*>').encode('ascii') for base in 'ACGT'}


def to_vcf_ref(ref):
    """Map every character of a REF allele to the VCF REF alphabet."""
    return ''.join(base if base in _VCF_REF_BASES else 'N'
                   for base in ref.upper())


def has_ambiguous_bases(allele):
    """True when a sequence allele contains an IUPAC ambiguity code.

    Symbolic alleles (``<...>``), the spanning-deletion placeholder ``*`` and
    the no-ALT placeholder ``.`` carry no sequence and are not flagged.
    """
    if allele in ('.', '*') or allele.startswith('<'):
        return False
    return any(base not in _UNAMBIGUOUS_BASES for base in allele.upper())


class BedIntervalError(ValueError):
    """A BED interval cannot be read from the reference."""


def validate_bed_interval(chrom, start, end, contig_lengths):
    """Reject a BED interval that cannot be read from the reference.

    Args:
        chrom: Contig name of the interval.
        start: 0-based interval start.
        end: 0-based, exclusive interval end.
        contig_lengths: Mapping of contig name to reference length.

    Raises:
        BedIntervalError: when the contig is absent from the reference,
            the coordinates are negative or reversed, or the interval ends
            past the contig length.
    """
    length = contig_lengths.get(chrom)
    if length is None:
        raise BedIntervalError(
            "BED contig '{}' is not present in the reference".format(chrom))
    if start < 0 or end < start:
        raise BedIntervalError(
            'invalid BED interval {}:{}-{}: start must be >= 0 and '
            'end >= start'.format(chrom, start, end))
    if end > length:
        raise BedIntervalError(
            "BED interval {}:{}-{} extends past the end of contig "
            "'{}' ({} bp); clip the interval to the contig".format(
                chrom, start, end, chrom, length))


def batched_intervals(bed_file, max_gap, contig_lengths):
    """Group BED intervals into spans that share one pileup iterator.

    Args:
        bed_file: pysam.TabixFile of the (bgzipped) coverage BED.
        max_gap: Maximum gap between intervals of the same span.
        contig_lengths: Mapping of contig name to reference length.

    Yields:
        Tuples of (contig, [(start, end), ...]) where start/end are the 0-based
        half-open BED coordinates.

    Raises:
        BedIntervalError: when an interval cannot be read from the reference.
    """
    contig = None
    span = []
    span_end = None
    for line in bed_file.fetch():
        chrom, start, end = line.split('\t')[:3]
        start = int(start)
        end = int(end)
        validate_bed_interval(chrom, start, end, contig_lengths)
        if span and (chrom != contig or start - span_end > max_gap):
            yield contig, span
            span = []
        span.append((start, end))
        contig = chrom
        span_end = end
    if span:
        yield contig, span


class CallerCursor:
    """Stream the records of one caller BCF over a span.

    Keeping only the current head record instead of materialising
    ``list(fetch(...))`` bounds the memory use by the number of caller records
    at one position rather than by the size of the span.
    """

    __slots__ = ('_head', '_iterator')

    def __init__(self, iterator):
        self._iterator = iterator
        self._head = next(iterator, None)

    def position(self):
        """1-based position of the next record, or ``_NO_RECORD``."""
        if self._head is None:
            return _NO_RECORD
        return self._head.pos

    def skip_before(self, pos1, on_skip=None):
        """Discard records that start before 1-based ``pos1``.

        ``on_skip`` is called with every discarded record.  A deletion
        anchored before ``pos1`` can still cover it, so the caller inspects
        the discarded records before they are lost.
        """
        head = self._head
        iterator = self._iterator
        while head is not None and head.pos < pos1:
            if on_skip is not None:
                on_skip(head)
            head = next(iterator, None)
        self._head = head

    def take_at(self, pos1):
        """Consume and return all records at 1-based position ``pos1``."""
        records = []
        head = self._head
        iterator = self._iterator
        while head is not None and head.pos == pos1:
            records.append(head)
            head = next(iterator, None)
        self._head = head
        return records


def genotype_key(records, ref_base, filtered_vote='abstain'):
    """Consensus voting key of one caller at a position.

    ``records`` are all records of the caller at the position: normally one,
    but several after multi-allelic splitting or allele atomisation.  Each
    record contributes the copies of its own ``(REF, ALT)`` alleles, so
    atomised records with different REF alleles at the same POS stay distinct
    instead of forcing the caller to abstain.  Hom-reference calls (any ploidy)
    become the reference vote; records with a missing genotype are ignored.
    Records whose FILTER is not ``PASS``/``.`` are ignored; when every record
    is filtered the caller abstains (``filtered_vote='abstain'``) or votes
    reference (``filtered_vote='ref'``).

    Split records must be consistent with one ploidy: when the ALT copies of
    all records exceed the ploidy (for example a split ``1/1`` plus ``0/1`` in
    a diploid caller), the caller abstains instead of inventing a genotype.

    Args:
        records: Non-empty list of pysam.VariantRecord objects at one position.
        ref_base: Reference base at the position, upper case.
        filtered_vote: Vote of a caller whose records are all filtered.

    Returns:
        ``(reference allele, {(ref, alt): copies}, ploidy)``, or ``None`` when
        the records do not form a usable vote.  A reference vote has an empty
        copies dictionary.
    """
    usable = [record for record in records
              if not record.filter.keys()
              or set(record.filter.keys()) == {'PASS'}]
    if not usable:
        if filtered_vote == 'ref':
            return (ref_base, {}, 0)
        return None
    alt_copies = {}
    ploidy = 0
    for record in usable:
        gt = record.samples[0].get('GT')
        if gt is None or any(allele is None for allele in gt):
            continue
        ploidy = max(ploidy, len(gt))
        alts = tuple(allele.upper() for allele in (record.alts or ()))
        if not alts or alts == _REFERENCE_ALTS:
            continue
        ref = to_vcf_ref(record.ref)
        for allele_index, alt in enumerate(alts, 1):
            copies = sum(1 for allele in gt if allele == allele_index)
            if copies:
                allele = (ref, alt)
                alt_copies[allele] = max(alt_copies.get(allele, 0), copies)
    total_copies = sum(alt_copies.values())
    if total_copies > ploidy:
        # inconsistent split genotypes (for example 1/1 plus 0/1 in a diploid
        # caller): abstain instead of inventing a genotype
        return None
    if not alt_copies:
        if ploidy:
            return (ref_base, {}, ploidy)
        return None
    return (ref_base, alt_copies, ploidy)


def supported_deletion_spans(record, minimum_pos, filtered_vote='abstain'):
    """REF spans of the deletion alleles a caller record carries.

    A record whose POS lies before a BED chunk is still fetched when its REF
    span overlaps the chunk.  The previous chunk already voted on it, but the
    cursor would discard it before the walk reaches the covered positions.
    Returning the spans lets the caller carry the deletion suppression across
    the chunk boundary.

    Args:
        record: pysam.VariantRecord.
        minimum_pos: 1-based position that a REF span must reach.
        filtered_vote: Vote policy for filtered records.  A filtered record
            supplies no ALT copies under either policy, because ``ref`` turns
            it into a reference vote.

    Returns:
        List of ``(REF, ALT, span_end)`` tuples for the deletion alleles the
        genotype carries.  ``span_end`` is the 1-based inclusive end of the
        record's REF span.  Empty when the record does not reach
        ``minimum_pos`` or carries no usable deletion allele.
    """
    if record.stop < minimum_pos:
        return []
    key = genotype_key([record], '', filtered_vote)
    if key is None:
        return []
    _, copies, _ = key
    return [(ref, alt, record.stop)
            for ref, alt in copies if len(ref) > len(alt)]


def note_deletion_support(record, minimum_pos, caller_index, support,
                          filtered_vote):
    """Add one caller's boundary-crossing deletions to the support map.

    ``support`` maps ``(POS, REF, ALT, span_end)`` to the set of caller
    indexes carrying the allele, so the site-level consensus threshold can be
    applied without retaining the records.
    """
    for ref, alt, span_end in supported_deletion_spans(
            record, minimum_pos, filtered_vote):
        support.setdefault((record.pos, ref, alt, span_end), set()).add(
            caller_index)


def _zygosity(key, allele):
    """Classify a vote supporting ``allele`` as ``het`` or ``hom`` ALT.

    ``allele`` is a ``(ref, alt)`` pair, and copies are keyed by that pair, so
    deletions whose anchor ALT is the same but whose REF differs are distinct.
    Copies equal to the caller ploidy are homozygous ALT; fewer copies are
    heterozygous.  A haploid ALT call is therefore homozygous (one ALT copy
    out of one).
    """
    _, copies, ploidy = key
    count = copies.get(allele, 0)
    if count == 0:
        return None
    if count >= ploidy:
        return 'hom'
    return 'het'


def _genotype_for(zygosity, ploidy):
    """Output GT for a zygosity and ploidy (ALT index 1)."""
    if zygosity == 'hom':
        return (1,) * ploidy
    return (0,) * (ploidy - 1) + (1,)


def consensus_calls(keys, consensus_threshold, tie_policy='nocall'):
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

    Every qualifying ALT produces its own split consensus record; the output
    stage re-joins them unless split output was requested.  The split records
    share one downstream genotype, so their combined ALT copies must fit the
    site ploidy.  When the qualifying ALTs would use more copies than the
    ploidy (for example two hom-ALT alleles at a diploid site), the position is
    a no-call and the list is empty.  When no ALT reaches the threshold, the
    position is hom-reference if at least ``consensus_threshold`` reference
    votes remain, otherwise the list is empty.  A caller key of ``None`` (an
    abstention) does not support any allele; a caller whose split records carry
    more ALT copies than its ploidy abstains in ``genotype_key`` and therefore
    supports nothing here.

    Args:
        keys: List of ``(ref, {(allele_ref, allele_alt): copies}, ploidy)``
            keys and ``None`` for abstentions, one entry per caller.  An ALT is
            identified by its ``(ref, alt)`` pair, not by the ALT string alone,
            so deletions of different lengths are distinct alleles and a caller
            may support several alleles at one position.
        consensus_threshold: Minimum number of callers supporting an ALT (or
            the reference) at the site stage.
        tie_policy: ``nocall`` or ``first``.

    Returns:
        A list of (ref, alt tuple, gt tuple) records; normally zero, one or one
        split record per qualifying ALT, ordered by ``(ref, alt)``.
    """
    reference_key = None
    reference_count = 0
    alt_support = {}
    for index, key in enumerate(keys):
        if key is None:
            continue
        ref, copies, ploidy = key
        if not copies:
            reference_count += 1
            if reference_key is None:
                reference_key = (ref, _REFERENCE_ALTS, (0, 0))
            continue
        for allele in copies:
            alt_support.setdefault(allele, []).append((index, key))

    candidates = sorted(
        (allele, supporters) for allele, supporters in alt_support.items()
        if len(supporters) >= consensus_threshold)
    if not candidates:
        if reference_count >= consensus_threshold:
            return [reference_key]
        return []

    calls = []
    for allele, supporters in candidates:
        classes = {}
        for index, key in supporters:
            zygosity = _zygosity(key, allele)
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
            if len(winners) > 1:
                # unresolved genotype tie: this ALT gets no consensus record
                continue
        zygosity, members = winners[0]
        ploidy = max(key[2] for _, key in members)
        allele_ref, allele_alt = allele
        calls.append((allele_ref, (allele_alt,),
                      _genotype_for(zygosity, ploidy)))
    # The split calls are re-joined into a single multiallelic genotype by the
    # output stage, so their ALT copies must fit the site ploidy.  For example,
    # two hom-ALT alleles at a diploid site would contribute four copies and
    # the joined genotype could keep only one of them.
    site_ploidy = max((len(gt) for _, _, gt in calls), default=0)
    alt_copies = sum(sum(1 for allele in gt if allele)
                     for _, _, gt in calls)
    if alt_copies > site_ploidy:
        return []
    return calls


def site_from_entries(names, seqs, quals):
    """Per-allele counters from materialised pileup column entries.

    Overlapping PE mates are deduplicated (one observation per query name,
    keeping the highest base quality), and DP is the number of collapsed
    observations, matching the overlap correction of the mosdepth coverage BED.
    Reference skips are not observations; a deletion is an observation without
    an allele to count.

    Args:
        names: Query names at the column.
        seqs: Bases with samtools-style indel annotations ('A+2TT', 'A-2',
            '*', '>', '<').
        quals: Base qualities (255 means missing).

    Returns:
        PileupSite(snp dict, insertion dict, deletion-length dict, depth).
    """
    best = {}
    observations = set()
    for i, name in enumerate(names):
        seq = seqs[i]
        if seq == '>' or seq == '<':
            # reference skip: no observation at this position
            continue
        observations.add(name)
        if seq == '*':
            # deleted base: an observation, but no allele to count
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
    return PileupSite(snp=snp, ins=ins, dele=dele, dp=len(observations))


def build_pileup_site(column):
    """Summarise a single pileup column into per-allele counters.

    DP is the number of collapsed observations, matching the mosdepth coverage
    BED, which corrects for overlapping mates by default.  The per-read data
    comes from pysam's bulk accessors, which are much cheaper than iterating
    PileupRead/AlignedSegment objects.

    Args:
        column: pysam.PileupColumn.

    Returns:
        PileupSite(snp dict, insertion dict, deletion-length dict, depth).
    """
    names = column.get_query_names()
    seqs = column.get_query_sequences(add_indels=True)
    return site_from_entries(names, seqs, column.get_query_qualities())


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
    """Reference-base support, depth and observation count of a pileup column.

    Equivalent to ``build_pileup_site`` followed by
    ``count_allele_depth(site, ref_base, ref_base)``, but avoids building the
    insertion/deletion counters and the allele tuples that the monomorphic fast
    path never needs.

    Args:
        column: pysam.PileupColumn.
        ref_base: Reference base at the position (upper case).

    Returns:
        (depth, number of collapsed reads supporting the reference base, number
        of collapsed observations).  Depth is the number of collapsed
        observations, matching the mosdepth coverage BED; the last value feeds
        the ``<*>`` AD slot.  Reference skips are not observations.
    """
    names = column.get_query_names()
    if not names:
        return 0, 0, 0
    seqs = column.get_query_sequences(add_indels=True)
    if len(set(names)) == len(names):
        # no overlapping mate pair at this column: nothing to collapse; every
        # read that is not a reference skip is one observation
        observations = len(seqs) - seqs.count('>') - seqs.count('<')
        support = seqs.count(ref_base) + seqs.count(ref_base.lower())
        return observations, support, observations
    support, observations = _collapsed_support(
        names, seqs, column.get_query_qualities(), ref_base)
    return observations, support, observations


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
        position are one observation each; reference skips are ignored.  A
        deletion is not an allele, so an overlapping mate's base or insertion
        is always preferred over it.
    """
    best = {}
    observations = set()
    for i, name in enumerate(names):
        seq = seqs[i]
        if seq == '>' or seq == '<':
            continue
        observations.add(name)
        if seq == '*':
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
    return support, len(observations)


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
            names, seqs, quals = self.tally.entries()
            support, observations = _collapsed_support(
                names, seqs, quals, ref_base)
            return observations, support, observations
        return counts

    def site(self):
        return site_from_entries(*self.tally.entries())


def write_monomorphic_record(output_file, contig, start0, end0, ref_base,
                            ad, dp):
    """Write a gVCF-style reference record or reference block.

    The record carries the symbolic ALT ``<*>``; ``end0`` is the 1-based
    inclusive end position (equal to ``start0 + 1`` for a single base).  ``ad``
    is a ``(reference depth, other-allele depth)`` tuple, or ``None`` for a
    multi-base block, whose per-position depth cannot be represented by one AD
    value.
    """
    # pysam keeps INFO/END for symbolic alleles even when stop == start + 1,
    # matching BcfLineWriter, which sets it for every reference record.
    if has_ambiguous_bases(ref_base):
        return
    new_record = output_file.new_record(contig=contig, start=start0, stop=end0,
                                        alleles=(ref_base, '<*>'), qual=None)
    sample = new_record.samples[0]
    sample['GT'] = (0, 0)
    if ad is not None:
        sample['AD'] = ad
    sample['DP'] = dp
    output_file.write(new_record)


def write_consensus_record(output_file, contig, start0, stop0, alleles, gt, ad, dp):
    """Write a consensus record from winning caller genotypes."""
    if any(has_ambiguous_bases(allele) for allele in alleles):
        return
    new_record = output_file.new_record(contig=contig, start=start0, stop=stop0,
                                        alleles=alleles, qual=None)
    sample = new_record.samples[0]
    sample['GT'] = gt
    sample['AD'] = ad
    sample['DP'] = dp
    output_file.write(new_record)


def write_reference_record(output_file, writer, use_accel, rid, contig,
                           start0, end0, ref_base, ad, dp):
    """Write a single-base reference record or a multi-base block.

    ``end0`` is the 1-based inclusive end.  ``ad`` is ``None`` for multi-base
    blocks, whose FORMAT/AD is deliberately omitted.
    """
    if has_ambiguous_bases(ref_base):
        return
    if use_accel:
        alleles = _MONOMORPHIC_ALLELES.get(ref_base)
        if alleles is None:
            alleles = (ref_base + ',<*>').encode('ascii')
        if ad is None:
            writer.write_reference_block(rid, start0, alleles, dp, end0)
        else:
            writer.write_monomorphic(rid, start0, alleles, ad[0], ad[1], dp,
                                     end0)
    else:
        write_monomorphic_record(output_file, contig, start0, end0, ref_base,
                                 ad, dp)


def flush_reference_block(block, output_file, writer, use_accel):
    """Write and forget a pending reference block."""
    if block is None:
        return
    contig, rid, start0, end0, ref_base, ad, dp = block[:7]
    if end0 != start0 + 1:
        # One AD value cannot describe every position of a multi-base block.
        ad = None
    write_reference_record(output_file, writer, use_accel, rid, contig,
                           start0, end0, ref_base, ad, dp)


def generate_minimal_header(sample_name,
                            reference_file,
                            consensus_threshold,
                            version,
                            callers,
                            callers_failed,
                            tie_policy,
                            absent_vote,
                            filtered_vote,
                            min_coverage,
                            compress_reference_blocks=False,
                            use_accel=False):
    """
    Generate a minimal VCF header for the consensus file.

    Args:
        sample_name: String sample name.
        reference_file: pysam.FastaFile object.
        consensus_threshold: Integer threshold for the site-level vote.
        version: String version of the ChoCallate pipeline.
        callers: Callers that voted for this sample, in voting order.
        callers_failed: Configured callers without a usable record set.
        tie_policy: Tie policy applied by the voting model.
        absent_vote: Policy for callers without a record at a position.
        filtered_vote: Policy for caller records with a non-PASS FILTER.
        min_coverage: Coverage threshold used by ``absent_vote=abstain``.
        compress_reference_blocks: Whether reference blocks were requested.
        use_accel: Whether the Cython accelerator wrote the records.

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
    meta.append(f'FilteredVote={filtered_vote}')
    meta.append(f'ReferenceBlocks={compress_reference_blocks}')
    meta.append(f'Accel={str(bool(use_accel)).lower()}')
    if callers_failed:
        meta.append(f'CallersFailed="{",".join(callers_failed)}"')
    header.add_line('##ChoCallate=<' + ','.join(meta) + '>')
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


def _bounded_int(value, minimum):
    """argparse type for an integer of at least ``minimum``."""
    try:
        number = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError(
            'invalid integer: {!r}'.format(value)) from None
    if number < minimum:
        raise argparse.ArgumentTypeError(
            'must be at least {}, got {}'.format(minimum, number))
    return number


def positive_int(value):
    """argparse type for an integer of at least 1."""
    return _bounded_int(value, 1)


def non_negative_int(value):
    """argparse type for an integer of at least 0."""
    return _bounded_int(value, 0)


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
    parser.add_argument('--consensus_threshold', required=True,
                        type=positive_int,
                        help='Minimum caller support at the site level '
                             '(at least 1)')
    parser.add_argument('--version', required=True, type=str, help='Version of ChoCallate pipeline')
    parser.add_argument('--absent_vote', choices=('ref', 'abstain'), default='ref',
                        help='Vote of a caller without a record at a position')
    parser.add_argument('--filtered_vote', choices=('abstain', 'ref'),
                        default='abstain',
                        help='Vote of a caller record whose FILTER is not '
                             'PASS or .')
    parser.add_argument('--tie_policy', choices=('nocall', 'first'),
                        default='nocall', help='How to break voting ties')
    parser.add_argument('--min_coverage', type=non_negative_int, default=0,
                        help='Coverage required for an absent caller to vote '
                             'reference when --absent_vote is abstain '
                             '(non-negative)')
    parser.add_argument('--compress_reference_blocks', '--reference_blocks',
                        action='store_true',
                        dest='compress_reference_blocks',
                        help='Compress runs of hom-reference positions into '
                             'multi-base gVCF blocks')
    parser.add_argument('--variants_only', action='store_true',
                        help='Write only consensus variant records, omitting '
                             'reference records')
    args = parser.parse_args()

    threshold = args.consensus_threshold
    absent_vote = args.absent_vote
    filtered_vote = args.filtered_vote
    tie_policy = args.tie_policy
    min_coverage = args.min_coverage
    use_accel = BcfLineWriter is not None
    variants_only = args.variants_only
    compress_blocks = args.compress_reference_blocks and not variants_only
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
        contig_lengths = dict(zip(reference_file.references,
                                  reference_file.lengths))
        bam_file = stack.enter_context(pysam.AlignmentFile(args.bam, 'rb'))
        header = generate_minimal_header(args.sample_name,
                                         reference_file,
                                         threshold,
                                         args.version,
                                         received,
                                         failed,
                                         tie_policy,
                                         absent_vote,
                                         filtered_vote,
                                         min_coverage,
                                         compress_blocks,
                                         use_accel)
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

        pending_block = None
        # 1-based inclusive end of the furthest called deletion.  Reference
        # records must not be written over that span: the deletion record's
        # REF already covers those positions.
        deletion_end = 0
        current_contig = None
        try:
            for contig, intervals in batched_intervals(
                    bed_file, SPAN_MAX_GAP, contig_lengths):
                if contig != current_contig:
                    current_contig = contig
                    deletion_end = 0
                span_start = intervals[0][0]
                span_end = intervals[-1][1]
                rid = contig_ids[contig]
                # Caller records are streamed once per span (not per BED
                # interval) through a per-file cursor: only the records at the
                # current position are held in memory.
                cursors = [CallerCursor(variant_file.fetch(
                               contig, span_start, span_end))
                           for variant_file in variant_files]
                min_next = min(cursor.position() for cursor in cursors)
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
                    ref_seq = to_vcf_ref(
                        reference_file.fetch(contig, start0, end0))
                    rel = 0
                    # BED 0-based half-open [start0, end0) covers 1-based positions
                    # start0+1 .. end0, so the walk stops after end0.
                    vcf_start = start0 + 1
                    while vcf_start <= end0:
                        # keep the pileup in step with the walk
                        has_column = sync(vcf_start)

                        if min_next < vcf_start:
                            # the start of the next interval, or a record that
                            # overlaps the span start, moved past some caller
                            # records
                            min_next = _NO_RECORD
                            # A deletion anchored in the previous BED chunk is
                            # fetched here because its REF span overlaps this
                            # chunk, but its POS is before the chunk start.
                            # Seed its span before the cursors discard the
                            # record, so reference records are not written
                            # over bases the previous chunk called deleted.
                            # Records skipped inside a BED gap are not seeded:
                            # no earlier interval wrote their deletion record.
                            seed_deletions = vcf_start == span_start + 1
                            deletion_support = {}
                            for fi, cursor in enumerate(cursors):
                                if seed_deletions:
                                    cursor.skip_before(
                                        vcf_start,
                                        lambda record, fi=fi: (
                                            note_deletion_support(
                                                record, vcf_start, fi,
                                                deletion_support,
                                                filtered_vote)))
                                else:
                                    cursor.skip_before(vcf_start)
                                position = cursor.position()
                                if position < min_next:
                                    min_next = position
                            for key, supporters in deletion_support.items():
                                if len(supporters) >= threshold:
                                    deletion_end = max(deletion_end, key[3])

                        if min_next != vcf_start:
                            # fast path: no caller record here.  With
                            # variants_only such positions are not written and
                            # need no reference support count.
                            if variants_only:
                                vcf_start += 1
                                rel += 1
                                continue
                            if vcf_start <= deletion_end:
                                # The position is inside a called deletion.
                                # Its REF span covers the site, so a
                                # hom-reference record would contradict the
                                # deletion genotype.
                                flush_reference_block(
                                    pending_block, output_file, writer,
                                    use_accel)
                                pending_block = None
                                vcf_start += 1
                                rel += 1
                                continue
                            ref_base = ref_seq[rel]
                            if has_ambiguous_bases(ref_base):
                                # An IUPAC code (N included) cannot be
                                # represented in the output; end the block and
                                # omit the position.
                                flush_reference_block(
                                    pending_block, output_file, writer,
                                    use_accel)
                                pending_block = None
                                vcf_start += 1
                                rel += 1
                                continue
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
                                ad = (ad_ref, ad_other)
                                if compress_blocks:
                                    if (pending_block is not None
                                            and pending_block[0] == contig
                                            and vcf_start - 1
                                            == pending_block[3]
                                            and (dp >= min_coverage)
                                            == pending_block[7]):
                                        pending_block[3] = vcf_start
                                        if dp < pending_block[6]:
                                            pending_block[6] = dp
                                    else:
                                        flush_reference_block(
                                            pending_block, output_file,
                                            writer, use_accel)
                                        pending_block = [
                                            contig, rid, vcf_start - 1,
                                            vcf_start, ref_base, ad, dp,
                                            dp >= min_coverage]
                                else:
                                    write_reference_record(
                                        output_file, writer, use_accel,
                                        rid, contig, vcf_start - 1,
                                        vcf_start, ref_base, ad, dp)
                            else:
                                flush_reference_block(
                                    pending_block, output_file, writer,
                                    use_accel)
                                pending_block = None
                            vcf_start += 1
                            rel += 1
                            continue

                        # A caller record bounds a reference block: the
                        # position is evaluated individually below.
                        flush_reference_block(pending_block, output_file,
                                              writer, use_accel)
                        pending_block = None
                        # some caller has a record here: collect all of its
                        # records at this position (one per split ALT after
                        # normalisation) and build one vote per caller
                        ref_base = ref_seq[rel]
                        keys = [None] * n_files
                        present = [False] * n_files
                        for fi, cursor in enumerate(cursors):
                            caller_records = cursor.take_at(vcf_start)
                            if caller_records:
                                key = genotype_key(caller_records, ref_base,
                                                   filtered_vote)
                                keys[fi] = key
                                present[fi] = True
                        min_next = _NO_RECORD
                        for cursor in cursors:
                            position = cursor.position()
                            if position < min_next:
                                min_next = position

                        reference_detail = None
                        if absent_vote == 'ref':
                            absent_key = (ref_base, {}, 0)
                        elif any(not seen for seen in present):
                            if has_column:
                                reference_detail = reference_support_of(ref_base)
                                absent_key = (
                                    (ref_base, {}, 0)
                                    if reference_detail[0] >= min_coverage
                                    else None)
                            else:
                                absent_key = None
                        else:
                            absent_key = None
                        for fi in range(n_files):
                            if not present[fi]:
                                keys[fi] = absent_key

                        calls = consensus_calls(keys, threshold, tie_policy)
                        if not calls:
                            vcf_start += 1
                            rel += 1
                            continue

                        if len(calls) == 1 and calls[0][1] == _REFERENCE_ALTS \
                                and all(allele == 0 for allele in calls[0][2]):
                            # hom-reference consensus: write the same gVCF-style
                            # record as the fast path, with honest <*> depth
                            if (not variants_only
                                    and vcf_start > deletion_end
                                    and not has_ambiguous_bases(ref_base)):
                                if reference_detail is None:
                                    if has_column:
                                        reference_detail = \
                                            reference_support_of(ref_base)
                                    else:
                                        reference_detail = (0, 0, 0)
                                dp, ad_ref, observations = reference_detail
                                ad_other = observations - ad_ref
                                write_reference_record(
                                    output_file, writer, use_accel, rid,
                                    contig, vcf_start - 1, vcf_start, ref_base,
                                    (ad_ref, ad_other), dp)
                            vcf_start += 1
                            rel += 1
                            continue

                        # one split record per qualifying ALT; every BED
                        # position is evaluated, including positions covered
                        # by a deletion, so caller records inside the deletion
                        # span are not discarded.  Reference records inside
                        # that span are suppressed above.
                        site = site_of() if has_column else None
                        for ref, alts, gt in calls:
                            alleles = [ref] + list(alts)
                            if any(has_ambiguous_bases(allele)
                                   for allele in alleles):
                                # IUPAC codes cannot be written to the output.
                                continue
                            ad, dp = recover_ad_dp(site, alleles)
                            if use_accel:
                                writer.write_record(rid, vcf_start - 1, alleles,
                                                    gt, ad, dp)
                            else:
                                write_consensus_record(output_file, contig,
                                                       vcf_start - 1,
                                                       vcf_start, alleles, gt,
                                                       ad, dp)
                            if (len(ref) > len(alts[0])
                                    and any(allele > 0 for allele in gt)):
                                deletion_end = max(
                                    deletion_end, vcf_start + len(ref) - 1)
                        vcf_start += 1
                        rel += 1

                flush_reference_block(pending_block, output_file, writer,
                                      use_accel)
                pending_block = None
        except BedIntervalError as error:
            sys.exit('ERROR: {}'.format(error))
        finally:
            if writer is not None:
                writer.close()


if __name__ == '__main__':
    main()
