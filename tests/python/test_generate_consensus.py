"""End-to-end regression tests for bin/generate_consensus.py.

Covers interval-end handling, missing genotypes, haploid calls, the voting
model, provenance headers, multi-allelic output and the pure/Cython equivalence
paths.  Every test runs twice, on the pure-pysam path and (when built) on the
Cython path, and asserts the documented all-sites behaviour: one output record
per evaluated BED base.
"""

import subprocess

import pytest
from conftest import (  # noqa: I100,I202  (tests run from this directory)
    ACCEL_AVAILABLE,
    BCFTOOLS,
    CHR,
    REF_SEQ,
    alt_base,
    ref_base,
    run_consensus,
    write_bam,
    write_bed,
    write_caller_bcf,
    write_reference,
    write_sparse_reference,
)

# 0-based indexes of two reference bases ('ACGT' repeats, so index % 4 picks it)
POS_20_0B = 19   # 1-based position 20, reference 'A'
POS_24_0B = 23   # 1-based position 24, covered by all five test reads
POS_30_0B = 29   # 1-based position 30, reference 'C'


def read_records(path):
    """(pos, alleles, GT) triples of a BCF written by the tests."""
    import pysam
    with pysam.VariantFile(str(path), index_filename=None) as variant_file:
        return [(record.pos, (record.ref,) + tuple(record.alts),
                 tuple(record.samples[0]['GT'])) for record in variant_file]


def read_records_full(path):
    """(pos, alleles, GT, AD, DP) tuples of a BCF written by the tests."""
    import pysam
    with pysam.VariantFile(str(path), index_filename=None) as variant_file:
        return [(record.pos, (record.ref,) + tuple(record.alts),
                 tuple(record.samples[0]['GT']), tuple(record.samples[0]['AD']),
                 record.samples[0]['DP']) for record in variant_file]


def read_spans(path):
    """(pos, inclusive end, alleles, GT, AD or None, DP) of every record."""
    import pysam
    spans = []
    with pysam.VariantFile(str(path), index_filename=None) as variant_file:
        for record in variant_file:
            sample = record.samples[0]
            ad = sample.get('AD')
            spans.append((record.pos, record.stop,
                          (record.ref,) + tuple(record.alts),
                          tuple(sample['GT']),
                          tuple(ad) if ad is not None else None,
                          sample['DP']))
    return spans


def span_positions(spans):
    """Expand POS..END of each span into the covered 1-based positions."""
    return [position for start, end, *_ in spans
            for position in range(start, end + 1)]


def test_interval_ends_are_written(consensus_runner):
    """Every base of every BED interval, last base included."""
    # chr1 9 10 is a single base (1-based 10); chr1 19 25 covers 20..25.
    output = consensus_runner(
        bed_intervals=[('chr1', 9, 10), ('chr1', 19, 25)],
        caller_records=[[], []],
        threshold=2,
    )
    assert [pos for pos, _, _ in read_records(output)] == [10, 20, 21, 22, 23, 24, 25]


def test_adjacent_intervals_are_not_truncated(consensus_runner):
    """Adjacent intervals keep their last bases (1-based 10 and 20)."""
    output = consensus_runner(
        bed_intervals=[('chr1', 0, 10), ('chr1', 10, 20)],
        caller_records=[[], []],
        threshold=2,
    )
    assert [pos for pos, _, _ in read_records(output)] == list(range(1, 21))


def test_all_sites_written_once(consensus_runner):
    """One record per BED base, in BED order, monomorphic sites included."""
    output = consensus_runner(
        bed_intervals=[('chr1', 0, 10), ('chr1', 12, 20)],
        caller_records=[[], []],
        threshold=2,
    )
    records = read_records(output)
    assert [pos for pos, _, _ in records] == list(range(1, 11)) + list(range(13, 21))
    for pos, alleles, gt in records:
        assert alleles == (ref_base(pos), '<*>')
        assert gt == (0, 0)


def test_reference_blocks_are_optional(tmp_path, use_accel):
    """Blocks compress hom-reference runs and omit multi-base AD."""
    import pysam

    reference = write_reference(tmp_path / 'ref.fa')
    reads = [(f'read{i}', i, REF_SEQ[i:i + 6]) for i in range(19, 24)]
    bam = write_bam(tmp_path / 'sample.bam', reads)
    bed_gz = write_bed(tmp_path / 'coverage.bed', [('chr1', 19, 25)])
    callers = [write_caller_bcf(tmp_path / 'c0.bcf', 'sample1', []),
               write_caller_bcf(tmp_path / 'c1.bcf', 'sample1', [])]

    plain = run_consensus(tmp_path, bam, bed_gz, callers, threshold=2,
                          use_accel=use_accel, reference=reference,
                          output_name='plain.bcf')
    assert [span[0] for span in read_spans(plain)] == list(range(20, 26))

    blocked = run_consensus(
        tmp_path, bam, bed_gz, callers, threshold=2,
        use_accel=use_accel, reference=reference, output_name='blocked.bcf',
        compress_reference_blocks=True)
    # 20..25 are one block; the minimum depth is 1 at position 20 and one AD
    # value cannot describe the whole block, so FORMAT/AD is omitted.
    assert read_spans(blocked) == [
        (20, 25, (ref_base(20), '<*>'), (0, 0), None, 1)]

    with pysam.VariantFile(str(blocked), index_filename=None) as variant_file:
        header = str(variant_file.header)
    assert 'ReferenceBlocks=True' in header

    # The block rlen must let an indexed fetch find a position inside it.
    subprocess.run([BCFTOOLS, 'index', '--csi', '-f', str(blocked)],
                   check=True)
    with pysam.VariantFile(str(blocked)) as variant_file:
        assert [record.pos for record in variant_file.fetch('chr1', 23, 24)] \
            == [20]


def test_reference_blocks_split_at_variant(tmp_path, use_accel):
    """A caller record ends the current block and is written separately."""
    reference = write_reference(tmp_path / 'ref.fa')
    reads = [(f'read{i}', i, REF_SEQ[i:i + 6]) for i in range(19, 24)]
    bam = write_bam(tmp_path / 'sample.bam', reads)
    bed_gz = write_bed(tmp_path / 'coverage.bed', [('chr1', 19, 25)])
    caller = write_caller_bcf(tmp_path / 'caller.bcf', 'sample1', [
        (21, ref_base(22), [alt_base(22)], (0, 1))])

    output = run_consensus(
        tmp_path, bam, bed_gz, [caller], threshold=1,
        use_accel=use_accel, reference=reference,
        compress_reference_blocks=True)
    spans = read_spans(output)
    assert [(start, end) for start, end, *_ in spans] == [
        (20, 21), (22, 22), (23, 25)]
    assert span_positions(spans) == list(range(20, 26))
    assert spans[1][2] == (ref_base(22), alt_base(22))


def test_reference_blocks_split_at_reference_caller_record(tmp_path,
                                                          use_accel):
    """Any caller record splits a block, even a hom-reference call."""
    reference = write_reference(tmp_path / 'ref.fa')
    reads = [(f'read{i}', i, REF_SEQ[i:i + 6]) for i in range(19, 24)]
    bam = write_bam(tmp_path / 'sample.bam', reads)
    bed_gz = write_bed(tmp_path / 'coverage.bed', [('chr1', 19, 25)])
    caller = write_caller_bcf(tmp_path / 'caller.bcf', 'sample1', [
        (21, ref_base(22), [], (0, 0))])

    output = run_consensus(
        tmp_path, bam, bed_gz, [caller], threshold=1,
        use_accel=use_accel, reference=reference,
        compress_reference_blocks=True)
    assert [(start, end) for start, end, *_ in read_spans(output)] == [
        (20, 21), (22, 22), (23, 25)]


def test_reference_blocks_do_not_span_bed_gaps(tmp_path, use_accel):
    """Discontiguous BED intervals become separate blocks."""
    reference = write_reference(tmp_path / 'ref.fa')
    reads = [(f'read{i}', i, REF_SEQ[i:i + 6]) for i in range(19, 24)]
    bam = write_bam(tmp_path / 'sample.bam', reads)
    bed_gz = write_bed(tmp_path / 'coverage.bed',
                       [('chr1', 19, 21), ('chr1', 23, 25)])
    callers = [write_caller_bcf(tmp_path / 'c0.bcf', 'sample1', [])]

    output = run_consensus(
        tmp_path, bam, bed_gz, callers, threshold=1,
        use_accel=use_accel, reference=reference,
        compress_reference_blocks=True)
    spans = read_spans(output)
    assert [(start, end) for start, end, *_ in spans] == [(20, 21), (24, 25)]
    assert span_positions(spans) == [20, 21, 24, 25]


def test_reference_blocks_split_when_depth_crosses_min_coverage(
        tmp_path, use_accel):
    """A block never spans positions on both sides of min_coverage."""
    reference = write_reference(tmp_path / 'ref.fa')
    reads = [(f'read{i}', i, REF_SEQ[i:i + 6]) for i in range(19, 24)]
    bam = write_bam(tmp_path / 'sample.bam', reads)
    bed_gz = write_bed(tmp_path / 'coverage.bed', [('chr1', 19, 25)])
    callers = [write_caller_bcf(tmp_path / 'c0.bcf', 'sample1', [])]

    output = run_consensus(
        tmp_path, bam, bed_gz, callers, threshold=1, min_coverage=3,
        use_accel=use_accel, reference=reference,
        compress_reference_blocks=True)
    spans = read_spans(output)
    # depths at 1-based 20..25 are 1, 2, 3, 4, 5, 5
    assert [(start, end, dp) for start, end, *_, dp in spans] == [
        (20, 21, 1), (22, 25, 3)]
    assert span_positions(spans) == list(range(20, 26))


def test_accel_and_pure_paths_agree_on_blocks(tmp_path):
    """Both implementations write byte-identical reference blocks."""
    if not ACCEL_AVAILABLE:
        pytest.skip('consensus_accel is not built for this interpreter')

    reference = write_reference(tmp_path / 'ref.fa')
    reads = [('r1', 19, REF_SEQ[19:25]), ('r2', 29, REF_SEQ[29:35])]
    bam = write_bam(tmp_path / 'sample.bam', reads)
    bed_gz = write_bed(tmp_path / 'coverage.bed', [('chr1', 19, 36)])
    caller_bcfs = [
        write_caller_bcf(tmp_path / 'caller0.bcf', 'sample1',
                         [(19, REF_SEQ[19], [], (0, 0)),
                          (POS_30_0B, ref_base(30), [alt_base(30)], (0, 1))]),
        write_caller_bcf(tmp_path / 'caller1.bcf', 'sample1', []),
    ]
    pure = run_consensus(
        tmp_path, bam, bed_gz, caller_bcfs, threshold=1, use_accel=False,
        output_name='pure.bcf', reference=reference,
        compress_reference_blocks=True)
    accel = run_consensus(
        tmp_path, bam, bed_gz, caller_bcfs, threshold=1, use_accel=True,
        output_name='accel.bcf', reference=reference,
        compress_reference_blocks=True)

    def text(path):
        return subprocess.run([BCFTOOLS, 'view', '-H', str(path)],
                              check=True, capture_output=True,
                              text=True).stdout

    assert text(pure) == text(accel)
    assert text(pure)


def test_missing_genotype_does_not_crash(consensus_runner):
    """A ``./.`` caller record is handled instead of raising TypeError."""
    output = consensus_runner(
        bed_intervals=[('chr1', POS_20_0B, 25)],
        caller_records=[[(POS_20_0B, REF_SEQ[POS_20_0B], [], (None, None))], [], []],
        threshold=2,
    )
    # The no-call abstains; the two other callers' placeholder votes reach the
    # threshold, so the position is written as hom-reference.
    assert [pos for pos, _, _ in read_records(output)] == [20, 21, 22, 23, 24, 25]


def test_missing_genotype_abstains_from_the_vote(consensus_runner):
    """A no-call record is not a hom-reference vote.

    Two callers, threshold 2: where one caller emitted ``./.`` only the other
    caller's placeholder vote is left, so the position stays below the threshold
    and is not reported as hom-reference (it shows up as missing after the merge).
    """
    output = consensus_runner(
        bed_intervals=[('chr1', POS_20_0B, POS_20_0B + 2)],
        caller_records=[[(POS_20_0B, REF_SEQ[POS_20_0B], [], (None, None))], []],
        threshold=2,
    )
    assert [pos for pos, _, _ in read_records(output)] == [21]


def test_missing_genotype_does_not_mask_variant(consensus_runner):
    """A no-call must not out-vote a real call from another caller."""
    output = consensus_runner(
        bed_intervals=[('chr1', POS_30_0B, POS_30_0B + 1)],
        caller_records=[
            [(POS_30_0B, ref_base(30), [], (None, None))],
            [(POS_30_0B, ref_base(30), [alt_base(30)], (0, 1))],
        ],
        threshold=1,
    )
    assert read_records(output) == [(30, (ref_base(30), alt_base(30)), (0, 1))]


def test_partial_missing_genotype_abstains(consensus_runner):
    """A ``0/.`` record abstains, it is not a hom-reference vote."""
    output = consensus_runner(
        bed_intervals=[('chr1', POS_30_0B, POS_30_0B + 1)],
        caller_records=[
            [(POS_30_0B, ref_base(30), [], (0, None))],
            [(POS_30_0B, ref_base(30), [alt_base(30)], (0, 1))],
        ],
        threshold=1,
    )
    assert read_records(output) == [(30, (ref_base(30), alt_base(30)), (0, 1))]


def test_haploid_hom_ref_is_reference(consensus_runner):
    """A haploid ``0`` is a hom-reference vote, not a variant."""
    output = consensus_runner(
        bed_intervals=[('chr1', POS_30_0B, 31)],
        caller_records=[[(POS_30_0B, ref_base(30), [], (0,))]],
        threshold=1,
    )
    assert read_records(output) == [
        (30, (ref_base(30), '<*>'), (0, 0)),
        (31, (ref_base(31), '<*>'), (0, 0)),
    ]


def test_haploid_alt_is_written(consensus_runner):
    """A haploid ``1`` does not crash and is reported as a variant."""
    output = consensus_runner(
        bed_intervals=[('chr1', POS_30_0B, 31)],
        caller_records=[[(POS_30_0B, ref_base(30), [alt_base(30)], (1,))]],
        threshold=1,
    )
    assert read_records(output) == [
        (30, (ref_base(30), alt_base(30)), (1,)),
        (31, (ref_base(31), '<*>'), (0, 0)),
    ]


def test_consensus_below_threshold_leaves_placeholder(consensus_runner):
    """A caller 0/0 record loses to the placeholders of the other callers."""
    output = consensus_runner(
        bed_intervals=[('chr1', POS_20_0B, POS_20_0B + 1)],
        caller_records=[[(POS_20_0B, REF_SEQ[POS_20_0B], [], (0, 0))], []],
        threshold=2,
    )
    assert read_records(output) == [(20, (ref_base(20), '<*>'), (0, 0))]


def test_monomorphic_record_carries_other_allele_depth(tmp_path, use_accel):
    """The <*> AD slot holds the non-reference read depth."""
    reference = write_reference(tmp_path / 'ref.fa')
    ref = ref_base(30)
    alt = alt_base(30)
    reads = [(f'ref{start}', start, REF_SEQ[start:start + 10])
             for start in range(19, 30)]
    for i in range(2):
        seq = list(REF_SEQ[POS_30_0B:POS_30_0B + 10])
        seq[0] = alt
        reads.append((f'alt{i}', POS_30_0B, ''.join(seq)))
    reads.sort(key=lambda item: item[1])
    bam = write_bam(tmp_path / 'sample.bam', reads)
    bed_gz = write_bed(tmp_path / 'coverage.bed', [('chr1', POS_30_0B, POS_30_0B + 1)])
    caller_bcfs = [write_caller_bcf(tmp_path / 'c0.bcf', 'sample1', []),
                   write_caller_bcf(tmp_path / 'c1.bcf', 'sample1', [])]

    output = run_consensus(tmp_path, bam, bed_gz, caller_bcfs, threshold=2,
                           use_accel=use_accel, reference=reference)
    # 10 reads support the reference base, 2 support the ALT base; the <*> AD
    # slot must report the non-reference depth instead of a fabricated 0.
    assert read_records_full(output) == [(30, (ref, '<*>'), (0, 0), (10, 2), 12)]


def test_concordant_callers_are_merged(consensus_runner):
    """Two matching ALT calls win over one placeholder vote.

    The genotype classes (0/1 and 1/1) tie, so the explicit ``first`` policy
    picks the class of the earliest caller instead of relying on the default
    no-call policy.
    """
    output = consensus_runner(
        bed_intervals=[('chr1', POS_30_0B, 31), ('chr1', 40, 42)],
        caller_records=[
            [(POS_30_0B, ref_base(30), [alt_base(30)], (0, 1), 50)],
            [(POS_30_0B, ref_base(30), [alt_base(30)], (1, 1), 20)],
        ],
        threshold=1,
        tie_policy='first',
    )
    records = read_records(output)
    # The two callers agree on the ALT allele, so the site passes the site
    # stage; the first caller's 0/1 class breaks the zygosity tie.  Placeholder
    # votes must not win, and the interval ends are written.
    assert records[0] == (30, (ref_base(30), alt_base(30)), (0, 1))
    assert [pos for pos, _, _ in records] == [30, 31, 41, 42]


@pytest.mark.parametrize('order', [(0, 1), (1, 0)])
@pytest.mark.parametrize('threshold,expected_count', [(2, 0), (1, 2)])
def test_deletions_of_different_lengths_are_distinct_alleles(
        consensus_runner, order, threshold, expected_count):
    """Normalised deletion ALTs with the same anchor base do not pool.

    ``REF=CG ALT=C`` (1-bp deletion) and ``REF=CGT ALT=C`` (2-bp deletion)
    share the ALT string but are different alleles.  With threshold 2 neither
    has enough supporters and no variant is written; with threshold 1 both
    are written and the order is independent of the caller order.
    """
    anchor = REF_SEQ[POS_30_0B].upper()
    short = (POS_30_0B, REF_SEQ[POS_30_0B:POS_30_0B + 2].upper(),
             [anchor], (0, 1))
    long = (POS_30_0B, REF_SEQ[POS_30_0B:POS_30_0B + 3].upper(),
            [anchor], (0, 1))
    records = [short, long]
    output = consensus_runner(
        bed_intervals=[('chr1', POS_30_0B, POS_30_0B + 1)],
        caller_records=[[records[index]] for index in order],
        threshold=threshold,
    )
    output_records = read_records(output)
    assert len(output_records) == expected_count
    if expected_count:
        assert output_records == [
            (30, (short[1], short[2][0]), (0, 1)),
            (30, (long[1], long[2][0]), (0, 1)),
        ]


def test_deletion_does_not_skip_caller_records_inside_span(tmp_path,
                                                           use_accel):
    """A het deletion at POS 100 must not hide a SNP at POS 102.

    Two callers support the heterozygous deletion and two callers (one of them
    also carries the deletion) support the SNP on the other copy.  The old
    jump by the deletion REF length skipped POS 102 and discarded its caller
    records, so the SNP never reached the vote.
    """
    seq = REF_SEQ * 2
    reference = write_reference(tmp_path / 'ref.fa', seq=seq)
    reads = [(f'read{i}', i, seq[i:i + 6]) for i in range(99, 109)]
    bam = write_bam(tmp_path / 'sample.bam', reads, length=len(seq))
    bed_gz = write_bed(tmp_path / 'coverage.bed', [('chr1', 99, 110)])
    deletion_ref = seq[99:103].upper()
    deletion_alt = seq[99].upper()
    snp_ref = seq[101].upper()
    snp_alt = {'A': 'G', 'C': 'T', 'G': 'C', 'T': 'A'}[snp_ref]
    caller_bcfs = [
        write_caller_bcf(tmp_path / 'c0.bcf', 'sample1', [
            (99, deletion_ref, [deletion_alt], (0, 1))],
            length=len(seq)),
        write_caller_bcf(tmp_path / 'c1.bcf', 'sample1', [
            (99, deletion_ref, [deletion_alt], (0, 1)),
            (101, snp_ref, [snp_alt], (0, 1))], length=len(seq)),
        write_caller_bcf(tmp_path / 'c2.bcf', 'sample1', [
            (101, snp_ref, [snp_alt], (0, 1))], length=len(seq)),
    ]
    output = run_consensus(
        tmp_path, bam, bed_gz, caller_bcfs, threshold=2,
        use_accel=use_accel, reference=reference)
    records = read_records(output)
    assert (100, (deletion_ref, deletion_alt), (0, 1)) in records
    assert (102, (snp_ref, snp_alt), (0, 1)) in records


def test_deletion_keeps_all_sites(consensus_runner):
    """Every BED position is present, including those covered by a deletion."""
    deletion_ref = REF_SEQ[POS_20_0B:POS_20_0B + 4].upper()
    deletion_alt = REF_SEQ[POS_20_0B].upper()
    output = consensus_runner(
        bed_intervals=[('chr1', POS_20_0B, 30)],
        caller_records=[[
            (POS_20_0B, deletion_ref, [deletion_alt], (0, 1))]],
        threshold=1,
    )
    spans = read_spans(output)
    assert spans[0][:3] == (20, 23, (deletion_ref, deletion_alt))
    # POS 21..23 are evaluated individually after the deletion record, so the
    # output has one record per BED position and POS..END covers every base.
    assert [start for start, *_ in spans] == list(range(20, 31))
    assert set(span_positions(spans)) == set(range(20, 31))


def test_deletion_does_not_duplicate_or_skip_interval_boundaries(
        consensus_runner):
    """A deletion spanning two BED intervals writes each position once."""
    deletion_ref = REF_SEQ[POS_20_0B:POS_20_0B + 4].upper()
    deletion_alt = REF_SEQ[POS_20_0B].upper()
    output = consensus_runner(
        bed_intervals=[('chr1', POS_20_0B, POS_20_0B + 2),
                       ('chr1', POS_20_0B + 2, 24)],
        caller_records=[[
            (POS_20_0B, deletion_ref, [deletion_alt], (0, 1))]],
        threshold=1,
    )
    assert [pos for pos, _, _ in read_records(output)] == [20, 21, 22, 23, 24]


def _caller_record(ref, alts, gt):
    """Minimal pysam.VariantRecord stand-in for genotype_key tests."""
    from types import SimpleNamespace

    return SimpleNamespace(ref=ref, alts=alts, samples=[{'GT': gt}],
                           filter=SimpleNamespace(keys=lambda: ()))


def test_consensus_calls_keys_alleles_by_ref_and_alt():
    """The documented reproduction: same anchor ALT, different REF."""
    import generate_consensus

    short = generate_consensus.genotype_key(
        [_caller_record('AT', ('A',), (0, 1))], 'A')
    long = generate_consensus.genotype_key(
        [_caller_record('ATT', ('A',), (0, 1))], 'A')
    reference = generate_consensus.genotype_key(
        [_caller_record('A', (), (0, 0))], 'A')
    # Each deletion has one supporter, so threshold 2 leaves only one
    # reference vote: below the threshold, no consensus call at all.
    assert generate_consensus.consensus_calls(
        [short, long, reference], 2) == []
    assert generate_consensus.consensus_calls(
        [long, short, reference], 2) == []
    expected_variants = [
        ('AT', ('A',), (0, 1)),
        ('ATT', ('A',), (0, 1)),
    ]
    assert generate_consensus.consensus_calls(
        [short, long, reference], 1) == expected_variants
    assert generate_consensus.consensus_calls(
        [long, short, reference], 1) == expected_variants


def test_genotype_key_votes_per_ref_alt_allele():
    """Atomised SNP and deletion at one POS stay separate alleles.

    A caller reporting ``A>G 0/1`` and ``AT>A 0/1`` at the same POS supports
    both alleles with one copy each; the ploidy is 2.  Two such callers reach
    threshold 2 and produce one consensus record per allele.
    """
    import generate_consensus

    key = generate_consensus.genotype_key([
        _caller_record('A', ('G',), (0, 1)),
        _caller_record('AT', ('A',), (0, 1)),
    ], 'A')
    assert key == ('A', {('A', 'G'): 1, ('AT', 'A'): 1}, 2)
    assert generate_consensus.consensus_calls([key, key], 2) == [
        ('A', ('G',), (0, 1)),
        ('AT', ('A',), (0, 1)),
    ]


def test_genotype_key_abstains_on_inconsistent_split_copies():
    """Split records whose ALT copies exceed the ploidy abstain.

    A diploid caller cannot be ``1/1`` for the SNP and ``0/1`` for the
    deletion at the same time, so the combination is a no-call instead of a
    fabricated genotype.  Two such callers produce no consensus record.
    """
    import generate_consensus

    key = generate_consensus.genotype_key([
        _caller_record('A', ('G',), (1, 1)),
        _caller_record('AT', ('A',), (0, 1)),
    ], 'A')
    assert key is None
    assert generate_consensus.consensus_calls([key, key], 2) == []


def test_atomized_snp_and_deletion_with_different_refs(tmp_path, use_accel):
    """Both callers support A>G and AT>A at one POS; both reach consensus.

    Atomisation emits the SNP and the deletion as separate records with
    different REF alleles at the same POS.  Keying votes by ``(REF, ALT)``
    keeps the caller from abstaining, so threshold 2 writes one record per
    allele.
    """
    seq = 'AT' * 30
    reference = write_reference(tmp_path / 'ref.fa', seq=seq)
    reads = [(f'read{i}', 18, seq[18:24]) for i in range(2)]
    bam = write_bam(tmp_path / 'sample.bam', reads, length=len(seq))
    bed_gz = write_bed(tmp_path / 'coverage.bed', [('chr1', 20, 21)])
    callers = [
        write_caller_bcf(tmp_path / f'c{index}.bcf', 'sample1', [
            (20, 'A', ['G'], (0, 1)),
            (20, 'AT', ['A'], (0, 1)),
        ], length=len(seq))
        for index in range(2)
    ]
    output = run_consensus(
        tmp_path, bam, bed_gz, callers, threshold=2,
        use_accel=use_accel, reference=reference)
    assert read_records(output) == [
        (21, ('A', 'G'), (0, 1)),
        (21, ('AT', 'A'), (0, 1)),
    ]


# Reference for the pileup edge-case dataset: upper-case bases, then a
# lower-case run and an N run, so both implementations must agree on case and
# unknown bases.
EDGE_REF_SEQ = (REF_SEQ[:30] + REF_SEQ[30:40].lower()
                + 'N' * 10 + REF_SEQ[50:])
EDGE_START0 = 19
EDGE_END0 = 60


def _bases(start, end):
    """Upper-case reference bases at 0-based [start, end)."""
    return EDGE_REF_SEQ[start:end].upper()


def write_edge_dataset(tmp_path):
    """Write a BAM covering the pileup edge cases of the equivalence tests.

    The reads use htslib's 255 missing base quality, overlapping mates whose
    qualities decide the observation, reference skips on both strands,
    insertions and deletions at and after their anchor, and a reference that
    mixes lower-case and N bases.
    """
    reference = write_reference(tmp_path / 'ref.fa', seq=EDGE_REF_SEQ)
    reads = [
        # full-length unique reads; the second has 255 qualities and an ALT
        # base at the first caller record
        ('normal', 19, _bases(19, 60), ((0, 41),), 0, 'I' * 41),
        ('missingq', 19, _bases(19, 29) + 'A' + _bases(30, 60),
         ((0, 41),), 0, None),
        # same-name mates: the higher quality is the observation.  The
        # reference base wins the first pair, the ALT base the second.
        ('mate1', 20, _bases(20, 26), ((0, 6),), 0, 'I' * 6),
        ('mate1', 20, _bases(20, 21) + 'A' + _bases(22, 26),
         ((0, 6),), 0, None),
        ('mate2', 20, _bases(20, 26), ((0, 6),), 0, None),
        ('mate2', 20, _bases(20, 23) + 'A' + _bases(24, 26),
         ((0, 6),), 0, 'I' * 6),
        # reference skips on the forward and reverse strand
        ('skipf', 19, _bases(19, 24) + _bases(34, 40),
         ((0, 5), (3, 10), (0, 6)), 0, 'I' * 11),
        ('skipr', 20, _bases(20, 25) + _bases(35, 41),
         ((0, 5), (3, 10), (0, 6)), 16, 'I' * 11),
        # insertions and deletions at a caller record, and again away from
        # any caller record so the fast path sees the affected columns
        ('ins', 28, _bases(28, 34) + 'AA' + _bases(34, 38),
         ((0, 6), (1, 2), (0, 4)), 0, 'I' * 12),
        ('del', 35, _bases(35, 39) + _bases(41, 46),
         ((0, 4), (2, 2), (0, 5)), 0, 'I' * 9),
        ('del2', 44, _bases(44, 47) + _bases(48, 52),
         ((0, 3), (2, 1), (0, 4)), 0, 'I' * 7),
        ('ins2', 48, _bases(48, 51) + 'TT' + _bases(51, 54),
         ((0, 3), (1, 2), (0, 3)), 0, 'I' * 8),
    ]
    reads.sort(key=lambda entry: entry[1])
    bam = write_bam(tmp_path / 'sample.bam', reads)
    bed_gz = write_bed(
        tmp_path / 'coverage.bed',
        [(CHR, 9, 10), (CHR, EDGE_START0, EDGE_END0), (CHR, 70, 72)])
    callers = [
        write_caller_bcf(tmp_path / 'caller0.bcf', 'sample1', [
            (29, ref_base(30), [alt_base(30)], (0, 1)),
            (33, ref_base(34), [ref_base(34) + 'AA'], (0, 1)),
            (38, EDGE_REF_SEQ[38:41].upper(),
             [EDGE_REF_SEQ[38].upper()], (0, 1)),
            (44, 'N', ['A'], (0, 1)),
            # two deletions with the same anchor ALT but different REF lengths,
            # plus a SNP inside the longer deletion span
            (51, _bases(51, 54), [_bases(51, 52)], (0, 1)),
            (53, ref_base(54), [alt_base(54)], (0, 1)),
        ]),
        write_caller_bcf(tmp_path / 'caller1.bcf', 'sample1', [
            (51, _bases(51, 53), [_bases(51, 52)], (0, 1)),
        ]),
    ]
    return {'reference': reference, 'bam': bam, 'bed': bed_gz,
            'callers': callers}


def test_accel_and_pure_paths_agree(tmp_path):
    """The Cython and pure-pysam implementations produce identical text.

    The dataset exercises 255 base qualities, reference skips, deletions and
    insertions after their anchor, deletions of different lengths with the
    same anchor ALT, a SNP inside a deletion span, and a reference that mixes
    lower-case and N bases.
    """
    if not ACCEL_AVAILABLE:
        pytest.skip('consensus_accel is not built for this interpreter')

    dataset = write_edge_dataset(tmp_path)
    pure = run_consensus(
        tmp_path, dataset['bam'], dataset['bed'], dataset['callers'],
        threshold=1, use_accel=False, output_name='pure.bcf',
        reference=dataset['reference'])
    accel = run_consensus(
        tmp_path, dataset['bam'], dataset['bed'], dataset['callers'],
        threshold=1, use_accel=True, output_name='accel.bcf',
        reference=dataset['reference'])

    def text(path):
        return subprocess.run([BCFTOOLS, 'view', '-H', str(path)],
                              check=True, capture_output=True, text=True).stdout

    assert text(pure) == text(accel)
    assert text(pure)  # not trivially empty
    records = read_records(pure)
    assert (30, (ref_base(30), alt_base(30)), (0, 1)) in records
    assert (34, (ref_base(34), ref_base(34) + 'AA'), (0, 1)) in records
    assert (39, (EDGE_REF_SEQ[38:41].upper(),
                 EDGE_REF_SEQ[38].upper()), (0, 1)) in records
    assert (45, ('N', 'A'), (0, 1)) in records
    # The 1-bp and 2-bp deletions are distinct alleles at POS 52, and the SNP
    # inside the longer deletion span at POS 54 is still written.
    assert (52, (_bases(51, 53), _bases(51, 52)), (0, 1)) in records
    assert (52, (_bases(51, 54), _bases(51, 52)), (0, 1)) in records
    assert (54, (ref_base(54), alt_base(54)), (0, 1)) in records


class _LimitedTally:
    """PileupTally facade with a deliberately small scratch table."""

    def __init__(self, tally, max_slots):
        self.tally = tally
        self.max_slots = max_slots
        self.fallbacks = 0

    def sync(self, pos1):
        return self.tally.sync(pos1)

    def simple(self, ref_base):
        counts = self.tally.simple(ref_base, self.max_slots)
        if counts is None:
            self.fallbacks += 1
        return counts

    def entries(self):
        return self.tally.entries()


def test_accel_and_pure_paths_agree_beyond_scratch_table(tmp_path):
    """A column too deep for the C scratch table falls back to Python."""
    if not ACCEL_AVAILABLE:
        pytest.skip('consensus_accel is not built for this interpreter')

    import consensus_accel
    import generate_consensus
    import pysam

    dataset = write_edge_dataset(tmp_path)
    reference = EDGE_REF_SEQ.upper()
    # Each walker owns an independent pileup iterator: pysam's columns point
    # into the read buffer of their file handle.
    with pysam.AlignmentFile(str(dataset['bam']), 'rb') as pure_file, \
            pysam.AlignmentFile(str(dataset['bam']), 'rb') as accel_file:
        def pileup(handle):
            return handle.pileup(
                CHR, EDGE_START0, EDGE_END0, truncate=True, stepper='all',
                min_base_quality=0, max_depth=1000000)

        pure = generate_consensus.ColumnWalker(pileup(pure_file))
        limited = _LimitedTally(
            consensus_accel.PileupTally(pileup(accel_file)), max_slots=4)
        accel = generate_consensus.AccelWalker(limited)
        for pos1 in range(EDGE_START0 + 1, EDGE_END0 + 1):
            has_column = accel.sync(pos1)
            assert has_column == pure.sync(pos1)
            if not has_column:
                continue
            ref = reference[pos1 - 1]
            assert accel.reference_support(ref) == pure.reference_support(ref)
            assert accel.site() == pure.site()
    # At least one column was deeper than the reduced table and took the
    # Python fallback instead of the C tally.
    assert limited.fallbacks


MAX_POSITION = (1 << 31) - 1


def _sparse_files_supported(directory):
    """True when a seek-created hole is stored compactly."""
    probe = directory / '.sparse-probe'
    with open(probe, 'wb') as handle:
        handle.seek(1 << 20)
        handle.write(b'x')
    try:
        return probe.stat().st_blocks * 512 < (1 << 19)
    finally:
        probe.unlink()


def test_accel_and_pure_paths_agree_near_position_limit(tmp_path):
    """Positions and columns just below 2**31 agree byte for byte."""
    if not ACCEL_AVAILABLE:
        pytest.skip('consensus_accel is not built for this interpreter')
    if not _sparse_files_supported(tmp_path):
        pytest.skip('filesystem does not support sparse files')

    contig = 'chrBig'
    length = MAX_POSITION
    seq = 'ACGTACGTACGTAC'
    start0 = length - len(seq)
    reference = write_sparse_reference(
        tmp_path / 'ref.fa', contig, length, start0, seq)

    # Three staggered reads cover all but the last two BED bases, so the walk
    # also asks for columns past the end of the pileup near the limit.
    reads = []
    for index in range(3):
        reads.append((f'read{index}', start0 + index, seq[index:index + 10],
                      ((0, 10),), 0, 'I' * 10))
    bam = write_bam(tmp_path / 'sample.bam', reads, contig=contig,
                    length=length, csi=True)
    bed_gz = write_bed(tmp_path / 'coverage.bed',
                       [(contig, start0, length)])
    alt = 'A' if seq[5] != 'A' else 'C'
    callers = [
        write_caller_bcf(tmp_path / 'caller0.bcf', 'sample1',
                         [(start0 + 5, seq[5], [alt], (0, 1))],
                         contig=contig, length=length),
        write_caller_bcf(tmp_path / 'caller1.bcf', 'sample1', [],
                         contig=contig, length=length),
    ]

    pure = run_consensus(
        tmp_path, bam, bed_gz, callers, threshold=1, use_accel=False,
        output_name='pure.bcf', reference=reference)
    accel = run_consensus(
        tmp_path, bam, bed_gz, callers, threshold=1, use_accel=True,
        output_name='accel.bcf', reference=reference)

    def text(path):
        return subprocess.run([BCFTOOLS, 'view', '-H', str(path)],
                              check=True, capture_output=True, text=True).stdout

    assert text(pure) == text(accel)
    assert f'{contig}\t{MAX_POSITION}\t' in text(pure)


def test_provenance_headers_list_used_and_missing_callers(tmp_path, use_accel):
    """The structured header describes the vote and its inputs."""
    import pysam

    reference = write_reference(tmp_path / 'ref.fa')
    reads = [(f'read{i}', i, REF_SEQ[i:i + 6]) for i in range(19, 24)]
    bam = write_bam(tmp_path / 'sample.bam', reads)
    bed_gz = write_bed(tmp_path / 'coverage.bed',
                       [('chr1', POS_30_0B, POS_30_0B + 1)])
    caller_bcfs = [
        write_caller_bcf(tmp_path / 'bcftools.bcf', 'sample1',
                         [(POS_30_0B, ref_base(30), [alt_base(30)], (0, 1))]),
        write_caller_bcf(tmp_path / 'gatk.bcf', 'sample1',
                         [(POS_30_0B, ref_base(30), [alt_base(30)], (0, 1))]),
    ]
    output = run_consensus(
        tmp_path, bam, bed_gz, caller_bcfs, threshold=2,
        callers=['bcftools', 'gatk'],
        configured_callers=['bcftools', 'freebayes', 'gatk'],
        use_accel=use_accel, reference=reference)
    with pysam.VariantFile(str(output), index_filename=None) as variant_file:
        header = str(variant_file.header)
        records = list(variant_file)
    expected = ('##ChoCallate=<Version=test,Threshold=2,Callers="bcftools,gatk",'
                'Rule=two-stage,Tie=nocall,AbsentVote=ref,'
                'FilteredVote=abstain,ReferenceBlocks=False,'
                'CallersFailed="freebayes">')
    assert expected in header
    assert [(record.pos, (record.ref,) + tuple(record.alts),
             tuple(record.samples[0]['GT']))
            for record in records] == [
        (30, (ref_base(30), alt_base(30)), (0, 1))]



def test_threshold_above_surviving_callers_fails(tmp_path, use_accel):
    """Degraded mode refuses a threshold no surviving caller can reach."""
    reference = write_reference(tmp_path / 'ref.fa')
    reads = [(f'read{i}', i, REF_SEQ[i:i + 6]) for i in range(19, 24)]
    bam = write_bam(tmp_path / 'sample.bam', reads)
    bed_gz = write_bed(tmp_path / 'coverage.bed',
                       [('chr1', POS_30_0B, POS_30_0B + 1)])
    caller = write_caller_bcf(tmp_path / 'bcftools.bcf', 'sample1',
                              [(POS_30_0B, ref_base(30), [alt_base(30)], (0, 1))])
    with pytest.raises(subprocess.CalledProcessError) as excinfo:
        run_consensus(tmp_path, bam, bed_gz, [caller], threshold=2,
                      callers=['bcftools'],
                      configured_callers=['bcftools', 'freebayes'],
                      use_accel=use_accel, reference=reference)
    assert 'surviving caller' in excinfo.value.stderr


def test_multiple_alts_are_emitted_regardless_of_caller_order(tmp_path, use_accel):
    """Distinct supported ALTs become separate order-independent records."""
    reference = write_reference(tmp_path / 'ref.fa')
    reads = [(f'read{i}', i, REF_SEQ[i:i + 6]) for i in range(19, 24)]
    bam = write_bam(tmp_path / 'sample.bam', reads)
    bed_gz = write_bed(tmp_path / 'coverage.bed',
                       [('chr1', POS_30_0B, POS_30_0B + 1)])
    alt1 = alt_base(30)
    alt2 = next(base for base in 'ACGT'
                if base not in (ref_base(30), alt1))
    caller_bcfs = [
        write_caller_bcf(tmp_path / 'bcftools.bcf', 'sample1',
                         [(POS_30_0B, ref_base(30), [alt1], (0, 1))]),
        write_caller_bcf(tmp_path / 'freebayes.bcf', 'sample1',
                         [(POS_30_0B, ref_base(30), [alt2], (0, 1))]),
    ]

    def consensus_records(order, output_name):
        output = run_consensus(
            tmp_path, bam, bed_gz, caller_bcfs, threshold=1,
            callers=['bcftools', 'freebayes'], configured_callers=order,
            use_accel=use_accel, reference=reference,
            output_name=output_name)
        return sorted(read_records(output))

    expected = sorted([
        (30, (ref_base(30), alt1), (0, 1)),
        (30, (ref_base(30), alt2), (0, 1)),
    ])
    assert consensus_records(['bcftools', 'freebayes'], 'first.bcf') == expected
    assert consensus_records(['freebayes', 'bcftools'], 'second.bcf') == expected


@pytest.mark.parametrize('callers', [
    ('bcftools',),
    ('freebayes',),
    ('gatk',),
    ('bcftools', 'freebayes'),
    ('bcftools', 'gatk'),
    ('freebayes', 'gatk'),
    ('bcftools', 'freebayes', 'gatk'),
])
def test_any_caller_subset(tmp_path, use_accel, callers):
    """Any non-empty caller subset produces a consensus and a header."""
    import pysam

    reference = write_reference(tmp_path / 'ref.fa')
    reads = [(f'read{i}', i, REF_SEQ[i:i + 6]) for i in range(19, 24)]
    bam = write_bam(tmp_path / 'sample.bam', reads)
    bed_gz = write_bed(tmp_path / 'coverage.bed',
                       [('chr1', POS_30_0B, POS_30_0B + 1)])
    caller_bcfs = [
        write_caller_bcf(tmp_path / f'{caller}.bcf', 'sample1',
                         [(POS_30_0B, ref_base(30), [alt_base(30)], (0, 1))])
        for caller in callers]
    output = run_consensus(tmp_path, bam, bed_gz, caller_bcfs,
                           threshold=len(callers), callers=list(callers),
                           use_accel=use_accel, reference=reference)
    with pysam.VariantFile(str(output), index_filename=None) as variant_file:
        header = str(variant_file.header)
        records = list(variant_file)
    assert f'Callers="{",".join(callers)}"' in header
    assert 'Rule=two-stage' in header
    assert [(record.pos, (record.ref,) + tuple(record.alts),
             tuple(record.samples[0]['GT']))
            for record in records] == [
        (30, (ref_base(30), alt_base(30)), (0, 1))]


# Two-stage decision table with three configured callers.  ``h`` is a 0/1 call,
# ``m`` a 1/1 call and ``-`` a caller without a record.  The expected values
# use the defaults: absent_vote=ref and tie_policy=nocall.
VOTING_TABLE = {
    ('--h', 1): 'variant',
    ('--h', 2): 'reference',
    ('--h', 3): None,
    ('h--', 1): 'variant',
    ('h--', 2): 'reference',
    ('h--', 3): None,
    ('hm-', 1): None,
    ('hm-', 2): None,
    ('hm-', 3): None,
    ('-hm', 1): None,
    ('-hm', 2): None,
    ('-hm', 3): None,
}


def _voting_caller_records(placement, alt):
    """Caller records for a row of the decision table."""
    records = []
    for marker in placement:
        if marker == 'h':
            records.append([(POS_30_0B, ref_base(30), [alt], (0, 1))])
        elif marker == 'm':
            records.append([(POS_30_0B, ref_base(30), [alt], (1, 1))])
        else:
            records.append([])
    return records


@pytest.mark.parametrize('placement,threshold', VOTING_TABLE)
def test_voting_table(consensus_runner, placement, threshold):
    """The two-stage vote reproduces the documented table."""
    alt = alt_base(30)
    output = consensus_runner(
        bed_intervals=[('chr1', POS_30_0B, POS_30_0B + 1)],
        caller_records=_voting_caller_records(placement, alt),
        threshold=threshold,
    )
    expected = VOTING_TABLE[(placement, threshold)]
    if expected == 'variant':
        assert read_records(output) == [
            (30, (ref_base(30), alt), (0, 1))]
    elif expected == 'reference':
        assert read_records(output) == [
            (30, (ref_base(30), '<*>'), (0, 0))]
    else:
        assert read_records(output) == []


@pytest.mark.parametrize('placement', ['hm-', '-hm'])
def test_genotype_tie_with_first_policy(consensus_runner, placement):
    """``first`` resolves a zygosity tie by configured caller order."""
    alt = alt_base(30)
    output = consensus_runner(
        bed_intervals=[('chr1', POS_30_0B, POS_30_0B + 1)],
        caller_records=_voting_caller_records(placement, alt),
        threshold=1,
        tie_policy='first',
    )
    assert read_records(output) == [(30, (ref_base(30), alt), (0, 1))]


def test_absent_vote_abstain_requires_coverage(consensus_runner):
    """An absent caller only votes ref when the position is covered."""
    alt = alt_base(24)
    covered_bed = [('chr1', POS_24_0B, POS_24_0B + 1)]
    covered_variant = [(POS_24_0B, ref_base(24), [alt], (0, 1))]

    reference_vote = consensus_runner(
        bed_intervals=covered_bed,
        caller_records=[covered_variant, [], []],
        threshold=2, absent_vote='ref')
    assert read_records(reference_vote) == [
        (24, (ref_base(24), '<*>'), (0, 0))]

    abstain_with_coverage = consensus_runner(
        bed_intervals=covered_bed,
        caller_records=[covered_variant, [], []],
        threshold=2, absent_vote='abstain', min_coverage=1,
        output_name='abstain_covered.bcf')
    assert read_records(abstain_with_coverage) == [
        (24, (ref_base(24), '<*>'), (0, 0))]

    abstain_too_strict = consensus_runner(
        bed_intervals=covered_bed,
        caller_records=[covered_variant, [], []],
        threshold=2, absent_vote='abstain', min_coverage=100,
        output_name='abstain_strict.bcf')
    assert read_records(abstain_too_strict) == []

    # A BED position without a pileup column has no coverage evidence, so an
    # absent caller abstains even with min_coverage=0.
    uncovered_bed = [('chr1', POS_30_0B, POS_30_0B + 1)]
    uncovered_variant = [(POS_30_0B, ref_base(30), [alt_base(30)], (0, 1))]
    abstain_uncovered = consensus_runner(
        bed_intervals=uncovered_bed,
        caller_records=[uncovered_variant, [], []],
        threshold=2, absent_vote='abstain', min_coverage=0,
        output_name='abstain_uncovered.bcf')
    assert read_records(abstain_uncovered) == []

    ref_uncovered = consensus_runner(
        bed_intervals=uncovered_bed,
        caller_records=[uncovered_variant, [], []],
        threshold=2, absent_vote='ref', output_name='ref_uncovered.bcf')
    assert read_records(ref_uncovered) == [
        (30, (ref_base(30), '<*>'), (0, 0))]


def test_filtered_record_abstains_from_the_vote(consensus_runner):
    """A soft-filtered record is not a reference vote.

    Caller A calls 1/1 with FILTER=PASS, caller B has the same variant with
    FILTER=LowQual and caller C has no record.  With threshold 2 the filtered
    record abstains (the default), so the position has neither a variant nor a
    hom-reference record.  Counting it as reference would produce a confident
    hom-reference call.
    """
    alt = alt_base(30)
    output = consensus_runner(
        bed_intervals=[('chr1', POS_30_0B, POS_30_0B + 1)],
        caller_records=[
            [(POS_30_0B, ref_base(30), [alt], (1, 1), None, 'PASS')],
            [(POS_30_0B, ref_base(30), [alt], (1, 1), None, 'LowQual')],
            [],
        ],
        threshold=2,
    )
    assert read_records(output) == []


def test_filtered_vote_ref_counts_as_reference(consensus_runner):
    """``filtered_vote=ref`` turns a filtered record into a reference vote."""
    alt = alt_base(30)
    output = consensus_runner(
        bed_intervals=[('chr1', POS_30_0B, POS_30_0B + 1)],
        caller_records=[
            [(POS_30_0B, ref_base(30), [alt], (1, 1), None, 'PASS')],
            [(POS_30_0B, ref_base(30), [alt], (1, 1), None, 'LowQual')],
            [],
        ],
        threshold=2,
        filtered_vote='ref',
    )
    assert read_records(output) == [
        (30, (ref_base(30), '<*>'), (0, 0))]


def test_header_describes_abstain_and_tie_policies(tmp_path, use_accel):
    """The single header line records the applied voting rule."""
    import pysam

    reference = write_reference(tmp_path / 'ref.fa')
    reads = [(f'read{i}', i, REF_SEQ[i:i + 6]) for i in range(19, 24)]
    bam = write_bam(tmp_path / 'sample.bam', reads)
    bed_gz = write_bed(tmp_path / 'coverage.bed',
                       [('chr1', POS_30_0B, POS_30_0B + 1)])
    caller_bcfs = [
        write_caller_bcf(tmp_path / 'bcftools.bcf', 'sample1',
                         [(POS_30_0B, ref_base(30), [alt_base(30)], (0, 1))]),
        write_caller_bcf(tmp_path / 'gatk.bcf', 'sample1',
                         [(POS_30_0B, ref_base(30), [alt_base(30)], (0, 1))]),
    ]
    output = run_consensus(
        tmp_path, bam, bed_gz, caller_bcfs, threshold=2,
        callers=['bcftools', 'gatk'],
        configured_callers=['bcftools', 'freebayes', 'gatk'],
        absent_vote='abstain', tie_policy='first', min_coverage=3,
        use_accel=use_accel, reference=reference)
    with pysam.VariantFile(str(output), index_filename=None) as variant_file:
        header = str(variant_file.header)
    expected = ('##ChoCallate=<Version=test,Threshold=2,'
                'Callers="bcftools,gatk",Rule=two-stage,Tie=first,'
                'AbsentVote=abstain,MinCoverage=3,'
                'FilteredVote=abstain,ReferenceBlocks=False,'
                'CallersFailed="freebayes">')
    assert expected in header


def test_split_multiallelic_records_are_aggregated(tmp_path, use_accel):
    """All records of a caller at one position form one vote."""
    reference = write_reference(tmp_path / 'ref.fa')
    reads = [(f'read{i}', i, REF_SEQ[i:i + 6]) for i in range(19, 24)]
    bam = write_bam(tmp_path / 'sample.bam', reads)
    bed_gz = write_bed(tmp_path / 'coverage.bed',
                       [('chr1', POS_30_0B, POS_30_0B + 1)])
    alt1 = alt_base(30)
    alt2 = next(base for base in 'ACGT'
                if base not in (ref_base(30), alt1))
    # Caller 0 carries two split records.  The second one (alt1) must not be
    # skipped: together with caller 1 it reaches threshold 2.
    caller0 = write_caller_bcf(tmp_path / 'caller0.bcf', 'sample1', [
        (POS_30_0B, ref_base(30), [alt2], (0, 1)),
        (POS_30_0B, ref_base(30), [alt1], (0, 1)),
    ])
    caller1 = write_caller_bcf(tmp_path / 'caller1.bcf', 'sample1', [
        (POS_30_0B, ref_base(30), [alt1], (0, 1)),
    ])
    output = run_consensus(tmp_path, bam, bed_gz, [caller0, caller1],
                           threshold=2, use_accel=use_accel,
                           reference=reference)
    assert read_records(output) == [
        (30, (ref_base(30), alt1), (0, 1))]


@pytest.mark.parametrize('split_records', [True, False])
def test_triallelic_caller_produces_split_consensus(tmp_path, use_accel,
                                                    split_records):
    """A 1/2 call produces one consensus record per ALT."""
    reference = write_reference(tmp_path / 'ref.fa')
    reads = [(f'read{i}', i, REF_SEQ[i:i + 6]) for i in range(19, 24)]
    bam = write_bam(tmp_path / 'sample.bam', reads)
    bed_gz = write_bed(tmp_path / 'coverage.bed',
                       [('chr1', POS_30_0B, POS_30_0B + 1)])
    alt1 = alt_base(30)
    alt2 = next(base for base in 'ACGT'
                if base not in (ref_base(30), alt1))
    if split_records:
        records = [(POS_30_0B, ref_base(30), [alt1], (0, 1)),
                   (POS_30_0B, ref_base(30), [alt2], (0, 1))]
    else:
        records = [(POS_30_0B, ref_base(30), [alt1, alt2], (1, 2))]
    caller = write_caller_bcf(tmp_path / 'caller.bcf', 'sample1', records)
    output = run_consensus(tmp_path, bam, bed_gz, [caller], threshold=1,
                           use_accel=use_accel, reference=reference)
    assert sorted(read_records(output)) == sorted([
        (30, (ref_base(30), alt1), (0, 1)),
        (30, (ref_base(30), alt2), (0, 1)),
    ])


def test_atomized_mnp_matches_snp(tmp_path, use_accel):
    """An atomised freebayes-style MNP votes with the matching SNP."""
    reference = write_reference(tmp_path / 'ref.fa')
    reads = [(f'read{i}', i, REF_SEQ[i:i + 6]) for i in range(19, 24)]
    bam = write_bam(tmp_path / 'sample.bam', reads)
    bed_gz = write_bed(tmp_path / 'coverage.bed',
                       [('chr1', POS_30_0B, POS_30_0B + 3)])
    mnp_ref = REF_SEQ[POS_30_0B:POS_30_0B + 3].upper()
    mnp_alt = ''.join({'A': 'T', 'C': 'A', 'G': 'T', 'T': 'G'}[base]
                      for base in mnp_ref)
    mnp = write_caller_bcf(tmp_path / 'mnp.bcf', 'sample1', [
        (POS_30_0B, mnp_ref, [mnp_alt], (0, 1)),
    ])
    normalized = tmp_path / 'mnp.normalized.bcf'
    subprocess.run([BCFTOOLS, 'norm', '-f', str(reference), '-m', '-any',
                    '--atomize', '-Ob', '-o', str(normalized), str(mnp)],
                   check=True)
    subprocess.run([BCFTOOLS, 'index', '--csi', str(normalized)], check=True)
    snp = write_caller_bcf(tmp_path / 'snp.bcf', 'sample1', [
        (POS_30_0B, ref_base(30), [alt_base(30)], (0, 1)),
    ])
    output = run_consensus(tmp_path, bam, bed_gz, [normalized, snp],
                           threshold=1, use_accel=use_accel,
                           reference=reference)
    assert (30, (ref_base(30), alt_base(30)), (0, 1)) in read_records(output)


def test_norm_dedup_removes_chunk_boundary_duplicates(tmp_path):
    """Exact duplicates from chunk boundaries are removed."""
    import pysam

    reference = write_reference(tmp_path / 'ref.fa')
    duplicated = write_caller_bcf(tmp_path / 'duplicated.bcf', 'sample1', [
        (POS_30_0B, ref_base(30), [alt_base(30)], (0, 1)),
        (POS_30_0B, ref_base(30), [alt_base(30)], (0, 1)),
    ])
    output = tmp_path / 'deduplicated.bcf'
    chain = (
        f'{BCFTOOLS} norm -f {reference} -m -any --atomize -Ou {duplicated}'
        f' | {BCFTOOLS} sort -Ou'
        f' | {BCFTOOLS} norm -d exact -Ob -o {output}'
    )
    subprocess.run(['bash', '-c', chain], check=True)
    with pysam.VariantFile(str(output), index_filename=None) as variant_file:
        records = list(variant_file)
    assert len(records) == 1
    assert (records[0].pos, records[0].ref,
            tuple(records[0].alts)) == (30, ref_base(30), (alt_base(30),))


def test_split_consensus_records_join_for_output(tmp_path, use_accel):
    """Per-ALT consensus records re-join when split output is off."""
    import pysam

    reference = write_reference(tmp_path / 'ref.fa')
    reads = [(f'read{i}', i, REF_SEQ[i:i + 6]) for i in range(19, 24)]
    bam = write_bam(tmp_path / 'sample.bam', reads)
    bed_gz = write_bed(tmp_path / 'coverage.bed',
                       [('chr1', POS_30_0B, POS_30_0B + 1)])
    alt1 = alt_base(30)
    alt2 = next(base for base in 'ACGT'
                if base not in (ref_base(30), alt1))
    caller = write_caller_bcf(tmp_path / 'caller.bcf', 'sample1', [
        (POS_30_0B, ref_base(30), [alt1], (0, 1)),
        (POS_30_0B, ref_base(30), [alt2], (0, 1)),
    ])
    split = run_consensus(tmp_path, bam, bed_gz, [caller], threshold=1,
                          use_accel=use_accel, reference=reference,
                          output_name='split.bcf')
    joined = tmp_path / 'joined.bcf'
    subprocess.run([BCFTOOLS, 'norm', '-m', '+any', '-Ob', '-o', str(joined),
                    str(split)], check=True)
    with pysam.VariantFile(str(joined), index_filename=None) as variant_file:
        records = list(variant_file)
    assert len(records) == 1
    record = records[0]
    assert record.pos == 30
    assert record.ref == ref_base(30)
    assert tuple(record.alts) == tuple(sorted((alt1, alt2)))
    assert tuple(sorted(record.samples[0]['GT'])) == (1, 2)
