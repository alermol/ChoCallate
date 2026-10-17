"""End-to-end regression tests for bin/generate_consensus.py.

Covers the correctness bugs from the code review (ChoCallate_review.md):

* C-01: the last base of every BED interval was dropped (1-bp intervals produced
  no output at all).
* C-02: a caller record with a missing genotype (``./.``, ``0/.``) crashed with
  ``TypeError``; haploid hom-reference calls were mis-keyed as variants.

Every test runs twice, on the pure-pysam path and (when built) on the Cython
path, and asserts the documented all-sites behaviour: one output record per
evaluated BED base.
"""

import subprocess

import pytest

from conftest import (  # noqa: I100,I202  (tests run from this directory)
    ACCEL_AVAILABLE,
    BCFTOOLS,
    REF_SEQ,
    alt_base,
    ref_base,
    run_consensus,
    write_bam,
    write_bed,
    write_caller_bcf,
    write_reference,
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


def test_interval_ends_are_written(consensus_runner):
    """C-01: every base of every BED interval, last base included."""
    # chr1 9 10 is a single base (1-based 10); chr1 19 25 covers 20..25.
    output = consensus_runner(
        bed_intervals=[('chr1', 9, 10), ('chr1', 19, 25)],
        caller_records=[[], []],
        threshold=2,
    )
    assert [pos for pos, _, _ in read_records(output)] == [10, 20, 21, 22, 23, 24, 25]


def test_adjacent_intervals_are_not_truncated(consensus_runner):
    """C-01: adjacent intervals keep their last bases (1-based 10 and 20)."""
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


def test_missing_genotype_does_not_crash(consensus_runner):
    """C-02: a ``./.`` caller record is handled instead of raising TypeError."""
    output = consensus_runner(
        bed_intervals=[('chr1', POS_20_0B, 25)],
        caller_records=[[(POS_20_0B, REF_SEQ[POS_20_0B], [], (None, None))], [], []],
        threshold=2,
    )
    # The no-call abstains; the two other callers' placeholder votes reach the
    # threshold, so the position is written as hom-reference.
    assert [pos for pos, _, _ in read_records(output)] == [20, 21, 22, 23, 24, 25]


def test_missing_genotype_abstains_from_the_vote(consensus_runner):
    """C-02: a no-call record is not a hom-reference vote.

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
    """C-02: a no-call must not out-vote a real call from another caller."""
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
    """C-02: a ``0/.`` record abstains, it is not a hom-reference vote."""
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
    """C-02: a haploid ``0`` is a hom-reference vote, not a variant."""
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
    """C-02: a haploid ``1`` does not crash and is reported as a variant."""
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
    """C-07 option 1: the <*> AD slot holds non-reference read depth."""
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

    The genotype classes (0/1 and 1/1) tie, so the explicit ``qual`` policy
    picks the class with the larger summed QUAL instead of relying on caller
    order.
    """
    output = consensus_runner(
        bed_intervals=[('chr1', POS_30_0B, 31), ('chr1', 40, 42)],
        caller_records=[
            [(POS_30_0B, ref_base(30), [alt_base(30)], (0, 1), 50)],
            [(POS_30_0B, ref_base(30), [alt_base(30)], (1, 1), 20)],
        ],
        threshold=1,
        tie_policy='qual',
    )
    records = read_records(output)
    # The two callers agree on the ALT allele, so the site passes the site
    # stage; the 0/1 class has the higher summed QUAL.  Placeholder votes must
    # not win, and the interval ends are written.
    assert records[0] == (30, (ref_base(30), alt_base(30)), (0, 1))
    assert [pos for pos, _, _ in records] == [30, 31, 41, 42]


def test_accel_and_pure_paths_agree(tmp_path):
    """T-02: the Cython and pure-pysam implementations produce identical text."""
    if not ACCEL_AVAILABLE:
        pytest.skip('consensus_accel is not built for this interpreter')

    reference = write_reference(tmp_path / 'ref.fa')
    # overlapping mates (same query name) and a position without coverage
    reads = [('r1', 19, REF_SEQ[19:25]), ('r1', 20, REF_SEQ[20:26]),
             ('r2', 29, REF_SEQ[29:35])]
    bam = write_bam(tmp_path / 'sample.bam', reads)
    bed_gz = write_bed(tmp_path / 'coverage.bed',
                       [('chr1', 9, 10), ('chr1', 19, 25), ('chr1', 29, 36)])
    caller_bcfs = [
        write_caller_bcf(tmp_path / 'caller0.bcf', 'sample1',
                         [(20, REF_SEQ[20], [], (0, 0)),
                          (POS_30_0B, ref_base(30), [alt_base(30)], (0, 1))]),
        write_caller_bcf(tmp_path / 'caller1.bcf', 'sample1', []),
    ]
    pure = run_consensus(tmp_path, bam, bed_gz, caller_bcfs, threshold=1,
                         use_accel=False, output_name='pure.bcf',
                         reference=reference)
    accel = run_consensus(tmp_path, bam, bed_gz, caller_bcfs, threshold=1,
                          use_accel=True, output_name='accel.bcf',
                          reference=reference)

    def text(path):
        return subprocess.run([BCFTOOLS, 'view', '-H', str(path)],
                              check=True, capture_output=True, text=True).stdout

    assert text(pure) == text(accel)
    assert text(pure)  # not trivially empty


def test_provenance_headers_list_used_and_missing_callers(tmp_path, use_accel):
    """D.2/C-04: the structured header describes the vote and its inputs."""
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
                'Rule=two-stage,Tie=nocall,AbsentVote=ref,RemoveInvariant=False,'
                'SplitMultiallelic=False,CallersFailed="freebayes">')
    assert expected in header
    assert [(record.pos, (record.ref,) + tuple(record.alts),
             tuple(record.samples[0]['GT']))
            for record in records] == [
        (30, (ref_base(30), alt_base(30)), (0, 1))]


def test_threshold_above_surviving_callers_fails(tmp_path, use_accel):
    """D.2: degraded mode refuses a threshold no surviving caller can reach."""
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


def test_caller_order_only_matters_with_first_policy(tmp_path, use_accel):
    """D.3: the default is order-independent; ``first`` opts into caller order."""
    reference = write_reference(tmp_path / 'ref.fa')
    reads = [(f'read{i}', i, REF_SEQ[i:i + 6]) for i in range(19, 24)]
    bam = write_bam(tmp_path / 'sample.bam', reads)
    bed_gz = write_bed(tmp_path / 'coverage.bed',
                       [('chr1', POS_30_0B, POS_30_0B + 1)])
    alt1 = alt_base(30)
    alt2 = {'A': 'C', 'C': 'A', 'G': 'A', 'T': 'A'}[ref_base(30)]
    caller_bcfs = [
        write_caller_bcf(tmp_path / 'bcftools.bcf', 'sample1',
                         [(POS_30_0B, ref_base(30), [alt1], (0, 1))]),
        write_caller_bcf(tmp_path / 'freebayes.bcf', 'sample1',
                         [(POS_30_0B, ref_base(30), [alt2], (0, 1))]),
    ]

    def consensus_alleles(order, policy, output_name):
        output = run_consensus(
            tmp_path, bam, bed_gz, caller_bcfs, threshold=1,
            callers=['bcftools', 'freebayes'], configured_callers=order,
            tie_policy=policy, use_accel=use_accel, reference=reference,
            output_name=output_name)
        records = read_records(output)
        return records[0][1] if records else None

    # Two distinct ALT alleles with equal support: the default no-call policy
    # is independent of the configured caller order.
    assert consensus_alleles(['bcftools', 'freebayes'], 'nocall', 'nocall1.bcf') is None
    assert consensus_alleles(['freebayes', 'bcftools'], 'nocall', 'nocall2.bcf') is None
    # The explicit ``first`` policy is the only one that uses caller order.
    assert consensus_alleles(['bcftools', 'freebayes'], 'first', 'first1.bcf') == \
        (ref_base(30), alt1)
    assert consensus_alleles(['freebayes', 'bcftools'], 'first', 'first2.bcf') == \
        (ref_base(30), alt2)


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
    """D.2: any non-empty caller subset produces a consensus and a header."""
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


# C-03 decision table with three configured callers.  ``h`` is a 0/1 call,
# ``m`` a 1/1 call and ``-`` a caller without a record.  The expected values
# use the D.3 defaults: absent_vote=ref and tie_policy=nocall.
C03_TABLE = {
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


def _c03_caller_records(placement, alt):
    """Caller records for a row of the C-03 table."""
    records = []
    for marker in placement:
        if marker == 'h':
            records.append([(POS_30_0B, ref_base(30), [alt], (0, 1))])
        elif marker == 'm':
            records.append([(POS_30_0B, ref_base(30), [alt], (1, 1))])
        else:
            records.append([])
    return records


@pytest.mark.parametrize('placement,threshold', C03_TABLE)
def test_c03_voting_table(consensus_runner, placement, threshold):
    """D.3/C-03: the two-stage vote reproduces the documented table."""
    alt = alt_base(30)
    output = consensus_runner(
        bed_intervals=[('chr1', POS_30_0B, POS_30_0B + 1)],
        caller_records=_c03_caller_records(placement, alt),
        threshold=threshold,
    )
    expected = C03_TABLE[(placement, threshold)]
    if expected == 'variant':
        assert read_records(output) == [
            (30, (ref_base(30), alt), (0, 1))]
    elif expected == 'reference':
        assert read_records(output) == [
            (30, (ref_base(30), '<*>'), (0, 0))]
    else:
        assert read_records(output) == []


@pytest.mark.parametrize('placement', ['hm-', '-hm'])
def test_c03_genotype_tie_with_first_policy(consensus_runner, placement):
    """D.3: ``first`` resolves a zygosity tie by configured caller order."""
    alt = alt_base(30)
    output = consensus_runner(
        bed_intervals=[('chr1', POS_30_0B, POS_30_0B + 1)],
        caller_records=_c03_caller_records(placement, alt),
        threshold=1,
        tie_policy='first',
    )
    assert read_records(output) == [(30, (ref_base(30), alt), (0, 1))]


def test_tie_policy_qual_uses_summed_quality(consensus_runner):
    """D.3: ``qual`` picks the genotype class with the larger summed QUAL."""
    alt = alt_base(30)
    bed = [('chr1', POS_30_0B, POS_30_0B + 1)]

    het_wins = consensus_runner(
        bed_intervals=bed,
        caller_records=[
            [(POS_30_0B, ref_base(30), [alt], (0, 1), 90)],
            [(POS_30_0B, ref_base(30), [alt], (1, 1), 20)],
        ],
        threshold=1, tie_policy='qual')
    assert read_records(het_wins) == [(30, (ref_base(30), alt), (0, 1))]

    hom_wins = consensus_runner(
        bed_intervals=bed,
        caller_records=[
            [(POS_30_0B, ref_base(30), [alt], (0, 1), 20)],
            [(POS_30_0B, ref_base(30), [alt], (1, 1), 90)],
        ],
        threshold=1, tie_policy='qual', output_name='hom_wins.bcf')
    assert read_records(hom_wins) == [(30, (ref_base(30), alt), (1, 1))]

    tied = consensus_runner(
        bed_intervals=bed,
        caller_records=[
            [(POS_30_0B, ref_base(30), [alt], (0, 1), 50)],
            [(POS_30_0B, ref_base(30), [alt], (1, 1), 50)],
        ],
        threshold=1, tie_policy='qual', output_name='tied.bcf')
    assert read_records(tied) == []


def test_absent_vote_abstain_requires_coverage(consensus_runner):
    """D.3: an absent caller only votes ref when the position is covered."""
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


def test_header_describes_abstain_and_qual_policies(tmp_path, use_accel):
    """D.3: the single header line records the applied voting rule."""
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
        absent_vote='abstain', tie_policy='qual', min_coverage=3,
        use_accel=use_accel, reference=reference)
    with pysam.VariantFile(str(output), index_filename=None) as variant_file:
        header = str(variant_file.header)
    expected = ('##ChoCallate=<Version=test,Threshold=2,Callers="bcftools,gatk",'
                'Rule=two-stage,Tie=qual,AbsentVote=abstain,MinCoverage=3,'
                'RemoveInvariant=False,SplitMultiallelic=False,'
                'CallersFailed="freebayes">')
    assert expected in header
