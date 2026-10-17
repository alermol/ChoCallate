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
    """The Cython and pure-pysam implementations produce identical text."""
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
                'Rule=two-stage,Tie=nocall,AbsentVote=ref,RemoveInvariant=False,'
                'SplitMultiallelic=False,CallersFailed="freebayes">')
    assert expected in header
    assert [(record.pos, (record.ref,) + tuple(record.alts),
             tuple(record.samples[0]['GT']))
            for record in records] == [
        (30, (ref_base(30), alt_base(30)), (0, 1))]


def test_command_and_caller_version_header_lines(tmp_path, use_accel):
    """The command line and caller versions are recorded in the header."""
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
        callers=['bcftools', 'gatk'], use_accel=use_accel,
        reference=reference,
        command_line='nextflow run main.nf -params-file run.yaml',
        caller_versions='bcftools=1.24,gatk=4.6.2.0')
    with pysam.VariantFile(str(output), index_filename=None) as variant_file:
        header = str(variant_file.header)
    assert ('##ChoCallateCommand="nextflow run main.nf '
            '-params-file run.yaml"' in header)
    assert ('##ChoCallateCallerVersions="bcftools=1.24,gatk=4.6.2.0"'
            in header)


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


def test_tie_policy_qual_uses_summed_quality(consensus_runner):
    """``qual`` picks the genotype class with the larger summed QUAL."""
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


def test_header_describes_abstain_and_qual_policies(tmp_path, use_accel):
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
        absent_vote='abstain', tie_policy='qual', min_coverage=3,
        use_accel=use_accel, reference=reference)
    with pysam.VariantFile(str(output), index_filename=None) as variant_file:
        header = str(variant_file.header)
    expected = ('##ChoCallate=<Version=test,Threshold=2,Callers="bcftools,gatk",'
                'Rule=two-stage,Tie=qual,AbsentVote=abstain,MinCoverage=3,'
                'RemoveInvariant=False,SplitMultiallelic=False,'
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
