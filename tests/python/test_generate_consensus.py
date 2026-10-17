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
    """Two matching variant calls win over one placeholder vote."""
    output = consensus_runner(
        bed_intervals=[('chr1', POS_30_0B, 31), ('chr1', 40, 42)],
        caller_records=[
            [(POS_30_0B, ref_base(30), [alt_base(30)], (0, 1))],
            [(POS_30_0B, ref_base(30), [alt_base(30)], (1, 1))],
        ],
        threshold=1,
    )
    records = read_records(output)
    # 0/1 and 1/1 are different keys (C-03, voting-model work); the tie goes to
    # the first caller and the threshold is met.  The point here is that the
    # placeholder votes do not win and the interval ends are written.
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
