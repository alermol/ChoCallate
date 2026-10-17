"""Fixtures for the Python unit tests of the ChoCallate helper scripts.

Tiny artificial datasets (reference, indexed BAM, bgzipped+tabixed coverage BED,
indexed caller BCFs) are built in code with pysam so that the tests run without
external test data.  The ChoCallate conda environment provides pysam and
bcftools; run the suite with:

    pytest -q tests/python

Both the pure-pysam and the Cython-accelerated consensus paths are exercised
(``use_accel`` parameter); the accelerated case is skipped when the optional
extension has not been built (``python bin/build_consensus_accel.py``).
"""

import os
import shutil
import subprocess
import sys
from pathlib import Path

import pysam
import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
BIN_DIR = REPO_ROOT / 'bin'
CONSENSUS_SCRIPT = BIN_DIR / 'generate_consensus.py'

CHR = 'chr1'
# 100 bp of "ACGT" repeats: position i (0-based) holds 'ACGT'[i % 4]
REF_SEQ = 'ACGT' * 25

# Prefer the bcftools shipped next to the running interpreter (the ChoCallate
# environment): other installations on PATH may lack the plugins the pipeline
# uses (bcftools +fill-tags).
_ENV_BCFTOOLS = Path(sys.executable).with_name('bcftools')
BCFTOOLS = str(_ENV_BCFTOOLS) if _ENV_BCFTOOLS.exists() else (shutil.which('bcftools') or 'bcftools')


def accel_available():
    """True when bin/consensus_accel can be imported by this interpreter."""
    if str(BIN_DIR) not in sys.path:
        sys.path.insert(0, str(BIN_DIR))
    try:
        import consensus_accel  # noqa: F401
    except ImportError:
        return False
    return True


ACCEL_AVAILABLE = accel_available()


def ref_base(pos1):
    """Reference base at a 1-based position of the test contig."""
    return REF_SEQ[(pos1 - 1) % len(REF_SEQ)]


def alt_base(pos1):
    """A different base than the reference at a 1-based position."""
    return {'A': 'G', 'C': 'T', 'G': 'C', 'T': 'A'}[ref_base(pos1)]


def write_reference(path, seq=REF_SEQ, contig=CHR):
    """Write a FASTA file plus .fai index and return its path."""
    with open(path, 'w', encoding='utf-8') as handle:
        handle.write(f'>{contig}\n')
        for start in range(0, len(seq), 60):
            handle.write(seq[start:start + 60] + '\n')
    pysam.faidx(str(path))
    return path


def write_bam(path, reads, contig=CHR, length=len(REF_SEQ)):
    """Write a coordinate-sorted indexed BAM.

    Args:
        path: Output BAM path.
        reads: Iterable of (query_name, 0-based start, sequence) tuples, in
            ascending start order.
        contig: Contig name of the single reference sequence.
        length: Contig length.

    Returns:
        The BAM path.
    """
    header = {'HD': {'VN': '1.6', 'SO': 'coordinate'},
              'SQ': [{'SN': contig, 'LN': length}]}
    with pysam.AlignmentFile(str(path), 'wb', header=header) as bam:
        for name, start, seq in reads:
            read = pysam.AlignedSegment()
            read.query_name = name
            read.query_sequence = seq
            read.flag = 0
            read.reference_id = 0
            read.reference_start = start
            read.mapping_quality = 60
            read.cigar = ((0, len(seq)),)
            read.query_qualities = pysam.qualitystring_to_array('I' * len(seq))
            bam.write(read)
    pysam.index(str(path))
    return path


def write_bed(path, intervals, preset='bed', csi=True):
    """Write, bgzip and index (0-based half-open) intervals.

    Args:
        path: Path of the plain text BED to create; the compressed file gets
            the '.gz' suffix.
        intervals: Iterable of (contig, start, end) tuples.
        preset: tabix preset.
        csi: Write a '.csi' index (instead of '.tbi') when true.

    Returns:
        Path of the bgzipped BED (the file the tools actually read).
    """
    plain = Path(path)
    with open(plain, 'w', encoding='utf-8') as handle:
        for contig, start, end in intervals:
            handle.write(f'{contig}\t{start}\t{end}\n')
    gz = f'{plain}.gz'
    pysam.tabix_compress(str(plain), gz, force=True)
    pysam.tabix_index(gz, preset=preset, csi=csi, force=True)
    return Path(gz)


def write_caller_bcf(path, sample_name, records, contig=CHR, length=len(REF_SEQ)):
    """Write an indexed single-sample caller BCF.

    Args:
        path: Output BCF path.
        sample_name: Name of the single sample.
        records: Iterable of (start0, ref, alts, gt) tuples; ``gt`` is a tuple
            such as (0, 1), (0, 0) or (None, None).  An empty ``alts`` is
            written as the VCF placeholder ALT ``.``, the same representation
            generate_consensus.py itself emits for hom-reference records.
        contig: Contig name.
        length: Contig length.

    Returns:
        The BCF path.
    """
    header = pysam.VariantHeader()
    header.add_line(f'##contig=<ID={contig},length={length}>')
    header.add_line('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    header.add_sample(sample_name)
    with pysam.VariantFile(str(path), 'wb', header=header) as out:
        for start0, ref, alts, gt in records:
            alleles = (ref,) + (tuple(alts) if alts else ('.',))
            record = out.new_record(contig=contig, start=start0,
                                    stop=start0 + len(ref),
                                    alleles=alleles)
            record.samples[0]['GT'] = gt
            out.write(record)
    subprocess.run([BCFTOOLS, 'index', '--csi', str(path)], check=True)
    return path


def run_consensus(tmp_path, bam, bed_gz, caller_bcfs, threshold, sample_name='sample1',
                  use_accel=False, output_name='out.bcf', reference=None):
    """Run bin/generate_consensus.py as a subprocess and return its output path."""
    output = tmp_path / output_name
    if reference is None:
        reference = tmp_path / 'ref.fa'
    env = os.environ.copy()
    if not use_accel:
        env['CHOCALLATE_NO_ACCEL'] = '1'
    cmd = [
        sys.executable, str(CONSENSUS_SCRIPT),
        '--input', *[str(path) for path in caller_bcfs],
        '--bed', str(bed_gz),
        '--bam', str(bam),
        '--output', str(output),
        '--sample_name', sample_name,
        '--reference', str(reference),
        '--consensus_threshold', str(threshold),
        '--version', 'test',
    ]
    subprocess.run(cmd, check=True, env=env, cwd=str(tmp_path))
    return output


def read_positions(bcf_path):
    """1-based positions (in file order) of a BCF written by the tests."""
    with pysam.VariantFile(str(bcf_path)) as variant_file:
        return [record.pos for record in variant_file]


@pytest.fixture(params=[False, True], ids=['pure', 'accel'])
def use_accel(request):
    """Consensus implementation under test (pure-pysam or Cython)."""
    if request.param and not ACCEL_AVAILABLE:
        pytest.skip('consensus_accel is not built for this interpreter')
    return request.param


@pytest.fixture
def basic_dataset(tmp_path, use_accel):
    """One reference/BAM and the paths needed to run generate_consensus.py."""
    reference = write_reference(tmp_path / 'ref.fa')
    # A single read per position 20..24 (1-based), i.e. 0-based starts 19..23.
    reads = [(f'read{i}', i, REF_SEQ[i:i + 6]) for i in range(19, 24)]
    bam = write_bam(tmp_path / 'sample.bam', reads)
    return {
        'reference': reference,
        'bam': bam,
        'tmp_path': tmp_path,
    }


@pytest.fixture
def consensus_runner(tmp_path, use_accel, basic_dataset):
    """Callable that writes caller BCFs and runs the consensus script."""
    def runner(bed_intervals, caller_records, threshold, sample_name='sample1'):
        bed_gz = write_bed(tmp_path / 'coverage.bed', bed_intervals)
        caller_bcfs = []
        for index, records in enumerate(caller_records):
            caller_bcfs.append(write_caller_bcf(
                tmp_path / f'caller{index}.bcf', 'sample1', records))
        output = run_consensus(tmp_path, basic_dataset['bam'], bed_gz,
                               caller_bcfs, threshold,
                               sample_name=sample_name,
                               use_accel=use_accel,
                               reference=basic_dataset['reference'])
        return output

    return runner
