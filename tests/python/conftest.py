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


def write_bam(path, reads, contig=CHR, length=len(REF_SEQ), csi=False):
    """Write a coordinate-sorted indexed BAM.

    Args:
        path: Output BAM path.
        reads: Iterable of tuples, in ascending start order:
            ``(query_name, 0-based start, sequence)`` optionally followed by a
            CIGAR tuple, a SAM flag and a base-quality string.  ``None``
            qualities are stored as the BAM missing-quality sentinel (255).
        contig: Contig name of the single reference sequence.
        length: Contig length.
        csi: Write a CSI index instead of a BAI index (needed for positions
            at or above 2**29).

    Returns:
        The BAM path.
    """
    header = {'HD': {'VN': '1.6', 'SO': 'coordinate'},
              'SQ': [{'SN': contig, 'LN': length}]}
    with pysam.AlignmentFile(str(path), 'wb', header=header) as bam:
        for entry in reads:
            name, start, seq = entry[:3]
            cigar = entry[3] if len(entry) > 3 else ((0, len(seq)),)
            flag = entry[4] if len(entry) > 4 else 0
            qualities = entry[5] if len(entry) > 5 else 'I' * len(seq)
            read = pysam.AlignedSegment()
            read.query_name = name
            read.query_sequence = seq
            read.flag = flag
            read.reference_id = 0
            read.reference_start = start
            read.mapping_quality = 60
            read.cigar = cigar
            if qualities is not None:
                read.query_qualities = pysam.qualitystring_to_array(
                    qualities)
            bam.write(read)
    if csi:
        pysam.index('-c', str(path))
    else:
        pysam.index(str(path))
    return path


def write_sparse_reference(path, contig, length, start0, seq,
                           line_bases=60):
    """Write a mostly-sparse FASTA with known bases and a hand-made .fai.

    Only ``seq``, placed at 0-based ``start0``, is written to the file; the
    gap before it is a filesystem hole, so a contig close to the 32-bit
    coordinate limit costs almost no disk space.  The ``.fai`` is written
    directly because indexing would read the whole sparse file.

    Args:
        path: Output FASTA path.
        contig: Contig name.
        length: Declared contig length.
        start0: 0-based position of the first base of ``seq``.
        seq: Known bases to place at ``start0``.
        line_bases: FASTA line width used by the .fai arithmetic.

    Returns:
        The FASTA path.
    """
    line_width = line_bases + 1
    header = f'>{contig}\n'
    offset = len(header)
    with open(path, 'w', encoding='utf-8') as handle:
        handle.write(header)
        for index, base in enumerate(seq, start=start0):
            line, column = divmod(index, line_bases)
            handle.seek(offset + line * line_width + column)
            handle.write(base)
    with open(f'{path}.fai', 'w', encoding='utf-8') as handle:
        handle.write(f'{contig}\t{length}\t{offset}\t{line_bases}\t'
                     f'{line_width}\n')
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
            such as (0, 1), (0, 0) or (None, None).  An optional fifth item sets
            the record QUAL.  An optional sixth item sets the FILTER: a name
            such as ``'PASS'``/``'LowQual'``, a sequence of names, or
            ``None``/``'.'`` for no filter.  An empty ``alts`` is written as the
            VCF placeholder ALT ``.``, the same representation
            generate_consensus.py itself emits for hom-reference records.
        contig: Contig name.
        length: Contig length.

    Returns:
        The BCF path.
    """
    header = pysam.VariantHeader()
    header.add_line(f'##contig=<ID={contig},length={length}>')
    header.add_line('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    header.add_line('##FILTER=<ID=LowQual,Description="Low quality">')
    header.add_sample(sample_name)
    with pysam.VariantFile(str(path), 'wb', header=header) as out:
        for entry in records:
            start0, ref, alts, gt = entry[:4]
            qual = entry[4] if len(entry) > 4 else None
            filters = entry[5] if len(entry) > 5 else None
            alleles = (ref,) + (tuple(alts) if alts else ('.',))
            record = out.new_record(contig=contig, start=start0,
                                    stop=start0 + len(ref),
                                    alleles=alleles)
            if qual is not None:
                record.qual = qual
            if isinstance(filters, str):
                filters = (filters,)
            for name in filters or ():
                if name != '.':
                    record.filter.add(name)
            record.samples[0]['GT'] = gt
            out.write(record)
    subprocess.run([BCFTOOLS, 'index', '--csi', '-f', str(path)], check=True)
    return path


def run_consensus(tmp_path, bam, bed_gz, caller_bcfs, threshold, sample_name='sample1',
                  use_accel=False, output_name='out.bcf', reference=None,
                  callers=None, configured_callers=None, check=True,
                  absent_vote='ref', tie_policy='nocall', min_coverage=0,
                  filtered_vote='abstain', compress_reference_blocks=False,
                  variants_only=False):
    """Run bin/generate_consensus.py as a subprocess and return its output path.

    Args:
        caller_bcfs: Caller BCF paths, in the order they should be passed.
        callers: Caller names matching ``caller_bcfs`` (defaults to
            ``caller0``, ``caller1``, ...).
        configured_callers: Caller list for ``--callers`` (defaults to
            ``callers``); extra names count as failed callers.
        check: Raise on a non-zero exit status when true.
        absent_vote: ``ref`` or ``abstain`` voting for absent callers.
        tie_policy: ``nocall`` or ``first`` tie-breaking.
        min_coverage: Coverage required by ``absent_vote=abstain``.
        filtered_vote: ``abstain`` or ``ref`` voting for filtered records.
        compress_reference_blocks: Compress hom-reference runs into blocks.
        variants_only: Write only consensus variant records, omitting reference
            records.
    """
    output = tmp_path / output_name
    if reference is None:
        reference = tmp_path / 'ref.fa'
    if callers is None:
        callers = [f'caller{index}' for index in range(len(caller_bcfs))]
    input_args = []
    for caller, path in zip(callers, caller_bcfs):
        input_args += ['--input', f'{caller}={path}']
    if configured_callers is None:
        configured_callers = callers
    env = os.environ.copy()
    if not use_accel:
        env['CHOCALLATE_NO_ACCEL'] = '1'
    cmd = [
        sys.executable, str(CONSENSUS_SCRIPT),
        *input_args,
        '--callers', ','.join(configured_callers),
        '--bed', str(bed_gz),
        '--bam', str(bam),
        '--output', str(output),
        '--sample_name', sample_name,
        '--reference', str(reference),
        '--consensus_threshold', str(threshold),
        '--absent_vote', absent_vote,
        '--filtered_vote', filtered_vote,
        '--tie_policy', tie_policy,
        '--min_coverage', str(min_coverage),
        '--version', 'test',
    ]
    if compress_reference_blocks:
        cmd.append('--compress_reference_blocks')
    if variants_only:
        cmd.append('--variants_only')
    subprocess.run(cmd, check=check, env=env, cwd=str(tmp_path),
                   capture_output=True, text=True)
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
    def runner(bed_intervals, caller_records, threshold, sample_name='sample1',
               **kwargs):
        bed_gz = write_bed(tmp_path / 'coverage.bed', bed_intervals)
        caller_bcfs = []
        for index, records in enumerate(caller_records):
            caller_bcfs.append(write_caller_bcf(
                tmp_path / f'caller{index}.bcf', 'sample1', records))
        output = run_consensus(tmp_path, basic_dataset['bam'], bed_gz,
                               caller_bcfs, threshold,
                               sample_name=sample_name,
                               use_accel=use_accel,
                               reference=basic_dataset['reference'],
                               **kwargs)
        return output

    return runner
