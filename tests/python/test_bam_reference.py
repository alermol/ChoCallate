"""Tests for the BAM header vs. reference .fai guard."""

import subprocess
import sys

import pysam
import pytest

# Tests run from this directory, so conftest is importable by name.
from conftest import BIN_DIR  # noqa: I100,I202

CHECK = BIN_DIR / 'check_bam_reference.py'


def write_bam(path, contigs, m5=None):
    """Write a header-only BAM with the given (name, length) @SQ entries."""
    entries = []
    for name, length in contigs:
        entry = {'SN': name, 'LN': length}
        if m5 and name in m5:
            entry['M5'] = m5[name]
        entries.append(entry)
    header = {'HD': {'VN': '1.6'}, 'SQ': entries}
    with pysam.AlignmentFile(str(path), 'wb', header=header):
        pass
    return path


def write_fai(path, contigs):
    """Write a minimal .fai with the given (name, length) pairs."""
    lines = [
        '{}\t{}\t6\t60\t61'.format(name, length)
        for name, length in contigs
    ]
    path.write_text('\n'.join(lines) + '\n')
    return path


def write_dict(path, contigs, m5=None):
    """Write a minimal SAM sequence dictionary."""
    lines = ['@HD\tVN:1.6']
    for name, length in contigs:
        line = '@SQ\tSN:{}\tLN:{}'.format(name, length)
        if m5 and name in m5:
            line += '\tM5:{}'.format(m5[name])
        lines.append(line)
    path.write_text('\n'.join(lines) + '\n')
    return path


def run_check(bam, fai, dict_path=None):
    command = [sys.executable, str(CHECK), str(bam), str(fai)]
    if dict_path is not None:
        command += ['--dict', str(dict_path)]
    return subprocess.run(
        command,
        capture_output=True,
        text=True,
    )


def test_matching_names_and_lengths_are_accepted(tmp_path):
    contigs = [('chr1', 100), ('chr2', 250)]
    bam = write_bam(tmp_path / 'in.bam', contigs)
    fai = write_fai(tmp_path / 'ref.fa.fai', contigs)

    result = run_check(bam, fai)

    assert result.returncode == 0
    assert result.stdout == ''
    assert result.stderr == ''


def test_length_mismatch_fails_with_a_diff(tmp_path):
    bam = write_bam(tmp_path / 'in.bam', [('chr1', 100), ('chr2', 200)])
    fai = write_fai(tmp_path / 'ref.fa.fai', [('chr1', 100), ('chr2', 201)])

    result = run_check(bam, fai)

    assert result.returncode != 0
    assert 'do not match' in result.stderr
    assert '-chr2\t200' in result.stderr
    assert '+chr2\t201' in result.stderr


def test_missing_and_extra_contigs_are_reported(tmp_path):
    bam = write_bam(tmp_path / 'in.bam', [('chr1', 100), ('chr2', 200)])
    fai = write_fai(tmp_path / 'ref.fa.fai', [('chr1', 100), ('chr3', 300)])

    result = run_check(bam, fai)

    assert result.returncode != 0
    assert '-chr2\t200' in result.stderr
    assert '+chr3\t300' in result.stderr


def test_matching_m5_checksums_are_accepted(tmp_path):
    contigs = [('chr1', 100), ('chr2', 250)]
    m5 = {'chr1': 'abcd', 'chr2': 'ef01'}
    bam = write_bam(tmp_path / 'in.bam', contigs, m5=m5)
    fai = write_fai(tmp_path / 'ref.fa.fai', contigs)
    dictionary = write_dict(tmp_path / 'ref.dict', contigs, m5=m5)

    result = run_check(bam, fai, dictionary)

    assert result.returncode == 0
    assert result.stderr == ''


def test_m5_mismatch_fails_for_same_names_and_lengths(tmp_path):
    contigs = [('chr1', 100)]
    bam = write_bam(tmp_path / 'in.bam', contigs, m5={'chr1': 'aaaa'})
    fai = write_fai(tmp_path / 'ref.fa.fai', contigs)
    dictionary = write_dict(
        tmp_path / 'ref.dict', contigs, m5={'chr1': 'bbbb'})

    result = run_check(bam, fai, dictionary)

    assert result.returncode != 0
    assert 'M5 checksums do not match' in result.stderr
    assert 'chr1\tBAM aaaa != reference bbbb' in result.stderr


def test_missing_bam_m5_is_not_an_error(tmp_path):
    contigs = [('chr1', 100)]
    bam = write_bam(tmp_path / 'in.bam', contigs)
    fai = write_fai(tmp_path / 'ref.fa.fai', contigs)
    dictionary = write_dict(
        tmp_path / 'ref.dict', contigs, m5={'chr1': 'bbbb'})

    result = run_check(bam, fai, dictionary)

    assert result.returncode == 0
    assert result.stderr == ''


def test_prepare_bam_invokes_the_guard():
    module = (BIN_DIR.parent / 'modules' / 'local' / 'prepare_bam'
              / 'main.nf').read_text()
    assert ('check_bam_reference.py tmp/read1 '
            'tmp/ref_genome.fasta.fai' in module)
    assert '--dict tmp/ref_genome.dict' in module


@pytest.mark.parametrize('mapper', [
    'mapping_bwa', 'mapping_bowtie2', 'mapping_minimap2'])
def test_mapping_passes_the_reference_index_to_the_guard(mapper):
    module = (BIN_DIR.parent / 'modules' / 'local' / mapper
              / 'main.nf').read_text()
    assert '--fai tmp/ref_genome.fasta.fai' in module
