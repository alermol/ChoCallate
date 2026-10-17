"""Tests for the reference contig length guard."""

import subprocess
import sys

from conftest import BIN_DIR  # noqa: I100,I202  (tests run from this directory)

CHECK = BIN_DIR / 'check_contig_lengths.py'
MAX_LENGTH = (1 << 31) - 1


def write_fai(path, contigs):
    """Write a minimal .fai with the given (name, length) pairs."""
    lines = [
        '{}\t{}\t6\t60\t61'.format(name, length)
        for name, length in contigs
    ]
    path.write_text('\n'.join(lines) + '\n')
    return path


def test_contigs_within_the_limit_are_accepted(tmp_path):
    fai = write_fai(tmp_path / 'ref.fa.fai',
                    [('chr1', 100), ('chr2', MAX_LENGTH)])
    result = subprocess.run([sys.executable, str(CHECK), str(fai)],
                            capture_output=True, text=True)
    assert result.returncode == 0
    assert result.stdout == ''
    assert result.stderr == ''


def test_contig_over_the_limit_fails_with_a_clear_message(tmp_path):
    fai = write_fai(tmp_path / 'ref.fa.fai',
                    [('chr1', 100), ('chrBig', MAX_LENGTH + 1)])
    result = subprocess.run([sys.executable, str(CHECK), str(fai)],
                            capture_output=True, text=True)
    assert result.returncode != 0
    assert "'chrBig' ({} bp)".format(MAX_LENGTH + 1) in result.stderr
    assert str(MAX_LENGTH) in result.stderr
    assert '2^31-1' in result.stderr
    assert '32-bit' in result.stderr


def test_fai_creation_process_invokes_the_guard():
    module = (BIN_DIR.parent / 'modules' / 'local' / 'create_fai_index'
              / 'main.nf').read_text()
    assert 'check_contig_lengths.py ref_genome.fasta.fai' in module
