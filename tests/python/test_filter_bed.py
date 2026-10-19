"""Tests for the BED header stripping and contig validation guard."""

import subprocess
import sys

# Tests run from this directory, so conftest is importable by name.
from conftest import BIN_DIR  # noqa: I100,I202

FILTER = BIN_DIR / 'filter_bed.py'


def run_filter(tmp_path, bed_text, contigs=('chr1',)):
    """Write a .fai and BED, then run filter_bed.py on them."""
    fai = tmp_path / 'ref.fa.fai'
    fai.write_text(''.join(
        '{}\t100\t6\t60\t61\n'.format(name) for name in contigs))
    bed = tmp_path / 'in.bed'
    bed.write_text(bed_text)
    return subprocess.run(
        [sys.executable, str(FILTER), str(fai), str(bed)],
        capture_output=True,
        text=True,
    )


def test_headers_and_comments_are_removed(tmp_path):
    bed = (
        'track name=test\n'
        'browser position chr1:1-100\n'
        '# leading comment\n'
        '\n'
        'chr1\t0\t50\n'
        '# trailing comment\n'
        'chr1\t60\t80\n'
    )

    result = run_filter(tmp_path, bed)

    assert result.returncode == 0
    assert result.stderr == ''
    assert result.stdout == 'chr1\t0\t50\nchr1\t60\t80\n'


def test_header_words_are_only_removed_as_standalone_lines(tmp_path):
    bed = 'track1\t0\t10\nbrowser2\t20\t30\n'

    result = run_filter(tmp_path, bed, contigs=('track1', 'browser2'))

    assert result.returncode == 0
    assert result.stdout == bed


def test_unknown_contig_fails_with_a_clear_message(tmp_path):
    bed = 'chr1\t0\t50\nchr2\t0\t10\nchr3\t0\t10\n'

    result = run_filter(tmp_path, bed)

    assert result.returncode != 0
    assert result.stdout == 'chr1\t0\t50\n'
    assert "'chr2'" in result.stderr
    assert "'chr3'" in result.stderr
    assert 'reference .fai' in result.stderr


def test_crlf_line_endings_are_normalised(tmp_path):
    bed = 'track name=test\r\nchr1\t0\t50\r\n'

    result = run_filter(tmp_path, bed)

    assert result.returncode == 0
    assert result.stdout == 'chr1\t0\t50\n'
