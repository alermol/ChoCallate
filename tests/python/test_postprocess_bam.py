"""Tests for the shared BAM post-processing helper.

The helper normalises aligner output before calling: it sorts, optionally
removes duplicates, filters MAPQ after duplicate marking, optionally left-aligns
indels, and writes the BAM index exactly once.
"""

import shutil
import subprocess

import pysam
import pytest

from conftest import BIN_DIR  # noqa: I100,I202  (tests run from this directory)

SAMTOOLS = shutil.which('samtools')
pytestmark = pytest.mark.skipif(SAMTOOLS is None, reason='samtools is not available')


def write_paired_bam(path):
    """Two duplicate pairs, one low-MAPQ pair and one unique pair."""
    header = {
        'HD': {'VN': '1.6', 'SO': 'coordinate'},
        'SQ': [{'SN': 'chr1', 'LN': 1000}],
        'RG': [{'ID': 's1', 'SM': 's1'}],
    }
    reads = []

    def add(name, start, flag, mapq, mate_start, template_length):
        read = pysam.AlignedSegment()
        read.query_name = name
        read.query_sequence = 'ACGTACGT'
        read.flag = flag
        read.reference_id = 0
        read.reference_start = start
        read.mapping_quality = mapq
        read.cigar = ((0, 8),)
        read.query_qualities = pysam.qualitystring_to_array('IIIIIIII')
        read.next_reference_id = 0
        read.next_reference_start = mate_start
        read.template_length = template_length
        read.set_tag('RG', 's1')
        reads.append(read)

    # duplicate pair 1
    add('dup1', 10, 99, 60, 30, 28)
    add('dup1', 30, 147, 60, 10, -28)
    # duplicate pair 2
    add('dup2', 10, 99, 60, 30, 28)
    add('dup2', 30, 147, 60, 10, -28)
    # low MAPQ pair, kept by markdup and removed by the later MAPQ filter
    add('low', 100, 99, 5, 120, 28)
    add('low', 120, 147, 5, 100, -28)
    # unique pair
    add('uniq', 200, 99, 60, 220, 28)
    add('uniq', 220, 147, 60, 200, -28)

    reads.sort(key=lambda read: (read.reference_start, read.query_name, read.flag))
    with pysam.AlignmentFile(str(path), 'wb', header=header) as out:
        for read in reads:
            out.write(read)
    pysam.index(str(path))
    return path


def run_helper(tmp_path, *args):
    output = tmp_path / 'mapping.bam'
    subprocess.run(
        [str(BIN_DIR / 'postprocess_bam.sh'),
         '--input', str(tmp_path / 'in.bam'),
         '--output', str(output),
         '--mapq', '10',
         '--threads', '1',
         *args],
        check=True, cwd=str(tmp_path))
    return output


def read_names(path):
    with pysam.AlignmentFile(str(path), 'rb') as bam:
        return [read.query_name for read in bam]


def test_postprocess_removes_duplicates_filters_mapq_and_indexes(tmp_path):
    write_paired_bam(tmp_path / 'in.bam')
    output = run_helper(tmp_path, '--rmdup')
    assert sorted(read_names(output)) == ['dup1', 'dup1', 'uniq', 'uniq']
    assert (tmp_path / 'mapping.bam.csi').is_file()
    subprocess.run([SAMTOOLS, 'quickcheck', str(output)], check=True)
    with pysam.AlignmentFile(str(output), 'rb') as bam:
        assert bam.header.to_dict()['HD']['SO'] == 'coordinate'


def test_postprocess_without_rmdup_keeps_duplicates(tmp_path):
    write_paired_bam(tmp_path / 'in.bam')
    output = run_helper(tmp_path)
    assert sorted(read_names(output)) == [
        'dup1', 'dup1', 'dup2', 'dup2', 'uniq', 'uniq']
    assert (tmp_path / 'mapping.bam.csi').is_file()
    subprocess.run([SAMTOOLS, 'quickcheck', str(output)], check=True)
