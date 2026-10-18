"""Tests for the shared BAM post-processing helper.

The helper normalises aligner output before calling: it sorts, optionally
removes duplicates, filters MAPQ after duplicate marking, optionally left-aligns
indels, and writes the BAM index exactly once.
"""

import os
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


def run_helper(tmp_path, *args, env=None, threads=1, memory_mb=None):
    output = tmp_path / 'mapping.bam'
    command = [
        str(BIN_DIR / 'postprocess_bam.sh'),
        '--input', str(tmp_path / 'in.bam'),
        '--output', str(output),
        '--mapq', '10',
        '--threads', str(threads),
    ]
    if memory_mb is not None:
        command += ['--memory-mb', str(memory_mb)]
    command += list(args)
    subprocess.run(command, check=True, cwd=str(tmp_path), env=env)
    return output


def write_stub_samtools(tmp_path):
    """Replace samtools with a wrapper that logs and forwards arguments."""
    real_samtools = shutil.which('samtools')
    stub_dir = tmp_path / 'stub_samtools_bin'
    stub_dir.mkdir()
    log = tmp_path / 'samtools_args.txt'
    stub = stub_dir / 'samtools'
    stub.write_text(
        '#!/bin/bash\n'
        'printf "%s\\n" "$*" >> "$SAMTOOLS_ARGS_LOG"\n'
        'exec "$REAL_SAMTOOLS" "$@"\n'
    )
    stub.chmod(0o755)
    env = os.environ.copy()
    env['PATH'] = str(stub_dir) + os.pathsep + env['PATH']
    env['SAMTOOLS_ARGS_LOG'] = str(log)
    env['REAL_SAMTOOLS'] = real_samtools
    return env, log


def write_stub_gatk(tmp_path):
    """Replace GATK with a stub that logs its arguments and copies the BAM."""
    stub_dir = tmp_path / 'stub_bin'
    stub_dir.mkdir()
    log = tmp_path / 'gatk_args.txt'
    stub = stub_dir / 'gatk'
    stub.write_text(
        '#!/bin/bash\n'
        'printf "%s\\n" "$@" > "$GATK_ARGS_LOG"\n'
        'input=""\n'
        'output=""\n'
        'while [ "$#" -gt 0 ]; do\n'
        '  case "$1" in\n'
        '    -I) input="$2"; shift 2 ;;\n'
        '    -O) output="$2"; shift 2 ;;\n'
        '    *) shift ;;\n'
        '  esac\n'
        'done\n'
        'cp "$input" "$output"\n'
    )
    stub.chmod(0o755)
    return stub_dir, log


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


@pytest.mark.parametrize('extra_args', [[], ['--rmdup']],
                         ids=['sort', 'markdup'])
def test_sort_memory_comes_from_task_memory(tmp_path, extra_args):
    write_paired_bam(tmp_path / 'in.bam')
    env, log = write_stub_samtools(tmp_path)

    run_helper(
        tmp_path, *extra_args, threads=2, memory_mb=4096, env=env)

    sort_calls = [
        line for line in log.read_text().splitlines()
        if line.startswith('sort ')
    ]
    assert sort_calls
    assert all('-m 1024M' in line for line in sort_calls)


def test_left_align_receives_task_java_options(tmp_path):
    write_paired_bam(tmp_path / 'in.bam')
    stub_dir, log = write_stub_gatk(tmp_path)
    reference = tmp_path / 'ref.fa'
    reference.write_text('>chr1\nACGT\n')
    dictionary = tmp_path / 'ref.dict'
    dictionary.write_text('')

    env = os.environ.copy()
    env['PATH'] = str(stub_dir) + os.pathsep + env['PATH']
    env['GATK_ARGS_LOG'] = str(log)
    java_options = '-Djava.io.tmpdir=/tmp/scratch -Xmx1024m'

    run_helper(
        tmp_path,
        '--left-align',
        '--reference', str(reference),
        '--dict', str(dictionary),
        '--java-options', java_options,
        env=env,
    )

    arguments = log.read_text().splitlines()
    assert arguments[0] == 'LeftAlignIndels'
    option_index = arguments.index('--java-options')
    assert arguments[option_index + 1] == java_options
