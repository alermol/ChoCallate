"""Tests for the shared BAM post-processing helper.

The helper normalises aligner output before calling: it sorts, optionally
marks or removes duplicates, filters MAPQ after duplicate marking, optionally
left-aligns indels, and writes the BAM index exactly once.
"""

import os
import shutil
import subprocess

import pysam
import pytest
from conftest import BIN_DIR  # noqa: I100,I202  (tests run from this directory)

SAMTOOLS = shutil.which('samtools')
pytestmark = pytest.mark.skipif(SAMTOOLS is None, reason='samtools is not available')


def write_paired_bam(path, m5=None):
    """Two duplicate pairs, one low-MAPQ pair and one unique pair."""
    sq = {'SN': 'chr1', 'LN': 1000}
    if m5:
        sq['M5'] = m5
    header = {
        'HD': {'VN': '1.6', 'SO': 'coordinate'},
        'SQ': [sq],
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


def run_helper(tmp_path, *args, env=None, threads=1, memory_mb=None,
               reserve_mb=None):
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
    if reserve_mb is not None:
        command += ['--reserve-mb', str(reserve_mb)]
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
    # --rmdup is the backwards-compatible alias for --rmdup-mode remove.
    output = run_helper(tmp_path, '--rmdup')
    assert sorted(read_names(output)) == ['dup1', 'dup1', 'uniq', 'uniq']
    assert (tmp_path / 'mapping.bam.csi').is_file()
    subprocess.run([SAMTOOLS, 'quickcheck', str(output)], check=True)
    with pysam.AlignmentFile(str(output), 'rb') as bam:
        assert bam.header.to_dict()['HD']['SO'] == 'coordinate'


def test_postprocess_marks_duplicates_without_removing_them(tmp_path):
    write_paired_bam(tmp_path / 'in.bam')
    output = run_helper(tmp_path, '--rmdup-mode', 'mark')
    assert sorted(read_names(output)) == [
        'dup1', 'dup1', 'dup2', 'dup2', 'uniq', 'uniq']
    assert (tmp_path / 'mapping.bam.csi').is_file()
    subprocess.run([SAMTOOLS, 'quickcheck', str(output)], check=True)
    with pysam.AlignmentFile(str(output), 'rb') as bam:
        flags = {}
        for read in bam:
            flags.setdefault(read.query_name, set()).add(read.is_duplicate)
    assert flags['uniq'] == {False}
    # Exactly one of the two identical pairs is flagged, both mates included.
    assert len([name for name in ('dup1', 'dup2')
                if flags[name] == {True}]) == 1


@pytest.mark.parametrize('extra_args', [[], ['--rmdup-mode', 'none']],
                         ids=['default', 'none'])
def test_postprocess_without_rmdup_keeps_duplicates(tmp_path, extra_args):
    write_paired_bam(tmp_path / 'in.bam')
    output = run_helper(tmp_path, *extra_args)
    assert sorted(read_names(output)) == [
        'dup1', 'dup1', 'dup2', 'dup2', 'uniq', 'uniq']
    assert (tmp_path / 'mapping.bam.csi').is_file()
    subprocess.run([SAMTOOLS, 'quickcheck', str(output)], check=True)
    with pysam.AlignmentFile(str(output), 'rb') as bam:
        assert not any(read.is_duplicate for read in bam)


def write_fai(path, length=1000):
    path.write_text('chr1\t{}\t6\t60\t61\n'.format(length))
    return path


def test_postprocess_accepts_a_matching_reference_index(tmp_path):
    write_paired_bam(tmp_path / 'in.bam')
    fai = write_fai(tmp_path / 'ref.fa.fai')

    output = run_helper(tmp_path, '--fai', str(fai))

    assert output.is_file()


def test_postprocess_rejects_a_bam_from_another_reference(tmp_path):
    write_paired_bam(tmp_path / 'in.bam')
    fai = write_fai(tmp_path / 'ref.fa.fai', length=999)

    with pytest.raises(subprocess.CalledProcessError):
        run_helper(tmp_path, '--fai', str(fai))


def test_postprocess_rejects_mismatched_m5_checksums(tmp_path):
    write_paired_bam(tmp_path / 'in.bam', m5='aaaa')
    fai = write_fai(tmp_path / 'ref.fa.fai')
    dictionary = tmp_path / 'ref.dict'
    dictionary.write_text(
        '@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\tM5:bbbb\n')

    with pytest.raises(subprocess.CalledProcessError):
        run_helper(tmp_path, '--fai', str(fai),
                   '--dict', str(dictionary))


def test_postprocess_rejects_unknown_rmdup_mode(tmp_path):
    write_paired_bam(tmp_path / 'in.bam')
    with pytest.raises(subprocess.CalledProcessError):
        run_helper(tmp_path, '--rmdup-mode', 'dedup')


@pytest.mark.parametrize(
    'extra_args,expected',
    [([], '292M'),
     (['--rmdup-mode', 'mark'], '512M'),
     (['--rmdup-mode', 'remove'], '512M')],
    ids=['sort', 'mark', 'remove'])
def test_sort_memory_comes_from_task_memory(tmp_path, extra_args, expected):
    write_paired_bam(tmp_path / 'in.bam')
    env, log = write_stub_samtools(tmp_path)

    run_helper(
        tmp_path, *extra_args, threads=8, memory_mb=4096, env=env)

    sort_calls = [
        line for line in log.read_text().splitlines()
        if line.startswith('sort ')
    ]
    assert sort_calls
    assert all(f'-m {expected}' in line for line in sort_calls)


@pytest.mark.parametrize(
    'extra_args,expected',
    [([], '146M'),
     (['--rmdup-mode', 'mark'], '256M'),
     (['--rmdup-mode', 'remove'], '256M')],
    ids=['sort', 'mark', 'remove'])
def test_sort_memory_excludes_the_aligner_reserve(tmp_path, extra_args,
                                                  expected):
    write_paired_bam(tmp_path / 'in.bam')
    env, log = write_stub_samtools(tmp_path)

    run_helper(
        tmp_path, *extra_args, threads=8, memory_mb=4096,
        reserve_mb=2048, env=env)

    sort_calls = [
        line for line in log.read_text().splitlines()
        if line.startswith('sort ')
    ]
    assert sort_calls
    assert all(f'-m {expected}' in line for line in sort_calls)


def test_postprocess_rejects_reserve_not_below_task_memory(tmp_path):
    write_paired_bam(tmp_path / 'in.bam')
    with pytest.raises(subprocess.CalledProcessError):
        run_helper(tmp_path, threads=1, memory_mb=1024, reserve_mb=1024)


@pytest.mark.parametrize(
    'extra_args,stages',
    [([], 2), (['--rmdup-mode', 'mark'], 5),
     (['--rmdup-mode', 'remove'], 5)],
    ids=['sort', 'mark', 'remove'])
def test_threads_are_split_across_concurrent_stages(
        tmp_path, extra_args, stages):
    write_paired_bam(tmp_path / 'in.bam')
    env, log = write_stub_samtools(tmp_path)

    run_helper(
        tmp_path, *extra_args, threads=8, memory_mb=8192, env=env)

    subcommands = {'collate', 'fixmate', 'sort', 'markdup', 'view'}
    calls = [
        line.split() for line in log.read_text().splitlines()
        if line.split() and line.split()[0] in subcommands
    ]
    assert len(calls) == stages
    allocations = [int(call[call.index('--threads') + 1]) for call in calls]
    assert sum(allocations) == 8
    # Sort is the only stage that receives more than one thread.
    sort_call = next(call for call in calls if call[0] == 'sort')
    assert int(sort_call[sort_call.index('--threads') + 1]) == 8 - stages + 1
    assert sorted(allocations) == [1] * (stages - 1) + [8 - stages + 1]


def test_markdup_is_streamed_into_mapq_filter(tmp_path):
    write_paired_bam(tmp_path / 'in.bam')
    env, log = write_stub_samtools(tmp_path)

    run_helper(tmp_path, '--rmdup-mode', 'mark', env=env)

    text = log.read_text()
    assert 'dedup.bam' not in text
    calls = [line.split() for line in text.splitlines()]
    markdup = next(call for call in calls if call[0] == 'markdup')
    view = next(call for call in calls if call[0] == 'view')
    # markdup reads and writes the stream; view filters that stream.
    assert markdup[-2:] == ['-', '-']
    assert view[-1] == '-'


def test_filter_runs_before_sort_without_rmdup(tmp_path):
    write_paired_bam(tmp_path / 'in.bam')
    env, log = write_stub_samtools(tmp_path)

    run_helper(tmp_path, env=env)

    text = log.read_text()
    assert 'sorted.bam' not in text
    calls = [line.split() for line in text.splitlines()]
    view = next(call for call in calls if call[0] == 'view')
    sort = next(call for call in calls if call[0] == 'sort')
    # The filter emits uncompressed BAM and sort consumes it from stdin.
    assert '-u' in view
    assert sort[-1] == '-'
    assert '-o' in sort


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
