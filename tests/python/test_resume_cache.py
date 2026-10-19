"""Acceptance test: adding -resume must not re-run completed tasks."""

import json
import shutil
import subprocess
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[2]
NEXTFLOW = shutil.which('nextflow')
BCFTOOLS = shutil.which('bcftools')
BOWTIE2_BUILD = shutil.which('bowtie2-build')
TEST_RUN = PROJECT_ROOT / 'test_run'
REFERENCE = TEST_RUN / 'arth_chr1.fasta.gz'
READ1 = TEST_RUN / 'test_reads_R1.fq.gz'
READ2 = TEST_RUN / 'test_reads_R2.fq.gz'

pytestmark = pytest.mark.skipif(
    not all([NEXTFLOW, BCFTOOLS, BOWTIE2_BUILD, REFERENCE.is_file(),
             READ1.is_file(), READ2.is_file()]),
    reason='nextflow, bcftools, bowtie2 or the test dataset is unavailable',
)


def _reference_index(tmp_path):
    index_dir = tmp_path / 'reference'
    index_dir.mkdir()
    subprocess.run(
        [BOWTIE2_BUILD, '--quiet', '--threads', '2', str(REFERENCE),
         str(index_dir / REFERENCE.name)],
        check=True, capture_output=True, text=True,
    )
    return index_dir


def _params_file(tmp_path, output_type, index_dir):
    params = {
        'input': {
            'samples_tsv': str(tmp_path / 'samples.tsv'),
            'reference_genome': str(REFERENCE),
            'reference_index_dir': str(index_dir),
            'reference_index_prefix': REFERENCE.name,
            'format': 'fastq',
            'reads_type': 'pe',
        },
        'output': {
            'directory': str(tmp_path / 'results'),
            'type': output_type,
            'format': 'bcf',
            'remove_invariant': True,
        },
        'mapping': {'mapper': 'bowtie2'},
        'left_align_indels': {'enabled': False},
        'coverage': {'min_coverage': 1},
        'calling': {'callers': 'bcftools'},
        'consensus': {'threshold': 1},
    }
    path = tmp_path / 'params.json'
    path.write_text(json.dumps(params))
    return path


def _run(tmp_path, params_file, resume):
    # tests/nextflow.config disables nf-boost cleanup: the plugin deletes work
    # files once no downstream process needs them, so cleanup and -resume
    # cannot be combined. The pipeline exposes the observer as
    # boost.cleanup, but this test isolates cache-key stability from cleanup.
    command = [
        NEXTFLOW, '-log', str(tmp_path / 'nextflow.log'),
        'run', str(PROJECT_ROOT / 'main.nf'),
        '-c', str(PROJECT_ROOT / 'tests' / 'nextflow.config'),
        '-params-file', str(params_file),
        '-work-dir', str(tmp_path / 'work'),
        '-ansi-log', 'false',
    ]
    if resume:
        command.append('-resume')
    return subprocess.run(command, cwd=str(tmp_path), text=True,
                          capture_output=True, timeout=1200)


def _trace_rows(trace):
    lines = trace.read_text().splitlines()
    header = lines[0].split('\t')
    return [dict(zip(header, line.split('\t'))) for line in lines[1:]]


@pytest.mark.parametrize('output_type', ['single', 'sample'])
def test_resume_keeps_all_tasks_cached(tmp_path, output_type):
    samples = tmp_path / 'samples.tsv'
    samples.write_text(f'sample1\t{READ1}\t{READ2}\n')
    params_file = _params_file(tmp_path, output_type,
                               _reference_index(tmp_path))

    first = _run(tmp_path, params_file, resume=False)
    assert first.returncode == 0, first.stdout + first.stderr

    second = _run(tmp_path, params_file, resume=True)
    assert second.returncode == 0, second.stdout + second.stderr

    rows = _trace_rows(tmp_path / 'ChoCallate_output' / 'trace.txt')
    assert rows
    statuses = {row['name']: row['status'] for row in rows}
    assert all(row['status'] == 'CACHED' for row in rows), statuses
    final_step = (
        'MERGE_OUTPUTS' if output_type == 'single' else 'EXPORT_OUTPUT'
    )
    assert any(row['name'].startswith(final_step) for row in rows)

    if output_type == 'single':
        output = tmp_path / 'results' / 'consensus.bcf'
    else:
        output = tmp_path / 'results' / 'per_sample' / 'sample1.bcf'
    header = subprocess.run(
        [BCFTOOLS, 'view', '-h', str(output)],
        check=True, capture_output=True, text=True,
    ).stdout
    assert '##ChoCallate=<' in header
