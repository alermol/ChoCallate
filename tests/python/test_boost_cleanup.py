"""End-to-end test for the enabled nf-boost cleanup observer.

The nf-test suite disables ``boost.cleanup`` so that work directories survive
its assertions.  This test runs a tiny BAM-only pipeline with cleanup enabled
and checks that the published result survives while a consumed task output is
removed from the work directory.
"""

import json
import os
import shutil
import subprocess
from pathlib import Path

import pytest
from conftest import (  # noqa: I100,I202  (tests run from this directory)
    REF_SEQ,
    write_bam,
    write_reference,
)

PROJECT_ROOT = Path(__file__).resolve().parents[2]
NEXTFLOW = shutil.which('nextflow')

pytestmark = pytest.mark.skipif(
    NEXTFLOW is None,
    reason='nextflow is not installed',
)


def _environment():
    env = os.environ.copy()
    if 'NXF_OFFLINE' not in env:
        plugins = Path.home() / '.nextflow' / 'plugins'
        has_plugins = (
            plugins.is_dir()
            and any(plugins.glob('nf-schema-*'))
            and any(plugins.glob('nf-boost-*'))
        )
        if has_plugins:
            env['NXF_OFFLINE'] = 'true'
    return env


def _params_file(tmp_path, reference, samplesheet, output_dir):
    params = {
        'input': {
            'samples_tsv': str(samplesheet),
            'reference_genome': str(reference),
            'format': 'bam',
        },
        'output': {
            'directory': str(output_dir),
            'type': 'single',
            'format': 'bcf',
        },
        'mapping': {'cpu': 1},
        'bam_filter': {'min_map_qual': 0},
        'rmdup': {'mode': 'none'},
        'left_align_indels': {'enabled': False},
        'coverage': {'min_coverage': 1},
        'calling': {'callers': 'bcftools', 'cpu': 1},
        'consensus': {'threshold': 1, 'cpu': 1},
        'boost': {'cleanup': 'v2'},
        'resources': {
            'max_retries': 0,
            'low': {'cpus': 1, 'memory_mb': 1024},
            'medium': {'cpus': 1, 'memory_mb': 1024},
            'high': {'cpus': 1, 'memory_mb': 1024},
        },
    }
    path = tmp_path / 'params.json'
    path.write_text(json.dumps(params))
    return path


def _trace_rows(path):
    lines = path.read_text().splitlines()
    header = lines[0].split('\t')
    return [dict(zip(header, line.split('\t'))) for line in lines[1:]]


def test_cleanup_enabled_removes_consumed_task_output(tmp_path):
    """Cleanup removes the per-sample BCF but keeps the published output."""
    reference = write_reference(tmp_path / 'ref.fa')
    reads = [(f'read{i}', i, REF_SEQ[i:i + 12]) for i in range(10, 30)]
    bam = write_bam(tmp_path / 'sample.bam', reads)
    samplesheet = tmp_path / 'samples.tsv'
    samplesheet.write_text(f'sample1\t{bam}\n')
    output_dir = tmp_path / 'results'
    params = _params_file(tmp_path, reference, samplesheet, output_dir)
    work_dir = tmp_path / 'work'

    result = subprocess.run(
        [NEXTFLOW, '-log', str(tmp_path / 'nextflow.log'), 'run',
         str(PROJECT_ROOT / 'main.nf'), '-params-file', str(params),
         '-work-dir', str(work_dir), '-ansi-log', 'false'],
        cwd=str(tmp_path), env=_environment(), capture_output=True,
        text=True, timeout=900,
    )

    assert result.returncode == 0, result.stdout + result.stderr
    assert (output_dir / 'consensus.bcf').is_file()

    # The trace always lands under the default output directory: it is opened
    # while nextflow.config is parsed, before the params file is applied.
    trace = tmp_path / 'ChoCallate_output' / 'trace.txt'
    consensus_rows = [
        row for row in _trace_rows(trace)
        if row['name'].startswith('GENERATE_CONSENSUS')]
    assert len(consensus_rows) == 1
    assert consensus_rows[0]['status'] == 'COMPLETED'
    # The trace records a shortened hash; the work directory name is full.
    shard, prefix = consensus_rows[0]['hash'].split('/')
    task_dirs = list((work_dir / shard).glob(prefix + '*'))
    assert len(task_dirs) <= 1
    # The downstream merge consumed sample1.bcf, so the observer deletes it
    # (and may remove the now-empty task directory as well).
    assert all(not list(task_dir.glob('*.bcf')) for task_dir in task_dirs)
