"""Tests for task labels, memory escalation and transient-failure retries.

A process that exits with status 137 (SIGKILL, typically the scheduler
out-of-memory killer) must be retried with more memory.  The fixture records
the task memory of the successful attempt, so the test also proves that the
retry escalated the request.
"""

import json
import os
import shutil
import subprocess
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[2]
BASE_CONFIG = PROJECT_ROOT / 'conf' / 'base.config'
FIXTURE = Path(__file__).parent / 'fixtures' / 'retry_flaky' / 'main.nf'
PRECEDENCE_FIXTURE = (
    Path(__file__).parent / 'fixtures' / 'resource_precedence' / 'main.nf'
)
CALLER_FIXTURE = (
    Path(__file__).parent / 'fixtures' / 'caller_strategy' / 'main.nf'
)
NEXTFLOW = shutil.which('nextflow')

pytestmark = pytest.mark.skipif(
    NEXTFLOW is None,
    reason='nextflow is not installed',
)


def _merge(dst, src):
    for key, value in src.items():
        if isinstance(value, dict) and isinstance(dst.get(key), dict):
            _merge(dst[key], value)
        else:
            dst[key] = value
    return dst


def _params(tmp_path, overrides=None):
    params = {
        'outdir': str(tmp_path / 'results'),
        'resources': {
            'max_retries': 3,
            'low': {'cpus': 1, 'memory_mb': 256, 'time': '5m'},
            'medium': {'cpus': 1, 'memory_mb': 256, 'time': '5m'},
            'high': {'cpus': 1, 'memory_mb': 256, 'time': '5m'},
        },
    }
    _merge(params, overrides or {})
    path = tmp_path / 'params.json'
    path.write_text(json.dumps(params))
    return path


def _run_fixture(fixture, params, tmp_path):
    return subprocess.run(
        [
            NEXTFLOW, 'run', str(fixture),
            '-c', str(BASE_CONFIG),
            '-params-file', str(params),
            '-ansi-log', 'false',
        ],
        cwd=tmp_path,
        env=os.environ.copy(),
        capture_output=True,
        text=True,
        timeout=180,
    )


def test_exit_137_is_retried_with_escalated_memory(tmp_path):
    result = _run_fixture(FIXTURE, _params(tmp_path), tmp_path)

    output = result.stdout + result.stderr
    assert result.returncode == 0, output
    assert 'error exit status (137)' in output
    assert 'retried' in output
    # The label requests 256 MiB; attempt 2 multiplies it by two.
    assert (tmp_path / 'results' / 'memory.txt').read_text().strip() == '512'


def test_stage_cpu_overrides_label_and_task_defaults(tmp_path):
    params = _params(
        tmp_path,
        {
            'mapping': {'cpu': 3},
            'resources': {'medium': {'cpus': 6}},
        },
    )
    result = _run_fixture(PRECEDENCE_FIXTURE, params, tmp_path)

    output = result.stdout + result.stderr
    assert result.returncode == 0, output
    resources = (tmp_path / 'results' / 'resources.txt').read_text()
    assert 'cpus=3' in resources
    assert 'memory=256' in resources


def test_caller_failure_is_ignored_when_enabled(tmp_path):
    params = _params(
        tmp_path,
        {
            'attempt_log': str(tmp_path / 'attempts.txt'),
            'calling': {'allow_caller_failure': True},
        },
    )
    result = _run_fixture(CALLER_FIXTURE, params, tmp_path)

    output = result.stdout + result.stderr
    assert result.returncode == 0, output
    assert (tmp_path / 'attempts.txt').read_text().split() == ['1']


@pytest.mark.parametrize(
    'max_retries, expected_attempts',
    [(3, ['1', '2', '3', '4']), (1, ['1', '2'])],
)
def test_caller_failure_retries_max_retries_times(
        tmp_path, max_retries, expected_attempts):
    params = _params(
        tmp_path,
        {
            'attempt_log': str(tmp_path / 'attempts.txt'),
            'calling': {'allow_caller_failure': False},
            'resources': {'max_retries': max_retries},
        },
    )
    result = _run_fixture(CALLER_FIXTURE, params, tmp_path)

    assert result.returncode != 0
    attempts = (tmp_path / 'attempts.txt').read_text().split()
    assert attempts == expected_attempts
