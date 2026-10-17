"""Tests for the per-task cleanup helper.

The helper must remove staged files only after a successful task and must
never delete a TMPDIR that was not created inside the task work directory.
"""

import os
import subprocess

from conftest import BIN_DIR  # noqa: I100,I202  (tests run from this directory)

CLEANUP = BIN_DIR / 'stage_cleanup.sh'


def _workdir(tmp_path):
    work = tmp_path / 'work'
    (work / 'tmp').mkdir(parents=True)
    (work / 'tmp' / 'payload.txt').write_text('payload')
    task_tmp = work / 'task_tmp'
    task_tmp.mkdir()
    (task_tmp / 'scratch.txt').write_text('scratch')
    return work, task_tmp


def _run(work, *args, tmpdir=None):
    env = os.environ.copy()
    env.pop('TMPDIR', None)
    if tmpdir is not None:
        env['TMPDIR'] = str(tmpdir)
    return subprocess.run(
        [str(CLEANUP), *args],
        cwd=work,
        env=env,
        capture_output=True,
        text=True,
    )


def test_failure_preserves_staged_and_temporary_files(tmp_path):
    work, task_tmp = _workdir(tmp_path)

    result = _run(work, '137', tmpdir=task_tmp)

    assert result.returncode == 0
    assert (work / 'tmp' / 'payload.txt').is_file()
    assert (task_tmp / 'scratch.txt').is_file()


def test_success_removes_staged_and_temporary_files(tmp_path):
    work, task_tmp = _workdir(tmp_path)

    result = _run(work, '0', tmpdir=task_tmp)

    assert result.returncode == 0
    assert not (work / 'tmp').exists()
    assert not task_tmp.exists()


def test_empty_tmpdir_is_ignored(tmp_path):
    work, _ = _workdir(tmp_path)

    result = _run(work, '0', tmpdir='')

    assert result.returncode == 0
    assert not (work / 'tmp').exists()


def test_tmpdir_outside_the_workdir_is_preserved(tmp_path):
    work, _ = _workdir(tmp_path)
    outside = tmp_path / 'outside'
    outside.mkdir()
    (outside / 'keep.txt').write_text('keep')

    result = _run(work, '0', tmpdir=outside)

    assert result.returncode == 0
    assert not (work / 'tmp').exists()
    assert (outside / 'keep.txt').is_file()


def test_missing_status_is_treated_as_success(tmp_path):
    work, task_tmp = _workdir(tmp_path)

    result = _run(work, tmpdir=task_tmp)

    assert result.returncode == 0
    assert not (work / 'tmp').exists()
    assert not task_tmp.exists()
