"""Tests for the version-recording helper."""

import shutil
import subprocess

import pytest

from conftest import BIN_DIR  # noqa: I100,I202  (tests run from this directory)

WRITE_VERSIONS = BIN_DIR / 'write_versions.sh'

pytestmark = pytest.mark.skipif(
    shutil.which('bcftools') is None or shutil.which('samtools') is None,
    reason='bcftools and samtools are required',
)


def test_write_versions_yaml_is_sorted(tmp_path):
    output = tmp_path / 'versions.yml'
    subprocess.run(
        [str(WRITE_VERSIONS), str(output), 'samtools', 'bcftools'],
        check=True,
    )

    lines = output.read_text().splitlines()
    assert any(line.startswith('bcftools: "') for line in lines)
    assert any(line.startswith('samtools: "') for line in lines)
    assert lines == sorted(lines)


def test_write_versions_records_pipeline_version(tmp_path):
    output = tmp_path / 'versions.yml'
    subprocess.run(
        [str(WRITE_VERSIONS), '--pipeline-version', '4.0.0',
         str(output), 'bcftools'],
        check=True,
    )

    lines = output.read_text().splitlines()
    assert 'choCallate: "4.0.0"' in lines
    assert lines == sorted(lines)


def test_write_versions_flat_format(tmp_path):
    output = tmp_path / 'versions.txt'
    subprocess.run(
        [str(WRITE_VERSIONS), '--format', 'flat', str(output),
         'bcftools', 'samtools'],
        check=True,
    )

    lines = output.read_text().splitlines()
    assert any(line.startswith('bcftools=') for line in lines)
    assert any(line.startswith('samtools=') for line in lines)
    assert all(': ' not in line for line in lines)
