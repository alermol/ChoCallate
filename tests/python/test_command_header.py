"""Tests for the final ##ChoCallateCommand header line."""

import os
import subprocess
import sys

from conftest import BIN_DIR  # noqa: I100,I202  (tests run from this dir)

SCRIPT = BIN_DIR / 'write_command_header.py'


def _write(tmp_path, command_line):
    env = os.environ.copy()
    if command_line is None:
        env.pop('NXF_CLI', None)
    else:
        env['NXF_CLI'] = command_line
    output = tmp_path / 'command.hdr'
    subprocess.run([sys.executable, str(SCRIPT), str(output)],
                   check=True, env=env)
    return output.read_text()


def test_header_records_command_line(tmp_path):
    assert _write(tmp_path, 'nextflow run main.nf -resume') == (
        '##ChoCallateCommand="nextflow run main.nf -resume"\n')


def test_header_value_is_escaped(tmp_path):
    text = _write(tmp_path, 'nextflow run "a b"\\c')
    assert text == '##ChoCallateCommand="nextflow run \\"a b\\"\\\\c"\n'


def test_header_defaults_to_unknown_without_nxf_cli(tmp_path):
    assert _write(tmp_path, None) == '##ChoCallateCommand="unknown"\n'
