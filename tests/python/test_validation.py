"""Early-validation tests for the Nextflow entry point.

Every case here is expected to fail before the first process starts, so the
suite does not need mapping or calling tools.  A Nextflow installation and the
pipeline plugins are required; the tests are skipped otherwise.
"""

import json
import os
import shlex
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[2]
NEXTFLOW = shutil.which("nextflow")

pytestmark = pytest.mark.skipif(
    NEXTFLOW is None,
    reason="nextflow is not installed",
)


def _merge(dst, src):
    for key, value in src.items():
        if isinstance(value, dict) and isinstance(dst.get(key), dict):
            _merge(dst[key], value)
        else:
            dst[key] = value
    return dst


def _prepare_inputs(tmp_path, rows=("sample1\tr1.fq\tr2.fq",)):
    (tmp_path / "idx").mkdir(exist_ok=True)
    (tmp_path / "ref.fa").write_text(">chr1\nACGTACGTACGT\n")
    (tmp_path / "r1.fq").write_text("@r1\nACGT\n+\nIIII\n")
    (tmp_path / "r2.fq").write_text("@r2\nACGT\n+\nIIII\n")
    (tmp_path / "sample.bam").write_text("not a real bam")
    samplesheet = tmp_path / "samples.tsv"
    samplesheet.write_text("\n".join(rows) + "\n")
    return samplesheet


def _config(tmp_path, samplesheet=None, overrides=None):
    config = {
        "input": {
            "samples_tsv": str(samplesheet or tmp_path / "samples.tsv"),
            "reference_genome": str(tmp_path / "ref.fa"),
            "reference_index_dir": str(tmp_path / "idx"),
        },
        "output": {"directory": str(tmp_path / "out")},
        "mapping": {"cpu": 1},
        "calling": {"cpu": 1},
        "consensus": {"cpu": 1, "threshold": 1},
    }
    _merge(config, overrides or {})
    params = tmp_path / "params.json"
    params.write_text(json.dumps(config))
    return params


def _environment():
    env = os.environ.copy()
    if "NXF_OFFLINE" not in env:
        plugins = Path.home() / ".nextflow" / "plugins"
        has_plugins = (
            plugins.is_dir()
            and any(plugins.glob("nf-schema-*"))
            and any(plugins.glob("nf-boost-*"))
        )
        if has_plugins:
            env["NXF_OFFLINE"] = "true"
    return env


def _run(tmp_path, params=None, ulimit=None):
    command = [NEXTFLOW, "run", str(PROJECT_ROOT / "main.nf")]
    if params is not None:
        command += ["-params-file", str(params)]
    command += ["-ansi-log", "false"]
    if ulimit is not None:
        shell_command = "ulimit -n {}; exec {}".format(
            ulimit, " ".join(shlex.quote(part) for part in command)
        )
        command = ["bash", "-c", shell_command]
    return subprocess.run(
        command,
        cwd=tmp_path,
        env=_environment(),
        capture_output=True,
        text=True,
        timeout=180,
    )


def _output(result):
    return result.stdout + result.stderr


def test_pipeline_defaults_for_depth_and_mapq(tmp_path):
    result = subprocess.run(
        [NEXTFLOW, "config", str(PROJECT_ROOT / "main.nf")],
        cwd=tmp_path,
        env=_environment(),
        capture_output=True,
        text=True,
        timeout=60,
    )

    assert result.returncode == 0
    output = _output(result)
    assert "min_map_qual = 10" in output
    assert "min_coverage = 5" in output
    assert "compress_reference_blocks = false" in output


def test_missing_required_parameter_is_named(tmp_path):
    _prepare_inputs(tmp_path)
    params = _config(tmp_path)
    config = json.loads(params.read_text())
    del config["input"]["samples_tsv"]
    params.write_text(json.dumps(config))

    result = _run(tmp_path, params)

    assert result.returncode != 0
    output = _output(result)
    assert "Missing required parameter(s): samples_tsv" in output
    assert "--input" in output


def test_enum_reports_key_and_accepted_values(tmp_path):
    _prepare_inputs(tmp_path)
    params = _config(tmp_path, overrides={"input": {"format": "cram"}})

    result = _run(tmp_path, params)

    assert result.returncode != 0
    output = _output(result)
    assert "--input.format" in output
    assert "Expected any of [fastq, bam]" in output


def test_negative_integer_is_rejected(tmp_path):
    _prepare_inputs(tmp_path)
    params = _config(tmp_path, overrides={"coverage": {"min_coverage": -1}})

    result = _run(tmp_path, params)

    assert result.returncode != 0
    output = _output(result)
    assert "--coverage.min_coverage" in output
    assert "-1 is less than 0" in output


def test_unknown_nested_parameter_is_rejected(tmp_path):
    _prepare_inputs(tmp_path)
    params = _config(tmp_path, overrides={"input": {"samples_tvs": "typo"}})

    result = _run(tmp_path, params)

    assert result.returncode != 0
    assert "--input.samples_tvs" in _output(result)


def test_duplicate_sample_ids_are_rejected(tmp_path):
    samplesheet = _prepare_inputs(
        tmp_path,
        rows=("sample1\tr1.fq\tr2.fq", "sample1\tr1.fq\tr2.fq"),
    )
    params = _config(tmp_path, samplesheet=samplesheet)

    result = _run(tmp_path, params)

    assert result.returncode != 0
    assert "duplicate sample id 'sample1'" in _output(result)


def test_unsafe_sample_id_is_rejected(tmp_path):
    samplesheet = _prepare_inputs(tmp_path, rows=("bad id!\tr1.fq\tr2.fq",))
    params = _config(tmp_path, samplesheet=samplesheet)

    result = _run(tmp_path, params)

    assert result.returncode != 0
    assert "not shell-safe" in _output(result)


def test_wrong_samplesheet_width_is_rejected(tmp_path):
    samplesheet = _prepare_inputs(
        tmp_path,
        rows=("sample1\tr1.fq\tr2.fq\textra.fq",),
    )
    params = _config(tmp_path, samplesheet=samplesheet)

    result = _run(tmp_path, params)

    assert result.returncode != 0
    output = _output(result)
    assert "expected 3 tab-separated columns" in output
    assert "found 4" in output


def test_missing_read_file_is_rejected(tmp_path):
    samplesheet = _prepare_inputs(
        tmp_path,
        rows=("sample1\tmissing.fq\tr2.fq",),
    )
    params = _config(tmp_path, samplesheet=samplesheet)

    result = _run(tmp_path, params)

    assert result.returncode != 0
    output = _output(result)
    assert "missing.fq" in output
    assert "file not found" in output


def test_header_and_comment_lines_are_skipped(tmp_path):
    samplesheet = _prepare_inputs(
        tmp_path,
        rows=(
            "# samples",
            "sample_id\tr1\tr2",
            "",
            "sample1\tr1.fq\tr2.fq",
            "sample1\tr1.fq\tr2.fq",
        ),
    )
    params = _config(tmp_path, samplesheet=samplesheet)

    result = _run(tmp_path, params)

    assert result.returncode != 0
    # If the header or comment were parsed as data, this would be a width
    # error instead of the duplicate-id error from the two data rows.
    assert "duplicate sample id 'sample1'" in _output(result)


def test_bam_samplesheet_uses_two_columns(tmp_path):
    samplesheet = _prepare_inputs(
        tmp_path,
        rows=("sample1\tsample.bam", "sample1\tsample.bam"),
    )
    params = _config(
        tmp_path,
        samplesheet=samplesheet,
        overrides={"input": {"format": "bam"}},
    )

    result = _run(tmp_path, params)

    assert result.returncode != 0
    output = _output(result)
    assert "duplicate sample id 'sample1'" in output
    assert "expected 2 tab-separated columns" not in output


@pytest.mark.parametrize(
    "caller, override, option",
    [
        (
            "bcftools",
            {"call": {"extra_args": "--ploidy 3"}},
            "--ploidy",
        ),
        ("freebayes", {"extra_args": "-p 1"}, "-p"),
        ("freebayes", {"extra_args": "--ploidy=4"}, "--ploidy"),
        ("gatk", {"extra_args": "--sample-ploidy 1"}, "--sample-ploidy"),
        ("gatk", {"extra_args": "--ploidy=1"}, "--ploidy"),
    ],
)
def test_non_diploid_ploidy_is_rejected(tmp_path, caller, override, option):
    _prepare_inputs(tmp_path)
    params = _config(
        tmp_path,
        overrides={
            "calling": {
                "callers": caller,
                caller: override,
            },
        },
    )

    result = _run(tmp_path, params)

    assert result.returncode != 0
    output = _output(result)
    assert option in output
    assert "ploidy is not supported yet except diploid" in output


@pytest.mark.skipif(
    not sys.platform.startswith("linux"),
    reason="uses a POSIX descriptor limit",
)
def test_merge_file_descriptor_guard(tmp_path):
    rows = ["sample{}\tr1.fq\tr2.fq".format(index) for index in range(300)]
    samplesheet = _prepare_inputs(tmp_path, rows=tuple(rows))
    params = _config(
        tmp_path,
        samplesheet=samplesheet,
        overrides={"consensus": {"cpu": 1}},
    )

    result = _run(tmp_path, params, ulimit=256)

    assert result.returncode != 0
    output = _output(result)
    assert "merge step needs at least" in output
    assert "soft limit is 256" in output
