"""Reproducibility guards for the environment specification."""

from conftest import (  # noqa: I100,I202  (tests run from this directory)
    REPO_ROOT,
)


def test_environment_uses_floating_minimum_versions():
    environment = (REPO_ROOT / 'environment.yaml').read_text()
    dependencies = environment.split('dependencies:', 1)[1]
    for line in dependencies.splitlines():
        stripped = line.strip()
        if not stripped.startswith('- '):
            continue
        assert '>=' in stripped
        assert stripped.count('=') == stripped.count('>=')
