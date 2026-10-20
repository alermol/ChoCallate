"""Reproducibility guards for the environment specification."""

from conftest import (  # noqa: I100,I202  (tests run from this directory)
    REPO_ROOT,
)


def test_environment_pins_exact_dependency_versions():
    environment = (REPO_ROOT / 'environment.yaml').read_text()
    dependencies = environment.split('dependencies:', 1)[1]
    pins = [
        line.strip()[2:]
        for line in dependencies.splitlines()
        if line.strip().startswith('- ')
    ]
    assert pins
    for pin in pins:
        assert '>=' not in pin
        name, version, *build = pin.split('=')
        assert name
        assert version[0].isdigit()
        assert len(build) <= 1
