"""Reproducibility guards for the pinned environment."""

from conftest import (  # noqa: I100,I202  (tests run from this directory)
    REPO_ROOT,
)


def test_environment_uses_exact_pins_and_has_a_lock_file():
    environment = (REPO_ROOT / 'environment.yaml').read_text()
    dependencies = environment.split('dependencies:', 1)[1]
    for line in dependencies.splitlines():
        stripped = line.strip()
        if not stripped.startswith('- '):
            continue
        assert '>=' not in stripped
        assert '=' in stripped

    lock = REPO_ROOT / 'environment.lock.yaml'
    assert lock.is_file()
    lock_text = lock.read_text()
    assert 'name: ChoCallate' in lock_text
    assert 'prefix:' not in lock_text
    assert 'bcftools=1.24=' in lock_text
    assert 'gatk4=4.6.2.0=' in lock_text
