#!/usr/bin/env python3
"""Build the optional Cython accelerator for generate_consensus.py in place.

Usage:
    python bin/build_consensus_accel.py

Produces ``bin/consensus_accel.<abi>.so``.  generate_consensus.py imports it if
present and otherwise falls back to an equivalent pure-pysam implementation.
"""

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    try:
        import pysam
    except ImportError:
        sys.exit("pysam is required to build the consensus accelerator")
    try:
        from Cython.Build import cythonize
    except ImportError:
        sys.exit("Cython is required to build the consensus accelerator")

    try:
        from setuptools import Distribution, Extension
    except ImportError:
        sys.exit("setuptools is required to build the consensus accelerator")

    ext = Extension(
        "consensus_accel",
        [os.path.join(HERE, "consensus_accel.pyx")],
        include_dirs=pysam.get_include(),
        define_macros=pysam.get_defines(),
        # resolve htslib symbols against the copies bundled in pysam
        extra_objects=pysam.get_libraries(),
    )
    dist = Distribution({
        "name": "consensus_accel",
        "ext_modules": cythonize([ext], language_level=3, quiet=True),
    })
    build_ext = dist.get_command_obj("build_ext")
    build_ext.inplace = True
    build_ext.build_lib = HERE
    build_ext.build_temp = os.path.join(HERE, ".build_accel")
    build_ext.ensure_finalized()
    # build from bin/ so the in-place copy lands next to the .pyx
    old_cwd = os.getcwd()
    os.chdir(HERE)
    try:
        build_ext.run()
    finally:
        os.chdir(old_cwd)
    built = [f for f in os.listdir(HERE) if f.startswith("consensus_accel")
             and f.endswith(".so")]
    print("built:", ", ".join(built))


if __name__ == "__main__":
    main()
