#!/usr/bin/env python3
"""Check that a BAM header matches a reference .fai index.

A BAM aligned to a different assembly, or to a reference whose contigs have
different lengths, would silently produce calls at wrong coordinates.
``PREPARE_BAM`` therefore compares the ``@SQ`` names and lengths of the input
BAM with the generated ``.fai`` and reports every difference as a diff.
"""

import argparse
import difflib
import sys

import pysam


def read_fai(path):
    """Return the ``(contig, length)`` pairs of a ``.fai`` file."""
    entries = []
    with open(path, encoding='utf-8') as handle:
        for line in handle:
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 2 or not fields[0]:
                continue
            try:
                length = int(fields[1])
            except ValueError:
                continue
            entries.append((fields[0], length))
    return entries


def read_bam_sq(path):
    """Return the ``@SQ`` ``(name, length)`` pairs of a BAM header."""
    with pysam.AlignmentFile(path, 'rb') as bam:
        return list(zip(bam.references, bam.lengths))


def format_entries(entries):
    """Format ``(name, length)`` pairs for a line-based diff."""
    return ['{}\t{}'.format(name, length) for name, length in entries]


def diff_entries(bam_entries, fai_entries):
    """Return a unified diff between the BAM and .fai contig lists."""
    return list(difflib.unified_diff(
        format_entries(bam_entries),
        format_entries(fai_entries),
        fromfile='BAM @SQ',
        tofile='reference .fai',
        lineterm='',
    ))


def main():
    parser = argparse.ArgumentParser(
        description='Fail when the BAM @SQ entries differ from the .fai.')
    parser.add_argument('bam', help='Input BAM')
    parser.add_argument('fai', help='Reference FASTA .fai index')
    args = parser.parse_args()

    bam_entries = read_bam_sq(args.bam)
    fai_entries = read_fai(args.fai)
    if bam_entries == fai_entries:
        return

    diff = diff_entries(bam_entries, fai_entries)
    sys.exit(
        ("ERROR: BAM header @SQ entries do not match reference '{}':\n{}\n"
         'The BAM must be aligned to the reference passed as '
         'input.reference_genome.').format(args.fai, '\n'.join(diff)))


if __name__ == '__main__':
    main()
