#!/usr/bin/env python3
"""Reject reference contigs that exceed the 32-bit coordinate range.

BAM ``pos``, BCF ``POS`` and GATK/Picard use signed 32-bit positions, so a
contig longer than ``2**31 - 1`` bp cannot be processed by this pipeline.
Checking the generated ``.fai`` before any mapping task gives an immediate,
actionable error instead of a failure deep inside a caller.
"""

import argparse
import sys

MAX_CONTIG_LENGTH = (1 << 31) - 1


def read_fai_lengths(path):
    """Yield ``(contig, length)`` pairs from a ``.fai`` file."""
    with open(path, encoding='utf-8') as handle:
        for line in handle:
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 2 or not fields[0]:
                continue
            try:
                length = int(fields[1])
            except ValueError:
                continue
            yield fields[0], length


def oversized_contigs(path, maximum=MAX_CONTIG_LENGTH):
    """Return the ``(contig, length)`` pairs longer than ``maximum``."""
    return [(contig, length)
            for contig, length in read_fai_lengths(path)
            if length > maximum]


def main():
    parser = argparse.ArgumentParser(
        description='Fail when a reference contig exceeds the 32-bit limit.')
    parser.add_argument('fai', help='Reference FASTA .fai index')
    args = parser.parse_args()

    oversized = oversized_contigs(args.fai)
    if not oversized:
        return
    details = ', '.join(
        "'{}' ({} bp)".format(contig, length)
        for contig, length in oversized)
    sys.exit(
        'ERROR: reference contig(s) exceed the maximum supported length of '
        '{} bp (2^31-1): {}. BAM/BCF positions and GATK use 32-bit '
        'coordinates; split these contigs or use a long-contig-capable '
        'workflow.'.format(MAX_CONTIG_LENGTH, details))


if __name__ == '__main__':
    main()
