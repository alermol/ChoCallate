#!/usr/bin/env python3
"""Strip BED header lines and validate contigs against a reference .fai.

BED inputs may contain UCSC ``track`` and ``browser`` lines and ``#``
comments.  ``sort-bed`` treats them as malformed records, so they must be
removed before sorting.  Every remaining record is checked against the
reference ``.fai``: a BED built for another assembly used to produce an
empty intersection and a generic "empty coverage" error, while this script
names the offending contigs and the BED file.
"""

import argparse
import sys

HEADER_WORDS = ('track', 'browser')


def is_header(line):
    """Return True for blank lines, comments and UCSC header lines."""
    text = line.rstrip('\r\n')
    if not text.strip():
        return True
    if text.startswith('#'):
        return True
    return any(
        text == word or text.startswith(word + ' ')
        or text.startswith(word + '\t')
        for word in HEADER_WORDS
    )


def read_contigs(path):
    """Return the contig names of a reference ``.fai`` index."""
    contigs = set()
    with open(path, encoding='utf-8') as handle:
        for line in handle:
            fields = line.rstrip('\r\n').split('\t')
            if len(fields) < 2 or not fields[0]:
                continue
            contigs.add(fields[0])
    return contigs


def main():
    parser = argparse.ArgumentParser(
        description='Strip BED headers and validate contigs against a .fai.')
    parser.add_argument('fai', help='Reference FASTA .fai index')
    parser.add_argument('bed', help='BED input to filter')
    args = parser.parse_args()

    contigs = read_contigs(args.fai)
    unknown = set()
    with open(args.bed, encoding='utf-8') as handle:
        for line in handle:
            if is_header(line):
                continue
            text = line.rstrip('\r\n')
            contig = text.split('\t', 1)[0]
            if contig not in contigs:
                unknown.add(contig)
                continue
            sys.stdout.write(text + '\n')

    if unknown:
        names = ', '.join(repr(name) for name in sorted(unknown))
        sys.exit(
            "ERROR: BED input '{}' uses contig(s) absent from the "
            'reference .fai: {}. Use the contig names of the reference '
            'FASTA.'.format(args.bed, names))


if __name__ == '__main__':
    main()
