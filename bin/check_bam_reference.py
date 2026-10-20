#!/usr/bin/env python3
"""Check that a BAM header matches a reference .fai index.

A BAM aligned to a different assembly, or to a reference whose contigs have
different lengths, would silently produce calls at wrong coordinates.
``PREPARE_BAM`` and ``postprocess_bam.sh`` therefore compare the ``@SQ``
names and lengths of the BAM with the generated ``.fai`` and report every
difference as a diff.  When a sequence dictionary is supplied, matching
``M5`` checksums are compared as well, which catches a different assembly
with identically named and sized contigs.
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


def read_bam_m5(path):
    """Return the ``@SQ`` ``M5`` checksums of a BAM, keyed by contig."""
    with pysam.AlignmentFile(path, 'rb') as bam:
        entries = bam.header.to_dict().get('SQ', [])
    return {
        entry['SN']: entry['M5'].lower()
        for entry in entries
        if 'SN' in entry and 'M5' in entry
    }


def read_dict_m5(path):
    """Return the ``@SQ`` ``M5`` checksums of a .dict, keyed by contig."""
    checksums = {}
    with open(path, encoding='utf-8') as handle:
        for line in handle:
            if not line.startswith('@SQ'):
                continue
            fields = {}
            for field in line.rstrip('\n').split('\t')[1:]:
                key, separator, value = field.partition(':')
                if separator:
                    fields[key] = value
            name = fields.get('SN')
            checksum = fields.get('M5')
            if name and checksum:
                checksums[name] = checksum.lower()
    return checksums


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


def diff_m5(bam_m5, dict_m5):
    """Return ``M5`` mismatches for contigs present in both headers."""
    return [
        '{}\tBAM {} != reference {}'.format(
            name, bam_m5[name], dict_m5[name])
        for name in sorted(bam_m5.keys() & dict_m5.keys())
        if bam_m5[name] != dict_m5[name]
    ]


def main():
    parser = argparse.ArgumentParser(
        description='Fail when the BAM @SQ entries differ from the .fai.')
    parser.add_argument('bam', help='Input BAM')
    parser.add_argument('fai', help='Reference FASTA .fai index')
    parser.add_argument('--dict', dest='dict_path',
                        help='Reference sequence dictionary (.dict)')
    args = parser.parse_args()

    errors = []
    bam_entries = read_bam_sq(args.bam)
    fai_entries = read_fai(args.fai)
    if bam_entries != fai_entries:
        diff = diff_entries(bam_entries, fai_entries)
        errors.append(
            ("ERROR: BAM header @SQ entries do not match reference "
             "'{}':\n{}\nThe BAM must be aligned to the reference passed "
             'as input.reference_genome.').format(
                 args.fai, '\n'.join(diff)))

    if args.dict_path:
        mismatches = diff_m5(read_bam_m5(args.bam),
                             read_dict_m5(args.dict_path))
        if mismatches:
            errors.append(
                ("ERROR: BAM @SQ M5 checksums do not match reference "
                 "dictionary '{}':\n{}\nThe BAM must be aligned to the "
                 'reference passed as input.reference_genome.').format(
                     args.dict_path, '\n'.join(mismatches)))

    if errors:
        sys.exit('\n'.join(errors))


if __name__ == '__main__':
    main()
