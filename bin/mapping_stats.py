#!/usr/bin/env python3
"""Collect per-sample mapping statistics into one TSV summary.

The pipeline reports ``samtools flagstat`` of the aligner output before MAPQ
filtering and duplicate handling (``<sample>.raw.flagstat``), plus
``samtools flagstat``, ``stats`` and ``idxstats`` of the final BAM
(``<sample>.final.flagstat``, ``<sample>.stats``, ``<sample>.idxstats``).
This script reads those files in samplesheet order and writes one row per
sample.  Counts are alignment records; a paired-end read pair contributes
two records.
"""

import argparse
import csv
import re
import sys
from pathlib import Path

COLUMNS = [
    'sample',
    'raw_total',
    'raw_primary',
    'raw_secondary',
    'raw_supplementary',
    'raw_mapped',
    'raw_mapped_pct',
    'raw_paired',
    'raw_properly_paired',
    'raw_properly_paired_pct',
    'raw_duplicates',
    'raw_duplicate_pct',
    'raw_singletons',
    'raw_mate_different_chr',
    'final_total',
    'final_secondary',
    'final_supplementary',
    'final_duplicates',
    'final_duplicate_pct',
    'final_properly_paired',
    'final_properly_paired_pct',
    'final_mq0',
    'final_mq0_pct',
    'average_read_length',
    'average_mapq',
    'insert_size_mean',
    'insert_size_sd',
    'error_rate',
    'mapped_bases',
    'contigs_with_reads',
    'contigs_without_reads',
]

# Flagstat labels that carry a read count.  The parenthetical suffix (a
# percentage or the QC explanation) is removed before the lookup.
_FLAGSTAT_LABELS = {
    'in total': 'total',
    'primary': 'primary',
    'secondary': 'secondary',
    'supplementary': 'supplementary',
    'duplicates': 'duplicates',
    'mapped': 'mapped',
    'paired in sequencing': 'paired',
    'properly paired': 'properly_paired',
    'singletons': 'singletons',
    'with mate mapped to a different chr': 'mate_different_chr',
}
_FLAGSTAT_LINE = re.compile(r'^(\d+) \+ (\d+) (.+)$')
_PARENTHETICAL = re.compile(r'\s+\(.*\)$')


def parse_flagstat(path):
    """Return flagstat read counts for one BAM.

    Args:
        path: ``samtools flagstat`` output.

    Returns:
        A dict with the labels from ``_FLAGSTAT_LABELS`` that are present.
    """
    counts = {}
    for line in path.read_text(encoding='utf-8').splitlines():
        match = _FLAGSTAT_LINE.match(line)
        if match is None:
            continue
        passed, failed, rest = match.groups()
        # The mapQ>=5 mate line has the same prefix as the line above it.
        if 'mapQ>=5' in rest:
            continue
        label = _PARENTHETICAL.sub('', rest)
        key = _FLAGSTAT_LABELS.get(label)
        if key is not None:
            counts[key] = int(passed) + int(failed)
    missing = [key for key in ('total', 'mapped') if key not in counts]
    if missing:
        raise ValueError(
            f"{path}: missing flagstat line(s): {', '.join(missing)}"
        )
    return counts


def parse_stats(path):
    """Return ``samtools stats`` SN values keyed by a normalised name."""
    values = {}
    for line in path.read_text(encoding='utf-8').splitlines():
        if not line.startswith('SN\t'):
            continue
        fields = line.split('\t')
        if len(fields) < 3:
            continue
        key = fields[1].rstrip(':').strip().lower()
        key = re.sub(r'[^a-z0-9]+', '_', key).strip('_')
        values[key] = fields[2].strip()
    return values


def parse_idxstats(path):
    """Return the number of contigs and the number with mapped reads.

    The ``*`` line reports unmapped reads and is not a reference sequence.
    """
    contigs = 0
    with_reads = 0
    for line in path.read_text(encoding='utf-8').splitlines():
        fields = line.split('\t')
        if len(fields) < 4 or fields[0] == '*':
            continue
        contigs += 1
        if int(fields[2]) > 0:
            with_reads += 1
    return contigs, with_reads


def percent(numerator, denominator):
    """Format a percentage with two decimals, or ``NA`` without reads."""
    if not denominator:
        return 'NA'
    return f'{100.0 * numerator / denominator:.2f}'


def integer(value):
    """Convert a stats value to int, tolerating floats and missing values."""
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return 0


def build_row(sample, raw, final, stats, contigs):
    """Build one summary row from the parsed per-sample reports."""
    raw_total = raw.get('total', 0)
    raw_mapped = raw.get('mapped', 0)
    raw_paired = raw.get('paired', 0)
    raw_proper = raw.get('properly_paired', 0)
    raw_dup = raw.get('duplicates', 0)

    final_total = final.get('total', 0)
    final_dup = final.get('duplicates', 0)
    final_proper = final.get('properly_paired', 0)
    sequences = integer(stats.get('sequences'))
    mq0 = integer(stats.get('reads_mq0'))
    contigs_total, contigs_with = contigs

    return {
        'sample': sample,
        'raw_total': raw_total,
        'raw_primary': raw.get('primary', 0),
        'raw_secondary': raw.get('secondary', 0),
        'raw_supplementary': raw.get('supplementary', 0),
        'raw_mapped': raw_mapped,
        'raw_mapped_pct': percent(raw_mapped, raw_total),
        'raw_paired': raw_paired,
        'raw_properly_paired': raw_proper,
        'raw_properly_paired_pct': percent(raw_proper, raw_total),
        'raw_duplicates': raw_dup,
        'raw_duplicate_pct': percent(raw_dup, raw_total),
        'raw_singletons': raw.get('singletons', 0),
        'raw_mate_different_chr': raw.get('mate_different_chr', 0),
        'final_total': final_total,
        'final_secondary': final.get('secondary', 0),
        'final_supplementary': final.get('supplementary', 0),
        'final_duplicates': final_dup,
        'final_duplicate_pct': percent(final_dup, final_total),
        'final_properly_paired': final_proper,
        'final_properly_paired_pct': percent(final_proper, final_total),
        'final_mq0': mq0,
        'final_mq0_pct': percent(mq0, sequences),
        'average_read_length': stats.get('average_length', 'NA'),
        'average_mapq': stats.get('average_quality', 'NA'),
        'insert_size_mean': stats.get('insert_size_average', 'NA'),
        'insert_size_sd': stats.get(
            'insert_size_standard_deviation', 'NA'
        ),
        'error_rate': stats.get('error_rate', 'NA'),
        'mapped_bases': stats.get('bases_mapped_cigar', 'NA'),
        'contigs_with_reads': contigs_with,
        'contigs_without_reads': contigs_total - contigs_with,
    }


def main(argv=None):
    """Parse the report directories and write ``mapping_summary.tsv``."""
    parser = argparse.ArgumentParser(
        description='Summarise per-sample mapping statistics.'
    )
    parser.add_argument('--sample-order', required=True, type=Path)
    parser.add_argument('--raw-dir', required=True, type=Path)
    parser.add_argument('--flagstat-dir', required=True, type=Path)
    parser.add_argument('--stats-dir', required=True, type=Path)
    parser.add_argument('--idxstats-dir', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args(argv)

    sample_ids = [
        line.strip()
        for line in args.sample_order.read_text(
            encoding='utf-8'
        ).splitlines()
        if line.strip()
    ]
    if not sample_ids:
        parser.error('sample order is empty')

    rows = []
    for sample in sample_ids:
        raw_path = args.raw_dir / f'{sample}.raw.flagstat'
        final_path = args.flagstat_dir / f'{sample}.final.flagstat'
        stats_path = args.stats_dir / f'{sample}.stats'
        idxstats_path = args.idxstats_dir / f'{sample}.idxstats'
        missing = [
            str(path)
            for path in (raw_path, final_path, stats_path, idxstats_path)
            if not path.is_file()
        ]
        if missing:
            raise SystemExit(
                'missing mapping statistics: ' + ', '.join(missing)
            )
        rows.append(build_row(
            sample,
            parse_flagstat(raw_path),
            parse_flagstat(final_path),
            parse_stats(stats_path),
            parse_idxstats(idxstats_path),
        ))

    with args.output.open('w', encoding='utf-8', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(COLUMNS)
        for row in rows:
            writer.writerow([row[column] for column in COLUMNS])
    return 0


if __name__ == '__main__':
    sys.exit(main())
