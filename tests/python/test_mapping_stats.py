"""Tests for the mapping statistics summary helper."""

import subprocess
import sys

from conftest import BIN_DIR  # noqa: I100,I202  (tests run from this directory)

SCRIPT = BIN_DIR / 'mapping_stats.py'

RAW_FLAGSTAT = """\
10 + 0 in total (QC-passed reads + QC-failed reads)
10 + 0 primary
0 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
8 + 0 mapped (80.00% : N/A)
10 + 0 paired in sequencing
5 + 0 read1
5 + 0 read2
6 + 0 properly paired (60.00% : N/A)
6 + 0 with itself and mate mapped
1 + 0 singletons (10.00% : N/A)
1 + 0 with mate mapped to a different chr
1 + 0 with mate mapped to a different chr (mapQ>=5)
"""

FINAL_FLAGSTAT = """\
8 + 0 in total (QC-passed reads + QC-failed reads)
8 + 0 primary
0 + 0 secondary
0 + 0 supplementary
1 + 0 duplicates
8 + 0 mapped (100.00% : N/A)
8 + 0 paired in sequencing
4 + 0 read1
4 + 0 read2
6 + 0 properly paired (75.00% : N/A)
6 + 0 with itself and mate mapped
0 + 0 singletons (0.00% : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)
"""

STATS = """\
SN\traw total sequences:\t8\t# excluding supplementary and secondary reads
SN\tfiltered sequences:\t0
SN\tsequences:\t8
SN\treads mapped:\t8
SN\treads MQ0:\t2
SN\taverage length:\t150
SN\taverage quality:\t36.5
SN\tinsert size average:\t250.5
SN\tinsert size standard deviation:\t20.1
SN\terror rate:\t1.25e-03
SN\tbases mapped (cigar):\t1200
"""

IDXSTATS = """\
chr1\t1000\t8\t0
chr2\t500\t0\t0
*\t0\t0\t0
"""


def write_reports(root, sample, raw=RAW_FLAGSTAT,
                  final=FINAL_FLAGSTAT, stats=STATS, idxstats=IDXSTATS):
    """Write the four reports used by the script for one sample."""
    dirs = {
        'raw': root / 'raw',
        'flagstat': root / 'flagstat',
        'stats': root / 'stats',
        'idxstats': root / 'idxstats',
    }
    for path in dirs.values():
        path.mkdir(parents=True, exist_ok=True)
    (dirs['raw'] / f'{sample}.raw.flagstat').write_text(raw)
    (dirs['flagstat'] / f'{sample}.final.flagstat').write_text(final)
    (dirs['stats'] / f'{sample}.stats').write_text(stats)
    (dirs['idxstats'] / f'{sample}.idxstats').write_text(idxstats)
    return dirs


def run_summary(root, samples, dirs):
    """Run the summary script and return its output path."""
    order = root / 'sample_order.txt'
    order.write_text('\n'.join(samples) + '\n')
    output = root / 'mapping_summary.tsv'
    result = subprocess.run(
        [
            sys.executable, str(SCRIPT),
            '--sample-order', str(order),
            '--raw-dir', str(dirs['raw']),
            '--flagstat-dir', str(dirs['flagstat']),
            '--stats-dir', str(dirs['stats']),
            '--idxstats-dir', str(dirs['idxstats']),
            '--output', str(output),
        ],
        capture_output=True,
        text=True,
    )
    return result, output


def table(path):
    """Return the summary header and rows split on tabs."""
    lines = path.read_text().splitlines()
    return lines[0].split('\t'), [
        dict(zip(lines[0].split('\t'), line.split('\t')))
        for line in lines[1:]
    ]


def test_summary_reports_counts_and_rates(tmp_path):
    dirs = write_reports(tmp_path, 'sample1')

    result, output = run_summary(tmp_path, ['sample1'], dirs)

    assert result.returncode == 0, result.stderr
    header, rows = table(output)
    assert header[0] == 'sample'
    assert 'raw_mapped_pct' in header
    assert len(rows) == 1
    row = rows[0]
    assert row['sample'] == 'sample1'
    assert row['raw_total'] == '10'
    assert row['raw_mapped'] == '8'
    assert row['raw_mapped_pct'] == '80.00'
    assert row['raw_properly_paired_pct'] == '60.00'
    assert row['raw_singletons'] == '1'
    assert row['raw_mate_different_chr'] == '1'
    assert row['final_total'] == '8'
    assert row['final_duplicates'] == '1'
    assert row['final_duplicate_pct'] == '12.50'
    assert row['final_properly_paired_pct'] == '75.00'
    assert row['final_mq0'] == '2'
    assert row['final_mq0_pct'] == '25.00'
    assert row['average_read_length'] == '150'
    assert row['average_mapq'] == '36.5'
    assert row['insert_size_mean'] == '250.5'
    assert row['error_rate'] == '1.25e-03'
    assert row['mapped_bases'] == '1200'
    assert row['contigs_with_reads'] == '1'
    assert row['contigs_without_reads'] == '1'


def test_summary_follows_samplesheet_order(tmp_path):
    dirs = write_reports(tmp_path, 'sample_b')
    write_reports(tmp_path, 'sample_a')

    result, output = run_summary(
        tmp_path, ['sample_b', 'sample_a'], dirs
    )

    assert result.returncode == 0, result.stderr
    _, rows = table(output)
    assert [row['sample'] for row in rows] == ['sample_b', 'sample_a']


def test_missing_report_is_an_error(tmp_path):
    dirs = write_reports(tmp_path, 'sample1')
    (dirs['stats'] / 'sample1.stats').unlink()

    result, _ = run_summary(tmp_path, ['sample1'], dirs)

    assert result.returncode != 0
    assert 'missing mapping statistics' in result.stderr
    assert 'sample1.stats' in result.stderr
