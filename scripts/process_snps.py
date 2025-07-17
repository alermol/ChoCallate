#!/usr/bin/env python3

import argparse
import sqlite3

def sort_gt(s):
    parts = s.split('/')
    numbers = [int(part) for part in parts]
    numbers_sorted = sorted(numbers)
    return '/'.join(map(str, numbers_sorted))


# def extract_contigs_list(vcf):
#     with open()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--vcf1', required=True, help='Path to VCF file 1 (bcftools)')
    parser.add_argument('--vcf2', required=True, help='Path to VCF file 2')
    parser.add_argument('--vcf3', required=True, help='Path to VCF file 3')
    parser.add_argument('--vcf4', required=True, help='Path to VCF file 4')
    parser.add_argument('--vcf5', required=True, help='Path to VCF file 5')
    parser.add_argument('--sample', required=True, help='Sample name')
    parser.add_argument('--chr', required=True, help='Chromosome name')
    args = parser.parse_args()

    conn = sqlite3.connect(':memory:')
    cursor = conn.cursor()
    
    cursor.execute('''
        CREATE TABLE variants (
            chrom TEXT NOT NULL,
            pos INTEGER NOT NULL,
            ref TEXT NOT NULL,
            alt TEXT NOT NULL,
            gt TEXT NOT NULL,
            source INTEGER NOT NULL,
            PRIMARY KEY (chrom, pos, source)
        )
    ''')
    cursor.execute('''
        CREATE TABLE base_positions (
            chrom TEXT NOT NULL,
            pos INTEGER NOT NULL,
            PRIMARY KEY (chrom, pos)
        )
    ''')
    conn.commit()

    base_positions_set = set()
    
    with open(args.vcf1) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if '.' in fields[9]:
                continue
            chrom, pos, ref, alt = fields[0], int(fields[1]), fields[3], fields[4]
            gt = sort_gt(fields[9])
            base_positions_set.add((chrom, pos))
            cursor.execute('''
                INSERT OR REPLACE INTO variants (chrom, pos, ref, alt, gt, source)
                VALUES (?, ?, ?, ?, ?, 1)
            ''', (chrom, pos, ref, alt, gt))
            cursor.execute('''
                INSERT OR REPLACE INTO base_positions (chrom, pos)
                VALUES (?, ?)
            ''', (chrom, pos))
    conn.commit()

    other_vcfs = [
        (2, args.vcf2),
        (3, args.vcf3),
        (4, args.vcf4),
        (5, args.vcf5)
    ]
    for source_id, file_path in other_vcfs:
        with open(file_path) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                fields = line.strip().split('\t')
                if '.' in fields[9]:
                    continue
                chrom, pos, ref, alt = fields[0], int(fields[1]), fields[3], fields[4]
                if (chrom, pos) not in base_positions_set:
                    continue
                gt = sort_gt(fields[9])
                cursor.execute('''
                    INSERT OR REPLACE INTO variants (chrom, pos, ref, alt, gt, source)
                    VALUES (?, ?, ?, ?, ?, ?)
                ''', (chrom, pos, ref, alt, gt, source_id))
        conn.commit()

    cursor.execute('''
        CREATE TABLE genotypes AS
        SELECT base.chrom, base.pos,
               COALESCE(v.ref, base.ref) AS ref,
               COALESCE(v.alt, base.alt) AS alt,
               COALESCE(v.gt, base.gt) AS gt,
               s.source
        FROM (SELECT chrom, pos, ref, alt, gt FROM variants WHERE source=1) base
        CROSS JOIN (SELECT 1 AS source UNION SELECT 2 UNION SELECT 3 UNION SELECT 4 UNION SELECT 5) s
        LEFT JOIN variants v ON base.chrom = v.chrom AND base.pos = v.pos AND v.source = s.source
    ''')

    cursor.execute('''
        SELECT chrom, pos, ref, alt, gt
        FROM (
            SELECT chrom, pos, ref, alt, gt, COUNT(*) AS cnt,
                   ROW_NUMBER() OVER (PARTITION BY chrom, pos, ref, alt, gt ORDER BY COUNT(*) DESC) AS rn
            FROM genotypes
            GROUP BY chrom, pos, ref, alt, gt
        ) t
        WHERE rn = 1 AND cnt >= 3 AND gt NOT LIKE '%.%'
        ORDER BY chrom, pos, ref, alt, gt
    ''')
    results = cursor.fetchall()

    with open(f"all_chrs_snps/{args.chr}.vcf", 'w') as out:
        out.write('##fileformat=VCFv4.3\n')
        out.write('##FORMAT=<ID=GT,Number=1,Type=String>\n')
        out.write(f'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{args.sample}\n')
        for chrom, pos, ref, alt, gt in results:
            out.write('\t'.join([chrom, str(pos), '.', ref.upper(), alt, '.', '.', '.', 'GT', f'{gt}\n']))

if __name__ == "__main__":
    main()