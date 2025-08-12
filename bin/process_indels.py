#!/usr/bin/env python3

import argparse
import sqlite3

def sort_gt(s):
    parts = s.split('/')
    numbers = sorted(int(part) for part in parts)
    return '/'.join(map(str, numbers))

def parse_vcf(vcf_file):
    with open(vcf_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if '.' in fields[9]:
                continue
            chrom = fields[0]
            pos = int(fields[1])
            ref = fields[3]
            alt = fields[4]
            gt = sort_gt(fields[9])
            yield (chrom, pos, ref, alt, gt)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--vcf1', required=True, help='Path to VCF file 1 (bcftools)')
    parser.add_argument('--vcf2', required=True, help='Path to VCF file 2')
    parser.add_argument('--vcf3', required=True, help='Path to VCF file 3')
    parser.add_argument('--vcf4', required=True, help='Path to VCF file 4')
    parser.add_argument('--vcf5', required=True, help='Path to VCF file 5')
    parser.add_argument('--sample', required=True, help='Sample name')
    parser.add_argument('--chr', required=True, help='Chromosome name')
    parser.add_argument('--cons_threshold', type=int, default=3, help='Consensus threshold')
    args = parser.parse_args()

    conn = sqlite3.connect(':memory:')
    c = conn.cursor()
    c.execute('''
        CREATE TABLE variants (
            chrom TEXT NOT NULL,
            pos INTEGER NOT NULL,
            ref TEXT NOT NULL,
            alt TEXT NOT NULL,
            gt TEXT NOT NULL,
            source INTEGER NOT NULL,
            PRIMARY KEY (chrom, pos, source))
    ''')

    vcf_files = [
        (1, args.vcf1),
        (2, args.vcf2),
        (3, args.vcf3),
        (4, args.vcf4),
        (5, args.vcf5)
    ]

    for source_id, vcf_file in vcf_files:
        for record in parse_vcf(vcf_file):
            chrom, pos, ref, alt, gt = record
            try:
                c.execute('''
                    INSERT INTO variants (chrom, pos, ref, alt, gt, source)
                    VALUES (?, ?, ?, ?, ?, ?)
                ''', (chrom, pos, ref, alt, gt, source_id))
            except sqlite3.IntegrityError:
                continue
    conn.commit()

    query = '''
        WITH base_coords AS (
        	SELECT DISTINCT chrom, pos, ref, alt, gt FROM variants
        ),
        aggregated AS (
        	SELECT 
        		v.chrom, 
        		v.pos, 
        		v.ref, 
        		v.alt, 
        		v.gt, 
        		COUNT(*) AS cnt
        	FROM variants v
        	JOIN base_coords b 
        		ON v.chrom = b.chrom AND v.pos = b.pos AND v.ref = b.ref AND v.alt = b.alt AND v.gt = b.gt
        	GROUP BY v.chrom, v.pos, v.ref, v.alt, v.gt
        ),
        ranked AS (
        	SELECT 
        		*,
        		ROW_NUMBER() OVER (
        			PARTITION BY chrom, pos, ref, alt, gt
        			ORDER BY cnt DESC
        		) AS rn
        	FROM aggregated
        )
        SELECT chrom, pos, ref, alt, gt, cnt
        FROM ranked
        WHERE rn = 1
    '''
    c.execute(query)
    results = c.fetchall()

    with open(f"all_chrs/{args.chr}.vcf", 'w') as out:
        out.write('##fileformat=VCFv4.3\n')
        out.write('##FORMAT=<ID=GT,Number=1,Type=String>\n')
        out.write(f'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{args.sample}\n')
        
        for row in results:
            chrom, pos, ref, alt, gt, cnt = row
            if cnt >= args.cons_threshold and '.' not in gt:
                out.write('\t'.join([chrom, str(pos), '.', ref.upper(), alt, '.', '.', '.', 'GT', f"{gt}\n"]))

if __name__ == "__main__":
    main()
