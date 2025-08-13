#!/usr/bin/env python3

import argparse
import sqlite3

def sort_gt(s):
    parts = s.split('/')
    numbers = [int(part) for part in parts]
    numbers_sorted = sorted(numbers)
    return '/'.join(map(str, numbers_sorted))

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--vcfs', required=True, help='Comma-separated paths to VCF files')
    parser.add_argument('--sample', required=True, help='Sample name')
    parser.add_argument('--chr', required=True, help='Chromosome name')
    parser.add_argument('--cons_threshold', type=int, help='Consensus threshold')
    args = parser.parse_args()

    # Split the comma-separated VCF paths into a list
    vcf_paths = [path.strip() for path in args.vcfs.split(',') if path.strip()]

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
    conn.commit()

    # Process all VCFs - these will be used for majority rule consensus
    # Create a list of (source_id, file_path) tuples dynamically
    vcf_sources = [(i+1, vcf_path) for i, vcf_path in enumerate(vcf_paths)]
    
    for source_id, file_path in vcf_sources:
        with open(file_path) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                fields = line.strip().split('\t')
                if '.' in fields[9]:
                    continue
                chrom, pos, ref, alt = fields[0], int(fields[1]), fields[3], fields[4]
                gt = sort_gt(fields[9])
                cursor.execute('''
                    INSERT OR REPLACE INTO variants (chrom, pos, ref, alt, gt, source)
                    VALUES (?, ?, ?, ?, ?, ?)
                ''', (chrom, pos, ref, alt, gt, source_id))
        conn.commit()

    # Create a table for majority rule consensus calculation
    # Build the source list dynamically based on the number of VCFs
    source_list = ','.join(map(str, range(1, len(vcf_paths) + 1)))
    
    cursor.execute(f'''
        CREATE TABLE consensus_candidates AS
        SELECT chrom, pos, ref, alt, gt, COUNT(*) AS cnt
        FROM variants
        WHERE source IN ({source_list}) AND gt NOT LIKE '%.%'
        GROUP BY chrom, pos, ref, alt, gt
        HAVING cnt >= {args.cons_threshold}
    ''')

    # Get the majority consensus for each position
    cursor.execute('''
        SELECT chrom, pos, ref, alt, gt, cnt
        FROM (
            SELECT chrom, pos, ref, alt, gt, cnt,
                   ROW_NUMBER() OVER (PARTITION BY chrom, pos ORDER BY cnt DESC, ref, alt, gt) AS rn
            FROM consensus_candidates
        ) t
        WHERE rn = 1
        ORDER BY chrom, pos
    ''')
    consensus_results = cursor.fetchall()

    with open(f"all_chrs/{args.chr}.vcf", 'w') as out:
        out.write('##fileformat=VCFv4.3\n')
        out.write('##FORMAT=<ID=GT,Number=1,Type=String>\n')
        out.write(f'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{args.sample}\n')
        
        for chrom, pos, ref, alt, gt, cnt in consensus_results:
            out.write('\t'.join([chrom, str(pos), '.', ref.upper(), alt, '.', '.', '.', 'GT', f'{gt}\n']))

if __name__ == "__main__":
    main()
