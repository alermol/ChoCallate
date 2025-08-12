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
    parser.add_argument('--zero_vcf', required=True, help='Path to zero VCF')
    parser.add_argument('--vcf1', required=True, help='Path to VCF file 1 (gatk)')
    parser.add_argument('--vcf2', required=True, help='Path to VCF file 2')
    parser.add_argument('--vcf3', required=True, help='Path to VCF file 3')
    parser.add_argument('--sample', required=True, help='Sample name')
    parser.add_argument('--chr', required=True, help='Chromosome name')
    parser.add_argument('--cons_threshold', type=int, default=3, help='Consensus threshold')
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
        CREATE TABLE zero_vcf_data (
            chrom TEXT NOT NULL,
            pos INTEGER NOT NULL,
            ref TEXT NOT NULL,
            alt TEXT NOT NULL,
            gt TEXT NOT NULL,
            PRIMARY KEY (chrom, pos)
        )
    ''')
    conn.commit()

    # Store zero VCF data separately
    zero_vcf_positions = set()
    
    with open(args.zero_vcf) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if '.' in fields[9]:
                continue
            chrom, pos, ref, alt = fields[0], int(fields[1]), fields[3], fields[4]
            gt = sort_gt(fields[9])
            zero_vcf_positions.add((chrom, pos))
            cursor.execute('''
                INSERT OR REPLACE INTO zero_vcf_data (chrom, pos, ref, alt, gt)
                VALUES (?, ?, ?, ?, ?)
            ''', (chrom, pos, ref, alt, gt))
    conn.commit()

    # Process other VCFs (VCF1-VCF3) - these will be used for majority rule consensus
    other_vcfs = [
        (1, args.vcf1),
        (2, args.vcf2),
        (3, args.vcf3)
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
                if (chrom, pos) not in zero_vcf_positions:
                    continue
                gt = sort_gt(fields[9])
                cursor.execute('''
                    INSERT OR REPLACE INTO variants (chrom, pos, ref, alt, gt, source)
                    VALUES (?, ?, ?, ?, ?, ?)
                ''', (chrom, pos, ref, alt, gt, source_id))
        conn.commit()

    # Create a table for majority rule consensus calculation (excluding zero VCF)
    cursor.execute(f'''
        CREATE TABLE consensus_candidates AS
        SELECT chrom, pos, ref, alt, gt, COUNT(*) AS cnt
        FROM variants
        WHERE source IN (1, 2, 3) AND gt NOT LIKE '%.%'
        GROUP BY chrom, pos, ref, alt, gt
        HAVING cnt >= {args.cons_threshold}
    ''')

    # Get the majority consensus for each position - select the most common allele combination per position
    cursor.execute('''
        SELECT chrom, pos, ref, alt, gt, cnt
        FROM (
            SELECT chrom, pos, ref, alt, gt, cnt,
                   ROW_NUMBER() OVER (PARTITION BY chrom, pos ORDER BY cnt DESC, ref, alt, gt) AS rn
            FROM consensus_candidates
        ) t
        WHERE rn = 1
    ''')
    consensus_results = cursor.fetchall()
    
    # Create a dictionary for quick lookup of consensus results
    consensus_dict = {(chrom, pos): (ref, alt, gt) for chrom, pos, ref, alt, gt, cnt in consensus_results}

    # Prepare final output: all positions from zero VCF must be present
    cursor.execute('''
        SELECT chrom, pos, ref, alt, gt
        FROM zero_vcf_data
        ORDER BY chrom, pos
    ''')
    zero_vcf_results = cursor.fetchall()

    with open(f"all_chrs/{args.chr}.vcf", 'w') as out:
        out.write('##fileformat=VCFv4.3\n')
        out.write('##FORMAT=<ID=GT,Number=1,Type=String>\n')
        out.write(f'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{args.sample}\n')
        
        for chrom, pos, ref, alt, gt in zero_vcf_results:
            # Check if we have a consensus for this position from other VCFs
            if (chrom, pos) in consensus_dict:
                # Use consensus from other VCFs (majority rule, excluding zero VCF)
                consensus_ref, consensus_alt, consensus_gt = consensus_dict[(chrom, pos)]
                out.write('\t'.join([chrom, str(pos), '.', consensus_ref.upper(), consensus_alt, '.', '.', '.', 'GT', f'{consensus_gt}\n']))
            else:
                # Position only exists in zero VCF, use zero VCF alleles
                out.write('\t'.join([chrom, str(pos), '.', ref.upper(), alt, '.', '.', '.', 'GT', f'{gt}\n']))

if __name__ == "__main__":
    main()
