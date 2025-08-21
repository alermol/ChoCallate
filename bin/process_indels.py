#!/usr/bin/env python3

import argparse
import sqlite3
import os
import sys
from typing import List, Tuple, Optional, Any
import pysam


def gt_tuple_to_sorted_str(gt: Optional[Tuple[Optional[int], ...]]) -> Optional[str]:
    if gt is None:
        return None
    if any(allele is None for allele in gt):
        return None
    sorted_gt: List[int] = sorted(int(a) for a in gt)
    return '/'.join(map(str, sorted_gt))


def gt_str_to_tuple(gt_str: str) -> Tuple[int, ...]:
    return tuple(int(x) for x in gt_str.replace('|', '/').split('/'))


def main() -> None:
    parser: argparse.ArgumentParser = argparse.ArgumentParser()
    parser.add_argument('--zero_bcf', required=True, help='Path to zero BCF')
    parser.add_argument('--bcfs', required=True, help='Comma-separated paths to BCF files')
    parser.add_argument('--sample', required=True, help='Sample name')
    parser.add_argument('--chr', required=True, help='Chromosome name')
    parser.add_argument('--cons_threshold', type=int, help='Consensus threshold')
    args: argparse.Namespace = parser.parse_args()

    # Split the comma-separated BCF paths into a list
    bcf_paths: List[str] = [path.strip() for path in args.bcfs.split(',') if path.strip()]

    conn: sqlite3.Connection = sqlite3.connect(':memory:')
    cursor: sqlite3.Cursor = conn.cursor()
    
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

    # Process all BCFs - these will be used for majority rule consensus
    # Create a list of (source_id, file_path) tuples dynamically
    bcf_sources: List[Tuple[int, str]] = [(i+1, bcf_path) for i, bcf_path in enumerate(bcf_paths)]
    
    total_variants = 0
    
    for source_id, file_path in bcf_sources:
        variants_count = 0
        with pysam.VariantFile(file_path, 'rb') as bcf_in:
            samples: List[str] = list(bcf_in.header.samples)
            sample_name: Optional[str] = samples[0] if len(samples) > 0 else None
            for rec in bcf_in:
                if sample_name is None:
                    continue
                chrom: str = rec.chrom
                pos: int = rec.pos
                ref: str = rec.ref
                alt: str = rec.alts[0] if rec.alts and len(rec.alts) > 0 else '.'
                gt_tuple: Optional[Tuple[Optional[int], ...]] = rec.samples[sample_name].get('GT', None)
                gt_str_opt: Optional[str] = gt_tuple_to_sorted_str(gt_tuple)
                if gt_str_opt is None:
                    continue
                cursor.execute('''
                    INSERT OR REPLACE INTO variants (chrom, pos, ref, alt, gt, source)
                    VALUES (?, ?, ?, ?, ?, ?)
                ''', (chrom, pos, ref, alt, gt_str_opt, source_id))
                variants_count += 1
        conn.commit()
        total_variants += variants_count
    
    # Create a table for majority rule consensus calculation
    # Build the source list dynamically based on the number of BCFs
    source_list: str = ','.join(map(str, range(1, len(bcf_paths) + 1)))
    
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
    consensus_results: List[Tuple[str, int, str, str, str, int]] = cursor.fetchall()

    # Generate final output
    # Create header and output VCF to stdout
    print('##fileformat=VCFv4.3')
    
    # Add all chromosome names to the header
    for chrom_name in args.chr.strip().split('\n'):
        print(f'##contig=<ID={chrom_name}>')
    
    print('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    print(f'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{args.sample}')

    # Output VCF records to stdout
    for chrom, pos, ref_val, alt_val, gt_val, cnt in consensus_results:
        # Format VCF line
        alt_field = alt_val if alt_val != '.' else '.'
        print(f'{chrom}\t{pos}\t.\t{ref_val.upper()}\t{alt_field}\t100\t.\t.\tGT\t{gt_val}')


if __name__ == "__main__":
    main()
