#!/usr/bin/env python3

import argparse
import sqlite3
import os
import sys
from typing import List, Tuple, Set, Dict, Optional, Any
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
    
    cursor.execute('''
        CREATE TABLE zero_bcf_data (
            chrom TEXT NOT NULL,
            pos INTEGER NOT NULL,
            ref TEXT NOT NULL,
            alt TEXT NOT NULL,
            gt TEXT NOT NULL,
            PRIMARY KEY (chrom, pos)
        )
    ''')
    conn.commit()

    # Store zero BCF data separately
    zero_bcf_positions: Set[Tuple[str, int]] = set()
    zero_bcf_count = 0
    
    with pysam.VariantFile(args.zero_bcf, 'rb') as zero_bcf_in:
        zero_samples: List[str] = list(zero_bcf_in.header.samples)
        zero_sample_name: Optional[str] = zero_samples[0] if len(zero_samples) > 0 else None
        for rec in zero_bcf_in:
            chrom: str = rec.chrom
            pos: int = rec.pos
            ref: str = rec.ref
            alt: str = rec.alts[0] if rec.alts and len(rec.alts) > 0 else '.'
            gt_tuple: Optional[Tuple[Optional[int], ...]] = None
            if zero_sample_name is not None:
                gt_tuple = rec.samples[zero_sample_name].get('GT', None)
            gt_str_opt: Optional[str] = gt_tuple_to_sorted_str(gt_tuple)
            if gt_str_opt is None:
                continue
            zero_bcf_positions.add((chrom, pos))
            cursor.execute('''
                INSERT OR REPLACE INTO zero_bcf_data (chrom, pos, ref, alt, gt)
                VALUES (?, ?, ?, ?, ?)
            ''', (chrom, pos, ref, alt, gt_str_opt))
            zero_bcf_count += 1
    conn.commit()

    # Process other BCFs
    # Create a list of (source_id, file_path) tuples dynamically
    other_bcf_sources: List[Tuple[int, str]] = [(i+1, bcf_path) for i, bcf_path in enumerate(bcf_paths)]
    
    total_variants = 0
    for source_id, file_path in other_bcf_sources:
        variants_count = 0
        with pysam.VariantFile(file_path, 'rb') as bcf_in:
            samples: List[str] = list(bcf_in.header.samples)
            sample_name: Optional[str] = samples[0] if len(samples) > 0 else None
            for rec in bcf_in:
                if sample_name is None:
                    continue
                chrom: str = rec.chrom
                pos: int = rec.pos
                if (chrom, pos) not in zero_bcf_positions:
                    continue
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
    
    # Create a table for majority rule consensus calculation (excluding zero BCF)
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
    ''')
    consensus_results: List[Tuple[str, int, str, str, str, int]] = cursor.fetchall()
    
    # Create a dictionary for quick lookup of consensus results
    consensus_dict: Dict[Tuple[str, int], Tuple[str, str, str]] = {(chrom, pos): (ref, alt, gt) for chrom, pos, ref, alt, gt, cnt in consensus_results}

    # Prepare final output: all positions from zero BCF must be present
    cursor.execute('''
        SELECT chrom, pos, ref, alt, gt
        FROM zero_bcf_data
        ORDER BY chrom, pos
    ''')
    zero_bcf_results: List[Tuple[str, int, str, str, str]] = cursor.fetchall()

    consensus_used = 0
    zero_bcf_used = 0
    
    # Create header and output VCF to stdout
    print('##fileformat=VCFv4.3')
    
    for chrom_name in args.chr.split('\n'):
        print(f'##contig=<ID={chrom_name}>')
    
    print('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    print(f'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{args.sample}')

    # Output VCF records to stdout
    for chrom, pos, ref_val, alt_val, gt_val in zero_bcf_results:
        if (chrom, pos) in consensus_dict:
            consensus_ref: str
            consensus_alt: str
            consensus_gt: str
            consensus_ref, consensus_alt, consensus_gt = consensus_dict[(chrom, pos)]
            # Format VCF line
            alt_field = consensus_alt if consensus_alt != '.' else '.'
            print(f'{chrom}\t{pos}\t.\t{consensus_ref.upper()}\t{alt_field}\t100\t.\t.\tGT\t{consensus_gt}')
            consensus_used += 1
        else:
            # Format VCF line
            alt_field = alt_val if alt_val != '.' else '.'
            print(f'{chrom}\t{pos}\t.\t{ref_val.upper()}\t{alt_field}\t100\t.\t.\tGT\t{gt_val}')
            zero_bcf_used += 1


if __name__ == "__main__":
    main()
