#!/usr/bin/env python3

import argparse
import sqlite3
import os
import sys
from typing import List, Tuple, Set, Dict, Optional, Any

# Add the current directory to the path to import logging utilities
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from logging_utils import setup_logger


def sort_gt(s: str) -> str:
    parts: List[str] = s.split('/')
    numbers: List[int] = [int(part) for part in parts]
    numbers_sorted: List[int] = sorted(numbers)
    return '/'.join(map(str, numbers_sorted))


def main() -> None:
    parser: argparse.ArgumentParser = argparse.ArgumentParser()
    parser.add_argument('--zero_vcf', required=True, help='Path to zero VCF')
    parser.add_argument('--vcfs', required=True, help='Comma-separated paths to VCF files')
    parser.add_argument('--sample', required=True, help='Sample name')
    parser.add_argument('--chr', required=True, help='Chromosome name')
    parser.add_argument('--cons_threshold', type=int, help='Consensus threshold')
    args: argparse.Namespace = parser.parse_args()

    # Initialize logger
    logger = setup_logger("process_snps", args.sample, "INFO")
    logger.log_process_start("SNPs consensus generation", {
        "zero_vcf": args.zero_vcf,
        "vcfs_count": len(args.vcfs.split(',')),
        "sample": args.sample,
        "chromosome": args.chr,
        "consensus_threshold": args.cons_threshold
    })

    # Split the comma-separated VCF paths into a list
    vcf_paths: List[str] = [path.strip() for path in args.vcfs.split(',') if path.strip()]
    logger.info(f"Processing {len(vcf_paths)} VCF files")

    conn: sqlite3.Connection = sqlite3.connect(':memory:')
    cursor: sqlite3.Cursor = conn.cursor()
    
    logger.info("Creating database tables")
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
    logger.info("Database tables created successfully")

    # Store zero VCF data separately
    logger.info("Processing zero VCF data")
    zero_vcf_positions: Set[Tuple[str, int]] = set()
    zero_vcf_count = 0
    
    with open(args.zero_vcf) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields: List[str] = line.strip().split('\t')
            if '.' in fields[9]:
                continue
            chrom: str = fields[0]
            pos: int = int(fields[1])
            ref: str = fields[3]
            alt: str = fields[4]
            gt: str = sort_gt(fields[9])
            zero_vcf_positions.add((chrom, pos))
            cursor.execute('''
                INSERT OR REPLACE INTO zero_vcf_data (chrom, pos, ref, alt, gt)
                VALUES (?, ?, ?, ?, ?)
            ''', (chrom, pos, ref, alt, gt))
            zero_vcf_count += 1
    conn.commit()
    logger.info(f"Zero VCF processing completed: {zero_vcf_count} variants loaded")

    # Process other VCFs
    logger.info("Processing variant caller VCFs")
    # Create a list of (source_id, file_path) tuples dynamically
    other_vcfs: List[Tuple[int, str]] = [(i+1, vcf_path) for i, vcf_path in enumerate(vcf_paths)]
    
    total_variants = 0
    for source_id, file_path in other_vcfs:
        logger.info(f"Processing VCF {source_id}: {os.path.basename(file_path)}")
        variants_count = 0
        with open(file_path) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                fields: List[str] = line.strip().split('\t')
                if '.' in fields[9]:
                    continue
                chrom: str = fields[0]
                pos: int = int(fields[1])
                ref: str = fields[3]
                alt: str = fields[4]
                if (chrom, pos) not in zero_vcf_positions:
                    continue
                gt: str = sort_gt(fields[9])
                cursor.execute('''
                    INSERT OR REPLACE INTO variants (chrom, pos, ref, alt, gt, source)
                    VALUES (?, ?, ?, ?, ?, ?)
                ''', (chrom, pos, ref, alt, gt, source_id))
                variants_count += 1
        conn.commit()
        total_variants += variants_count
        logger.info(f"VCF {source_id} processed: {variants_count} variants loaded")
    
    logger.info(f"All VCFs processed: {total_variants} total variants loaded")

    # Create a table for majority rule consensus calculation (excluding zero VCF)
    logger.info("Calculating consensus variants")
    # Build the source list dynamically based on the number of VCFs
    source_list: str = ','.join(map(str, range(1, len(vcf_paths) + 1)))
    
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
    
    logger.info(f"Consensus calculation completed: {len(consensus_results)} consensus variants found")

    # Prepare final output: all positions from zero VCF must be present
    logger.info("Generating final consensus VCF")
    cursor.execute('''
        SELECT chrom, pos, ref, alt, gt
        FROM zero_vcf_data
        ORDER BY chrom, pos
    ''')
    zero_vcf_results: List[Tuple[str, int, str, str, str]] = cursor.fetchall()

    # Ensure output directory exists
    os.makedirs("all_chrs", exist_ok=True)
    output_file = f"all_chrs/{args.chr}.vcf"
    
    consensus_used = 0
    zero_vcf_used = 0
    
    with open(output_file, 'w') as out:
        out.write('##fileformat=VCFv4.3\n')
        out.write('##FORMAT=<ID=GT,Number=1,Type=String>\n')
        out.write(f'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{args.sample}\n')
        
        for chrom, pos, ref, alt, gt in zero_vcf_results:
            # Check if we have a consensus for this position from other VCFs
            if (chrom, pos) in consensus_dict:
                # Use consensus from other VCFs (majority rule, excluding zero VCF)
                consensus_ref: str
                consensus_alt: str
                consensus_gt: str
                consensus_ref, consensus_alt, consensus_gt = consensus_dict[(chrom, pos)]
                out.write('\t'.join([chrom, str(pos), '.', consensus_ref.upper(), consensus_alt, '.', '.', '.', 'GT', f'{consensus_gt}\n']))
                consensus_used += 1
            else:
                # Position only exists in zero VCF, use zero VCF alleles
                out.write('\t'.join([chrom, str(pos), '.', ref.upper(), alt, '.', '.', '.', 'GT', f'{gt}\n']))
                zero_vcf_used += 1
    
    logger.info(f"Final VCF generated: {output_file}")
    logger.info(f"Variants breakdown: {consensus_used} consensus, {zero_vcf_used} zero VCF")
    
    # Log performance metrics
    logger.log_performance({
        "zero_vcf_variants": zero_vcf_count,
        "total_variants_loaded": total_variants,
        "consensus_variants": len(consensus_results),
        "final_variants": len(zero_vcf_results),
        "consensus_used": consensus_used,
        "zero_vcf_used": zero_vcf_used
    })
    
    logger.log_process_complete("SNPs consensus generation", {
        "output_file": output_file,
        "total_variants": len(zero_vcf_results)
    })

if __name__ == "__main__":
    main()
