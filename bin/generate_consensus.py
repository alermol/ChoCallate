#!/usr/bin/env python3

import argparse
import time
from collections import Counter
from contextlib import ExitStack
from typing import List, Optional, Tuple

import pysam


def genotype_at_position(records: List[pysam.VariantRecord], index: int, pos: int, ploidy: int) -> Tuple[int, Tuple[int, int], Tuple[int, int]]:
    i = index
    while i < len(records) and records[i].pos < pos:
        i += 1
    if i < len(records) and records[i].pos == pos:
        gt = records[i].ref.upper(), records[i].alts[0].upper()
        return i, gt, convert_numeric_consensus(sum(records[i].samples[0]['GT']), ploidy)
    return i, (0, 0), convert_numeric_consensus(0, ploidy) # ref and alt are reference


def get_consensus_genotype(genotypes: List[Tuple[str, str]], consensus_threshold: int) -> Optional[Tuple[str, str]]:
    counter: Counter[Tuple[str, str]] = Counter[Tuple[str, str]](genotypes)
    value, count = counter.most_common(1)[0]
    if count < consensus_threshold:
        return None, 0
    return value, count


def generate_mininimal_header(sample_name: str, reference_file: pysam.FastaFile, consensus_threshold: int) -> pysam.VariantHeader:
    header = pysam.VariantHeader()
    header.add_sample(sample_name)
    header.add_line(f'##Consensus threshold={consensus_threshold}')
    for ref in reference_file.references:
        header.add_line(f'##contig=<ID={ref},length={reference_file.get_reference_length(ref)}>')
    header.add_line('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    header.add_line('##INFO=<ID=NM,Number=1,Type=Integer,Description="Number of matching calls">')
    return header


def convert_numeric_consensus(numeric_cons: int, ploidy: int):
    return (0,) * (ploidy - numeric_cons) + (1,) * numeric_cons


def variant_is_snp(variant: Tuple[str, str]) -> bool:
    return len(variant[0]) == 1 and len(variant[1]) == 1


def variant_is_indel(variant: Tuple[str, str]) -> bool:
    return len(variant[0]) > 1 or len(variant[1]) > 1


def main():
    start_time = time.time()
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True, nargs='+', help='Input BCF/VCF files')
    parser.add_argument('--bed', required=True, help='Input BED file (gzipped and tabixed)')
    parser.add_argument('--output', required=True, help='Output BCF file')
    parser.add_argument('--sample_name', required=True, help='Sample name')
    parser.add_argument('--ploidy', required=True, type=int, help='Sample ploidy')
    parser.add_argument('--reference', required=True, help='Reference genome assembly (fasta file and fai index)')
    parser.add_argument('--consensus_threshold', required=True, type=int, help='Number of matching calls for consensus')
    parser.add_argument('--variant_types', required=True, choices=['snp', 'indel', 'both'], help='Variant types to process')
    args = parser.parse_args()

    with pysam.TabixFile(args.bed, index=args.bed + '.csi') as bed_file, \
         pysam.FastaFile(args.reference, filepath_index=args.reference + '.fai') as reference_file:
        with pysam.VariantFile(args.output, mode='w', 
                               header=generate_mininimal_header(args.sample_name, reference_file, args.consensus_threshold)) as output_file:
            with ExitStack[bool | None]() as stack:
                variant_files = [stack.enter_context(pysam.VariantFile(bcf, index_filename=bcf + '.csi')) for bcf in args.input]
                position_progress = 0
                for line in bed_file.fetch():
                    contig, start, end = line.split('\t')
                    vcf_start = int(start)
                    vcf_end = int(end)
                    ref_seq = reference_file.fetch(contig, vcf_start, vcf_end).upper()
                    records_per_file = [list(vf.fetch(contig, vcf_start, vcf_end)) for vf in variant_files]
                    indices = [0] * len(records_per_file)
                    ref_seq_rel_pos = 0
                    while vcf_start < vcf_end:
                        # collect genotypes at current position
                        genotypes = []
                        for fi in range(len(variant_files)):
                            indices[fi], gt, num_gt = genotype_at_position(records_per_file[fi], indices[fi], vcf_start, args.ploidy)
                            genotypes.append((gt, num_gt))
                        genotypes = [i if i != ((0, 0), (0, 0)) else ((ref_seq[ref_seq_rel_pos].upper(), '.'), i[1]) for i in genotypes]

                        # get consensus genotype at current position
                        consensus_genotype, num_matching_calls = get_consensus_genotype(genotypes, args.consensus_threshold)

                        # write consensus genotype to output file
                        if consensus_genotype is not None:
                            if args.variant_types == 'both':
                                new_record = output_file.new_record(
                                    contig=contig,
                                    start=vcf_start - 1,
                                    stop=vcf_start + 1,
                                    alleles=consensus_genotype[0],
                                    qual=42,
                                )
                                new_record.info['NM'] = num_matching_calls
                                new_record.samples[0]['GT'] = consensus_genotype[1]
                                output_file.write(new_record)
                            else:
                                # write consensus genotype to output file if it is of the same type as the requested variant type
                                current_variant_type = 'snp' if variant_is_snp(consensus_genotype[0]) else 'indel'
                                if current_variant_type == args.variant_types:
                                    new_record = output_file.new_record(
                                        contig=contig,
                                        start=vcf_start - 1,
                                        stop=vcf_start + 1,
                                        alleles=consensus_genotype[0],
                                        qual=42,
                                    )
                                    new_record.info['NM'] = num_matching_calls
                                    new_record.samples[0]['GT'] = consensus_genotype[1]
                                    output_file.write(new_record)
                        vcf_start += len(consensus_genotype[0][0]) if consensus_genotype is not None else 1 # shift on referene genotype length
                        ref_seq_rel_pos += len(consensus_genotype[0][0]) if consensus_genotype is not None else 1 # shift on reference genotype length
                        
                        position_progress += len(consensus_genotype[0][0]) if consensus_genotype is not None else 1
                        if position_progress % 1000000 == 0:
                            print(f'Position progress: {position_progress}')

    end_time = time.time()
    print(f'Time taken: {end_time - start_time} seconds')


if __name__ == '__main__':
    main()