#!/usr/bin/env python3

import argparse
import time
from collections import Counter
from contextlib import ExitStack
from typing import NamedTuple

import pysam


class Genotypes(NamedTuple):
    index: int = 0
    ref: str = ""
    filter: str = ""
    alt: tuple[str, ...] = ("",)
    gt: tuple[int, ...] = (0, 0)


def genotype_at_position(records, index, pos):
    """
    Get the genotype at a given position in a list of VCF records.

    Args:
        records: List of pysam.VariantRecord objects.
        index: Integer index of the record to get the genotype from.
        pos: Integer position to get the genotype from.

    Returns:
        A Genotypes named tuple containing the index, reference allele, alternate alleles, and numeric representation of the genotype.
    """
    i = index
    while i < len(records) and records[i].pos < pos:
        i += 1
    if (i < len(records)) and (records[i].pos == pos):
        return Genotypes(index=i, 
                         ref=records[i].ref.upper(), 
                         filter=list(records[i].filter)[0],
                         alt=tuple(sorted(records[i].alts)), 
                         gt=tuple(sorted(records[i].samples[0]['GT'])))
    return Genotypes()


def get_consensus_genotype(genotypes, consensus_threshold):
    """
    Get the consensus genotype from a list of genotypes.

    Args:
        genotypes: List of Genotypes named tuples. Each tuple contains the index, reference allele, 
        alternate alleles, and numeric representation of the genotype.
        consensus_threshold: Integer threshold for the consensus genotype.

    Returns:
        A tuple containing the consensus genotype (reference allele, alternate alleles, numeric representation of the genotype) and the number of matching calls.
        If the consensus genotype is not found, returns None and 0.
    """
    genotypes = [(i.ref, i.alt, i.gt) for i in genotypes if i.filter != "LowQual"]
    counter = Counter(genotypes)
    value, count = counter.most_common(1)[0]
    if count < consensus_threshold:
        return None, 0
    return value, count


def generate_mininimal_header(sample_name, reference_file, consensus_threshold):
    """
    Generate a minimal VCF header for the consensus file.

    Args:
        sample_name: String sample name.
        reference_file: pysam.FastaFile object.
        consensus_threshold: Integer threshold for the consensus genotype.

    Returns:
        A pysam.VariantHeader object.
    """
    header = pysam.VariantHeader()
    header.add_sample(sample_name)
    header.add_line(f'##Consensus threshold={consensus_threshold}')
    for ref in reference_file.references:
        header.add_line(f'##contig=<ID={ref},length={reference_file.get_reference_length(ref)}>')
    header.add_line('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    header.add_line('##FORMAT=<ID=NM,Number=1,Type=Integer,Description="Number of matching calls">')
    return header



def main():
    start_time = time.time()
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True, nargs='+', help='Input BCF/VCF files')
    parser.add_argument('--bed', required=True, help='Input BED file (gzipped and tabixed)')
    parser.add_argument('--output', required=True, help='Output BCF file')
    parser.add_argument('--sample_name', required=True, help='Sample name')
    parser.add_argument('--reference', required=True, help='Reference genome assembly (fasta file and fai index)')
    parser.add_argument('--consensus_threshold', required=True, type=int, help='Number of matching calls for consensus')
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
                    vcf_start = int(start) + 1
                    vcf_end = int(end)
                    ref_seq = reference_file.fetch(contig, int(start), int(end)).upper()
                    records_per_file = [list(vf.fetch(contig, int(start), int(end))) for vf in variant_files]
                    indices = [0] * len(records_per_file)
                    ref_seq_rel_pos = 0
                    while vcf_start < vcf_end:
                        # collect genotypes at current position
                        genotypes = []
                        for fi in range(len(variant_files)):
                            res = genotype_at_position(records_per_file[fi], indices[fi], vcf_start)
                            indices[fi] = res.index
                            genotypes.append(res)
                        genotypes = [i if i.gt != (0, 0) else i._replace(ref=ref_seq[ref_seq_rel_pos].upper(), alt=('.',)) for i in genotypes]

                        # get consensus genotype at current position
                        consensus_genotype, num_matching_calls = get_consensus_genotype(genotypes, args.consensus_threshold)

                        # write consensus genotype to output file
                        if consensus_genotype is not None:
                            new_record = output_file.new_record(
                                contig=contig,
                                start=vcf_start - 1,
                                stop=vcf_start,
                                alleles=[consensus_genotype[0]] + list(consensus_genotype[1]),
                                qual=42,
                            )
                            new_record.samples[0]['NM'] = num_matching_calls
                            new_record.ref = consensus_genotype[0]
                            new_record.samples[0]['GT'] = consensus_genotype[2]
                            output_file.write(new_record)
                        shift = len(consensus_genotype[0]) if consensus_genotype is not None else 1
                        vcf_start += shift # shift on referene genotype length
                        ref_seq_rel_pos += shift # shift on reference genotype length
                        
                        position_progress += shift
                        if position_progress % 1000000 == 0:
                            print(f'Position progress: {position_progress}')

    end_time = time.time()
    print(f'Processing of {args.bed} took {end_time - start_time} seconds')


if __name__ == '__main__':
    main()