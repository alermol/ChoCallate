#!/usr/bin/env python3

import pysam
import argparse
from contextlib import ExitStack
from collections import Counter
from typing import List, Tuple, Optional


def get_genotype(variant_file: pysam.VariantFile, 
                 contig: str, 
                 start: int, 
                 stop: int) -> int:
    result: List[pysam.VariantRecord] = list(variant_file.fetch(contig, start, stop))
    if len(result) == 0:
        return 0
    return int(sum(result[0].samples[0]['GT']))


def get_alleles(variant_file: pysam.VariantFile, 
                reference: pysam.FastaFile, 
                contig: str, 
                start: int, 
                stop: int) -> Tuple[str]:
    result: List[pysam.VariantRecord] = list(variant_file.fetch(contig, start, stop))
    if len(result) == 0:
        return reference.fetch(contig, start, stop).upper(), '.'
    return result[0].ref, result[0].alts[0]


def get_consensus_genotype(genotypes: List[int], consensus_threshold: int) -> Optional[int]:
    counter: Counter[int] = Counter[int](genotypes)
    value, count = counter.most_common(1)[0]
    if count < consensus_threshold:
        return None
    return value


def generate_mininimal_header(sample_name: str, reference_file: pysam.FastaFile) -> pysam.VariantHeader:
    header = pysam.VariantHeader()
    header.add_sample(sample_name)
    for ref in reference_file.references:
        header.add_line(f'##contig=<ID={ref},length={reference_file.get_reference_length(ref)}>')
    header.add_line('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    return header


def convert_numeric_consensus(numeric_cons: int, ploidy: int):
    return (0,) * (ploidy - numeric_cons) + (1,) * numeric_cons


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True, nargs='+', help='Input BCF/VCF files')
    parser.add_argument('--bed', required=True, help='Input BED file')
    parser.add_argument('--output', required=True, help='Output BCF file')
    parser.add_argument('--sample_name', required=True, help='Sample name')
    parser.add_argument('--ploidy', required=True, type=int, help='Sample ploidy')
    parser.add_argument('--reference', required=True, help='Reference genome assembly')
    parser.add_argument('--consensus_threshold', required=True, type=int, help='Number of matching calls for consensus')
    args = parser.parse_args()

    with pysam.TabixFile(args.bed, index=args.bed + '.csi') as bed_file, \
         pysam.FastaFile(args.reference, 
                         filepath_index=args.reference + '.fai', 
                         filepath_index_compressed=args.reference + '.gzi') as reference_file:
        with pysam.VariantFile(args.output, mode='wb', 
                               header=generate_mininimal_header(args.sample_name, reference_file)) as output_file:
            with ExitStack[bool | None]() as stack:
                variant_files = [stack.enter_context(pysam.VariantFile(bcf, index_filename=bcf + '.csi')) 
                                 for bcf in args.input]
                for line in bed_file.fetch():
                    contig, start, end = line.split('\t')
                    for pos in range(int(start) + 1, int(end) + 1):
                        genotypes = [get_genotype(vf, contig, pos, pos + 1, args.ploidy) for vf in variant_files]
                        consensus_genotype = get_consensus_genotype(genotypes, args.consensus_threshold)
                        if consensus_genotype is None:
                            continue
                        new_record = output_file.new_record(
                            contig=contig,
                            start=pos - 1,
                            stop=pos,
                            alleles=get_alleles(variant_files[0], reference_file, contig, pos, pos + 1),
                            qual=42,
                        )
                        new_record.samples[args.sample_name]['GT'] = convert_numeric_consensus(consensus_genotype, args.ploidy)
                        output_file.write(new_record)


if __name__ == '__main__':
    main()