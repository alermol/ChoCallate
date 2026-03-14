#!/usr/bin/env python3

import pysam
import argparse
from contextlib import ExitStack
from collections import Counter
from typing import List, Tuple, Optional


def genotype_at_position(records: List, index: int, pos: int) -> Tuple[int, int]:
    i = index
    while i < len(records) and records[i].pos < pos:
        i += 1
    if i < len(records) and records[i].pos == pos:
        gt = int(sum(records[i].samples[0]['GT']))
        return i, gt
    return i, 0


def alleles_at_position(records: List[pysam.VariantRecord], index: int, pos: int, ref_seq: str, start: int) -> Tuple[str, str]:
    if index < len(records) and records[index].pos == pos:
        return records[index].ref.upper(), records[index].alts[0].upper()
    return ref_seq[pos - start - 1].upper(), '.'


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
    parser.add_argument('--type', choices=['snp', 'indel'], required=True, help='Type of mutation')
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
                    start = int(start)
                    end = int(end)
                    ref_seq = reference_file.fetch(contig, start, end)
                    records_per_file = [list(vf.fetch(contig, start, end)) for vf in variant_files]
                    indices = [0] * len(records_per_file)
                    for pos in range(start + 1, end + 1):
                        genotypes = []
                        for fi in range(len(variant_files)):
                            indices[fi], gt = genotype_at_position(records_per_file[fi], indices[fi], pos)
                            genotypes.append(gt)
                        if args.type == 'indel' and all([g == 0 for g in genotypes]):
                            continue
                        consensus_genotype = get_consensus_genotype(genotypes, args.consensus_threshold)
                        if args.type == 'indel' and consensus_genotype == 0:
                            continue
                        if consensus_genotype is None:
                            continue
                        new_record = output_file.new_record(
                            contig=contig,
                            start=pos - 1,
                            stop=pos,
                            alleles=alleles_at_position(records_per_file[0], indices[0], pos, ref_seq, start),
                            qual=42,
                        )
                        if args.type == 'indel' and new_record.alts[0] == '.':
                            continue
                        new_record.samples[args.sample_name]['GT'] = convert_numeric_consensus(consensus_genotype, args.ploidy)
                        output_file.write(new_record)



if __name__ == '__main__':
    main()