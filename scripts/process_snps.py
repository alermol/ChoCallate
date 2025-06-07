#!/usr/bin/env python3

import argparse
from collections import defaultdict, Counter

def sort_gt(s):
    parts = s.split('/')
    numbers = [int(part) for part in parts]
    numbers_sorted = sorted(numbers)
    result = '/'.join(map(str, numbers_sorted))
    return result


def get_most_frequent(items):
    counter = Counter(items)
    most_common = counter.most_common(1)
    return most_common


def parse_vcf(vcf_file):
    variants = defaultdict(dict)
    with open(vcf_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if '.' in fields[9]:
                continue
            else:
                chrom = fields[0]
                pos = int(fields[1])
                ref = fields[3]
                alt = fields[4]
                gt = sort_gt(fields[9])
                variants[(chrom, pos)] = [(ref, alt, gt)]
    return variants

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--vcf1', required=True, help='Path to VCF file 1 (bcftools)')
    parser.add_argument('--vcf2', required=True, help='Path to VCF file 2')
    parser.add_argument('--vcf3', required=True, help='Path to VCF file 3')
    parser.add_argument('--vcf4', required=True, help='Path to VCF file 4')
    parser.add_argument('--vcf5', required=True, help='Path to VCF file 5')
    parser.add_argument('--sample', required=True, help='Sample name')
    args = parser.parse_args()

    base_vars = parse_vcf(args.vcf1)

    for file in [args.vcf2, args.vcf3, args.vcf4, args.vcf5]:
        polymorphic_vcf = parse_vcf(file)
        for coord, gen in base_vars.items():
            if coord in polymorphic_vcf.keys():
                base_vars[coord].append(polymorphic_vcf[coord][0])
            else:
                base_vars[coord].append(base_vars[coord][0])

    with open(f"{args.sample}.vcf", 'w') as out:
        out.write('##fileformat=VCFv4.3\n')
        out.write('##FORMAT=<ID=GT,Number=1,Type=String>\n')
        out.write(f'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{args.sample}\n')
        for coord, gen in base_vars.items():
            most_freq_var = get_most_frequent(gen)
            if ((most_freq_var[0][1] < 3) | 
                ('.' in most_freq_var[0][0][2]) | 
                (not most_freq_var[0][0][0].isupper())):
                continue
            else:
                chrom = coord[0]
                pos = coord[1]
                ref = most_freq_var[0][0][0]
                alt = most_freq_var[0][0][1]
                gt = most_freq_var[0][0][2]
                out.write('\t'.join([chrom, str(pos), '.', ref, alt, '.', '.', '.', 'GT', f'{gt}\n']))
