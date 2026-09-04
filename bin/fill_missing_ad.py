#!/usr/bin/env python3
"""Replace missing FORMAT/AD allele depths with 0 for samples with a called GT.

After multi-sample merge, hom-ref per-sample records (ALT=.) keep REF AD but lose the
placeholder ALT slot, leaving values like 0/0:4,.:4. Missing AD entries for called
samples are filled with 0. Samples with missing GT are left unchanged.
"""

import pysam


def fill_record_ad(record):
    for sample in record.samples.values():
        gt = sample.get('GT')
        if gt is None or all(allele is None for allele in gt):
            continue
        ad = sample.get('AD')
        if ad is None:
            continue
        if any(depth is None for depth in ad):
            sample['AD'] = tuple(0 if depth is None else depth for depth in ad)


def main():
    try:
        with pysam.VariantFile('/dev/stdin') as inp:
            with pysam.VariantFile('/dev/stdout', 'wb0', header=inp.header) as out:
                for record in inp:
                    fill_record_ad(record)
                    out.write(record)
    except BrokenPipeError:
        pass


if __name__ == '__main__':
    main()
