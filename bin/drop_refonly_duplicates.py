#!/usr/bin/env python3
"""Drop reference-only records that duplicate a real variant at the same
position.

Per-sample consensus files contain a gVCF-style reference record (ALT="<*>",
GT 0/0) at every covered position.  When `bcftools merge --gvcf` finds a real
variant at the same position in another sample, it keeps the `<*>` allele in
the merged ALT list (e.g. ``C,T,<*>``); the following ``bcftools norm -m -any``
splits that into two records at the same position: the variant, and a redundant
reference-only record carrying the hom-reference calls.

This filter reads a position-sorted BCF stream (merged and split by ``norm
-m -any``) on stdin and writes it on stdout, keeping the variant records of an
affected position and dropping its reference-only duplicate.  Sample values
that ``bcftools merge --gvcf`` assigned only to the dropped record are copied
to the surviving records first, so every split ALT keeps the genotype of a
hom-reference sample.  Reference-only records at positions without any real
allele are monomorphic sites and are kept: the consensus output contains all
sites, monomorphic ones included.

Used by ``MERGE_OUTPUTS`` between ``bcftools norm -m -any`` and the invariant
filter.  Records are buffered one chromosome position at a time, so the memory
use is constant.
"""

import pysam

# ALT alleles that carry no sequence: the plain "." placeholder and the gVCF
# "any other allele" symbol.
_REFERENCE_ONLY_ALTS = ('.', '<*>')


def is_ref_only(record):
    """True for a record without a real alternate allele."""
    if record.alts is None:
        return True
    return all(alt in _REFERENCE_ONLY_ALTS for alt in record.alts)


def fill_missing_samples(variant, references):
    """Copy sample values from a reference-only record to a variant record.

    ``bcftools merge --gvcf`` can assign a hom-reference sample's ``<*>``
    record to only one of several split ALT records at a position.  The
    sample then has a missing genotype on the other split records, although
    its per-sample consensus called the position hom-reference.  The original
    ``<*>`` record is about to be dropped as a duplicate, so its values are
    the last copy of that sample's reference call.

    Args:
        variant: Variant record that keeps its real ALT.
        references: Reference-only records at the same position.
    """
    for sample_name in variant.header.samples:
        sample = variant.samples[sample_name]
        gt = sample.get('GT')
        if gt is not None and all(allele is not None for allele in gt):
            continue
        for reference in references:
            ref_sample = reference.samples[sample_name]
            ref_gt = ref_sample.get('GT')
            if ref_gt is None or any(allele is None for allele in ref_gt):
                continue
            sample['GT'] = ref_gt
            ref_ad = ref_sample.get('AD')
            if (ref_ad is not None
                    and len(ref_ad) == len(variant.alts) + 1):
                sample['AD'] = ref_ad
            ref_dp = ref_sample.get('DP')
            if ref_dp is not None:
                sample['DP'] = ref_dp
            break


def keep_records(records):
    """Records of one position that should be written downstream.

    Args:
        records: List of VariantRecord copies belonging to one position.

    Returns:
        The variant records when the position carries at least one real allele;
        otherwise the records themselves (a monomorphic position).
    """
    variants = [record for record in records if not is_ref_only(record)]
    if not variants:
        return records
    references = [record for record in records if is_ref_only(record)]
    for variant in variants:
        fill_missing_samples(variant, references)
    return variants


def main():
    try:
        with pysam.VariantFile('/dev/stdin') as inp, \
                pysam.VariantFile(
                    '/dev/stdout', 'wb0', header=inp.header) as out:
            locus = None
            pending = []
            for record in inp:
                record_locus = (record.chrom, record.pos)
                if locus is not None and record_locus != locus:
                    for kept in keep_records(pending):
                        out.write(kept)
                    pending = []
                locus = record_locus
                pending.append(record.copy())
            for kept in keep_records(pending):
                out.write(kept)
    except BrokenPipeError:
        pass


if __name__ == '__main__':
    main()
