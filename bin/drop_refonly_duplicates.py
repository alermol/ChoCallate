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
affected position and dropping its reference-only duplicate.  Reference-only
records at positions without any real allele are monomorphic sites and are
kept: the consensus output contains all sites, monomorphic ones included.

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


def keep_records(records):
    """Records of one position that should be written downstream.

    Args:
        records: List of VariantRecord copies belonging to one position.

    Returns:
        The variant records when the position carries at least one real allele;
        otherwise the records themselves (a monomorphic position).
    """
    variants = [record for record in records if not is_ref_only(record)]
    return variants or records


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
