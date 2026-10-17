#!/usr/bin/env python3
"""Drop ALT="." records that duplicate a real variant at the same position.

Per-sample consensus files contain a hom-reference record (ALT=".") at every
covered position.  When `bcftools merge` finds a real variant at the same
position in another sample, it keeps the ALT="." placeholder as a literal
allele (e.g. alleles ``C,T,.``), and the following ``bcftools norm -m -any``
splits that into two records at the same position: the variant, and a redundant
ALT="." record carrying the hom-reference calls.

This filter reads a position-sorted BCF stream (merged and split by ``norm
-m -any``) on stdin and writes it on stdout, keeping the variant records of an
affected position and dropping its ALT="." duplicate.  ALT="." records at
positions without any real allele are monomorphic sites and are kept: the
consensus output contains all sites, monomorphic ones included.

Used by ``MERGE_OUTPUTS`` between ``bcftools norm -m -any`` and the invariant
filter.  Records are buffered one position at a time, so the memory use is
constant.
"""

import pysam


def is_ref_only(record):
    """True for a record without alternate alleles (ALT is ``.``)."""
    return record.alts is None or tuple(record.alts) == ('.',)


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
                pysam.VariantFile('/dev/stdout', 'wb0', header=inp.header) as out:
            position = None
            pending = []
            for record in inp:
                if position is not None and record.pos != position:
                    for kept in keep_records(pending):
                        out.write(kept)
                    pending = []
                position = record.pos
                pending.append(record.copy())
            for kept in keep_records(pending):
                out.write(kept)
    except BrokenPipeError:
        pass


if __name__ == '__main__':
    main()
