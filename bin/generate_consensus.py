#!/usr/bin/env python3

import argparse
import time
from collections import Counter, namedtuple
from contextlib import ExitStack
from typing import NamedTuple

import pysam


PileupSite = namedtuple('PileupSite', ['snp', 'ins', 'dele', 'dp'])


class Genotypes(NamedTuple):
    index: int = 0
    ref: str = ""
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
        The consensus genotype (reference allele, alternate alleles, numeric representation of the genotype).
        If the consensus genotype is not found, returns None.
    """
    genotypes = [(i.ref, i.alt, i.gt) for i in genotypes]
    counter = Counter(genotypes)
    value, count = counter.most_common(1)[0]
    if count < consensus_threshold:
        return None
    return value


def collect_region_pileup(bam_file, contig, start0, end0):
    """
    Build per-position pileup summaries for a region.

    Args:
        bam_file: pysam.AlignmentFile object.
        contig: Contig/chromosome name.
        start0: 0-based inclusive start.
        end0: 0-based exclusive end.

    Returns:
        Dict mapping 1-based positions to PileupSite(snp Counter, insertion Counter,
        deletion-length Counter, depth).
    """
    sites = {}
    for column in bam_file.pileup(
        contig,
        start0,
        end0,
        truncate=True,
        stepper='all',
        min_base_quality=0,
        max_depth=1000000,
    ):
        if column.reference_pos < start0 or column.reference_pos >= end0:
            continue
        # Deduplicate overlapping mates for allele counts (one obs per query name).
        # DP counts every covering read so it matches samtools depth used for the BED.
        best = {}
        dp = 0
        for read in column.pileups:
            if read.is_refskip:
                continue
            dp += 1
            aln = read.alignment
            if read.is_del or read.query_position is None:
                continue
            bq = aln.query_qualities[read.query_position] if aln.query_qualities else 0
            if read.indel > 0:
                allele = (
                    'ins',
                    aln.query_sequence[
                        read.query_position + 1 : read.query_position + 1 + read.indel
                    ].upper(),
                )
            elif read.indel < 0:
                allele = ('del', -read.indel)
            else:
                allele = ('snp', aln.query_sequence[read.query_position].upper())
            qname = aln.query_name
            prev = best.get(qname)
            if prev is None or bq > prev[0]:
                best[qname] = (bq, allele)
        snp = Counter()
        ins = Counter()
        dele = Counter()
        for _, allele in best.values():
            kind, value = allele
            if kind == 'snp':
                snp[value] += 1
            elif kind == 'ins':
                ins[value] += 1
            else:
                dele[value] += 1
        sites[column.reference_pos + 1] = PileupSite(
            snp=snp, ins=ins, dele=dele, dp=dp
        )
    return sites


def count_allele_depth(site, ref, allele):
    """
    Count reads supporting a single allele at a pileup site.

    Args:
        site: PileupSite or None.
        ref: Reference allele string.
        allele: Allele string to count ('.' is a non-allele placeholder for hom-ref).

    Returns:
        Integer read count supporting the allele.
    """
    if site is None:
        return 0
    # ALT=. is not a real allele; keep a 0 placeholder so Number=R matches pysam's
    # alleles=[REF, '.'] representation. Merge drops this slot and fill_missing_ad.py
    # writes 0 for the real ALT(s) afterward.
    if allele == '.':
        return 0
    if allele == ref:
        return site.snp.get(ref[0], 0)
    if len(ref) == 1 and len(allele) == 1:
        return site.snp.get(allele, 0)
    if len(allele) > len(ref) and allele.startswith(ref):
        return site.ins.get(allele[len(ref):], 0)
    if len(allele) < len(ref) and ref.startswith(allele):
        return site.dele.get(len(ref) - len(allele), 0)
    if len(allele) == len(ref):
        return site.snp.get(allele[0], 0)
    return 0


def recover_ad_dp(pileup_sites, pos, alleles):
    """
    Recover FORMAT/AD and FORMAT/DP for consensus alleles from pileup.

    Args:
        pileup_sites: Dict from collect_region_pileup.
        pos: 1-based variant position.
        alleles: List/tuple of alleles (REF followed by ALTs).

    Returns:
        Tuple of (AD tuple, DP int).
    """
    site = pileup_sites.get(pos)
    ref = alleles[0]
    ad = tuple(count_allele_depth(site, ref, allele) for allele in alleles)
    dp = site.dp if site is not None else 0
    return ad, dp


def generate_mininimal_header(sample_name,
                              reference_file,
                              consensus_threshold,
                              split_multiallelic,
                              remove_invariant,
                              version):
    """
    Generate a minimal VCF header for the consensus file.

    Args:
        sample_name: String sample name.
        reference_file: pysam.FastaFile object.
        consensus_threshold: Integer threshold for the consensus genotype.
        split_multiallelic: Whether multiallelic sites were split.
        remove_invariant: Whether monomorphic sites were removed.
        version: String version of the ChoCallate pipeline.

    Returns:
        A pysam.VariantHeader object.
    """
    header = pysam.VariantHeader()
    header.add_sample(sample_name)
    header.add_line(f'##tool=ChoCallate {version}')
    header.add_line(f'##consensusThreshold={consensus_threshold}')
    header.add_line(f'##multialleleSplit={split_multiallelic}')
    header.add_line(f'##noMonomorphic={remove_invariant}')
    for ref in reference_file.references:
        header.add_line(f'##contig=<ID={ref},length={reference_file.get_reference_length(ref)}>')
    header.add_line('##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of samples with data">')
    header.add_line('##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">')
    header.add_line('##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes">')
    header.add_line('##INFO=<ID=AC_Hom,Number=A,Type=Integer,Description="Allele counts in homozygous genotypes">')
    header.add_line('##INFO=<ID=AC_Het,Number=A,Type=Integer,Description="Allele counts in heterozygous genotypes">')
    header.add_line('##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">')
    header.add_line('##INFO=<ID=MAF,Number=1,Type=Float,Description="Frequency of the second most common allele">')
    header.add_line('##INFO=<ID=TYPE,Number=.,Type=String,Description="Variant type">')
    header.add_line('##INFO=<ID=F_MISSING,Number=1,Type=Float,Description="Fraction of missing genotypes">')
    header.add_line('##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth summed across samples">')
    header.add_line('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    header.add_line('##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">')
    header.add_line('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">')
    return header



def main():
    start_time = time.time()
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True, nargs='+', help='Input BCF/VCF files')
    parser.add_argument('--bed', required=True, help='Input BED file (gzipped and tabixed)')
    parser.add_argument('--bam', required=True, help='Input BAM file (indexed)')
    parser.add_argument('--output', required=True, help='Output BCF file')
    parser.add_argument('--sample_name', required=True, help='Sample name')
    parser.add_argument('--reference', required=True, help='Reference genome assembly (fasta file and fai index)')
    parser.add_argument('--consensus_threshold', required=True, type=int, help='Number of matching calls for consensus')
    parser.add_argument('--version', required=True, type=str, help='Version of ChoCallate pipeline')
    parser.add_argument('--split_multiallelic', action='store_true')
    parser.add_argument('--remove_invariant', action='store_true')
    args = parser.parse_args()

    with pysam.TabixFile(args.bed, index=args.bed + '.csi') as bed_file, \
         pysam.FastaFile(args.reference, filepath_index=args.reference + '.fai') as reference_file, \
         pysam.AlignmentFile(args.bam, 'rb') as bam_file:
        with pysam.VariantFile(args.output, mode='wb', 
                               header=generate_mininimal_header(args.sample_name,
                                                                reference_file,
                                                                args.consensus_threshold,
                                                                args.split_multiallelic,
                                                                args.remove_invariant,
                                                                args.version)) as output_file:
            with ExitStack[bool | None]() as stack:
                variant_files = [stack.enter_context(pysam.VariantFile(bcf, index_filename=bcf + '.csi')) for bcf in args.input]
                position_progress = 0
                for line in bed_file.fetch():
                    contig, start, end = line.split('\t')
                    vcf_start = int(start) + 1
                    vcf_end = int(end)
                    ref_seq = reference_file.fetch(contig, int(start), int(end)).upper()
                    records_per_file = [list(vf.fetch(contig, int(start), int(end))) for vf in variant_files]
                    pileup_sites = collect_region_pileup(bam_file, contig, int(start), int(end))
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
                        consensus_genotype = get_consensus_genotype(genotypes, args.consensus_threshold)

                        # write consensus genotype to output file
                        if consensus_genotype is not None:
                            alleles = [consensus_genotype[0]] + list(consensus_genotype[1])
                            ad, dp = recover_ad_dp(pileup_sites, vcf_start, alleles)
                            new_record = output_file.new_record(
                                contig=contig,
                                start=vcf_start - 1,
                                stop=vcf_start,
                                alleles=alleles,
                                qual=None,
                            )
                            new_record.ref = consensus_genotype[0]
                            new_record.samples[0]['GT'] = consensus_genotype[2]
                            new_record.samples[0]['AD'] = ad
                            new_record.samples[0]['DP'] = dp
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
