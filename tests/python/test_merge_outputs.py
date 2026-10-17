"""Tests for the merged (multi-sample) output stage.

Covers the review findings on `MERGE_OUTPUTS`:

* C-08: with `remove_invariant: false` the merged output must contain all sites,
  monomorphic ones included.  Per-sample files use gVCF-style `ALT="<*>"`
  records now, so the chain drops the redundant reference-only record that
  `bcftools merge --gvcf` + `norm -m -any` split off at variant positions.
* C-07 (option 1): a hom-reference sample's non-reference read depth must reach
  the merged ALT AD through the `<*>` slot, instead of a fabricated 0.
* The merge path for N >= 2 samples as a whole (T-01): sample order, genotypes
  at the variant position, and the invariant filter when enabled.

The command chain below mirrors modules/local/merge_outputs/main.nf; keep both in
sync when the module changes.
"""

import subprocess
import sys

import pysam

from conftest import (  # noqa: I100,I202  (tests run from this directory)
    BCFTOOLS,
    BIN_DIR,
    REF_SEQ,
    REPO_ROOT,
    alt_base,
    ref_base,
    run_consensus,
    write_bam,
    write_bed,
    write_caller_bcf,
    write_reference,
)

FILL_TAGS = ("AN,AC,AF,NS,AC_Hom,AC_Het,MAF,TYPE,F_MISSING,"
             "'DP:1=int(sum(FORMAT/DP))'")

# 1-based position and 0-based start of the variant carried by sampleA
VARIANT_POS = 30
VARIANT_0B = 29

# The variant position carries 10 reference-supporting and 2 ALT-supporting
# reads; a hom-reference sample must report (10, 2), not (10, 0).
REF_DEPTH = 10
ALT_DEPTH = 2


def build_two_samples(tmp_path, use_accel):
    """Per-sample consensus BCFs: sampleA is 0/1 at position 30, sampleB is ref."""
    reference = write_reference(tmp_path / 'ref.fa')
    alt = alt_base(VARIANT_POS)
    # reference-matching reads covering the whole BED window, plus ALT reads at
    # the variant position
    reads = [(f'ref{start}', start, REF_SEQ[start:start + 10])
             for start in range(19, 30)]
    for i in range(ALT_DEPTH):
        seq = list(REF_SEQ[VARIANT_0B:VARIANT_0B + 10])
        seq[0] = alt
        reads.append((f'alt{i}', VARIANT_0B, ''.join(seq)))
    reads.sort(key=lambda item: item[1])
    bam = write_bam(tmp_path / 'sample.bam', reads)
    bed_gz = write_bed(tmp_path / 'coverage.bed', [('chr1', 19, 32)])

    variant_calls = [
        write_caller_bcf(tmp_path / f'sampleA_caller{i}.bcf', 'sampleA',
                         [(VARIANT_0B, ref_base(VARIANT_POS),
                           [alt], (0, 1))])
        for i in range(2)]
    reference_calls = [
        write_caller_bcf(tmp_path / f'sampleB_caller{i}.bcf', 'sampleB', [])
        for i in range(2)]

    sample_files = {}
    for sample, calls in (('sampleA', variant_calls), ('sampleB', reference_calls)):
        sample_files[sample] = run_consensus(
            tmp_path, bam, bed_gz, calls, threshold=2, sample_name=sample,
            use_accel=use_accel, output_name=f'{sample}.bcf', reference=reference)
        subprocess.run([BCFTOOLS, 'index', '--csi', str(sample_files[sample])],
                       check=True)
    return sample_files, reference


def merge_chain(tmp_path, sample_files, reference, remove_invariant,
                split_multiallelic=False):
    """Run the MERGE_OUTPUTS command chain and return the merged BCF path."""
    list_file = tmp_path / 'merge_list.txt'
    list_file.write_text('\n'.join(str(path) for path in sample_files.values()) + '\n')
    # Mirror MERGE_OUTPUTS: bcftools merge drops generic meta lines, so carry
    # the first sample's ##ChoCallate line over to the merged header.
    first_sample = next(iter(sample_files.values()))
    first_header = subprocess.run(
        [BCFTOOLS, 'view', '-h', str(first_sample)], check=True,
        capture_output=True, text=True).stdout
    meta_file = tmp_path / 'cho_meta.hdr'
    meta_file.write_text(''.join(
        line + '\n' for line in first_header.splitlines()
        if line.startswith('##ChoCallate=')))
    if remove_invariant:
        expression = 'ALT="." || ALT="<*>" || COUNT(GT="RR")=N_SAMPLES'
        invariant_filter = f'{BCFTOOLS} filter -e \'{expression}\' -Ou'
    else:
        invariant_filter = f'{BCFTOOLS} view -Ou'
    join_multiallelic = ('' if split_multiallelic
                         else f'| {BCFTOOLS} norm -m +any -Ou')
    output = tmp_path / 'consensus.bcf'
    chain = (
        f"{sys.executable} {BIN_DIR / 'merge_bcf_tree.py'}"
        f" --gvcf {reference} --cpus 1 --file-list {list_file}"
        f" | {BCFTOOLS} norm -m -any -Ou"
        f" | {sys.executable} {BIN_DIR / 'drop_refonly_duplicates.py'}"
        f" | {invariant_filter} {join_multiallelic}"
        f" | {BCFTOOLS} +fill-tags -Ou -- -t {FILL_TAGS}"
        f" | {BCFTOOLS} annotate -h {meta_file} -Ou"
        f" | {BCFTOOLS} sort -Ob -o {output}"
    )
    subprocess.run(['bash', '-c', chain], check=True, cwd=str(REPO_ROOT))
    return output


def read_merged(path):
    with pysam.VariantFile(str(path), index_filename=None) as variant_file:
        return [(record.pos, (record.ref,) + tuple(record.alts),
                 tuple(record.samples['sampleA']['GT']),
                 tuple(record.samples['sampleB']['GT'])) for record in variant_file]


def read_merged_ad(path):
    with pysam.VariantFile(str(path), index_filename=None) as variant_file:
        return [(record.pos,
                 tuple(record.samples['sampleA']['AD']),
                 tuple(record.samples['sampleB']['AD'])) for record in variant_file]


def test_all_sites_kept_without_remove_invariant(tmp_path, use_accel):
    """C-08: remove_invariant=false keeps every covered position exactly once."""
    sample_files, reference = build_two_samples(tmp_path, use_accel)
    merged = merge_chain(tmp_path, sample_files, reference, remove_invariant=False)
    records = read_merged(merged)
    assert [pos for pos, _, _, _ in records] == list(range(20, 33))
    variant = [record for record in records if record[0] == VARIANT_POS]
    assert variant == [(VARIANT_POS, (ref_base(VARIANT_POS), alt_base(VARIANT_POS)),
                        (0, 1), (0, 0))]
    for pos, alleles, gt_a, gt_b in records:
        if pos != VARIANT_POS:
            assert alleles[1] == '<*>'
            assert (gt_a, gt_b) == ((0, 0), (0, 0))


def test_other_allele_depth_reaches_the_merged_alt_slot(tmp_path, use_accel):
    """C-07 option 1: sampleB's <*> depth becomes its AD for the merged ALT."""
    sample_files, reference = build_two_samples(tmp_path, use_accel)
    merged = merge_chain(tmp_path, sample_files, reference, remove_invariant=False)
    variant_ad = [entry for entry in read_merged_ad(merged) if entry[0] == VARIANT_POS]
    assert variant_ad == [(VARIANT_POS, (REF_DEPTH, ALT_DEPTH),
                           (REF_DEPTH, ALT_DEPTH))]


def test_all_sites_kept_without_remove_invariant_split(tmp_path, use_accel):
    """Same, with output.split_multiallelic: true (no `norm -m +any` re-join)."""
    sample_files, reference = build_two_samples(tmp_path, use_accel)
    merged = merge_chain(tmp_path, sample_files, reference,
                         remove_invariant=False, split_multiallelic=True)
    assert [pos for pos, _, _, _ in read_merged(merged)] == list(range(20, 33))


def test_remove_invariant_drops_monomorphic_sites(tmp_path, use_accel):
    """remove_invariant=true still removes the monomorphic records."""
    sample_files, reference = build_two_samples(tmp_path, use_accel)
    merged = merge_chain(tmp_path, sample_files, reference, remove_invariant=True)
    records = read_merged(merged)
    assert [pos for pos, _, _, _ in records] == [VARIANT_POS]
    assert records[0][1] == (ref_base(VARIANT_POS), alt_base(VARIANT_POS))
    assert records[0][2:] == ((0, 1), (0, 0))


def test_merged_output_has_metadata_line(tmp_path, use_accel):
    """The merged consensus keeps the ##ChoCallate metadata line."""
    sample_files, reference = build_two_samples(tmp_path, use_accel)
    merged = merge_chain(tmp_path, sample_files, reference, remove_invariant=False)
    with pysam.VariantFile(str(merged), index_filename=None) as variant_file:
        header = str(variant_file.header)
    assert '##ChoCallate=' in header
    assert 'Rule=two-stage' in header


def test_drop_refonly_duplicates_keeps_monomorphic_positions(tmp_path):
    """The dedup helper drops '.'/'<*>' duplicates but keeps monomorphic sites."""
    header = pysam.VariantHeader()
    header.add_line('##contig=<ID=chr1,length=100>')
    header.add_line('##ALT=<ID=*,Description="Represents allele(s) other than observed.">')
    header.add_line('##INFO=<ID=END,Number=1,Type=Integer,Description="Stop position of the interval">')
    header.add_line('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    header.add_sample('sample1')
    source = tmp_path / 'in.bcf'
    with pysam.VariantFile(str(source), 'wb', header=header) as out:
        # monomorphic positions with both reference-only symbols: kept
        for start, alt in ((9, '.'), (14, '<*>')):
            record = out.new_record(contig='chr1', start=start, stop=start + 1,
                                    alleles=('A', alt))
            record.samples[0]['GT'] = (0, 0)
            out.write(record)
        # variant position with the duplicated reference-only record: only the
        # variant record survives
        for start, alleles, gt in ((29, ('C', '<*>'), (0, 0)),
                                   (29, ('C', 'T'), (0, 1))):
            record = out.new_record(contig='chr1', start=start, stop=start + 1,
                                    alleles=alleles)
            record.samples[0]['GT'] = gt
            out.write(record)

    result = subprocess.run([sys.executable, str(BIN_DIR / 'drop_refonly_duplicates.py')],
                            check=True, input=source.read_bytes(),
                            capture_output=True)
    out_path = tmp_path / 'out.bcf'
    out_path.write_bytes(result.stdout)
    with pysam.VariantFile(str(out_path), index_filename=None) as variant_file:
        assert [(record.pos, tuple(record.alts)) for record in variant_file] == \
            [(10, ('.',)), (15, ('<*>',)), (30, ('T',))]
