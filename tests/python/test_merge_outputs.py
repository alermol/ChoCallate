"""Tests for the merged (multi-sample) output stage.

Covers the review findings on `MERGE_OUTPUTS`:

* C-08: with `remove_invariant: false` the merged output must contain all sites,
  monomorphic ones included.  `bcftools merge` turns the per-sample ALT="."
  placeholder into a literal "." allele when another sample has a variant at the
  same position, so the chain also has to drop that duplicate record
  (bin/drop_refonly_duplicates.py).
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

# The same expression as modules/local/merge_outputs/main.nf; `\.` matches the
# literal ALT "." of monomorphic records.
INVARIANT_EXPRESSION = 'ALT~"\\." || COUNT(GT="RR")=N_SAMPLES'

# 1-based position and 0-based start of the variant carried by sampleA
VARIANT_POS = 30
VARIANT_0B = 29


def build_two_samples(tmp_path, use_accel):
    """Per-sample consensus BCFs: sampleA is 0/1 at position 30, sampleB is ref."""
    reference = write_reference(tmp_path / 'ref.fa')
    reads = [(f'read{i}', i, 'ACGT' * 4) for i in range(19, 31)]
    bam = write_bam(tmp_path / 'sample.bam', reads)
    bed_gz = write_bed(tmp_path / 'coverage.bed', [('chr1', 19, 32)])

    variant_calls = [
        write_caller_bcf(tmp_path / f'sampleA_caller{i}.bcf', 'sampleA',
                         [(VARIANT_0B, ref_base(VARIANT_POS),
                           [alt_base(VARIANT_POS)], (0, 1))])
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
    return sample_files


def merge_chain(tmp_path, sample_files, remove_invariant, split_multiallelic=False):
    """Run the MERGE_OUTPUTS command chain and return the merged BCF path."""
    import shlex
    list_file = tmp_path / 'merge_list.txt'
    list_file.write_text('\n'.join(str(path) for path in sample_files.values()) + '\n')
    if remove_invariant:
        invariant_filter = (f'{BCFTOOLS} filter -e '
                            f'{shlex.quote(INVARIANT_EXPRESSION)} -Ou')
    else:
        invariant_filter = f'{BCFTOOLS} view -Ou'
    join_multiallelic = ('' if split_multiallelic
                         else f'| {BCFTOOLS} norm -m +any -Ou')
    output = tmp_path / 'consensus.bcf'
    chain = (
        f"{sys.executable} {BIN_DIR / 'merge_bcf_tree.py'} --cpus 1"
        f" --file-list {list_file}"
        f" | {BCFTOOLS} norm -m -any -Ou"
        f" | {sys.executable} {BIN_DIR / 'drop_refonly_duplicates.py'}"
        f" | {invariant_filter} {join_multiallelic}"
        f" | {sys.executable} {BIN_DIR / 'fill_missing_ad.py'}"
        f" | {BCFTOOLS} +fill-tags -Ou -- -t {FILL_TAGS}"
        f" | {BCFTOOLS} sort -Ob -o {output}"
    )
    subprocess.run(['bash', '-c', chain], check=True, cwd=str(REPO_ROOT))
    return output


def read_merged(path):
    with pysam.VariantFile(str(path), index_filename=None) as variant_file:
        return [(record.pos, (record.ref,) + tuple(record.alts),
                 tuple(record.samples['sampleA']['GT']),
                 tuple(record.samples['sampleB']['GT'])) for record in variant_file]


def test_all_sites_kept_without_remove_invariant(tmp_path, use_accel):
    """C-08: remove_invariant=false keeps every covered position exactly once."""
    sample_files = build_two_samples(tmp_path, use_accel)
    merged = merge_chain(tmp_path, sample_files, remove_invariant=False)
    records = read_merged(merged)
    assert [pos for pos, _, _, _ in records] == list(range(20, 33))
    variant = [record for record in records if record[0] == VARIANT_POS]
    assert variant == [(VARIANT_POS, (ref_base(VARIANT_POS), alt_base(VARIANT_POS)),
                        (0, 1), (0, 0))]
    for pos, alleles, gt_a, gt_b in records:
        if pos != VARIANT_POS:
            assert alleles[1] == '.'
            assert (gt_a, gt_b) == ((0, 0), (0, 0))


def test_all_sites_kept_without_remove_invariant_split(tmp_path, use_accel):
    """Same, with output.split_multiallelic: true (no `norm -m +any` re-join)."""
    sample_files = build_two_samples(tmp_path, use_accel)
    merged = merge_chain(tmp_path, sample_files, remove_invariant=False,
                         split_multiallelic=True)
    assert [pos for pos, _, _, _ in read_merged(merged)] == list(range(20, 33))


def test_remove_invariant_drops_monomorphic_sites(tmp_path, use_accel):
    """remove_invariant=true still removes the monomorphic records."""
    sample_files = build_two_samples(tmp_path, use_accel)
    merged = merge_chain(tmp_path, sample_files, remove_invariant=True)
    records = read_merged(merged)
    assert [pos for pos, _, _, _ in records] == [VARIANT_POS]
    assert records[0][1] == (ref_base(VARIANT_POS), alt_base(VARIANT_POS))
    assert records[0][2:] == ((0, 1), (0, 0))


def test_drop_refonly_duplicates_keeps_monomorphic_positions(tmp_path):
    """The dedup helper drops ALT="." duplicates but keeps monomorphic sites."""
    header = pysam.VariantHeader()
    header.add_line('##contig=<ID=chr1,length=100>')
    header.add_line('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    header.add_sample('sample1')
    source = tmp_path / 'in.bcf'
    with pysam.VariantFile(str(source), 'wb', header=header) as out:
        # monomorphic position: kept
        record = out.new_record(contig='chr1', start=9, stop=10, alleles=('A', '.'))
        record.samples[0]['GT'] = (0, 0)
        out.write(record)
        # variant position with the duplicated ALT='.' record: only the variant
        for start, alleles, gt in ((29, ('C', '.'), (0, 0)), (29, ('C', 'T'), (0, 1))):
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
            [(10, ('.',)), (30, ('T',))]
