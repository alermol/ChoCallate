"""End-to-end regression tests for the merged (multi-sample) output stage.

Covers:

* all sites are kept when ``remove_invariant`` is false, including monomorphic
  records;
* a hom-reference sample's non-reference read depth reaches the merged ALT AD;
* the merge path for N >= 2 samples as a whole: sample order, genotypes at the
  variant position, and the invariant filter when enabled.

The command chain is exercised through `bin/merge_outputs.sh`, the same script
used by `MERGE_OUTPUTS`; changes to the production merge path are therefore
tested directly.
"""

import os
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
            use_accel=use_accel, output_name=f'{sample}.bcf',
            reference=reference,
            command_line='nextflow run main.nf -params-file run.yaml',
            caller_versions='bcftools=1.24,freebayes=1.3.10')
        subprocess.run([BCFTOOLS, 'index', '--csi', str(sample_files[sample])],
                       check=True)
    return sample_files, reference


def merge_chain(tmp_path, sample_files, reference, remove_invariant,
                split_multiallelic=False, output_format='bcf'):
    """Run bin/merge_outputs.sh and return the merged output path."""
    list_file = tmp_path / 'merge_list.txt'
    list_file.write_text('\n'.join(str(path) for path in sample_files.values()) + '\n')
    output = tmp_path / ('consensus.vcf.gz' if output_format == 'vcf'
                         else 'consensus.bcf')
    cmd = [
        str(BIN_DIR / 'merge_outputs.sh'),
        '--file-list', str(list_file),
        '--reference', str(reference),
        '--output', str(output),
        '--format', output_format,
        '--cpus', '1',
        '--bcftools', BCFTOOLS,
        '--python', sys.executable,
    ]
    if remove_invariant:
        cmd.append('--remove-invariant')
    if split_multiallelic:
        cmd.append('--split-multiallelic')
    env = os.environ.copy()
    env['TMPDIR'] = str(tmp_path)
    subprocess.run(cmd, check=True, cwd=str(REPO_ROOT), env=env)
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
    """With remove_invariant=false every covered position is kept exactly once."""
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
    """A hom-reference sample's <*> depth becomes its AD for the merged ALT."""
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
    assert '##ChoCallateCommand="nextflow run main.nf' in header
    assert '##ChoCallateCallerVersions="bcftools=1.24' in header


def test_merged_vcf_output_is_readable(tmp_path, use_accel):
    """The production merge script writes a valid VCF.gz as well as BCF."""
    sample_files, reference = build_two_samples(tmp_path, use_accel)
    merged = merge_chain(tmp_path, sample_files, reference, remove_invariant=False,
                         output_format='vcf')
    with pysam.VariantFile(str(merged), index_filename=None) as variant_file:
        records = list(variant_file)
        header = str(variant_file.header)
    assert '##ChoCallate=' in header
    assert [record.pos for record in records] == list(range(20, 33))


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
