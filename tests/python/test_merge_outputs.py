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
import time

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


def build_two_samples(tmp_path, use_accel, compress_reference_blocks=False):
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
            compress_reference_blocks=compress_reference_blocks)
        subprocess.run([BCFTOOLS, 'index', '--csi', str(sample_files[sample])],
                       check=True)
    return sample_files, reference


def merge_chain(tmp_path, sample_files, reference, remove_invariant,
                split_multiallelic=False, output_format='bcf', cpus=1,
                samples_dropped=None):
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
        '--cpus', str(cpus),
        '--bcftools', BCFTOOLS,
        '--python', sys.executable,
    ]
    if remove_invariant:
        cmd.append('--remove-invariant')
    if split_multiallelic:
        cmd.append('--split-multiallelic')
    if samples_dropped:
        cmd += ['--samples-dropped', samples_dropped]
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
        entries = []
        for record in variant_file:
            ad_a = record.samples['sampleA'].get('AD')
            ad_b = record.samples['sampleB'].get('AD')
            entries.append((record.pos,
                            tuple(ad_a) if ad_a is not None else None,
                            tuple(ad_b) if ad_b is not None else None))
        return entries


def read_merged_spans(path):
    """(pos, inclusive end, alleles, GT A, GT B, DP A, DP B) tuples."""
    with pysam.VariantFile(str(path), index_filename=None) as variant_file:
        return [(record.pos, record.stop,
                 (record.ref,) + tuple(record.alts),
                 tuple(record.samples['sampleA']['GT']),
                 tuple(record.samples['sampleB']['GT']),
                 record.samples['sampleA']['DP'],
                 record.samples['sampleB']['DP'])
                for record in variant_file]


def expanded_positions(spans):
    return [position for start, end, *_ in spans
            for position in range(start, end + 1)]


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


def test_reference_blocks_survive_the_merge(tmp_path, use_accel):
    """Compressed per-sample input stays compressed and covers all sites."""
    sample_files, reference = build_two_samples(
        tmp_path, use_accel, compress_reference_blocks=True)
    merged = merge_chain(tmp_path, sample_files, reference,
                         remove_invariant=False)
    spans = read_merged_spans(merged)
    assert expanded_positions(spans) == list(range(20, 33))
    # 20..29 is one block, 30 the only variant, 31..32 a trailing block.
    assert [(start, end) for start, end, *_ in spans] == [
        (20, 29), (30, 30), (31, 32)]
    assert any(end > start for start, end, *_ in spans)
    with pysam.VariantFile(str(merged), index_filename=None) as variant_file:
        assert 'ReferenceBlocks=True' in str(variant_file.header)
    variant = [span for span in spans if span[0] == VARIANT_POS][0]
    assert variant[2] == (ref_base(VARIANT_POS), alt_base(VARIANT_POS))
    assert variant[3:5] == ((0, 1), (0, 0))
    # Block DP is the minimum over the block (1 at position 20).
    assert variant[6] == 1


def test_reference_block_alt_depth_is_missing_not_zero(tmp_path, use_accel):
    """A compressed hom-reference sample reports missing, not fake, ALT AD."""
    sample_files, reference = build_two_samples(
        tmp_path, use_accel, compress_reference_blocks=True)
    merged = merge_chain(tmp_path, sample_files, reference,
                         remove_invariant=False)
    variant_ad = [entry for entry in read_merged_ad(merged)
                  if entry[0] == VARIANT_POS]
    assert len(variant_ad) == 1
    _, ad_a, ad_b = variant_ad[0]
    assert ad_a == (REF_DEPTH, ALT_DEPTH)
    assert ad_b == (None, None)


def test_reference_blocks_survive_a_chunked_merge(tmp_path, use_accel):
    """Blocks pass through the chunked named-pipe merge as well."""
    sample_files, reference = build_two_samples(
        tmp_path, use_accel, compress_reference_blocks=True)
    merged = merge_chain(tmp_path, sample_files, reference,
                         remove_invariant=False, cpus=2)
    spans = read_merged_spans(merged)
    assert expanded_positions(spans) == list(range(20, 33))
    variant = [span for span in spans if span[0] == VARIANT_POS][0]
    assert variant[3:5] == ((0, 1), (0, 0))


def test_reference_blocks_survive_remove_invariant(tmp_path, use_accel):
    """The invariant filter sees block records and keeps only the variant."""
    sample_files, reference = build_two_samples(
        tmp_path, use_accel, compress_reference_blocks=True)
    merged = merge_chain(tmp_path, sample_files, reference,
                         remove_invariant=True)
    spans = read_merged_spans(merged)
    assert [(start, end) for start, end, *_ in spans] == [(VARIANT_POS,
                                                          VARIANT_POS)]


def test_merged_output_has_metadata_line(tmp_path, use_accel):
    """The merged consensus keeps the ##ChoCallate metadata line."""
    sample_files, reference = build_two_samples(tmp_path, use_accel)
    merged = merge_chain(tmp_path, sample_files, reference, remove_invariant=False)
    with pysam.VariantFile(str(merged), index_filename=None) as variant_file:
        header = str(variant_file.header)
    assert '##ChoCallate=' in header
    assert 'Rule=two-stage' in header


def test_merged_header_records_dropped_samples(tmp_path, use_accel):
    """Dropped sample ids are appended to the merged ##ChoCallate line."""
    sample_files, reference = build_two_samples(tmp_path, use_accel)
    merged = merge_chain(tmp_path, sample_files, reference,
                         remove_invariant=False,
                         samples_dropped='sampleC,sampleD')
    with pysam.VariantFile(str(merged), index_filename=None) as variant_file:
        header = str(variant_file.header)
    assert 'SamplesDropped="sampleC,sampleD"' in header
    assert header.count('##ChoCallate=') == 1


def test_single_file_merge_records_dropped_samples(tmp_path, use_accel):
    """The single-input merge path patches the header as well."""
    sample_files, reference = build_two_samples(tmp_path, use_accel)
    one_sample = {'sampleA': sample_files['sampleA']}
    merged = merge_chain(tmp_path, one_sample, reference,
                         remove_invariant=False,
                         samples_dropped='sampleB')
    with pysam.VariantFile(str(merged), index_filename=None) as variant_file:
        header = str(variant_file.header)
    assert 'SamplesDropped="sampleB"' in header
    assert header.count('##ChoCallate=') == 1


def test_dropped_samples_patches_only_chocallate_line(tmp_path, use_accel):
    """Other structured header lines must not gain SamplesDropped."""
    sample_files, reference = build_two_samples(tmp_path, use_accel)
    with pysam.VariantFile(
            str(sample_files['sampleA']), index_filename=None) as source:
        source_header = str(source.header).splitlines()
    merged = merge_chain(tmp_path, sample_files, reference,
                         remove_invariant=False,
                         samples_dropped='sampleC,sampleD')
    with pysam.VariantFile(str(merged), index_filename=None) as variant_file:
        header = str(variant_file.header)
    merged_header = header.splitlines()
    cho_callate = [line for line in merged_header
                   if line.startswith('##ChoCallate=<')]
    assert len(cho_callate) == 1
    assert cho_callate[0].endswith(',SamplesDropped="sampleC,sampleD">')
    for prefix in ('##INFO=<', '##FORMAT=<', '##ALT=<', '##contig=<'):
        expected = [line for line in source_header
                    if line.startswith(prefix)]
        actual = [line for line in merged_header if line.startswith(prefix)]
        assert expected, f'test fixture lacks {prefix} lines'
        assert actual == expected


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


def test_drop_refonly_duplicates_does_not_cross_contigs(tmp_path):
    """Equal positions on different contigs are separate loci.

    The last record of chr1 (``chr1:5 <*>``) and the first record of chr2
    (``chr2:5 A>T``) share POS 5 and must both reach the output.
    """
    header = pysam.VariantHeader()
    header.add_line('##contig=<ID=chr1,length=100>')
    header.add_line('##contig=<ID=chr2,length=100>')
    header.add_line('##ALT=<ID=*,Description="Represents allele(s) other ' +
                    'than observed.">')
    header.add_line('##INFO=<ID=END,Number=1,Type=Integer,' +
                    'Description="Stop position of the interval">')
    header.add_line('##FORMAT=<ID=GT,Number=1,Type=String,' +
                    'Description="Genotype">')
    header.add_sample('sample1')
    source = tmp_path / 'in.bcf'
    with pysam.VariantFile(str(source), 'wb', header=header) as out:
        record = out.new_record(contig='chr1', start=4, stop=5,
                                alleles=('A', '<*>'))
        record.samples[0]['GT'] = (0, 0)
        out.write(record)
        record = out.new_record(contig='chr2', start=4, stop=5,
                                alleles=('A', 'T'))
        record.samples[0]['GT'] = (0, 1)
        out.write(record)

    result = subprocess.run(
        [sys.executable, str(BIN_DIR / 'drop_refonly_duplicates.py')],
        check=True, input=source.read_bytes(), capture_output=True)
    out_path = tmp_path / 'out.bcf'
    out_path.write_bytes(result.stdout)
    with pysam.VariantFile(str(out_path), index_filename=None) as variant_file:
        assert [(record.chrom, record.pos, tuple(record.alts))
                for record in variant_file] == [
            ('chr1', 5, ('<*>',)), ('chr2', 5, ('T',))]


def test_failed_final_merge_does_not_hang_on_fifos(tmp_path):
    """A failing final merge must not leave chunk writers blocked on FIFOs."""
    inputs = [
        write_caller_bcf(tmp_path / f'sample{index}.bcf',
                         f'sample{index}', [])
        for index in range(2)]
    file_list = tmp_path / 'merge_list.txt'
    file_list.write_text('\n'.join(str(path) for path in inputs) + '\n')
    env = os.environ.copy()
    env['TMPDIR'] = str(tmp_path)

    started = time.monotonic()
    result = subprocess.run(
        [sys.executable, str(BIN_DIR / 'merge_bcf_tree.py'),
         '--file-list', str(file_list),
         '--gvcf', str(tmp_path / 'missing.fa'),
         '--cpus', '2'],
        capture_output=True, text=True, timeout=10,
        cwd=str(REPO_ROOT), env=env)
    elapsed = time.monotonic() - started

    assert result.returncode != 0
    assert elapsed < 5


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
