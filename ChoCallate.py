#!/usr/bin/env python3

from subprocess import run
import argparse
from pathlib import Path
import tempfile
from itertools import islice
import sys


def parser_resolve_path(path):
    try:
        full_path = Path(path).resolve(strict=True)
    except FileNotFoundError:
        sys.exit(f'ERROR! The file {path} does not exist. Emergency termination.')
    return full_path


def split_samplesheet(file, n):
    chunk_files = []
    with tempfile.TemporaryDirectory(delete=False) as tmpdir, open(file) as samplesheet:
        tmpdir = Path(tmpdir)
        lines = samplesheet.readlines()
        it = iter(lines)
        lines = [list(islice(it, n)) for _ in range((len(lines) + n - 1) // n)]
        for cid, chunk in enumerate(lines, start=1):
            filename = tmpdir.joinpath(f'{Path(file).name}.{cid}')
            chunk_files.append(filename)
            with open(filename, 'w') as chunk_file:
                chunk_file.write("".join(chunk))
    return chunk_files


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--samples_tsv', type=parser_resolve_path, default='samples.tsv')
    parser.add_argument('--outdir', type=parser_resolve_path, default='ChoCallate_output')
    parser.add_argument('--reference_genome', type=parser_resolve_path, required=True)
    parser.add_argument('--reference_index', type=parser_resolve_path, required=True)
    parser.add_argument('--min_coverage', type=int, default=5)
    parser.add_argument('--min_base_quality', type=int, default=20)
    parser.add_argument('--samtools_min_map_qual', type=int, default=10)
    parser.add_argument('--min_snp_qual', type=int, default=20)
    parser.add_argument('--reads_type', type=str, default='pe', choices=['pe', 'se'])
    parser.add_argument('--reads_source', type=str, default='gbs', choices=['gbs', 'wgs'])
    parser.add_argument('--bowtie2_cpu', type=int, default=10)
    parser.add_argument('--bowtie2_forks', type=int, default=1)
    parser.add_argument('--freebayes_forks', type=int, default=1)
    parser.add_argument('--bcftools_cpu', type=int, default=1)
    parser.add_argument('--bcftools_forks', type=int, default=1)
    parser.add_argument('--gatk4_forks', type=int, default=1)
    parser.add_argument('--vardict_cpu', type=int, default=1)
    parser.add_argument('--vardict_forks', type=int, default=1)
    parser.add_argument('--snver_forks', type=int, default=1)
    parser.add_argument('--cons_forks', type=int, default=1)
    parser.add_argument('--ploidy', type=int, default=2)
    parser.add_argument('--debug', action='store_true', default=False)
    parser.add_argument('--chunk_size', type=int, default=0)
    args = parser.parse_args()
    
    # add validation of arguments


    diploid_pipeline_path = parser_resolve_path('modules/local/nextflow/run/diploid_calling.nf')
    polyploid_pipeline_path = parser_resolve_path('modules/local/nextflow/run/polyploid_calling.nf')

    chunks = [args.samples_tsv] if args.chunk_size == 0 else split_samplesheet(args.samples_tsv, args.chunk_size)
    if args.chunk_size == 0:
        print('The input file will be processed as a single file.')
    else:
        print(f'The input file was split into {len(chunks)} files.')


    for c in chunks:
        run(
            f"""
            {polyploid_pipeline_path if args.ploidy > 2 else diploid_pipeline_path} \
                --samples_tsv {str(c)} \
                --outdir {args.outdir} \
                --reference_genome {args.reference_genome} \
                --reference_index {args.reference_index} \
                --min_coverage {args.min_coverage} \
                --min_base_quality {args.min_base_quality} \
                --samtools_min_map_qual {args.samtools_min_map_qual} \
                --min_snp_qual {args.min_snp_qual} \
                --reads_type {args.reads_type} \
                --reads_source {args.reads_source} \
                --bowtie2_cpu {args.bowtie2_cpu} \
                --bowtie2_forks {args.bowtie2_forks} \
                --freebayes_forks {args.freebayes_forks} \
                {"" if args.ploidy > 2 else f"--bcftools_cpu {args.bcftools_cpu}"} \
                {"" if args.ploidy > 2 else f"--bcftools_forks {args.bcftools_forks}"} \
                --gatk4_forks {args.gatk4_forks} \
                {"" if args.ploidy > 2 else f"--vardict_cpu {args.vardict_cpu}"} \
                {"" if args.ploidy > 2 else f"--vardict_forks {args.vardict_forks}"} \
                --snver_forks {args.snver_forks} \
                --cons_forks {args.cons_forks} \
                {f"--ploidy {args.ploidy}" if args.ploidy > 2 else ""}
                {"--debug" if args.debug else ""}
                """,
                shell=True
            )
