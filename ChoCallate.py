#!/usr/bin/env python3

import argparse
import os
import sys
import tempfile
from itertools import islice
from pathlib import Path
from subprocess import run
from typing import List, Optional, Union, Any


class RangeCheckAction(argparse.Action):
    def __init__(self,
                 option_strings: List[str],
                 dest: str,
                 min_val: Optional[int] = None,
                 max_val: Optional[int] = None,
                 no_upper_limit: bool = False,
                 **kwargs: Any) -> None:
        super().__init__(option_strings, dest, **kwargs)
        self.min_val: Optional[int] = min_val
        self.max_val: Optional[int] = max_val
        self.no_upper_limit: bool = no_upper_limit

    def __call__(self, parser: argparse.ArgumentParser, namespace: argparse.Namespace, 
                 value: str, option_string: Optional[str] = None) -> None:
        ivalue: int = int(value)
        if self.no_upper_limit == True:
            self.max_val = None
            if ((self.min_val is not None) and (ivalue < self.min_val)):
                sys.exit(
                    f"Argument {self.dest}: Invalid value: {value} (must be greater or equal to {self.min_val}). Emergency termination."
                )
        else:
            if (self.min_val is not None and ivalue < self.min_val) or \
               (self.max_val is not None and ivalue > self.max_val):
                sys.exit(
                    f"Argument {self.dest}: Invalid value: {value} (must be in range {self.min_val}-{self.max_val}). Emergency termination."
                )
        setattr(namespace, self.dest, ivalue)


def parser_resolve_path(path: str) -> Path:
    try:
        full_path: Path = Path(path).resolve(strict=True)
    except FileNotFoundError:
        sys.exit(f'ERROR! The file {path} does not exist. Emergency termination.')
    return full_path


def split_samplesheet(file: Union[str, Path], n: int) -> List[Path]:
    chunk_files: List[Path] = []
    with tempfile.TemporaryDirectory(delete=False) as tmpdir, open(file) as samplesheet:
        tmpdir = Path(tmpdir)
        lines: List[str] = samplesheet.readlines()
        it = iter(lines)
        lines = [list(islice(it, n)) for _ in range((len(lines) + n - 1) // n)]
        for cid, chunk in enumerate(lines, start=1):
            filename: Path = tmpdir.joinpath(f'{Path(file).name}.{cid}')
            chunk_files.append(filename)
            with open(filename, 'w') as chunk_file:
                chunk_file.write("".join(chunk))
    return chunk_files


if __name__ == "__main__":
    parser: argparse.ArgumentParser = argparse.ArgumentParser()
    parser.add_argument('--samples_tsv', type=parser_resolve_path, default='samples.tsv')
    parser.add_argument('--outdir', type=str, default='ChoCallate_output')
    parser.add_argument('--reference_genome', type=parser_resolve_path, required=True)
    parser.add_argument('--reference_index', type=parser_resolve_path, required=True)
    parser.add_argument('--min_coverage', type=int, action=RangeCheckAction, no_upper_limit=True, min_val=0, default=5)
    parser.add_argument('--min_base_quality', type=int, action=RangeCheckAction, no_upper_limit=True, min_val=0, default=20)
    parser.add_argument('--min_map_qual', type=int, action=RangeCheckAction, no_upper_limit=True, min_val=0, default=10)
    parser.add_argument('--min_snp_qual', type=int, action=RangeCheckAction, no_upper_limit=True, min_val=0, default=20)
    parser.add_argument('--reads_type', type=str, default='pe', choices=['pe', 'se', 'mx'])
    parser.add_argument('--reads_source', type=str, default='gbs', choices=['gbs', 'wgs'])
    parser.add_argument('--bowtie2_cpu', type=int, action=RangeCheckAction, min_val=1, max_val=os.cpu_count(), default=10)
    parser.add_argument('--bowtie2_forks', type=int, action=RangeCheckAction, min_val=1, max_val=os.cpu_count(), default=1)
    parser.add_argument('--calling_forks', type=int, action=RangeCheckAction, min_val=1, max_val=os.cpu_count(), default=1)
    parser.add_argument('--zero_vcf_forks', type=int, action=RangeCheckAction, min_val=1, max_val=os.cpu_count(), default=1)
    parser.add_argument('--zero_vcf_cpu', type=int, action=RangeCheckAction, min_val=1, max_val=os.cpu_count(), default=1)
    parser.add_argument('--vardict_cpu', type=int, action=RangeCheckAction, min_val=1, max_val=os.cpu_count(), default=1)
    parser.add_argument('--bcftools_cpu', type=int, action=RangeCheckAction, min_val=1, max_val=os.cpu_count(), default=1)
    parser.add_argument('--cons_forks', type=int, action=RangeCheckAction, min_val=1, max_val=os.cpu_count(), default=1)
    parser.add_argument('--cons_cpus', type=int, action=RangeCheckAction, min_val=1, max_val=os.cpu_count(), default=5)
    parser.add_argument('--win_size', type=int, action=RangeCheckAction, no_upper_limit=True, min_val=1, default=1000000)
    parser.add_argument('--ploidy', type=int, action=RangeCheckAction, no_upper_limit=True, min_val=2, default=2)
    parser.add_argument('--debug', action='store_true', default=False)
    parser.add_argument('--chunk_size', type=int, action=RangeCheckAction, no_upper_limit=True, min_val=0, default=0)
    parser.add_argument('--cons_type', type=str, default='mj', choices=['mj', 'n1', 'fc'],
                        help='Consensus type: mj (majority rule), n1 (n-1 consensus), fc (full consensus)')
    parser.add_argument('--effective_callers', type=str, default='-',
                        help='Available callers: bcftools,freebayes,snver,vardict,gatk')
    args: argparse.Namespace = parser.parse_args()

    Path.mkdir(Path(args.outdir), parents=True, exist_ok=True)

    pipeline_path: Path = Path(__file__).parent.joinpath('ChoCallate.nf')

    chunks: List[Path] = [args.samples_tsv] if args.chunk_size == 0 else split_samplesheet(args.samples_tsv, args.chunk_size)
    if args.chunk_size == 0:
        print('The input file will be processed as a single file.')
    else:
        print(f'The input file was split into {len(chunks)} files.')


    for c in chunks:
        run(
            f"""
            {pipeline_path} \
                --samples_tsv {str(c)} \
                --outdir {parser_resolve_path(args.outdir)} \
                --reference_genome {args.reference_genome} \
                --reference_index {args.reference_index} \
                --min_coverage {args.min_coverage} \
                --min_base_quality {args.min_base_quality} \
                --min_map_qual {args.min_map_qual} \
                --min_snp_qual {args.min_snp_qual} \
                --reads_type {args.reads_type} \
                --reads_source {args.reads_source} \
                --bowtie2_cpu {args.bowtie2_cpu} \
                --bowtie2_forks {args.bowtie2_forks} \
                --zero_vcf_forks {args.zero_vcf_forks} \
                {"" if args.ploidy > 2 else f"--bcftools_cpu {args.bcftools_cpu}"} \
                {"" if args.ploidy > 2 else f"--vardict_cpu {args.vardict_cpu}"} \
                --cons_forks {args.cons_forks} \
                --cons_cpus {args.cons_cpus} \
                --win_size {args.win_size} \
                --effective_callers {args.effective_callers} \
                {f"--ploidy {args.ploidy}" if args.ploidy > 2 else ""} \
                {f"--cons_type {args.cons_type}" if args.cons_type else ""} \
                {"--debug" if args.debug else ""}
                """,
                shell=True
            )
