#!/usr/bin/env python3

import argparse
import atexit
import os
import shlex
import signal
import subprocess
import sys
import tempfile


class BcfChunkMerge:
    def __init__(self, cpus, pipe_dir=None):
        self.cpus = cpus
        self.pipe_dir = pipe_dir or os.environ.get("TMPDIR", "/tmp")
        self.pipes = []
        self.processes = []

    def cleanup(self):
        for proc in self.processes:
            if proc.poll() is None:
                proc.terminate()
                try:
                    proc.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    proc.kill()
        for pipe in self.pipes:
            try:
                os.unlink(pipe)
            except FileNotFoundError:
                pass

    def _make_pipe(self):
        fd, path = tempfile.mkstemp(prefix="bcfmerge_pipe.", suffix=".bcf", dir=self.pipe_dir)
        os.close(fd)
        os.unlink(path)
        os.mkfifo(path)
        self.pipes.append(path)
        return path

    def _spawn_shell(self, shell_cmd, stdout=None):
        proc = subprocess.Popen(["sh", "-c", shell_cmd], stdout=stdout)
        self.processes.append(proc)
        return proc

    def _split_chunks(self, files):
        n = len(files)
        k = min(self.cpus, n)
        base, extra = divmod(n, k)
        chunks = []
        idx = 0
        for i in range(k):
            size = base + (1 if i < extra else 0)
            chunks.append(files[idx : idx + size])
            idx += size
        return chunks

    def _chunk_cmd(self, files, out_pipe):
        out = shlex.quote(out_pipe)
        if len(files) == 1:
            return f"bcftools view -Ou {shlex.quote(files[0])} > {out}"
        inputs = " ".join(shlex.quote(path) for path in files)
        return (
            "bcftools merge -Ou -m all --force-samples --force-single --threads 1 "
            f"{inputs} > {out}"
        )

    def _final_cmd(self, pipes):
        inputs = " ".join(shlex.quote(pipe) for pipe in pipes)
        return (
            f"bcftools merge --no-index -Ou -m all --force-samples --threads {self.cpus} "
            f"{inputs}"
        )

    def _wait_all(self):
        errors = []
        for proc in self.processes:
            rc = proc.wait()
            if rc != 0:
                errors.append(rc)
        if errors:
            raise SystemExit(f"Merge failed with exit codes: {errors}")

    def run(self, files):
        for path in files:
            if not os.path.isfile(path):
                raise FileNotFoundError(f"Input file not found: {path}")
        os.makedirs(self.pipe_dir, exist_ok=True)

        groups = self._split_chunks(files)
        pipes = [self._make_pipe() for _ in groups]

        self._spawn_shell(self._final_cmd(pipes), stdout=sys.stdout.buffer)
        for group, pipe in zip(groups, pipes):
            self._spawn_shell(self._chunk_cmd(group, pipe))

        self._wait_all()


def parse_args():
    parser = argparse.ArgumentParser(
        description="Merge BCF files for ChoCallate MERGE_OUTPUTS (chunked via named pipes)."
    )
    parser.add_argument(
        "--file-list",
        required=True,
        help="Input BCF paths, one per line",
    )
    parser.add_argument(
        "--cpus",
        type=int,
        required=True,
        help="Chunk count and final-merge thread count (consensus.cpu / task.cpus)",
    )
    return parser.parse_args()


def collect_files(args):
    with open(args.file_list, encoding="utf-8") as handle:
        files = [line.strip() for line in handle if line.strip()]
    if not files:
        print("No input BCF files provided", file=sys.stderr)
        sys.exit(1)
    return files


def main():
    args = parse_args()
    if args.cpus < 1:
        print("--cpus must be >= 1", file=sys.stderr)
        sys.exit(1)

    merger = BcfChunkMerge(cpus=args.cpus)

    def _cleanup_handler(signum=None, frame=None):
        merger.cleanup()
        raise SystemExit(128 + signum if signum else 1)

    atexit.register(merger.cleanup)
    signal.signal(signal.SIGINT, _cleanup_handler)
    signal.signal(signal.SIGTERM, _cleanup_handler)
    merger.run(collect_files(args))


if __name__ == "__main__":
    main()
