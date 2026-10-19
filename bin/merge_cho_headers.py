#!/usr/bin/env python3
"""Aggregate per-sample ``##ChoCallate=`` metadata for the merged output.

``bcftools merge`` drops generic meta lines, so ``MERGE_OUTPUTS`` rebuilds the
structured line from the per-sample headers.  The run-level fields are
identical across samples, but ``Callers`` and ``CallersFailed`` can differ in
degraded mode (``calling.allow_caller_failure: true``).  The merged line
therefore records the union of the surviving and failed callers, and one
``##ChoCallateSamples=<...>`` line per sample preserves the per-sample values
whenever they differ from the aggregate.
"""

import argparse
import subprocess


def split_fields(line):
    """Split a ``##ChoCallate=<...>`` line into ``(key, raw_value)`` pairs.

    Values may be quoted and contain commas, so splitting on commas alone is
    not enough.
    """
    body = line[line.index('<') + 1:line.rindex('>')]
    parts = []
    current = []
    quoted = False
    for char in body:
        if char == '"':
            quoted = not quoted
            current.append(char)
        elif char == ',' and not quoted:
            parts.append(''.join(current))
            current = []
        else:
            current.append(char)
    parts.append(''.join(current))

    fields = []
    for part in parts:
        key, separator, raw = part.partition('=')
        if separator:
            fields.append((key.strip(), raw.strip()))
    return fields


def unquote(raw):
    """Strip the surrounding quotes of a structured-field value."""
    if len(raw) >= 2 and raw[0] == raw[-1] == '"':
        return raw[1:-1]
    return raw


def caller_list(fields, key):
    """Names in a quoted comma-separated field, or an empty list."""
    for field_key, raw in fields:
        if field_key == key:
            return [name for name in unquote(raw).split(',') if name]
    return []


def ordered_union(lists):
    """Union of lists, keeping the first occurrence order."""
    seen = set()
    merged = []
    for values in lists:
        for value in values:
            if value not in seen:
                seen.add(value)
                merged.append(value)
    return merged


def set_field(fields, key, value):
    """Replace a field's value, appending it when the key is absent."""
    for index, (field_key, _raw) in enumerate(fields):
        if field_key == key:
            fields[index] = (key, f'"{value}"')
            return
    fields.append((key, f'"{value}"'))


def remove_field(fields, key):
    """Drop every occurrence of a field."""
    fields[:] = [field for field in fields if field[0] != key]


def read_header(bcftools, path):
    """Return the text header of a BCF/VCF file."""
    result = subprocess.run(
        [bcftools, 'view', '-h', path],
        check=True,
        capture_output=True,
        text=True,
    )
    return result.stdout


def sample_name(header_text, path):
    """Sample name from the ``#CHROM`` line, falling back to the file stem."""
    for line in header_text.splitlines():
        if line.startswith('#CHROM'):
            columns = line.split('\t')
            if len(columns) > 9:
                return columns[9]
    name = path.rsplit('/', 1)[-1]
    for suffix in ('.bcf', '.vcf.gz', '.vcf'):
        if name.endswith(suffix):
            return name[:-len(suffix)]
    return name


def chocallate_line(header_text):
    """The structured ChoCallate line, or None when the header lacks one."""
    for line in header_text.splitlines():
        if line.startswith('##ChoCallate=<'):
            return line
    return None


def collect(paths, bcftools):
    """Metadata of every input that carries a ``##ChoCallate=`` line."""
    entries = []
    for path in paths:
        header_text = read_header(bcftools, path)
        line = chocallate_line(header_text)
        if line is None:
            continue
        fields = split_fields(line)
        entries.append({
            'sample': sample_name(header_text, path),
            'fields': fields,
            'callers': caller_list(fields, 'Callers'),
            'failed': caller_list(fields, 'CallersFailed'),
        })
    return entries


def build_header(entries, samples_dropped):
    """Render the merged ``##ChoCallate=`` header lines."""
    if not entries:
        return []

    union_callers = ordered_union(
        entry['callers'] for entry in entries
    )
    union_failed = ordered_union(
        entry['failed'] for entry in entries
    )

    fields = list(entries[0]['fields'])
    set_field(fields, 'Callers', ','.join(union_callers))
    remove_field(fields, 'CallersFailed')
    if union_failed:
        fields.append(('CallersFailed', f'"{",".join(union_failed)}"'))
    remove_field(fields, 'SamplesDropped')
    if samples_dropped:
        fields.append(('SamplesDropped', f'"{samples_dropped}"'))

    lines = [
        '##ChoCallate=<' + ','.join(
            f'{key}={raw}' for key, raw in fields
        ) + '>'
    ]

    differs = union_failed or any(
        entry['callers'] != union_callers for entry in entries
    )
    if differs:
        for entry in entries:
            parts = [
                f'Sample={entry["sample"]}',
                'Callers="' + ','.join(entry['callers']) + '"',
            ]
            if entry['failed']:
                parts.append(
                    'CallersFailed="' + ','.join(entry['failed']) + '"'
                )
            lines.append('##ChoCallateSamples=<' + ','.join(parts) + '>')
    return lines


def parse_args():
    parser = argparse.ArgumentParser(
        description='Aggregate per-sample ##ChoCallate headers for a merge.'
    )
    parser.add_argument(
        '--file-list',
        required=True,
        help='Input BCF paths, one per line',
    )
    parser.add_argument(
        '--bcftools',
        default='bcftools',
        help='bcftools executable (default: bcftools)',
    )
    parser.add_argument(
        '--samples-dropped',
        default='',
        help='Comma-separated sample ids removed before the merge',
    )
    return parser.parse_args()


def main():
    args = parse_args()
    with open(args.file_list, encoding='utf-8') as handle:
        paths = [line.strip() for line in handle if line.strip()]
    entries = collect(paths, args.bcftools)
    for line in build_header(entries, args.samples_dropped):
        print(line)


if __name__ == '__main__':
    main()
