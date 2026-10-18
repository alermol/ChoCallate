#!/usr/bin/env python3
"""Write the ##ChoCallateCommand VCF header line.

The command line is read from the NXF_CLI environment variable at task
runtime, so it never becomes part of a process script or task hash.
"""

import os
import sys


def _header_value(value):
    """Escape a value for a quoted VCF meta line."""
    return (str(value).replace('\\', '\\\\')
            .replace('"', '\\"')
            .replace('\t', ' ')
            .replace('\r', ' ')
            .replace('\n', ' '))


def main():
    if len(sys.argv) != 2:
        sys.exit('usage: write_command_header.py OUTPUT')
    command_line = os.environ.get('NXF_CLI', '').strip() or 'unknown'
    with open(sys.argv[1], 'w', encoding='utf-8') as handle:
        handle.write('##ChoCallateCommand="{}"\n'.format(
            _header_value(command_line)))


if __name__ == '__main__':
    main()
