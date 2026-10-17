#!/bin/bash
set -u

# The exit status of the task script is passed by the afterScript directive.
# Nextflow only sets it on failure, so an absent argument means success.
exit_status="${1:-0}"
if [ "${exit_status}" -ne 0 ]; then
    exit 0
fi

rm -rf -- tmp/

# Only remove a temporary directory created inside this work directory.
if [ -n "${TMPDIR:-}" ] && [ "${TMPDIR}" != "/" ]; then
    case "${TMPDIR}" in
        "${PWD}"/*)
            rm -rf -- "${TMPDIR}"
            ;;
    esac
fi
