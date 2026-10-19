/*
 * Configuration and samplesheet validation used by the pipeline entry point.
 */
import java.nio.charset.StandardCharsets
import java.nio.file.Files
import java.nio.file.Path

import nextflow.exception.AbortOperationException

class CLIParamsValidation {

    static final List<String> AVAILABLE_CALLERS =
        ['bcftools', 'freebayes', 'gatk']
    static final String SAMPLE_ID_PATTERN = /^[A-Za-z0-9][A-Za-z0-9_.-]*$/
    static final String READ_GROUP_PATTERN = /^[A-Za-z0-9_.:-]+$/

    // Caller options that select ploidy or copy number.  A true value means
    // that an explicit diploid setting (2) is accepted; a false value means
    // the option is rejected because its contents are not inspected.
    static final Map<String, Map<String, Boolean>> PLOIDY_OPTIONS = [
        bcftools: [
            '--ploidy': true,
            '--ploidy-file': false
        ],
        freebayes: [
            '--ploidy': true,
            '-p': true,
            '--cnv-map': false,
            '-A': false
        ],
        gatk: [
            '--sample-ploidy': true,
            '--ploidy': true,
            '-ploidy': true,
            '--ploidy-regions': false
        ]
    ]

    private static void fail(String message) {
        throw new AbortOperationException("ERROR: ${message}")
    }

    static List<String> parse_callers(String callers) {
        if (callers == null) {
            return []
        }
        // trimmed, non-empty and deduplicated, preserving the configured order
        return callers.split(',')
            .collect { it.trim() }
            .findAll { it }
            .unique()
    }

    /**
     * Resolves the effective duplicate-handling mode.
     */
    static String rmdup_mode(Map rmdup) {
        def mode = rmdup?.mode ?: 'mark'
        if (!(mode in ['none', 'mark', 'remove'])) {
            fail(
                "rmdup.mode: unknown mode '${mode}'. Accepted values are: " +
                'none, mark, remove'
            )
        }
        return mode
    }

    static void effective_callers_validation(String callers) {
        def callersList = parse_callers(callers)
        if (!callersList) {
            fail(
                'calling.callers is required and must name at least one ' +
                "caller (available: ${AVAILABLE_CALLERS.join(', ')})"
            )
        }
        def unknown = callersList.findAll { !(it in AVAILABLE_CALLERS) }
        if (unknown) {
            fail(
                "calling.callers: unknown caller(s) ${unknown.join(', ')}. " +
                "Available callers are: ${AVAILABLE_CALLERS.join(', ')}"
            )
        }
    }

    static void mixed_reads_validation(
        String mapper,
        String reads_type,
        String format
    ) {
        if (format == 'bam') {
            return
        }
        if (reads_type == 'mx' && mapper != 'bowtie2') {
            fail(
                "input.reads_type: 'mx' is only supported with " +
                "mapping.mapper 'bowtie2' (got '${mapper}')"
            )
        }
    }

    /**
     * Rejects read-group values that are unsafe in the mapper command line
     * or in the SAM header.
     *
     * The schema applies the same pattern, but JSON Schema pattern matching
     * can accept a trailing line terminator.  A full-string match keeps
     * values such as `platform: "ILLUMINA\n"` out of the @RG tag.
     */
    static void read_group_validation(Map mapping) {
        def values = [
            'mapping.platform': mapping?.platform ?: 'ILLUMINA',
            'mapping.library': mapping?.library
        ]
        values.each { key, value ->
            if (value != null && !(value ==~ READ_GROUP_PATTERN)) {
                fail(
                    "${key}: '${value}' is not safe for a SAM read group; " +
                    'accepted characters are letters, digits, underscore, ' +
                    'dot, colon and hyphen'
                )
            }
        }
    }

    /**
     * Rejects an index prefix that selects a subdirectory.
     *
     * Index files are resolved from reference_index_dir before staging but
     * staged under a flat `index/` directory.  A prefix containing a path
     * separator would therefore point at the files during validation and at
     * a non-existent nested path in the mapping command.
     */
    static void reference_index_prefix_validation(String prefix) {
        if (prefix == null) {
            return
        }
        if (prefix.contains('/') || prefix.contains('\\')) {
            fail(
                "input.reference_index_prefix: '${prefix}' must not " +
                'contain a path separator'
            )
        }
    }

    /**
     * Rejects caller arguments that request a non-diploid genotype.
     *
     * Consensus reference records are always diploid (0/0), while a caller
     * that runs with another ploidy emits variants with a different GT arity.
     * Only an explicit diploid setting is accepted.  Custom ploidy,
     * copy-number and region files are rejected because their contents are
     * not inspected.
     */
    static void ploidy_validation(Map calling) {
        if (calling == null) {
            return
        }
        def bcftools = calling['bcftools'] ?: [:]
        validate_ploidy_args(
            'calling.bcftools.call.extra_args',
            bcftools['call']?.extra_args,
            PLOIDY_OPTIONS['bcftools']
        )
        def freebayes = calling['freebayes'] ?: [:]
        validate_ploidy_args(
            'calling.freebayes.extra_args',
            freebayes['extra_args'],
            PLOIDY_OPTIONS['freebayes']
        )
        def gatk = calling['gatk'] ?: [:]
        validate_ploidy_args(
            'calling.gatk.extra_args',
            gatk['extra_args'],
            PLOIDY_OPTIONS['gatk']
        )
    }

    private static void validate_ploidy_args(
        String key,
        String args,
        Map<String, Boolean> options
    ) {
        if (args == null || args.trim().isEmpty()) {
            return
        }
        options.each { option, diploid_allowed ->
            def value = option_value(args, option)
            if (value == null) {
                return
            }
            if (diploid_allowed && unquote(value) == '2') {
                return
            }
            def found = value ? "${option} ${value}" : option
            fail(
                "${key}: ploidy is not supported yet except diploid (2); " +
                "remove '${found}'"
            )
        }
    }

    /**
     * Returns the value of an option in a shell-like argument string.
     * Returns an empty string when the option has no value and null when the
     * option is absent.
     */
    private static String option_value(String args, String option) {
        def with_value = args =~ (
            /(?<!\S)/ + option +
            /(?:\s*=\s*|\s+)("[^"]*"|'[^']*'|[^\s"']+)/
        )
        if (with_value.find()) {
            return with_value.group(1)
        }
        if (option.length() == 2 && option.startsWith('-')) {
            def attached = args =~ /(?<!\S)/ + option + /(\S+)/
            if (attached.find()) {
                return attached.group(1)
            }
        }
        def present = args =~ /(?<!\S)/ + option + /(?=\s|=|$)/
        return present.find() ? '' : null
    }

    private static String unquote(String value) {
        return value.replaceAll(/^['"]|['"]$/, '')
    }

    static void cons_threshold_validation(
        Number cons_threshold,
        String callers
    ) {
        if (cons_threshold == null) {
            fail('consensus.threshold is required')
        }
        def callersList = parse_callers(callers)
        if (!callersList) {
            fail(
                'consensus.threshold: calling.callers must name at least ' +
                'one caller'
            )
        }
        def threshold
        try {
            threshold = cons_threshold as Integer
        } catch (NumberFormatException ignored) {
            fail(
                'consensus.threshold must be an integer ' +
                "(got '${cons_threshold}')"
            )
        }
        if (threshold < 1) {
            fail(
                'consensus.threshold must be at least 1 ' +
                "(got ${threshold})"
            )
        }
        if (threshold > callersList.size()) {
            fail(
                "consensus.threshold (${threshold}) must be less than or " +
                'equal to the number of configured callers ' +
                "(${callersList.size()})"
            )
        }
    }

    /**
     * Reads and validates the samplesheet before any task runs.
     *
     * Blank lines and lines starting with '#' are ignored.  A header line whose
     * first column is sample, sample_id, sampleid or id is skipped.  The row
     * width depends on input.format and input.reads_type.  Relative read paths
     * are resolved from the samplesheet directory and must exist.
     *
     * @return list of maps with sample id and resolved read paths
     */
    static List<Map<String, Object>> parse_samplesheet(
        Path samplesheet,
        String format,
        String reads_type
    ) {
        File sheet = samplesheet.toFile()
        if (!sheet.isFile()) {
            fail("input.samples_tsv: file not found: ${sheet}")
        }
        int expected = expected_columns(format, reads_type)
        File base_dir = sheet.absoluteFile.parentFile

        def samples = []
        def seen = new LinkedHashMap<String, Integer>()
        int line_number = 0

        Files.newBufferedReader(sheet.toPath(), StandardCharsets.UTF_8)
            .withCloseable { reader ->
                String line
                while ((line = reader.readLine()) != null) {
                    line_number++
                    if (line_number == 1 && line.startsWith('\uFEFF')) {
                        line = line.substring(1)
                    }
                    if (line.trim().isEmpty() ||
                        line.trim().startsWith('#')) {
                        continue
                    }
                    def columns = line.split('\t', -1)
                        .collect { it.replaceAll(/\r$/, '') }
                    if (is_header_row(columns)) {
                        continue
                    }
                    if (columns.size() != expected) {
                        fail(
                            "input.samples_tsv line ${line_number}: " +
                            "expected ${expected} tab-separated columns " +
                            "for format='${format}', " +
                            "reads_type='${reads_type}' but found " +
                            "${columns.size()}"
                        )
                    }

                    String sample_id = columns[0]
                    if (!sample_id) {
                        fail(
                            "input.samples_tsv line ${line_number}: " +
                            'sample id is empty'
                        )
                    }
                    if (!(sample_id ==~ SAMPLE_ID_PATTERN)) {
                        fail(
                            "input.samples_tsv line ${line_number}: " +
                            "sample id '${sample_id}' is not shell-safe. " +
                            "Use letters, digits, '.', '_' or '-' and " +
                            'start with a letter or digit.'
                        )
                    }
                    if (seen.containsKey(sample_id)) {
                        fail(
                            "input.samples_tsv line ${line_number}: " +
                            "duplicate sample id '${sample_id}' " +
                            "(first seen on line ${seen[sample_id]})"
                        )
                    }
                    seen[sample_id] = line_number

                    def reads = (1..<expected).collect { index ->
                        resolve_read_path(
                            base_dir,
                            columns[index],
                            line_number,
                            index + 1
                        )
                    }
                    samples << [id: sample_id, reads: reads]
                }
            }

        if (samples.isEmpty()) {
            fail("input.samples_tsv: no samples found in ${sheet}")
        }
        return samples
    }

    private static int expected_columns(String format, String reads_type) {
        if (format == 'bam') {
            return 2
        }
        switch (reads_type) {
            case 'se':
                return 2
            case 'pe':
                return 3
            case 'mx':
                return 4
            default:
                fail(
                    "input.reads_type: unsupported value '${reads_type}' " +
                    '(accepted: se, pe, mx)'
                )
        }
    }

    private static boolean is_header_row(List<String> columns) {
        def first = columns[0].trim().toLowerCase()
        return first in ['sample', 'sample_id', 'sampleid', 'id']
    }

    private static String resolve_read_path(
        File base_dir,
        String value,
        int line_number,
        int column_number
    ) {
        if (value == null || value.trim().isEmpty()) {
            fail(
                "input.samples_tsv line ${line_number}: " +
                "column ${column_number} is empty"
            )
        }
        def raw = value.trim()
        // Leave remote URIs to Nextflow staging; local paths are checked here.
        if (raw ==~ /^[A-Za-z][A-Za-z0-9+.-]*:\/\/.*/) {
            return raw
        }
        def candidate = new File(raw)
        if (!candidate.isAbsolute()) {
            candidate = new File(base_dir, raw)
        }
        def resolved = candidate.canonicalFile
        if (!resolved.isFile()) {
            fail(
                "input.samples_tsv line ${line_number}: " +
                "column ${column_number} file not found: ${resolved}"
            )
        }
        return resolved.absolutePath
    }

    /**
     * Detects BGZF from the gzip header.  A plain gzip stream has the FEXTRA
     * flag unset; BGZF carries an extra field with the 'BC' subfield.
     */
    static boolean is_bgzip(Path fasta) {
        try {
            return Files.newInputStream(fasta).withCloseable { stream ->
                byte[] header = new byte[12]
                int read = stream.readNBytes(header, 0, header.length)
                if (read < header.length) {
                    return false
                }
                if ((header[0] & 0xff) != 0x1f ||
                    (header[1] & 0xff) != 0x8b) {
                    return false
                }
                int flags = header[3] & 0xff
                if ((flags & 0x04) == 0) {
                    return false
                }
                int xlen = (header[10] & 0xff) |
                    ((header[11] & 0xff) << 8)
                if (xlen == 0) {
                    return false
                }
                byte[] extra = new byte[xlen]
                int extra_read = stream.readNBytes(extra, 0, xlen)
                int i = 0
                while (i + 3 < extra_read) {
                    if ((extra[i] & 0xff) == 0x42 &&
                        (extra[i + 1] & 0xff) == 0x43) {
                        int slen = (extra[i + 2] & 0xff) |
                            ((extra[i + 3] & 0xff) << 8)
                        return slen == 2
                    }
                    int slen = (extra[i + 2] & 0xff) |
                        ((extra[i + 3] & 0xff) << 8)
                    i += 4 + slen
                }
                return false
            }
        } catch (IOException exception) {
            fail(
                "input.reference_genome: cannot read '${fasta}': " +
                exception.message
            )
        }
        return false
    }

    /**
     * Detects plain gzip from the magic bytes, ignoring the extra fields.
     * A gzip-compressed reference that is not BGZF cannot be indexed by
     * samtools faidx, so the workflow rejects it before any task starts.
     */
    static boolean is_gzip(Path fasta) {
        try {
            return Files.newInputStream(fasta).withCloseable { stream ->
                byte[] magic = new byte[2]
                int read = stream.readNBytes(magic, 0, magic.length)
                return read == magic.length &&
                    (magic[0] & 0xff) == 0x1f &&
                    (magic[1] & 0xff) == 0x8b
            }
        } catch (IOException exception) {
            fail(
                "input.reference_genome: cannot read '${fasta}': " +
                exception.message
            )
        }
        return false
    }

    /**
     * Fails early when the merge step cannot open the descriptors it needs.
     * Each chunk opens its BCFs plus their indexes; the final merge opens one
     * named pipe per chunk.  A small fixed margin covers stdio and htslib
     * temporary handles.
     */
    static void merge_ulimit_validation(int sample_count, int cpus) {
        if (sample_count < 1) {
            return
        }
        int chunk_count = Math.min(Math.max(cpus, 1), sample_count)
        long largest_chunk = (long) Math.ceil(
            sample_count / (double) chunk_count
        )
        long required = Math.max(
            largest_chunk * 2L,
            (long) chunk_count
        ) + 64L
        long limit = open_file_limit()
        if (limit > 0 && limit < required) {
            fail(
                'The merge step needs at least ' +
                "${required} open file descriptors for " +
                "${sample_count} sample(s) with " +
                "consensus.cpu=${cpus}, but the current soft limit is " +
                "${limit}. Raise it with 'ulimit -n ${required}' or " +
                'lower consensus.cpu.'
            )
        }
    }

    /**
     * Returns the soft open-file limit of the current process, or -1 when it
     * cannot be determined.
     */
    static long open_file_limit() {
        def limits = new File('/proc/self/limits')
        if (limits.isFile()) {
            try {
                def line = limits.readLines('UTF-8').find {
                    it.trim().startsWith('Max open files')
                }
                if (line) {
                    def parts = line.trim().split(/\s+/)
                    int files_index = -1
                    for (int i = 0; i < parts.length; i++) {
                        if (parts[i] == 'files') {
                            files_index = i
                        }
                    }
                    if (files_index >= 3) {
                        return parts[files_index - 2] as long
                    }
                }
            } catch (Exception ignored) {
                // fall back to ulimit below
            }
        }
        try {
            def process = ['sh', '-c', 'ulimit -n'].execute()
            process.waitFor()
            def text = process.text.trim()
            if (process.exitValue() == 0 && text ==~ /\d+/) {
                return text as long
            }
        } catch (Exception ignored) {
            // unknown limit: do not block the run
        }
        return -1L
    }
}
