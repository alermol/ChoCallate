workflow GENERATE_SAMPLE_CHANNEL {
    take:
    samples_tsv
    input_format
    reads_type

    main:
    if (input_format == 'bam') {
        sample_run_ch = channel
            .fromPath(samples_tsv)
            .splitCsv(header: false, sep: '\t')
            .map { row ->
                def (sample, bam_file) = row
                tuple(sample, bam_file)
            }
    }
    if (input_format == 'fastq') {
        if (reads_type == 'pe') {
            sample_run_ch = channel
                .fromPath(samples_tsv)
                .splitCsv(header: false, sep: '\t')
                .map { row ->
                    def (sample, read1, read2) = row
                    tuple(sample, read1, read2)
                }
        }
        if (reads_type == 'se') {
            sample_run_ch = channel
                .fromPath(samples_tsv)
                .splitCsv(header: false, sep: '\t')
                .map { row ->
                    def (sample, read1) = row
                    tuple(sample, read1)
                }
        }
        if (reads_type == 'mx') {
            sample_run_ch = channel
                .fromPath(samples_tsv)
                .splitCsv(header: false, sep: '\t')
                .map { row ->
                    def (sample, read1, read2, read3) = row
                    tuple(sample, read1, read2, read3)
                }
        }
    }
    
    emit: sample_run_ch
}
