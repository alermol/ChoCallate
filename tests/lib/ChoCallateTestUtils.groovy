import java.nio.file.Path

/**
 * Shared fixtures and assertions for the ChoCallate nf-test suite.
 */
class ChoCallateTestUtils {

    static final String REFERENCE_GZ = 'test_run/arth_chr1.fasta.gz'
    static final String READ1 = 'test_run/test_reads_R1.fq.gz'
    static final String READ2 = 'test_run/test_reads_R2.fq.gz'
    static final String CONTIG = 'CP002684.1'
    static final int CONTIG_LENGTH = 200001

    static File projectFile(Object projectDir, String relative) {
        return new File(projectDir.toString(), relative)
    }

    static File outputFile(Object path) {
        return new File(path.toString())
    }

    /** Runs a shell command and fails with its output on a non-zero exit. */
    static String run(File workDir, String command) {
        def builder = new ProcessBuilder('bash', '-c', command)
        builder.directory(workDir)
        builder.redirectErrorStream(true)
        def process = builder.start()
        String output = process.inputStream.getText('UTF-8')
        int status = process.waitFor()
        if (status != 0) {
            throw new RuntimeException(
                "Command failed with exit ${status}: ${command}\n${output}"
            )
        }
        return output
    }

    /**
     * Decompresses the test reference and creates the fai index and sequence
     * dictionary used by the mapping and calling modules.
     */
    static Map referenceBundle(Object projectDir, File launchDir,
                               String name = 'reference') {
        File dir = new File(launchDir, name)
        dir.mkdirs()
        File fasta = new File(dir, 'genome.fasta')
        if (!fasta.isFile()) {
            run(
                dir,
                "gzip -dc '${projectFile(projectDir, REFERENCE_GZ)}' " +
                "> '${fasta}'"
            )
            run(dir, "samtools faidx '${fasta}'")
            run(
                dir,
                "gatk CreateSequenceDictionary -R '${fasta}' " +
                "-O '${dir}/genome.dict'"
            )
        }
        return [
            dir: dir,
            fasta: fasta,
            fai: new File(dir, 'genome.fasta.fai'),
            dict: new File(dir, 'genome.dict'),
            prefix: 'genome.fasta'
        ]
    }

    static Map bowtie2Index(Object projectDir, File launchDir,
                            String name = 'reference') {
        Map ref = referenceBundle(projectDir, launchDir, name)
        run(
            ref.dir,
            "bowtie2-build --quiet --threads 2 '${ref.fasta}' " +
            "'${ref.dir}/${ref.prefix}'"
        )
        return ref
    }

    static Map bwaIndex(Object projectDir, File launchDir,
                        String name = 'reference') {
        Map ref = referenceBundle(projectDir, launchDir, name)
        run(ref.dir, "bwa index '${ref.fasta}'")
        return ref
    }

    static Map minimap2Index(Object projectDir, File launchDir,
                             String name = 'reference') {
        Map ref = referenceBundle(projectDir, launchDir, name)
        run(
            ref.dir,
            "minimap2 -x sr -t 2 -d '${ref.dir}/genome.fasta.mmi' " +
            "'${ref.fasta}'"
        )
        return ref
    }

    /** Maps the test reads with bowtie2 and returns a sorted, indexed BAM. */
    static File bowtie2Bam(Object projectDir, File launchDir,
                           String sampleId = 'sample1',
                           boolean singleEnd = false) {
        Map ref = bowtie2Index(projectDir, launchDir)
        File bam = new File(launchDir, "${sampleId}.bam")
        File read1 = projectFile(projectDir, READ1)
        File read2 = projectFile(projectDir, READ2)
        String reads = singleEnd
            ? "-U '${read1}'"
            : "-1 '${read1}' -2 '${read2}'"
        run(
            launchDir,
            "set -o pipefail; bowtie2 --quiet --threads 2 " +
            "--rg-id '${sampleId}' --rg 'SM:${sampleId}' " +
            "--rg 'PL:ILLUMINA' --rg 'LB:${sampleId}' " +
            "-x '${ref.dir}/${ref.prefix}' ${reads} | " +
            "samtools sort --threads 2 -o '${bam}' - && " +
            "samtools index --csi '${bam}'"
        )
        return bam
    }

    /** Writes a coordinate-sorted, indexed BAM with no alignments. */
    static File emptyBam(File launchDir, String sampleId = 'empty') {
        File bam = new File(launchDir, "${sampleId}.bam")
        run(
            launchDir,
            "printf '@HD\\tVN:1.6\\tSO:coordinate\\n" +
            "@SQ\\tSN:${CONTIG}\\tLN:${CONTIG_LENGTH}\\n' | " +
            "samtools view -b -o '${bam}' - && " +
            "samtools index --csi '${bam}'"
        )
        return bam
    }

    static File samplesheet(File dir, List<List<String>> rows) {
        File sheet = new File(dir, 'samples.tsv')
        sheet.text = rows.collect { row -> row.join('\t') }.join('\n') + '\n'
        return sheet
    }

    private static File writeVcfText(File dir, String name, String sample,
                                     List<List<String>> records) {
        dir.mkdirs()
        def header = [
            '##fileformat=VCFv4.2',
            '##contig=<ID=CP002684.1,length=200001>',
            '##ALT=<ID=*,Description="Represents allele(s) other than ' +
                'observed.">',
            '##INFO=<ID=END,Number=1,Type=Integer,' +
                'Description="Stop position of the interval">',
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
            '##FORMAT=<ID=AD,Number=R,Type=Integer,' +
                'Description="Allelic depths">',
            '##FORMAT=<ID=DP,Number=1,Type=Integer,' +
                'Description="Read depth">',
            '##ChoCallate=<Version=test,Threshold=1,Callers="bcftools",' +
                'Rule=two-stage,Tie=nocall,AbsentVote=ref,' +
                'FilteredVote=abstain,ReferenceBlocks=False>',
            '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\t' +
                "FORMAT\t${sample}"
        ]
        File vcf = new File(dir, "${name}.vcf")
        vcf.text = (header + records.collect { it.join('\t') })
            .join('\n') + '\n'
        return vcf
    }

    /** Writes a minimal single-sample BCF against the test contig. */
    static File writeBcf(File dir, String name, String sample,
                         List<List<String>> records) {
        File vcf = writeVcfText(dir, name, sample, records)
        File bcf = new File(dir, "${name}.bcf")
        run(dir, "bcftools view -Ob -o '${bcf}' '${vcf}'")
        return bcf
    }

    /** Writes a minimal bgzipped single-sample VCF against the test contig. */
    static File writeVcfGz(File dir, String name, String sample,
                           List<List<String>> records) {
        File vcf = writeVcfText(dir, name, sample, records)
        File gz = new File(dir, "${name}.vcf.gz")
        run(dir, "bgzip -c '${vcf}' > '${gz}'")
        return gz
    }

    static String bamHeader(File bam) {
        return run(
            bam.absoluteFile.parentFile,
            "samtools view -H '${bam.absolutePath}'"
        )
    }

    static int bamRecords(File bam) {
        String count = run(
            bam.absoluteFile.parentFile,
            "samtools view -c -F 4 '${bam.absolutePath}'"
        )
        return count.trim() as int
    }

    static String vcfHeader(File vcf) {
        return run(
            vcf.absoluteFile.parentFile,
            "bcftools view -h '${vcf.absolutePath}'"
        )
    }

    static List<String> vcfSamples(File vcf) {
        String chrom = vcfHeader(vcf).readLines().find { line ->
            line.startsWith('#CHROM')
        }
        return chrom.split('\t', -1).drop(9).collect { it.trim() }
    }

    /** Returns CHROM, POS, REF and ALT for every record in the VCF/BCF. */
    static List<String> vcfRecords(File vcf) {
        String output = run(
            vcf.absoluteFile.parentFile,
            "bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\n' " +
            "'${vcf.absolutePath}'"
        )
        return output.readLines()
    }

    static List<String> statsLines(File vcf) {
        String output = run(
            vcf.absoluteFile.parentFile,
            "bcftools stats '${vcf.absolutePath}' | " +
            "grep -E '^(SN|TSTV)[[:space:]]'"
        )
        return output.readLines()
    }

    static int recordCount(File vcf) {
        String line = statsLines(vcf).find { entry ->
            entry.startsWith('SN\t0\tnumber of records:')
        }
        if (line == null) {
            throw new IllegalStateException("No record count for ${vcf}")
        }
        return line.split('\t', -1)[3] as int
    }

    static int statsValue(File vcf, String key) {
        String line = statsLines(vcf).find { entry ->
            entry.startsWith("SN\t0\t${key}:")
        }
        if (line == null) {
            throw new IllegalStateException("No ${key} for ${vcf}")
        }
        return line.split('\t', -1)[3] as int
    }
}
