process VERSIONS {
    publishDir "${params.output.directory}", mode: 'copy',
        pattern: 'versions.yml'

    input:
    val tools

    output:
    path("versions.yml"), emit: versions

    script:
    """
    write_versions.sh \
        --pipeline-version "${workflow.manifest.version}" \
        versions.yml ${tools.join(' ')}
    """
}
