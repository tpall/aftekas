process FIX_FASTA_HEADERS {
    tag "$meta.id"
    label 'process_small'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5a/5aae5977ff9de3e01ff962dc495bfa23f4304c676446b5fdf2de5c7edfa2dc4e/data' :
        'community.wave.seqera.io/library/bbmap_pigz:07416fe99b090fa9' }"

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("*fixed*"), emit: fixed
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def base = input.name.matches(/.*gz$/) ? input.getBaseName(2) : input.getBaseName()
    def ext = input.name.minus(base)
    output = "${base}_fixed${ext}"

    """
    reformat.sh \\
        -Xmx${task.memory.toGiga()}g \\
        in=${input} \\
        out=${output} \\
        trd=t \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bbmap: \$(bbversion.sh | grep -v "Duplicate cpuset")
        pigz: \$( pigz --version 2>&1 | sed 's/pigz //g' )
    END_VERSIONS
    """
}
