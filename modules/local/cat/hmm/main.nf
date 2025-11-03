process CAT_HMM {
    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/52/52ccce28d2ab928ab862e25aae26314d69c8e38bd41ca9431c67ef05221348aa/data'
        : 'community.wave.seqera.io/library/coreutils_grep_gzip_lbzip2_pruned:838ba80435a629f8'}"

    input:
    tuple val(meta), path(tigr), path(pfam)

    output:
    tuple val(meta), path('*.hmm'), emit: hmm

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    zcat $tigr | grep -v "^#" | awk '{print \$1"\t"\$3"\t"\$5}' > ${prefix}.tigr
    zcat $pfam | grep -v "^#" | awk '{print \$1"\t"\$4"\t"\$5}' > ${prefix}.pfam
    cat ${prefix}.pfam ${prefix}.tigr > ${prefix}.hmm
    """
}