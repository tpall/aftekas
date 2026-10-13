// Prepare one sample's refined MAGs for cross-sample clustering.
//
// Two jobs, both boring and both necessary:
//   1. bin names collide between samples (every sample has a "bin_1"), so every fasta is renamed
//      <sample>__<bin>.fa. Without this the genome set would silently merge two different genomes.
//   2. Binette's quality report is rewritten in CheckM2 column names, under the new fasta names, so
//      galah can pick representatives on completeness and contamination without re-scoring anything.
process PREP_MAGS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.10' :
        'quay.io/biocontainers/python:3.10' }"

    input:
    tuple val(meta), path(bins, stageAs: 'bins/*'), path(quality)

    output:
    path "prepped/*.fa"       , emit: mags
    path "quality_rows.tsv"   , emit: quality
    path "versions.yml"       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    prep_mags.py \\
        --sample ${meta.id} \\
        --bins bins \\
        --quality ${quality} \\
        --out prepped \\
        --quality-out quality_rows.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/^Python //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p prepped
    printf '>c1\\nACGT\\n' > prepped/${meta.id}__bin_1.fa
    printf 'Name\\tCompleteness\\tContamination\\n${meta.id}__bin_1\\t95.0\\t1.0\\n' > quality_rows.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.10
    END_VERSIONS
    """
}
