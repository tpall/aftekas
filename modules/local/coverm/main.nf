// Abundance of every genome in the client's genome set, for one sample, from mapping reads to the
// set. This is the only source of abundance in a series run: a sample's own assembly is never the
// denominator, because two samples' assemblies are different sequences.
//
// Three methods on purpose:
//   relative_abundance  what the report prints (share of mapped reads)
//   mean                coverage, for the "is this genome really here" judgement
//   count               read counts, which the monitoring noise floor needs to know when a route is
//                       too thinly covered to call a change (Poisson limit)
process COVERM_GENOME {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/coverm:0.8.0--h750ce8b_0' :
        'quay.io/biocontainers/coverm:0.8.0--h750ce8b_0' }"

    input:
    tuple val(meta), path(reads)
    path(representatives, stageAs: 'representatives/*')
    val(set_version)

    output:
    tuple val(meta), path("${meta.id}.coverm.tsv"), emit: abundance
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '--min-read-percent-identity 95 --min-read-aligned-percent 75'
    def input_reads = meta.single_end ? "--single ${reads}" : "--coupled ${reads}"
    """
    coverm genome \\
        ${args} \\
        ${input_reads} \\
        --genome-fasta-directory representatives/ \\
        --genome-fasta-extension fa \\
        --methods relative_abundance mean count \\
        --threads ${task.cpus} \\
        --output-file ${meta.id}.coverm.tsv

    # The set version travels with the numbers: an abundance table is only meaningful against the
    # genome set it was mapped to, and the report refuses to compare two versions (design doc).
    sed -i "1s/\$/\\tgenome_set_version/; 2,\$ s/\$/\\t${set_version}/" ${meta.id}.coverm.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coverm: \$(coverm --version 2>&1 | sed 's/^coverm //')
        genome_set_version: ${set_version}
    END_VERSIONS
    """

    stub:
    """
    printf 'Genome\\tRelative Abundance (%%)\\tMean\\tRead Count\\tgenome_set_version\\n' > ${meta.id}.coverm.tsv
    printf 'CGS-0001\\t12.5\\t8.1\\t120000\\t${set_version}\\n' >> ${meta.id}.coverm.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coverm: 0.8.0
        genome_set_version: ${set_version}
    END_VERSIONS
    """
}
