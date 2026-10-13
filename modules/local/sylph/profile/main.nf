// Profile every sketch against one database in a single call: sylph loads the
// database once (GTDB r226 -c 200 .syl2db is ~14 GB on disk) and writes one
// table for all samples, which is also what sylph-tax wants downstream.
process SYLPH_PROFILE {
    tag "${db.name}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "quay.io/biocontainers/sylph:1.0.0--hb42e459_0"

    input:
    path sketches
    path db

    output:
    path "sylph_profile.tsv", emit: profile
    path "versions.yml"     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    sylph profile \\
        ${args} \\
        -t ${task.cpus} \\
        ${db} \\
        ${sketches} \\
        -o sylph_profile.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sylph: \$(sylph --version 2>&1 | sed 's/^sylph //')
        sylph_db: ${db.name}
    END_VERSIONS
    """

    stub:
    """
    printf 'Sample_file\\tGenome_file\\tTaxonomic_abundance\\tSequence_abundance\\tAdjusted_ANI\\tEff_cov\\tANI_5-95_percentile\\tEff_lambda\\tLambda_5-95_percentile\\tMedian_cov\\tMean_cov_geq1\\tContainment_ind\\tNaive_ANI\\tContig_name\\n' > sylph_profile.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sylph: \$(sylph --version 2>/dev/null | sed 's/^sylph //')
        sylph_db: ${db.name}
    END_VERSIONS
    """
}
