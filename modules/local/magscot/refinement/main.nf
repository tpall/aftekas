process MAGSCOT {
    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "ikmb/magscot:v1.1"

    input:
    tuple val(meta), path(contigs_to_bin), path(hmm)

    output:
    path "magscot_bins/"      , emit: bins
    path "magscot_report.txt" , emit: report
    path "versions.yml"       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    Rscript scripts/MAGScoT.R -i ${contigs_to_bin} --hmm ${hmm} --out ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        magscot: v1.1
    END_VERSIONS
    """
}