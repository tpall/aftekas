// Sketch one sample's reads into a sylph .sylsp file. One task per sample so
// the cluster profile can run them as a SLURM array like the rest of PREPROCESS.
// Input must be the host/phiX-removed reads, never the BBNorm output: normalisation
// flattens coverage on purpose, and sylph's abundances are coverage estimates.
process SYLPH_SKETCH {
    tag "${meta.id}"
    label 'process_small'

    conda "${moduleDir}/environment.yml"
    container "quay.io/biocontainers/sylph:1.0.0--hb42e459_0"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.sylsp"), emit: sketch
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // paired sketches are written as <name>.paired.sylsp, single-end as <name>.sylsp
    def input  = meta.single_end ? "${reads}" : "-1 ${reads[0]} -2 ${reads[1]}"
    """
    sylph sketch \\
        ${args} \\
        -t ${task.cpus} \\
        -S ${prefix} \\
        -d . \\
        ${input}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sylph: \$(sylph --version 2>&1 | sed 's/^sylph //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.paired.sylsp
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sylph: \$(sylph --version 2>/dev/null | sed 's/^sylph //')
    END_VERSIONS
    """
}
