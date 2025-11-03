process MAGSCOT {
    tag "$contigs.baseName"

    conda "${moduleDir}/environment.yml"
    container "ikmb/magscot:v1.1"

    input:
    tuple val(meta), path(sets), path(hmm)

    output:
    path "magscot_bins/"      , emit: bins
    path "magscot_report.txt" , emit: report
    path "versions.yml"       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    magscot \\
        -i $contigs \\
        -d $depth_stats \\
        -o magscot_bins/ \\
        -t ${task.cpus} \\
        --report magscot_report.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        magscot: \$(magscot --version | awk '{print \$2}')
    END_VERSIONS
    """
}