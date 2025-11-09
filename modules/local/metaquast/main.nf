process METAQUAST {
    tag "${meta.assembler}_${meta.id}"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/quast:5.0.2--py37pl526hb5aa323_2' :
        'quay.io/biocontainers/quast:5.0.2--py37pl526hb5aa323_2' }"

    input:
    tuple val(meta), path(assembly)

    output:
    path "QUAST/*"                       , emit: qc
    path "QUAST/*report.tsv"             , emit: report
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.assembler}_${meta.id}"
    
    """
    metaquast.py ${args} -o QUAST --threads ${task.cpus} -l ${prefix} ${assembly}
    mv QUAST/report.tsv QUAST/${prefix}.report.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        metaquast: \$(metaquast.py --version | sed "s/QUAST v//; s/ (MetaQUAST mode)//")
    END_VERSIONS
    """
}
