// Build or extend the client's genome set from the cluster definition. All identity rules live in
// bin/genome_set.py: stable ids, no replacement of a representative, no silent merge of two ids.
process GENOME_SET_UPDATE {
    tag "${set_id} v+1"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.10' :
        'quay.io/biocontainers/python:3.10' }"

    input:
    val(set_id)
    val(sample_label)
    val(ani)
    path(clusters)
    path(mags, stageAs: 'mags/*')
    path(previous)          // previous version directory, or [] for a first build
    path(tool_versions)

    output:
    path "genome_set"                        , emit: set_dir
    path "genome_set/representatives/*.fa"   , emit: representatives
    path "genome_set/set.json"               , emit: set_json
    path "genome_set/members.tsv"            , emit: members
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prev = previous ? "--previous ${previous}" : ''
    """
    genome_set.py \\
        --clusters ${clusters} \\
        --set-id ${set_id} \\
        --ani ${ani} \\
        --sample ${sample_label} \\
        --fasta-dir mags \\
        ${prev} \\
        --tool-versions "\$(cat ${tool_versions} | tr '\\n' ';')" \\
        --out genome_set

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/^Python //')
        genome_set_version: \$(python3 -c "import json;print(json.load(open('genome_set/set.json'))['version'])")
    END_VERSIONS
    """

    stub:
    """
    mkdir -p genome_set/representatives
    printf '>c\\nACGT\\n' > genome_set/representatives/CGS-0001.fa
    printf 'genome_id\\tfasta\\tsample\\tversion_added\\tis_representative\\nCGS-0001\\tx.fa\\t${sample_label}\\t1\\tyes\\n' > genome_set/members.tsv
    printf '{"set_id": "${set_id}", "version": 1, "genomes": 1, "ani": ${ani}}' > genome_set/set.json
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.10
        genome_set_version: 1
    END_VERSIONS
    """
}
