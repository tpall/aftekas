process CONCOCT {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/concoct:1.1.0--py39h8907335_8'
        : 'quay.io/biocontainers/concoct:1.1.0--py39h8907335_8'}"

    input:
    tuple val(meta), path(contigs), path(bam), path(bai)

    output:
    tuple val(meta), path("${prefix}*.fa.gz"), emit: fasta
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    assert contigs.extension == "gz" : "Input contigs file must be gzipped (.gz)"
    contigs_fa = contigs.baseName

    """
    gunzip -c ${contigs} > ${contigs_fa}
    cut_up_fasta.py \\
        ${contigs_fa} \\
        -c 10000 -o 0 --merge_last \\
        -b ${prefix}.contigs_10K.bed \\
        > ${prefix}.contigs_10K.fa
    concoct_coverage_table.py \\
        ${prefix}.contigs_10K.bed \\
        ${bam} \\
        > ${prefix}.coverage_table.tsv
    concoct \\
        ${args} \\
        --threads ${task.cpus} \\
        --coverage_file ${prefix}.coverage_table.tsv \\
        --composition_file ${prefix}.contigs_10K.fa \\
        -b ${prefix}
    merge_cutup_clustering.py ${prefix}_clustering_gt1000.csv > ${prefix}_clustering_merged.csv
    mkdir -p ${prefix}
    extract_fasta_bins.py \\
        ${contigs_fa} \\
        ${prefix}_clustering_merged.csv \\
        --output_path ${prefix}
    
    for i in ${prefix}/*.fa; do
        mv \${i} \${i/\\///${prefix}.}
        gzip \${i/\\///${prefix}.}
    done

    mv ${prefix}/*.fa.gz ./
    rm -r ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        concoct: \$(echo \$(concoct --version 2> /dev/null) | sed 's/concoct //g' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_args.txt
    touch ${prefix}_clustering_gt1000.csv
    touch ${prefix}_log.txt
    touch ${prefix}_original_data_gt1000.csv
    touch ${prefix}_PCA_components_data_gt1000.csv
    touch ${prefix}_PCA_transformed_data_gt1000.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        concoct: \$(echo \$(concoct --version 2> /dev/null) | sed 's/concoct //g' )
    END_VERSIONS
    """
}
