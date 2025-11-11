
process BINETTE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/checkm2_diamond_pip_binette:d2744f87fb84aab7"

    input:
    tuple val(meta), path(contigs2bin), path(bins), path(contigs), path(proteins)
    path(checkm2_db)

    output:
    tuple val(meta), path("final_bins/*.fa.gz"), emit: bins, optional: true
    tuple val(meta), path("*_final_bins_quality_reports.tsv"), emit: quality_report
    tuple val(meta), path("*_final_contig_to_bin.tsv"), emit: contig_to_bin
    tuple val(meta), path("input_bins_quality_reports/*.tsv"), emit: input_quality_report


    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}" 
    def bin_defs = contigs2bin ? "-b ${contigs2bin.join(' ')}" : "-d ${bins.join(' ')}"
    if (bins) {
        assert bins[0].isDirectory() : "Bins must be directory with bin fasta files."
    }
    def prots = proteins ? "-p ${proteins}" : ""

    """
    binette \\
        $args \\
        -c $contigs \\
        $bin_defs \\
        $prots \\
        -t $task.cpus \\
        --prefix $prefix \\
        --checkm2_db $checkm2_db

    for i in results/final_bins/${prefix}*.fa; do
        gzip \${i}
    done

    for i in results/*.tsv; do
        suffix=\$(basename \${i})
        mv \${i} results/${prefix}_\${suffix}
    done

    mv results/* .
    rm -rf results

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Binette: \$(binette --version |& sed '1!d ; s/Binette //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p reults
    touch results/${prefix}_bin1.fa.gz
    touch results/${prefix}_bin2.fa.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Binette: \$(binette --version |& sed '1!d ; s/Binette //')
    END_VERSIONS
    """
}
