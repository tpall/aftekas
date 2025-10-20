
process DATASETS_DOWNLOAD {

    tag "${accession}"

    container "staphb/ncbi-datasets"

    input:
    tuple val(meta), val(accession)

    output:
    tuple val(meta), path("*.fna"), emit: fasta
    path "versions.yml", emit: versions

    script:
    """
    datasets download virus genome accession ${accession} --filename ${accession}.zip --include genome
    unzip ${accession}.zip -d ${accession}_unzipped
    mv ${accession}_unzipped/ncbi_dataset/data/genomic.fna ./${accession}.fna
    rm -rf ${accession}.zip ${accession}_unzipped
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ncbi-datasets: \$(datasets --version 2>&1 | sed -e "s/datasets version//g")
    END_VERSIONS
    """
}
