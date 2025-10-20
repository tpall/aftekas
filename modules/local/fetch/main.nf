
process FETCH_INDEX {

    tag "${url.baseName}"

    input:
    tuple val(meta), path(url)

    output:
    tuple val(meta), path('*.{bt2,bt2l}'), emit: index

    script:
    """ 
    mkdir bowtie2
    tar -xf ${url}
    rm *.tar
    """
}

