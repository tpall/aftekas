process CONTIG_TO_BIN {
    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/52/52ccce28d2ab928ab862e25aae26314d69c8e38bd41ca9431c67ef05221348aa/data'
        : 'community.wave.seqera.io/library/coreutils_grep_gzip_lbzip2_pruned:838ba80435a629f8'}"

    input:
    tuple val(meta), path(bins)
    val(tool)

    output:
    tuple val(meta), path('*.contigs_to_bin.tsv'), emit: contigs_to_bin

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    if (tool.toLowerCase() == "magscot")
        """
        for i in ${bins}; do
            binname=\$(echo \$(basename \$i) | sed "s/\\.fa\\.gz//g")
            binner=\$(echo \$binname | awk -F_ '{print \$2}')
            zcat \$i | grep ">" | perl -pe "s/\n/\t\${binname}\n/g" | perl -pe "s/>//g" | awk -v binner=\$binner '{print \$2"\t"\$1"\t"binner}' >> ${prefix}.contigs_to_bin.tsv
        done
        """
    else if (tool.toLowerCase() == "binette")
        """
        for i in ${bins}; do
            binname=\$(echo \$(basename \$i) | sed "s/\\.fa\\.gz//g")
            bins=\$(echo \$binname | awk -F. '{print \$1}')
            zcat \$i | grep ">" | perl -pe "s/\n/\t\${binname}\n/g" | perl -pe "s/>//g" | awk '{print \$1"\t"\$2}' >> \${bins}.contigs_to_bin.tsv
        done
        """
    else 
        println "Ups! Unkown bin refinement tool!\nTool can be either\n'magscot' with 'bin contig binner' cols in one file or \n'binette' with 'contig bin' cols for each binner in separate files" 

}