// Cluster MAGs into a species-level genome set. galah is dRep's fast equivalent and reads the
// CheckM2 report Binette already produced, so no genome is re-scored here.
//
// This process only clusters. It does not decide what enters the client's genome set: that decision
// keeps representative ids stable across versions and lives in GENOME_SET, which reads the cluster
// definition written here. See docs/design-common-reference-set.md.
process GALAH_CLUSTER {
    tag "${set_id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/galah:0.5.2--hfa8f182_0' :
        'quay.io/biocontainers/galah:0.5.2--hfa8f182_0' }"
    // galah probes for fastANI by running `/usr/bin/which --tty-only fastANI`, by absolute path.
    // The biocontainers image ships BusyBox which, which rejects that flag, so galah aborts although
    // fastANI is present. --bird_which_shim binds a `which` that tolerates it. PATH does not help:
    // the call is by absolute path.
    containerOptions { params.bird_which_shim && workflow.containerEngine in ['singularity', 'apptainer']
                       ? "--bind ${params.bird_which_shim}:/usr/bin/which" : '' }

    input:
    val(set_id)
    path(genomes, stageAs: 'genomes/*')
    path(checkm2_report)
    val(ani)

    output:
    path "clusters.tsv"        , emit: clusters
    path "representatives.txt" , emit: representatives
    path "versions.yml"        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def quality = checkm2_report ? "--checkm2-quality-report ${checkm2_report}" : ''
    """
    # One fasta per line: galah takes the list, so a long series does not overflow the command line.
    find genomes/ -name '*.fa' -o -name '*.fna' -o -name '*.fasta' | sort > genome_list.txt

    galah cluster \\
        ${args} \\
        --genome-fasta-list genome_list.txt \\
        ${quality} \\
        --ani ${ani} \\
        --threads ${task.cpus} \\
        --output-cluster-definition clusters.tsv \\
        --output-representative-list representatives.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        galah: \$(galah --version 2>&1 | sed 's/^galah //')
        ani_threshold: ${ani}
    END_VERSIONS
    """

    stub:
    """
    printf 'genomes/rep1.fa\\tgenomes/rep1.fa\\n' > clusters.tsv
    printf 'genomes/rep1.fa\\n' > representatives.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        galah: 0.5.2
        ani_threshold: ${ani}
    END_VERSIONS
    """
}
