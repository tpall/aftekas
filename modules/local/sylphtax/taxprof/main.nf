// Attach GTDB lineages to sylph's genome-level profile: one MetaPhlAn-style
// .sylphmpa per sample plus two clade x sample matrices (relative and sequence
// abundance). --no-config keeps sylph-tax from writing to $HOME, which is
// read-only or absent inside containers; the ~50 MB taxonomy metadata is either
// staged from params.sylph_tax_dir or downloaded into the task directory.
process SYLPHTAX_TAXPROF {
    tag "${tax_id}"
    label 'process_small'

    conda "${moduleDir}/environment.yml"
    container "quay.io/biocontainers/sylph-tax:1.9.1--pyhdfd78af_0"

    input:
    path profile
    path taxdir
    val  tax_id

    output:
    path "taxprof/*.sylphmpa"           , emit: taxprof
    path "sylph_relative_abundance.tsv" , emit: relative_abundance
    path "sylph_sequence_abundance.tsv" , emit: sequence_abundance
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def taxpath = taxdir ? "${taxdir}" : "taxdb"
    def fetch   = taxdir ? "" : "mkdir -p taxdb && sylph-tax --no-config --taxonomy-dir taxdb download"
    """
    ${fetch}
    mkdir -p taxprof
    sylph-tax --no-config --taxonomy-dir ${taxpath} \\
        taxprof ${args} ${profile} -t ${tax_id} -o taxprof/

    # sylph names each output after the Sample_file column, i.e. the sketch file;
    # strip the sketch suffix so the files carry plain sample ids.
    for f in taxprof/*.sylphmpa; do
        n=\$(basename "\$f" .sylphmpa | sed -E 's/\\.paired\\.sylsp\$|\\.sylsp\$//')
        [ "taxprof/\$n.sylphmpa" = "\$f" ] || mv "\$f" "taxprof/\$n.sylphmpa"
    done

    sylph-tax --no-config --taxonomy-dir ${taxpath} \\
        merge taxprof/*.sylphmpa --column relative_abundance -o sylph_relative_abundance.tsv
    sylph-tax --no-config --taxonomy-dir ${taxpath} \\
        merge taxprof/*.sylphmpa --column sequence_abundance -o sylph_sequence_abundance.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sylph-tax: \$(python -c 'from sylph_tax.version import __version__; print(__version__)')
        taxonomy: ${tax_id}
    END_VERSIONS
    """

    stub:
    """
    mkdir -p taxprof
    touch taxprof/stub.sylphmpa sylph_relative_abundance.tsv sylph_sequence_abundance.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sylph-tax: \$(python -c 'from sylph_tax.version import __version__; print(__version__)' 2>/dev/null)
        taxonomy: ${tax_id}
    END_VERSIONS
    """
}
