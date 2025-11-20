include { BINETTE } from '../modules/local/binette/main'
include { PRODIGAL } from '../modules/nf-core/prodigal/main' 

workflow BINREFINE {
    take:
    contigs
    contigs_to_bin_mappings

    main:
    ch_versions = channel.empty()

    // Predict genes with Prodigal
    ch_contigs = contigs
    PRODIGAL(ch_contigs, "gff")
    ch_prodigal_faa = PRODIGAL.out.amino_acid_fasta
    ch_versions = ch_versions.mix(PRODIGAL.out.versions)

    ch_binette_input = contigs_to_bin_mappings
    .map { meta, mappings ->  [ meta.id, meta, mappings ] } 
    .join( ch_contigs.map { meta, cont -> [ meta.id, meta, cont ] }, by: 0 )
    .map { id, meta, mappings, _meta2, cont -> [ id, meta, mappings, cont ] }
    .join( ch_prodigal_faa.map { meta, proteins -> [ meta.id, meta, proteins] } )
    .map { _id, meta, mappings, cont, _meta2, proteins -> [ meta, mappings, [], cont, proteins ] }
    
    BINETTE(ch_binette_input, params.checkm2_db)
    ch_refined_bins = BINETTE.out.bins
    ch_refined_bins_qc = BINETTE.out.quality_report
    ch_input_bins_qc = BINETTE.out.input_quality_report
    ch_versions = ch_versions.mix(BINETTE.out.versions)

    emit:
    refined_bins = ch_refined_bins
    refined_bins_qc = ch_refined_bins_qc
    contig_to_bin = BINETTE.out.contig_to_bin
    input_bins_qc = ch_input_bins_qc
    prodigal_faa = ch_prodigal_faa
    versions = ch_versions

}