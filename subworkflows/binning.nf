include { METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS as GET_ABUN } from '../modules/nf-core/metabat2/jgisummarizebamcontigdepths/main'
include {
    CONVERT_DEPTHS as MAXBIN2_DEPTHS
    CONVERT_DEPTHS as VAMB_DEPTHS
    } from '../modules/local/convert_depths/main'
include { MAXBIN2 } from '../modules/local/maxbin2/main' 
include { METABAT2_METABAT2 as METABAT2 } from '../modules/nf-core/metabat2/metabat2/main'
include { CONCOCT } from '../modules/local/concoct/main'
include { VAMB_BIN } from '../modules/nf-core/vamb/bin/main'

workflow BINNING {

    take:
    assembly
    bambais

    main:
    ch_versions = channel.empty()

    // Calculate contig abundance
    GET_ABUN(
        assembly
        .map { meta, _contigs, bam -> [ meta, bam, [] ] }
        )
    ch_depths = GET_ABUN.out.depth
    ch_versions = ch_versions.mix(GET_ABUN.out.versions)

    // METABAT2 binning
    assembly
        .map { meta, contigs, _bam -> [ meta.id, meta, contigs ] }
        .combine(ch_depths.map { meta, depth -> [ meta.id, meta, depth ] }, by: 0)
        .map { _id, meta, contigs, _meta2, depth -> [ meta, contigs, depth ] }
        .set { ch_metabat2_input }
    METABAT2(ch_metabat2_input)
    ch_metabat2_binned = METABAT2.out.fasta
    ch_versions = ch_versions.mix(METABAT2.out.versions)
    
    // MAXBIN2 binning
    MAXBIN2_DEPTHS(ch_metabat2_input, [])
    ch_versions = ch_versions.mix(MAXBIN2_DEPTHS.out.versions)
    ch_maxbin2_input = MAXBIN2_DEPTHS.out.output
    MAXBIN2(ch_maxbin2_input)
    ch_maxbin2_binned = MAXBIN2.out.binned_fastas
    ch_versions = ch_versions.mix(MAXBIN2.out.versions)

    // CONCOCT binning
    assembly
        .map { meta, contigs, bam -> [ meta.id, meta, contigs, bam ] }
        .combine(bambais.map { meta, bai -> [ meta.id, meta, bai ] }, by: 0)
        .map { _id, meta, contigs, bam, _meta2, bai -> [ meta, contigs, bam, bai ] }
        .set { ch_concoct_input }
    CONCOCT(ch_concoct_input)
    ch_concoct_binned = CONCOCT.out.fasta
    ch_versions = ch_versions.mix(CONCOCT.out.versions)

    // VAMB binning
    VAMB_DEPTHS(ch_metabat2_input, "contigname")
    VAMB_DEPTHS.out.output
    .map { meta, contigs, _reads, depths -> [ meta, contigs, depths, [], [] ] }
    .set { ch_vamb_input }
    VAMB_BIN(ch_vamb_input)
    ch_vamb_binned = VAMB_BIN.out.bins
    ch_versions = ch_versions.mix(VAMB_BIN.out.versions)

    emit:
    maxbin2_bins = ch_maxbin2_binned
    metabat2_bins = ch_metabat2_binned
    concoct_bins = ch_concoct_binned
    vamb_bins = ch_vamb_binned
    versions = ch_versions

}
