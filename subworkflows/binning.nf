include { METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS as GET_ABUN } from '../modules/nf-core/metabat2/jgisummarizebamcontigdepths/main'
include { CONVERT_DEPTHS } from '../modules/local/convert_depths/main'
include { MAXBIN2 } from '../modules/nf-core/maxbin2/main' 
include { METABAT2_METABAT2 as METABAT2 } from '../modules/nf-core/metabat2/metabat2/main'
include { CONCOCT_CONCOCT as CONCOCT } from '../modules/nf-core/concoct/concoct/main'
include { CONCOCT_CONCOCTCOVERAGETABLE as COVERAGETABLE } from '../modules/nf-core/concoct/concoctcoveragetable/main'
include { CONCOCT_CUTUPFASTA as CUTUPFASTA } from '../modules/nf-core/concoct/cutupfasta/main'
include { CONCOCT_EXTRACTFASTABINS as EXTRACTFASTABINS } from '../modules/nf-core/concoct/extractfastabins/main' 
include { CONCOCT_MERGECUTUPCLUSTERING as MERGECUTUPCLUSTERING } from '../modules/nf-core/concoct/mergecutupclustering/main'
include { GUNZIP } from '../modules/nf-core/gunzip/main'
include { SAMTOOLS_BEDCOV } from '../modules/nf-core/samtools/bedcov/main'

workflow BINNING {

    take:
    assembly
    bambais

    main:
    ch_versions = Channel.empty()

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
    
    // // MAXBIN2 binning
    CONVERT_DEPTHS(ch_metabat2_input)
    ch_versions = ch_versions.mix(CONVERT_DEPTHS.out.versions)
    ch_maxbin2_input = CONVERT_DEPTHS.out.output
    MAXBIN2(ch_maxbin2_input)
    // ch_maxbin2_binned = MAXBIN2.out.binned_fastas
    ch_versions = ch_versions.mix(MAXBIN2.out.versions)

    // CONCOCT binning
    GUNZIP(
        assembly
        .map { meta, contigs, _bam -> [ meta, contigs ] }
    )
    ch_versions = ch_versions.mix(GUNZIP.out.versions)
    contigs_gunzipped = GUNZIP.out.gunzip

    contigs_gunzipped
    .map { meta, contigs -> [ meta.id, meta, contigs ] }
    .combine(assembly.map { meta, _contigs, bam -> [ meta.id, meta, bam ] }, by: 0)
    .map { _id, meta, contigs, _meta2, bam -> [ meta, contigs, bam ] }
    .set { assembly_gunzipped }

    CUTUPFASTA( 
        assembly_gunzipped
        .map { meta, contigs, _bam -> [ meta, contigs ] }
        , true )
    ch_versions = ch_versions.mix(CUTUPFASTA.out.versions)
    ch_cutup_bed = CUTUPFASTA.out.bed
    ch_cutup_bed.map { 
        meta, bed -> [ meta.id, meta, bed ] }
        .combine(assembly_gunzipped.map { meta, _contigs, bam -> [ meta.id, meta, bam ] }, by: 0)
        .map { _id, meta, bed, _meta2, bam -> [ _id, meta, bed, bam ] }
        .combine(bambais.map { meta, bai -> [ meta.id, meta, bai ] }, by: 0)
        .map { _id, meta, bed, bam, _meta3, bai -> [ meta, bed, bam, bai ] }
        .set { ch_coveragetable_input }
    
    COVERAGETABLE(ch_coveragetable_input)
    ch_versions = ch_versions.mix(COVERAGETABLE.out.versions)
    ch_concoct_coverage = COVERAGETABLE.out.tsv
    
    ch_concoct_coverage
    .map { meta, coverage -> [ meta.id, meta, coverage ] }
    .combine(assembly_gunzipped.map { meta, contigs, _bam -> [ meta.id, meta, contigs ] }, by: 0)
    .map { _id, meta, coverage, _meta2, contigs -> [ meta, coverage, contigs ] }
    .set { ch_concoct_input }    
    CONCOCT(ch_concoct_input)
    ch_versions = ch_versions.mix(CONCOCT.out.versions)
    MERGECUTUPCLUSTERING(
        CONCOCT.out.clustering_csv
        .map { meta, clustering -> [ meta.id, meta, clustering ] }
        .combine(ch_cutup_bed.map { meta, bed -> [ meta.id, meta, bed ] }, by: 0)
        .map { _id, meta, clustering, _meta2, bed -> [ meta, clustering, bed ] }
    )
    ch_versions = ch_versions.mix(MERGECUTUPCLUSTERING.out.versions)
    ch_merged_csv = MERGECUTUPCLUSTERING.out.csv
    EXTRACTFASTABINS(ch_merged_csv)
    ch_concoct_binned = EXTRACTFASTABINS.out.fasta
    ch_versions = ch_versions.mix(EXTRACTFASTABINS.out.versions)

    emit:
    // maxbin2_bins = ch_maxbin2_binned
    metabat2_bins = ch_metabat2_binned
    concoct_bins = ch_concoct_binned
    versions = ch_versions

}
