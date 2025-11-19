
include { MEGAHIT } from '../modules/nf-core/megahit/main' 
include { BBMAP_ALIGN as ALIGN } from '../modules/local/bbmap/align/main'
include { SAMTOOLS_SORT } from '../modules/nf-core/samtools/sort/main'
include { BBMAP_PILEUP as PILEUP } from '../modules/nf-core/bbmap/pileup/main'
include { FIX_FASTA_HEADER } from '../modules/local/fix_fasta_header/main'

workflow ASSEMBLY {
    take:
    normed_reads
    processed_reads

    main:
    ch_versions = channel.empty()
    ch_multiqc = channel.empty()

    // Assemble reads with MEGAHIT
    normed_reads
    .map { meta, reads -> 
        if ( meta.single_end ) { 
            [ meta, reads, [] ] }
        else {
            [ meta, reads[0], reads[1] ]
        }}
    .set { ch_normed_reads }
    MEGAHIT(ch_normed_reads)
    ch_versions = ch_versions.mix(MEGAHIT.out.versions)
    MEGAHIT.out.contigs
    .map { meta, contigs -> 
        def fmeta = [:]
        fmeta.id = meta.id
        fmeta.assembler = "MEGAHIT"
        return [ fmeta, contigs ]
    }
    .set { ch_contigs_megahit }

    // Fix contig headers by keeping only contig id
    FIX_FASTA_HEADER(ch_contigs_megahit)
    ch_contigs = FIX_FASTA_HEADER.out.fixed

    // Map reads back to assembly to calculate coverage
    processed_reads.map { meta, reads -> [ meta.id, meta, reads ] }
    .join( ch_contigs.map { meta, contigs -> [ meta.id, meta, contigs ] }, by: 0 )
    .map { _id, meta, reads, _meta2, contigs -> [ meta, reads, contigs ] }
    .set { ch_bbmap_input }

    ALIGN(ch_bbmap_input)
    ch_bams = ALIGN.out.bam
    ch_versions = ch_versions.mix(ALIGN.out.versions)

    SAMTOOLS_SORT (ch_bams, tuple([:], []), 'bai')
    ch_sorted_bams = SAMTOOLS_SORT.out.bam
    ch_sorted_bambais = SAMTOOLS_SORT.out.bai
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions)

    PILEUP(ch_sorted_bams)
    PILEUP.out.covstats
    .map { meta, stats -> 
        def fmeta = [:]
        fmeta.id = meta.id
        fmeta.assembler = "MEGAHIT"
        return [ fmeta, stats ]
    }
    .set { ch_covstats }
    ch_versions = ch_versions.mix(PILEUP.out.versions)

    ch_contigs.map { meta, contigs -> [ meta.id, meta, contigs ] }
    .combine(ch_sorted_bams.map { meta, bam -> [ meta.id, meta, bam ] }, by: 0)
    .map { _id, _meta, contigs, _meta2, bam -> [ _meta, contigs, bam ] }
    .set { ch_assemblies }

    emit:
    assemblies = ch_assemblies
    covstats = ch_covstats
    bambais = ch_sorted_bambais
    versions = ch_versions
    multiqc = ch_multiqc
}
