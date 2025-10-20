
include { MEGAHIT } from '../modules/nf-core/megahit/main' 
include { QUAST as METAQUAST_ASSEMBLY } from '../modules/local/metaquast/assembly/main'  
include { BBMAP_ALIGN as ALIGN_MEGAHIT } from '../modules/local/bbmap/align/main'
include { BBMAP_PILEUP as PILEUP_MEGAHIT } from '../modules/nf-core/bbmap/pileup/main'

workflow ASSEMBLY {
    take:
    normed_reads
    trimmed_reads

    main:
    ch_versions = Channel.empty()
    ch_multiqc = Channel.empty()

    // Assemble reads with MEGAHIT
    normed_reads
    .map { meta, reads -> 
        return [ meta, reads[0], reads[1] ]
    }
    .set { ch_normed_reads }
    MEGAHIT(ch_normed_reads)
    ch_versions = ch_versions.mix(MEGAHIT.out.versions)
    ch_multiqc = ch_multiqc.mix(MEGAHIT.out.log)
    MEGAHIT.out.contigs
    .map { meta, contigs -> 
        def fmeta = [:]
        fmeta.id = meta.id
        fmeta.assembler = "MEGAHIT"
        return [ fmeta, contigs ]
    }
    .set { ch_contigs }

    // Assess assembly quality with QUAST
    METAQUAST_ASSEMBLY(ch_contigs)
    ch_multiqc = ch_multiqc.mix(METAQUAST_ASSEMBLY.out.qc)
    ch_versions = ch_versions.mix(METAQUAST_ASSEMBLY.out.versions)

    // Map reads back to assembly to calculate coverage
    ch_contigs
    .map { meta, contigs -> [meta.id, meta, contigs] }
    .set { ch_bbmap_contigs }     
    
    trimmed_reads.map { meta, sample_reads -> [meta.id, meta, sample_reads] }
    .join(ch_bbmap_contigs, by: 0)
    .map { _id, reads_meta, sample_reads, _contigs_meta, contigs -> [ reads_meta, sample_reads, contigs ]
    }
    .set { ch_bbmap_input }

    ALIGN_MEGAHIT (ch_bbmap_input)
    ch_versions = ch_versions.mix(ALIGN_MEGAHIT.out.versions)
    ch_multiqc = ch_multiqc.mix(ALIGN_MEGAHIT.out.log)
    
    PILEUP_MEGAHIT (ALIGN_MEGAHIT.out.bam)
    ch_covstats = PILEUP_MEGAHIT.out.covstats
    ch_versions = ch_versions.mix(PILEUP_MEGAHIT.out.versions)
    ch_multiqc = ch_multiqc.mix(PILEUP_MEGAHIT.out.hist)

    emit:
    contigs = ch_contigs
    covstats = ch_covstats
    versions = ch_versions
    multiqc = ch_multiqc
}
