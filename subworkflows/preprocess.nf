include { 
    FASTQC as QC_RAW
    FASTQC as QC_PREPROC
    } from '../modules/nf-core/fastqc/main'
include { FASTP } from '../modules/nf-core/fastp/main'
include { BBMAP_BBNORM as BBNORM } from '../modules/nf-core/bbmap/bbnorm/main'
include { 
    BOWTIE2_ALIGN as REMOVE_HOST 
    BOWTIE2_ALIGN as REMOVE_PHIX 
    } from '../modules/nf-core/bowtie2/align/main'

workflow PREPROCESS {

    take:
    raw_reads
    host_index 
    phix_index
    host_fasta
    phix_fasta

    main:
    ch_multiqc = channel.empty( )
    ch_versions = channel.empty( )  
    
    // Initial QC on raw reads
    QC_RAW(raw_reads)
    ch_multiqc = ch_multiqc.mix(QC_RAW.out.zip)
    ch_versions = ch_versions.mix(QC_RAW.out.versions)
    
    // Trim and filter reads
    FASTP(raw_reads, [], false, false, false)
    ch_trimmed_reads = FASTP.out.reads
    ch_versions = ch_versions.mix(FASTP.out.versions)
    ch_multiqc = ch_multiqc.mix(FASTP.out.json)
    
    // Build Bowtie2 host index
    REMOVE_HOST(ch_trimmed_reads, host_index.first(), host_fasta.first(), true, true)
    ch_reads_host_removed = REMOVE_HOST.out.fastq
    ch_versions = ch_versions.mix(REMOVE_HOST.out.versions)
    ch_multiqc = ch_multiqc.mix(REMOVE_HOST.out.log)

    // Build Bowtie2 phix index and remove phix reads
    REMOVE_PHIX(ch_reads_host_removed, phix_index.first(), phix_fasta.first(), true, true)
    ch_reads_phix_removed = REMOVE_PHIX.out.fastq
    ch_multiqc = ch_multiqc.mix(REMOVE_PHIX.out.log)

    // QC on trimmed reads
    QC_PREPROC(ch_reads_phix_removed)
    ch_multiqc = ch_multiqc.mix(QC_PREPROC.out.zip)
    ch_versions = ch_versions.mix(QC_PREPROC.out.versions)

    // Normalize reads
    BBNORM(ch_reads_phix_removed)
    ch_reads_norm = BBNORM.out.fastq
    ch_versions = ch_versions.mix(BBNORM.out.versions)

    emit:
    processed_reads = ch_reads_phix_removed
    normed_reads = ch_reads_norm
    multiqc = ch_multiqc
    versions = ch_versions  

}
