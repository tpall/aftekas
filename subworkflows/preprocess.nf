include { FASTQC as FASTQC_RAW_READS; FASTQC as FASTQC_TRIMMED_READS } from '../modules/nf-core/fastqc/main'
include { FASTP } from '../modules/nf-core/fastp/main'
include { BBMAP_BBNORM } from '../modules/nf-core/bbmap/bbnorm/main'
include { BOWTIE2_ALIGN as REMOVE_HOST; BOWTIE2_ALIGN as REMOVE_PHIX } from '../modules/nf-core/bowtie2/align/main'
include { MULTIQC } from '../modules/nf-core/multiqc/main'


workflow PREPROCESS {

    take:
    raw_reads
    host_index 
    phix_index
    host_fasta
    phix_fasta

    main:
    ch_multiqc = Channel.empty( )
    ch_versions = Channel.empty( )  
    
    // Initial QC on raw reads
    FASTQC_RAW_READS(raw_reads)
    ch_multiqc = ch_multiqc.mix(FASTQC_RAW_READS.out.zip)
    ch_versions = ch_versions.mix(FASTQC_RAW_READS.out.versions)
    
    // Trim and filter reads
    FASTP(raw_reads, [], false, params.fastp_save_trimmed_fail, false)
    ch_trimmed_reads = FASTP.out.reads
    ch_versions = ch_versions.mix(FASTP.out.versions)
    ch_multiqc = ch_multiqc.mix(FASTP.out.json)
    
    // QC on trimmed reads
    FASTQC_TRIMMED_READS(ch_trimmed_reads)
    ch_multiqc = ch_multiqc.mix(FASTQC_TRIMMED_READS.out.zip)
    ch_versions = ch_versions.mix(FASTQC_TRIMMED_READS.out.versions)
    
    // Build Bowtie2 host index
    REMOVE_HOST(ch_trimmed_reads, host_index.first(), host_fasta.first(), true, true)
    ch_reads_host_removed = REMOVE_HOST.out.fastq
    ch_versions = ch_versions.mix(REMOVE_HOST.out.versions)
    ch_multiqc = ch_multiqc.mix(REMOVE_HOST.out.log)

    // Build Bowtie2 phix index and remove phix reads
    REMOVE_PHIX(ch_reads_host_removed, phix_index.first(), phix_fasta.first(), true, true)
    ch_reads_phix_removed = REMOVE_PHIX.out.fastq
    ch_multiqc = ch_multiqc.mix(REMOVE_PHIX.out.log)

    // Normalize reads
    BBMAP_BBNORM(ch_reads_phix_removed)
    ch_reads_norm = BBMAP_BBNORM.out.fastq
    ch_versions = ch_versions.mix(BBMAP_BBNORM.out.versions)
    ch_multiqc = ch_multiqc.mix(BBMAP_BBNORM.out.log)

    emit:
    trimmed_reads = ch_reads_phix_removed
    normed_reads = ch_reads_norm
    multiqc_files = ch_multiqc
    versions = ch_versions  

}
