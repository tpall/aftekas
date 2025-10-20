#!/usr/bin/env nextflow
nextflow.enable.dsl=2
nextflow.preview.output = true

include { samplesheetToList } from 'plugin/nf-schema'
include { INPUT_CHECK } from './subworkflows/input_check.nf'
include { SETUP } from './subworkflows/setup.nf'
include { PREPROCESS } from './subworkflows/preprocess.nf'
include { MULTIQC } from './modules/nf-core/multiqc'
include { CAT_FASTQ } from './modules/nf-core/cat/fastq/main' 
include { ASSEMBLY } from './subworkflows/assembly.nf'

// Defaulting workflow
workflow {

    main:
    ch_input = file(params.input)
    ch_versions = Channel.empty()
    ch_multiqc = Channel.empty()
    // ch_multiqc_config = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
    // ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

    // Download and index host genome and phiX genome
    SETUP(params.phix_accession, params.host_fasta_url, params.host_index_url)
    ch_versions = ch_versions.mix(SETUP.out.versions)

    INPUT_CHECK (
        ch_input
    )
    .reads
    .map {
        meta, fastq ->
            meta.id = meta.id.split('_')[0..-2].join('_')
            [ meta, fastq ] }
    .groupTuple(by: [0])
    .branch {
        meta, fastq ->
            single  : fastq.size() == 1
                return [ meta, fastq.flatten() ]
            multiple: fastq.size() > 1
                return [ meta, fastq.flatten() ]
    }
    .set { ch_fastq }
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    // Concatenate FastQ files from same sample if required
    CAT_FASTQ (
        ch_fastq.multiple
    )
    .reads
    .mix(ch_fastq.single)
    .set { ch_cat_fastq }
    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions.first().ifEmpty(null))

    // Quality trimming, host and phiX removal
    PREPROCESS(ch_cat_fastq, SETUP.out.host_index, SETUP.out.phix_index, SETUP.out.host_fasta, SETUP.out.phix_fasta)
    ch_versions = ch_versions.mix(PREPROCESS.out.versions)
    ch_multiqc = ch_multiqc.mix(PREPROCESS.out.multiqc_files)
    ch_trimmed_reads = PREPROCESS.out.trimmed_reads
    ch_normed_reads = PREPROCESS.out.normed_reads


    // Assemble and assess assembly
    ASSEMBLY(
        ch_normed_reads,
        ch_trimmed_reads
    )
    ch_versions = ch_versions.mix(ASSEMBLY.out.versions)
    ch_multiqc = ch_multiqc.mix(ASSEMBLY.out.multiqc)
    ch_contigs = ASSEMBLY.out.contigs
    ch_covstats = ASSEMBLY.out.covstats

    // Generate MultiQC report


    //Publish workflow outputs
    publish:
    trimmed_reads = ch_trimmed_reads
    contigs = ch_contigs
    covstats = ch_covstats
    // multiqc_report = MULTIQC.out.report
    versions = ch_versions.collectFile(name: 'versions.yml')
    
}

output {
    trimmed_reads {
        path 'trimmed_reads/'
    }
    contigs {
        path 'assembly/contigs/'
    }
    covstats {
        path 'assembly/coverage_stats/'
    }
    // multiqc_report {
    //     path 'multiqc'
    // }
    versions {
        path 'versions.yml'
    }
}

