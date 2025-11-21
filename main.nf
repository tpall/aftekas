#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { validateParameters; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'
include { SETUP } from './subworkflows/setup.nf'
include { PREPROCESS } from './subworkflows/preprocess.nf'
include { CAT_FASTQ } from './modules/nf-core/cat/fastq/main' 
include { ASSEMBLY } from './subworkflows/assembly.nf'
include { METAQUAST } from './modules/local/metaquast/main'  
include { BINNING } from './subworkflows/binning.nf'
include { BINREFINE } from './subworkflows/binrefine.nf'
include { TAXONOMY } from './subworkflows/taxonomy.nf'
include { MULTIQC } from './modules/nf-core/multiqc'

// Defaulting workflow
workflow {

    main:
    // Validate input parameters
    validateParameters()

    // Print summary of supplied parameters
    log.info paramsSummaryLog(workflow)

    // Create a new channel of metadata from a sample sheet passed to the pipeline through the --input parameter
    ch_input = channel.fromList(samplesheetToList(params.input, "assets/schema_input.json"))
    ch_versions = channel.empty()
    ch_multiqc = channel.empty()

    // Download and index host genome and phiX genome
    SETUP(params.phix_accession, params.host_fasta_url, params.host_index_url)
    ch_versions = ch_versions.mix(SETUP.out.versions)

    ch_input
    .map { meta, fastq_1, fastq_2 -> 
            if ( fastq_2 == [] ) { 
                meta.single_end = true }
            else { 
                meta.single_end = false } 
            return [ meta, fastq_1, fastq_2 ]
        }
    .groupTuple()
    .branch {
        meta, fastq_1, fastq_2 ->
            single  : fastq_1.size() == 1
                return [ meta, fastq_1 + fastq_2 ]
            multiple: fastq_1.size() > 1
                return [ meta, fastq_1 + fastq_2 ]
    }
    .set { ch_fastq }
    
    CAT_FASTQ (ch_fastq.multiple).reads
    .mix(ch_fastq.single)
    .map { meta, reads -> [ meta, reads.findAll() ] }
    .set { ch_cat_fastq }
    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions.first().ifEmpty(null))

    // Quality trimming, host and phiX removal
    PREPROCESS(ch_cat_fastq, SETUP.out.host_index, SETUP.out.phix_index, SETUP.out.host_fasta, SETUP.out.phix_fasta)
    ch_versions = ch_versions.mix(PREPROCESS.out.versions)
    ch_multiqc = ch_multiqc.mix(PREPROCESS.out.multiqc.collect { it -> it[1] }.ifEmpty([]))
    ch_processed_reads = PREPROCESS.out.processed_reads
    ch_normed_reads = PREPROCESS.out.normed_reads
    
    // Assemble and assess assembly
    ASSEMBLY(
        ch_normed_reads,
        ch_processed_reads
    )
    ch_versions = ch_versions.mix(ASSEMBLY.out.versions)
    ch_multiqc = ch_multiqc.mix(ASSEMBLY.out.multiqc.collect { it -> it[1] }.ifEmpty([]))
    ch_assemblies = ASSEMBLY.out.assemblies
    ch_bambais = ASSEMBLY.out.bambais
    ch_assemblies
        .map { meta, contigs, _bam -> [ meta, contigs ] }
        .set { ch_contigs }
    ch_covstats = ASSEMBLY.out.covstats

    // Assess assembly quality with QUAST
    METAQUAST( ch_contigs.collect { _meta, contigs -> [ contigs ] } )
    ch_multiqc = ch_multiqc.mix(METAQUAST.out.report)
    ch_versions = ch_versions.mix(METAQUAST.out.versions)

    // Binning of contigs into MAGs
    BINNING(ch_assemblies, ch_bambais)
    ch_versions = ch_versions.mix(BINNING.out.versions)
    ch_contigs_to_bin = BINNING.out.contigs2bin
    ch_binning_results = BINNING.out.binning_results

    // Bin refinement
    BINREFINE(ch_contigs, ch_contigs_to_bin)
    ch_versions = ch_versions.mix(BINREFINE.out.versions)
    ch_final_bins = BINREFINE.out.refined_bins
    ch_multiqc = ch_multiqc.mix(BINREFINE.out.refined_bins_qc.collect { it -> it[1] }.ifEmpty([]))

    // Taxonomic classification of final bins
    TAXONOMY(ch_final_bins)
    ch_versions = ch_versions.mix(TAXONOMY.out.versions)

    // Generate MultiQC report
    ch_multiqc_config = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
    MULTIQC(ch_multiqc.collect(), ch_multiqc_config, [], [], [], [])
    ch_versions = ch_versions.mix(MULTIQC.out.versions)

    // Publish workflow outputs
    publish:
    processed_reads = ch_processed_reads
    .map { meta, reads -> 
        if ( meta.single_end ) { 
            [ id: meta.id, fq1: reads ] }
        else {
            [ id: meta.id, fq1: reads[0], fq2: reads[1] ] 
        }}
    contigs = ch_contigs.map { meta, contigs -> [ id: meta.id, assembler: meta.assembler, fa: contigs ] }
    covstats = ch_covstats.map { meta, stats -> [ id: meta.id, assembler: meta.assembler, tsv: stats ] }
    prodigal_faa = BINREFINE.out.prodigal_faa
    binning_results = ch_binning_results
    binning_results_qc = BINREFINE.out.input_bins_qc.map { _meta, file -> [ file ]}.flatten()
    final_bins = ch_final_bins
    final_bins_qc = BINREFINE.out.refined_bins_qc
    final_contig_to_bin = BINREFINE.out.contig_to_bin
    tax_summary = TAXONOMY.out.tax_summary
    tax_tree = TAXONOMY.out.tax_tree
    multiqc_report = MULTIQC.out.report
    versions = ch_versions.collectFile(name: 'versions.yml')
}

output {
    processed_reads {
        path { reads -> 
            if ( reads.fq2 == null ) {
                reads.fq1 >> "processed_reads/${reads.id}_processed.fq.gz"
            } else {
            reads.fq1 >> "processed_reads/${reads.id}_processed_1.fq.gz"
            reads.fq2 >> "processed_reads/${reads.id}_processed_2.fq.gz"
            }
    }}
    contigs {
        path { contigs -> contigs.fa >> "assembly/contigs/${contigs.assembler}_${contigs.id}_contigs.fa.gz" }
    }
    covstats {
        path { covstats -> covstats.tsv >> "assembly/covstats/${covstats.assembler}_${covstats.id}_contigs_covstats.tsv" }
    }
    prodigal_faa {
        path "assembly/prodigal/"
    }
    binning_results {
        path "binning/bins/"
    }
    binning_results_qc {
         path  { filename -> filename >> "binning/qc/${filename.baseName.tokenize('.')[1]}_quality_reports.tsv"}
    }
    final_contig_to_bin {
         path  { filename -> filename >> "binning/contig_to_bin/${filename.baseName.tokenize('_')[1]}_contig_to_bin.tsv"}
    }
    final_bins {
        path "binrefine/"
    }
    final_bins_qc {
        path "binrefine/qc/"
    }
    tax_summary {
        path "taxonomy/"
    }
    tax_tree {
        path "taxonomy"
    }
    multiqc_report {
        path "."
    }
    versions {
        path "."
    }
}
