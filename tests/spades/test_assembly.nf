#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Integration test for the --assembler spades path, exercising the ASSEMBLY
// subworkflow on its own. Preprocessing, binning and taxonomy are deliberately out of
// scope: they are unchanged, and GTDB-Tk/CheckM2 would pull in multi-GB databases that
// say nothing about this change.
//
//   nextflow run tests/spades/test_assembly.nf -profile cluster \
//       -c tests/spades/test.config --assembler spades \
//       --r1 <reads_1.fq.gz> --r2 <reads_2.fq.gz>

include { ASSEMBLY } from '../../subworkflows/assembly.nf'

params.r1 = null
params.r2 = null
params.sample = 'spades_test'

workflow {
    if (!params.r1 || !params.r2) {
        error "need --r1 and --r2"
    }

    // The same reads stand in for both inputs: ASSEMBLY takes normalised reads for the
    // assembler and processed reads for the mapping-back step, and this test is not
    // about BBNorm.
    channel.of([ [ id: params.sample, single_end: false ],
                 [ file(params.r1, checkIfExists: true),
                   file(params.r2, checkIfExists: true) ] ])
    .set { ch_reads }

    ASSEMBLY(ch_reads, ch_reads)

    ASSEMBLY.out.assemblies.view { meta, contigs, bam ->
        "ASSEMBLED ${meta.assembler} ${meta.id} contigs=${contigs} bam=${bam}" }
    ASSEMBLY.out.assembly_graph.view { meta, gfa ->
        "GRAPH ${meta.assembler} ${meta.id} gfa=${gfa}" }
    ASSEMBLY.out.bams.view { meta, bam -> "BAM ${bam}" }
    ASSEMBLY.out.covstats.view { meta, tsv -> "COVSTATS ${tsv}" }
}
