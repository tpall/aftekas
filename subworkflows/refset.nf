include { PREP_MAGS         } from '../modules/local/prep_mags/main'
include { GALAH_CLUSTER     } from '../modules/local/galah/main'
include { GENOME_SET_UPDATE } from '../modules/local/genome_set/main'
include { COVERM_GENOME     } from '../modules/local/coverm/main'

// One common reference set for a sample series, and abundances measured against it.
//
// Per-sample assembly and binning stay as they are. This subworkflow takes the refined MAGs of every
// sample in the run, clusters them together with the genomes the client's set already holds, and then
// maps every sample's reads to the set. A sample's own assembly is never the denominator: two
// samples' assemblies are different sequences, so their MAG shares are not comparable.
//
// See docs/design-common-reference-set.md for the identity rules and why this is not a co-assembly.
workflow REFSET {
    take:
    refined_bins        // [ meta, [ bins ] ]
    bins_quality        // [ meta, quality_report ]
    reads               // [ meta, [ reads ] ]  host- and phiX-removed, never the normalised ones
    previous_set        // path to the previous genome-set version directory, or []

    main:
    ch_versions = channel.empty()

    // Unique MAG names across samples, and one quality table in CheckM2 columns
    PREP_MAGS(refined_bins.join(bins_quality))
    ch_versions = ch_versions.mix(PREP_MAGS.out.versions.first())

    ch_new_mags = PREP_MAGS.out.mags.collect()
    ch_quality  = PREP_MAGS.out.quality
        .collectFile(name: 'quality_report.tsv', keepHeader: true, skip: 1)

    // The genomes the set already holds go into the same clustering, so an existing genome keeps its
    // id instead of being rediscovered under a new one.
    ch_previous = previous_set ? channel.fromPath(previous_set, checkIfExists: true) : channel.value([])
    ch_previous_reps = previous_set
        ? channel.fromPath("${previous_set}/representatives/*.fa").collect()
        : channel.value([])

    // new MAGs of this run plus the genomes the set already holds, as one flat list of fastas
    ch_all_genomes = ch_new_mags.mix(ch_previous_reps).collect().map { it.flatten() }

    GALAH_CLUSTER(
        params.genome_set_id,
        ch_all_genomes,
        ch_quality,
        params.genome_set_ani
    )
    ch_versions = ch_versions.mix(GALAH_CLUSTER.out.versions)

    // The label says which sample contributed the MAGs of this round; with several samples in one
    // run it is the run's sample list, kept for the history in members.tsv.
    ch_label = refined_bins.map { meta, _bins -> meta.id }.collect().map { it.sort().join(',') }

    GENOME_SET_UPDATE(
        params.genome_set_id,
        ch_label,
        params.genome_set_ani,
        GALAH_CLUSTER.out.clusters,
        ch_all_genomes,
        ch_previous,
        GALAH_CLUSTER.out.versions
    )
    ch_versions = ch_versions.mix(GENOME_SET_UPDATE.out.versions)

    // Abundance of every genome in the set, per sample, from mapping. When the set gained genomes,
    // earlier samples are remapped as well: that is what --remap_all carries into the run's input.
    ch_set_version = GENOME_SET_UPDATE.out.set_json
        .map { json -> new groovy.json.JsonSlurper().parseText(json.text).version }

    COVERM_GENOME(
        reads,
        GENOME_SET_UPDATE.out.representatives.collect(),
        ch_set_version
    )
    ch_versions = ch_versions.mix(COVERM_GENOME.out.versions.first())

    emit:
    set_dir         = GENOME_SET_UPDATE.out.set_dir
    set_json        = GENOME_SET_UPDATE.out.set_json
    members         = GENOME_SET_UPDATE.out.members
    representatives = GENOME_SET_UPDATE.out.representatives
    abundance       = COVERM_GENOME.out.abundance
    versions        = ch_versions
}
