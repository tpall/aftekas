include { SYLPH_SKETCH     } from '../modules/local/sylph/sketch/main'
include { SYLPH_PROFILE    } from '../modules/local/sylph/profile/main'
include { SYLPHTAX_TAXPROF } from '../modules/local/sylphtax/taxprof/main'

// Read-based taxonomic profiling with sylph against GTDB. Answers "which species
// are in the reads" independently of assembly and binning, so it works as the
// whole analysis (--skip_assembly, the amplicon-priced shotgun snapshot) and,
// next to GTDB-Tk on the bins, shows the part of the community binning lost.
workflow PROFILE {

    take:
    processed_reads   // [meta, reads]: host/phiX-removed, NOT normalised
    sylph_db          // path to a .syl2db / .syldb
    sylph_tax_dir     // directory from `sylph-tax download`, or [] to fetch in-task
    sylph_tax_id      // e.g. GTDB_r226 -- must match the release the db was built from

    main:
    ch_versions = channel.empty()

    SYLPH_SKETCH(processed_reads)
    ch_versions = ch_versions.mix(SYLPH_SKETCH.out.versions.first())

    SYLPH_PROFILE(SYLPH_SKETCH.out.sketch.map { _meta, sketch -> sketch }.collect(), sylph_db)
    ch_versions = ch_versions.mix(SYLPH_PROFILE.out.versions)

    SYLPHTAX_TAXPROF(SYLPH_PROFILE.out.profile, sylph_tax_dir, sylph_tax_id)
    ch_versions = ch_versions.mix(SYLPHTAX_TAXPROF.out.versions)

    emit:
    profile            = SYLPH_PROFILE.out.profile
    taxprof            = SYLPHTAX_TAXPROF.out.taxprof
    relative_abundance = SYLPHTAX_TAXPROF.out.relative_abundance
    sequence_abundance = SYLPHTAX_TAXPROF.out.sequence_abundance
    versions           = ch_versions
}
