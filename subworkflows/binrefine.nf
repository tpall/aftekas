include { 
    HMMER_HMMSEARCH as HMMSEARCH_PFAM
    HMMER_HMMSEARCH as HMMSEARCH_TIGRFAM
     } from '../modules/local/hmmer/hmmsearch/main'
include { PRODIGAL } from '../modules/nf-core/prodigal/main' 
include { CAT_HMM } from '../modules/local/cat/hmm/main'
include { MAGSCOT } from '../modules/local/magscot/refinement/main'

workflow BINREFINE {
    take:
    contigs
    contigs_to_bin

    main:
    ch_versions = channel.empty()

    // Predict genes with Prodigal
    ch_contigs = contigs
    PRODIGAL(ch_contigs, "gff")
    ch_prodigal_faa = PRODIGAL.out.amino_acid_fasta
    ch_versions = ch_versions.mix(PRODIGAL.out.versions)

    // Annotate genes with HMMER
    HMMSEARCH_PFAM(
        ch_prodigal_faa
        .combine( channel.of(params.hmmer_pfam_db) ) 
        .map ({ meta, faa, pfam_db -> [ meta, pfam_db, faa, false, true, false ]})
        )
    ch_tblouts_pfam = HMMSEARCH_PFAM.out.target_summary
    ch_versions = ch_versions.mix(HMMSEARCH_PFAM.out.versions)
    HMMSEARCH_TIGRFAM(
        ch_prodigal_faa
        .combine( channel.of(params.hmmer_tigrfam_db) )
        .map ( { meta, faa, tigrfam_db -> [ meta, tigrfam_db, faa, false, true, false ]})
        )
    ch_tblouts_tigrfam = HMMSEARCH_TIGRFAM.out.target_summary
    ch_tblouts_pfam
    .map ( { meta, pfam -> [ meta.id, meta, pfam ] } )
    .set { ch_tblouts_pfam_mapped }
    ch_tblouts_tigrfam
    .map( { meta, tigrfam -> [ meta.id, meta, tigrfam ] } )
    .set { ch_tblouts_tigrfam_mapped }

    ch_tblouts_pfam_mapped
    .combine ( ch_tblouts_tigrfam_mapped, by: 0 )
    .map( { _id, meta, pfam, _meta2, tigrfam -> [ meta, tigrfam, pfam ]} )
    .set { ch_tblouts }

    CAT_HMM(ch_tblouts)
    ch_hmm = CAT_HMM.out.hmm

    contigs_to_bin
    .map ( { meta, tsv -> [ meta.id, meta, tsv ] } )
    .join (
        ch_hmm
        .map ( { meta, hmm -> [ meta.id, meta, hmm ] } )
    )
    .map ( { _id, meta, tsv, _meta2, hmm -> [ meta, tsv, hmm ] } )
    .set { ch_magscot_input }
    MAGSCOT( ch_magscot_input )

    emit:
    // final_bins = refined_bins
    prodigal_faa = ch_prodigal_faa
    versions = ch_versions

}