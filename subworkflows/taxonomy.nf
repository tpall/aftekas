include { GTDBTK_CLASSIFYWF as CLASSIFY } from '../modules/nf-core/gtdbtk/classifywf/main'


workflow TAXONOMY {

    take:
    final_bins
    
    main:
    ch_versions = channel.empty()
    
    gtdbtk_db = tuple('release220', params.gtdbtk_db)
    CLASSIFY(final_bins, gtdbtk_db, params.use_pplacer_scratch_dir)
    ch_tax_summary = CLASSIFY.out.summary
    ch_tax_tree = CLASSIFY.out.tree
    ch_versions = CLASSIFY.out.versions

    emit:
    tax_summary = ch_tax_summary
    tax_tree = ch_tax_tree
    versions = ch_versions

}
