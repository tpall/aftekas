include { GTDBTK_CLASSIFYWF as CLASSIFY } from '../modules/local/gtdbtk/classifywf/main'


workflow TAXONOMY {

    take:
    final_bins
    
    main:
    ch_versions = channel.empty()
    
    final_bins
    .set { ch_bins }

    CLASSIFY(ch_bins, params.gtdbtk_db, false)
    ch_tax_summary = CLASSIFY.out.summary
    ch_tax_tree = CLASSIFY.out.tree
    ch_versions = CLASSIFY.out.versions

    emit:
    tax_summary = ch_tax_summary
    tax_tree = ch_tax_tree
    versions = ch_versions

}
