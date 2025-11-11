include { GTDBTK_CLASSIFYWF as CLASSIFY } from '../modules/local/gtdbtk/classifywf/main'


workflow TAXONOMY {

    take:
    final_bins
    
    main:
    ch_versions = channel.empty()
    
    final_bins
    .map { meta, bins -> [ meta, bins.parent ]}
    .set { ch_bins_dir }

    CLASSIFY(ch_bins_dir, params.gtdbtk_db, params.use_pplacer_scratch_dir)
    ch_tax_summary = CLASSIFY.out.summary
    ch_tax_tree = CLASSIFY.out.tree
    ch_versions = CLASSIFY.out.versions

    emit:
    tax_summary = ch_tax_summary
    tax_tree = ch_tax_tree
    versions = ch_versions

}
