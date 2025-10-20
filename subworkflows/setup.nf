include { DATASETS_DOWNLOAD as GET_PHIX_GENOME} from '../modules/local/datasets/main'
include { FETCH_INDEX as FETCH_HOST_INDEX} from '../modules/local/fetch/main'
include { 
    BOWTIE2_BUILD as BUILD_PHIX_INDEX 
    BOWTIE2_BUILD as BUILD_HOST_INDEX 
    } from '../modules/nf-core/bowtie2/build/main'

workflow SETUP {
    
    take:
    phix_accession   // [phix accession]
    host_fasta_url   // [host index path/url]
    host_index_url  // [host index path/url] (optional)

    main:
    ch_versions = Channel.empty()

    // Download and build PhiX index
    Channel.from(phix_accession)
        .map { v -> tuple( [ id: 'phix' ], v ) }
        .set { ch_phix_accession }   
    GET_PHIX_GENOME(ch_phix_accession)
    ch_phix_fasta = GET_PHIX_GENOME.out.fasta
    ch_versions = ch_versions.mix(GET_PHIX_GENOME.out.versions)

    BUILD_PHIX_INDEX(ch_phix_fasta)
    ch_phix_index = BUILD_PHIX_INDEX.out.index
    ch_versions = ch_versions.mix(BUILD_PHIX_INDEX.out.versions)

    // Download host fasta and build index
    Channel.fromPath(host_fasta_url)
        .map { url -> tuple( [ id : url.baseName ], url ) }
        .set { ch_host_fasta }
    if ( !params.host_index_url ) {
        ch_host_index = BUILD_HOST_INDEX( ch_host_fasta ).index
        ch_versions = ch_versions.mix( BUILD_HOST_INDEX.out.versions )
    } else {
        Channel.fromPath(host_index_url)
        .map { url -> tuple( [ id : url.baseName ], url ) }
        .set { ch_host_index_url }
        ch_host_index = FETCH_HOST_INDEX( ch_host_index_url ).index
    }

    emit:
    host_index = ch_host_index
    host_fasta = ch_host_fasta
    phix_fasta = ch_phix_fasta
    phix_index = ch_phix_index
    versions = ch_versions

}
