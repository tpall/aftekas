include { MEGAHIT } from '../modules/nf-core/megahit/main'
include { SPADES } from '../modules/nf-core/spades/main'
include { BBMAP_ALIGN as ALIGN } from '../modules/local/bbmap/align/main'
include { SAMTOOLS_SORT } from '../modules/nf-core/samtools/sort/main'
include { BBMAP_PILEUP as PILEUP } from '../modules/nf-core/bbmap/pileup/main'
include { FIX_FASTA_HEADER } from '../modules/local/fix_fasta_header/main'

workflow ASSEMBLY {
    take:
    normed_reads
    processed_reads

    main:
    ch_versions = channel.empty()
    ch_multiqc = channel.empty()
    // Only metaSPAdes produces an assembly graph, so this stays empty on the MEGAHIT path.
    ch_assembly_graph = channel.empty()

    def assembler = (params.assembler ?: 'megahit').toString().toLowerCase()
    if (!(assembler in ['megahit', 'spades'])) {
        error "--assembler must be 'megahit' or 'spades', got '${params.assembler}'"
    }
    def asm_name = assembler == 'spades' ? 'metaSPAdes' : 'MEGAHIT'

    if (assembler == 'megahit') {
        // Assemble reads with MEGAHIT
        normed_reads
        .map { meta, reads ->
            if ( meta.single_end ) {
                [ meta, reads, [] ] }
            else {
                [ meta, reads[0], reads[1] ]
            }}
        .set { ch_megahit_input }
        MEGAHIT(ch_megahit_input)
        ch_versions = ch_versions.mix(MEGAHIT.out.versions)
        ch_assembled = MEGAHIT.out.contigs
    }
    else {
        // metaSPAdes. Reads are already BBNorm-normalised, which is what keeps SPAdes'
        // memory tractable on a gut metagenome.
        normed_reads
        .map { meta, reads -> [ meta, reads, [], [] ] }
        .set { ch_spades_input }
        SPADES(ch_spades_input, [], [])
        ch_versions = ch_versions.mix(SPADES.out.versions)
        // Bin the SCAFFOLDS, not contigs.fasta. The graph's P (path) lines name
        // scaffolds, and scaffolds.fasta / contigs.fasta carry independent NODE_*
        // numbering, so binning contigs would leave every contig name unresolvable
        // against the graph. See accessory-annex docs/design.md.
        ch_assembled = SPADES.out.scaffolds
        SPADES.out.gfa
        .map { meta, gfa -> [ [ id: meta.id, assembler: asm_name ], gfa ] }
        .set { ch_assembly_graph }
    }

    ch_assembled
    .map { meta, contigs ->
        def fmeta = [:]
        fmeta.id = meta.id
        fmeta.assembler = asm_name
        return [ fmeta, contigs ]
    }
    .set { ch_contigs_raw }

    // Fix contig headers by keeping only contig id. reformat.sh trd=t trims at the
    // first whitespace; MEGAHIT headers carry flags/length after a space, while SPAdes
    // NODE_* names have none and pass through unchanged -- which is what keeps them
    // matching the GFA path lines.
    FIX_FASTA_HEADER(ch_contigs_raw)

    // Drop samples whose assembly is empty so they don't break BBMap indexing downstream
    FIX_FASTA_HEADER.out.fixed
    .join(FIX_FASTA_HEADER.out.check)
    .filter { meta, _contigs, check ->
        def passed = check.text.contains("passed")
        if (!passed) log.warn "ASSEMBLY: skipping ${meta.id} -- ${asm_name} produced empty assembly"
        passed
    }
    .map { meta, contigs, _check -> [ meta, contigs ] }
    .set { ch_contigs }

    // Map reads back to assembly to calculate coverage.
    // The meta carried forward is the reads' own (ALIGN needs single_end), with the
    // assembler added: the contigs' fmeta is dropped by this join, and without it
    // SAMTOOLS_SORT and PILEUP see no meta.assembler and name their files null_*.
    // Both consumers of these channels key on meta.id and discard this meta, so the
    // extra entry is local to the mapping-back steps.
    processed_reads.map { meta, reads -> [ meta.id, meta, reads ] }
    .join( ch_contigs.map { meta, contigs -> [ meta.id, meta, contigs ] }, by: 0 )
    .map { _id, meta, reads, _meta2, contigs ->
        [ meta + [ assembler: asm_name ], reads, contigs ] }
    .set { ch_bbmap_input }

    ALIGN(ch_bbmap_input)
    ch_bams = ALIGN.out.bam
    ch_versions = ch_versions.mix(ALIGN.out.versions)

    SAMTOOLS_SORT (ch_bams, tuple([:], []), 'bai')
    ch_sorted_bams = SAMTOOLS_SORT.out.bam
    ch_sorted_bambais = SAMTOOLS_SORT.out.bai
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions)

    PILEUP(ch_sorted_bams)
    PILEUP.out.covstats
    .map { meta, stats ->
        def fmeta = [:]
        fmeta.id = meta.id
        fmeta.assembler = asm_name
        return [ fmeta, stats ]
    }
    .set { ch_covstats }
    ch_versions = ch_versions.mix(PILEUP.out.versions)

    ch_contigs.map { meta, contigs -> [ meta.id, meta, contigs ] }
    .combine(ch_sorted_bams.map { meta, bam -> [ meta.id, meta, bam ] }, by: 0)
    .map { _id, _meta, contigs, _meta2, bam -> [ _meta, contigs, bam ] }
    .set { ch_assemblies }

    emit:
    assemblies = ch_assemblies
    covstats = ch_covstats
    bams = ch_sorted_bams
    bambais = ch_sorted_bambais
    assembly_graph = ch_assembly_graph
    versions = ch_versions
    multiqc = ch_multiqc
}
