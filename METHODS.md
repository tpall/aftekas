# Methods

## English

Sequencing reads were processed with *aftekas* (https://github.com/tpall/aftekas), a Nextflow (Di Tommaso et al., 2017) pipeline inspired by nf-core/mag and optimised for SLURM array execution. Read quality was assessed with FastQC v0.12.1 (Andrews, 2010), and paired-end reads were trimmed and filtered with fastp v1.0.1 (Chen et al., 2018). Cleaned reads were sequentially aligned with Bowtie2 v2.5.4 (Langmead & Salzberg, 2012) and samtools v1.21 (Danecek et al., 2021) against the human-T2T-HLA host reference masked with 150-mers from FDA-ARGOS bacterial, RefSeq viral and Millard Lab phage genomes (Constantinides et al., 2023) and the phiX174 genome, to remove host and contaminant reads. Surviving reads were normalised with BBNorm (BBMap v39.18; Bushnell, 2014) and assembled with MEGAHIT v1.2.9 (Li et al., 2015); assembly quality was evaluated with metaQUAST v5.3.0 (Mikheenko et al., 2016). Reads were re-mapped to contigs with BBMap, sorted with samtools, and per-contig depths summarised with `jgi_summarize_bam_contig_depths`. Contigs were independently binned with MetaBAT2 v2.17 (Kang et al., 2019), MaxBin2 v2.2.7 (Wu et al., 2016), CONCOCT v1.1.0 (Alneberg et al., 2014) and VAMB v5.0.4 (Nissen et al., 2021). Genes were predicted with Prodigal v2.6.3 (Hyatt et al., 2010), and bin sets were refined with Binette v1.2.0 (https://github.com/genotoul-bioinfo/Binette), using CheckM2 v1.1.0 (Chklovski et al., 2023) for bin quality scoring. Refined metagenome-assembled genomes were taxonomically classified with GTDB-Tk v2.4.1 (Chaumeil et al., 2022) against release 226. Run-level QC metrics were aggregated with MultiQC v1.31 (Ewels et al., 2016).

For a sample series, refined MAGs were dereplicated into one common, versioned reference set with galah v0.5.2 (https://github.com/wwood/galah). Genomes were preclustered at 90 % and clustered at 95 % average nucleotide identity with skani (Shaw & Yu, 2023), requiring an aligned fraction of at least 15 %, and the representative of each cluster was chosen from the CheckM2 estimates by the quality score completeness − 5 × contamination, with penalties for contig number and ambiguous bases (Parks et al., 2020, Supplementary Table 19). Representatives keep a stable identifier across versions of the set and are not replaced by later MAGs. The host- and phiX-removed reads of each sample were mapped to the representatives with CoverM v0.8.0 (Aroney et al., 2025) using strobealign (Sahlin, 2022), counting reads aligned at 95 % identity or more over at least 75 % of their length, and relative abundance, mean coverage and read count were reported per genome. Optionally, the same reads were profiled with sylph v1.0.0 (Shaw & Yu, 2024) against a GTDB database, with lineages assigned by sylph-tax v1.9.1.

## Eesti keeles

Sekveneerimislugemeid töödeldi Nextflow (Di Tommaso et al., 2017) põhise töövooga *aftekas* (https://github.com/tpall/aftekas). *aftekas* lähtub nf-core/mag töövoost, kuhu on alles jäetud ainult meie grupi jaoks vajalikud funktsioonid ja on optimeeritud SLURM-i massiivtööde (array jobs) jaoks. Sekveneeritud järjestuste kvaliteeti hinnati FastQC v0.12.1 (Andrews, 2010) tarkvaraga ning paaris-otsa lugemid trimmiti ja filtreeriti programmiga fastp v1.0.1 (Chen et al., 2018). Puhastatud lugemid joondati järjestikku Bowtie2 v2.5.4 (Langmead & Salzberg, 2012) ja samtools v1.21 (Danecek et al., 2021) abil esmalt inimese T2T-HLA referentsgenoomi vastu, mis on maskeeritud 150-meeride abil FDA-ARGOS bakteriaalsetest, RefSeq viraalsetest ning Millard Lab faagigenoomidest (Constantinides et al., 2023), ning seejärel phiX174 genoomi vastu, et eemaldada peremees- ja saastelugemid. Allesjäänud lugemid normaliseeriti BBNorm-iga (BBMap v39.18; Bushnell, 2014) ning kokku pandi MEGAHIT v1.2.9 (Li et al., 2015) abil; kokkupanekute kvaliteeti hinnati metaQUAST v5.3.0 (Mikheenko et al., 2016) abil. Lugemid joondati uuesti kontiigidele BBMap-iga, sorteeriti samtoolsiga ning kontiigi-tasandi katvused arvutati `jgi_summarize_bam_contig_depths` skriptiga. Kontiigid biniti sõltumatult tarkvaradega MetaBAT2 v2.17 (Kang et al., 2019), MaxBin2 v2.2.7 (Wu et al., 2016), CONCOCT v1.1.0 (Alneberg et al., 2014) ja VAMB v5.0.4 (Nissen et al., 2021). Geenid ennustati Prodigal v2.6.3 (Hyatt et al., 2010) abil ning binnide komplektid täpsustati Binette v1.2.0 (https://github.com/genotoul-bioinfo/Binette) abil, kasutades CheckM2 v1.1.0 (Chklovski et al., 2023) binnide kvaliteedi hindamiseks. Täpsustatud metagenoomidest kokku pandud genoomid (MAG-id) klassifitseeriti taksonoomiliselt GTDB-Tk v2.4.1 (Chaumeil et al., 2022) abil selle andmebaasi v226 vastu. Töövoo kvaliteedinäitajad koondati MultiQC v1.31 (Ewels et al., 2016) raportisse.

Proovide seeria puhul dereplitseeriti täpsustatud MAG-id üheks ühiseks versioneeritud referentskomplektiks galah v0.5.2 (https://github.com/wwood/galah) abil. Genoomid eelklasterdati 90 % ja klasterdati 95 % keskmise nukleotiidse identsuse (ANI) juures skani (Shaw & Yu, 2023) abil, nõudes vähemalt 15 % joondunud osa, ning iga klastri esindaja valiti CheckM2 hinnangute põhjal kvaliteediskoori järgi (täielikkus − 5 × saastatus, lisaks karistus kontiikide arvu ja määramata aluste eest; Parks et al., 2020, lisatabel 19). Esindajate tunnused püsivad komplekti versioonide vahel samad ning hilisemad MAG-id neid ei asenda. Iga proovi peremees- ja phiX-lugemitest puhastatud lugemid joondati esindajatele CoverM v0.8.0 (Aroney et al., 2025) ja strobealign (Sahlin, 2022) abil, arvestades lugemeid, mis joondusid vähemalt 95 % identsusega vähemalt 75 % ulatuses oma pikkusest, ning iga genoomi kohta raporteeriti suhteline arvukus, keskmine katvus ja lugemite arv. Soovi korral profileeriti samu lugemeid sylph v1.0.0 (Shaw & Yu, 2024) abil GTDB andmebaasi vastu ning taksonoomiline kuuluvus lisati sylph-tax v1.9.1 abil.

## References / Kasutatud kirjandus

Alneberg, J., Bjarnason, B. S., de Bruijn, I., Schirmer, M., Quick, J., Ijaz, U. Z., Lahti, L., Loman, N. J., Andersson, A. F., & Quince, C. (2014). Binning metagenomic contigs by coverage and composition. *Nature Methods*, 11(11), 1144–1146. https://doi.org/10.1038/nmeth.3103

Andrews, S. (2010). *FastQC: A quality control tool for high throughput sequence data*. http://www.bioinformatics.babraham.ac.uk/projects/fastqc/

Aroney, S. T. N., Newell, R. J. P., Nissen, J. N., Camargo, A. P., Tyson, G. W., & Woodcroft, B. J. (2025). CoverM: Read alignment statistics for metagenomics. *Bioinformatics*, 41(4), btaf147. https://doi.org/10.1093/bioinformatics/btaf147

Bushnell, B. (2014). *BBMap: A fast, accurate, splice-aware aligner*. Lawrence Berkeley National Laboratory. https://sourceforge.net/projects/bbmap/

Chaumeil, P.-A., Mussig, A. J., Hugenholtz, P., & Parks, D. H. (2022). GTDB-Tk v2: Memory friendly classification with the genome taxonomy database. *Bioinformatics*, 38(23), 5315–5316. https://doi.org/10.1093/bioinformatics/btac672

Chen, S., Zhou, Y., Chen, Y., & Gu, J. (2018). fastp: An ultra-fast all-in-one FASTQ preprocessor. *Bioinformatics*, 34(17), i884–i890. https://doi.org/10.1093/bioinformatics/bty560

Chklovski, A., Parks, D. H., Woodcroft, B. J., & Tyson, G. W. (2023). CheckM2: A rapid, scalable and accurate tool for assessing microbial genome quality using machine learning. *Nature Methods*, 20(8), 1203–1212. https://doi.org/10.1038/s41592-023-01940-w

Constantinides, B., Hunt, M., & Crook, D. W. (2023). Hostile: Accurate decontamination of microbial host sequences. *Bioinformatics*, 39(12), btad728. https://doi.org/10.1093/bioinformatics/btad728

Danecek, P., Bonfield, J. K., Liddle, J., Marshall, J., Ohan, V., Pollard, M. O., Whitwham, A., Keane, T., McCarthy, S. A., Davies, R. M., & Li, H. (2021). Twelve years of SAMtools and BCFtools. *GigaScience*, 10(2), giab008. https://doi.org/10.1093/gigascience/giab008

Di Tommaso, P., Chatzou, M., Floden, E. W., Barja, P. P., Palumbo, E., & Notredame, C. (2017). Nextflow enables reproducible computational workflows. *Nature Biotechnology*, 35(4), 316–319. https://doi.org/10.1038/nbt.3820

Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: Summarize analysis results for multiple tools and samples in a single report. *Bioinformatics*, 32(19), 3047–3048. https://doi.org/10.1093/bioinformatics/btw354

Hyatt, D., Chen, G.-L., LoCascio, P. F., Land, M. L., Larimer, F. W., & Hauser, L. J. (2010). Prodigal: Prokaryotic gene recognition and translation initiation site identification. *BMC Bioinformatics*, 11, 119. https://doi.org/10.1186/1471-2105-11-119

Kang, D. D., Li, F., Kirton, E., Thomas, A., Egan, R., An, H., & Wang, Z. (2019). MetaBAT 2: An adaptive binning algorithm for robust and efficient genome reconstruction from metagenome assemblies. *PeerJ*, 7, e7359. https://doi.org/10.7717/peerj.7359

Langmead, B., & Salzberg, S. L. (2012). Fast gapped-read alignment with Bowtie 2. *Nature Methods*, 9(4), 357–359. https://doi.org/10.1038/nmeth.1923

Li, D., Liu, C.-M., Luo, R., Sadakane, K., & Lam, T.-W. (2015). MEGAHIT: An ultra-fast single-node solution for large and complex metagenomics assembly via succinct de Bruijn graph. *Bioinformatics*, 31(10), 1674–1676. https://doi.org/10.1093/bioinformatics/btv033

Mikheenko, A., Saveliev, V., & Gurevich, A. (2016). MetaQUAST: Evaluation of metagenome assemblies. *Bioinformatics*, 32(7), 1088–1090. https://doi.org/10.1093/bioinformatics/btv697

Nissen, J. N., Johansen, J., Allesøe, R. L., Sønderby, C. K., Armenteros, J. J. A., Grønbech, C. H., Jensen, L. J., Nielsen, H. B., Petersen, T. N., Winther, O., & Rasmussen, S. (2021). Improved metagenome binning and assembly using deep variational autoencoders. *Nature Biotechnology*, 39(5), 555–560. https://doi.org/10.1038/s41587-020-00777-4

Parks, D. H., Chuvochina, M., Chaumeil, P.-A., Rinke, C., Mussig, A. J., & Hugenholtz, P. (2020). A complete domain-to-species taxonomy for Bacteria and Archaea. *Nature Biotechnology*, 38(9), 1079–1086. https://doi.org/10.1038/s41587-020-0501-8

Sahlin, K. (2022). Strobealign: Flexible seed size enables ultra-fast and accurate read alignment. *Genome Biology*, 23, 260. https://doi.org/10.1186/s13059-022-02831-7

Shaw, J., & Yu, Y. W. (2023). Fast and robust metagenomic sequence comparison through sparse chaining with skani. *Nature Methods*, 20(11), 1661–1665. https://doi.org/10.1038/s41592-023-02018-3

Shaw, J., & Yu, Y. W. (2024). Rapid species-level metagenome profiling and containment estimation with sylph. *Nature Biotechnology*, 43(8), 1348–1359. https://doi.org/10.1038/s41587-024-02412-y

Wu, Y.-W., Simmons, B. A., & Singer, S. W. (2016). MaxBin 2.0: An automated binning algorithm to recover genomes from multiple metagenomic datasets. *Bioinformatics*, 32(4), 605–607. https://doi.org/10.1093/bioinformatics/btv638
