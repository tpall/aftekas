# Aftekas (meta) -- a minimal metagenomic workflow

This [Nextflow](https://www.nextflow.io/docs/latest/index.html) workflow was conceived to implement more efficient execution of the [nf-core/mag](https://github.com/nf-core/mag) workflow in slurm hpc cluster.

1. The main feature is the [**use of array jobs**](https://www.nextflow.io/docs/latest/reference/process.html#array) when executed in slurm cluster.
2. Many data munging steps are moved into the customised local modules.
3. The workflow is simple and devoid of forking paths:  
      - uses **only double-ended short reads**,
      - uses **fastp** to trim raw reads,
      - always maps reads to host genome and phix
      - uses **MEGAHIT** for assembly
      - uses **MAXBIN2**, **METABAT2**, **VAMB**, **CONCOCT** for binning
      - uses **BINETTE** for bin refinement
      - assigns taxonomy to refined bins using **GTDBTk**.

By default, *aftekas* uses the "*human-t2t-hla masked with 150mers for 985 FDA-ARGOS bacterial, 18,719 RefSeq viral, and 26,928 Millard Lab phage genomes*" reference genome and bowtie2 index from the [bede/hostile](https://github.com/bede/hostile).

## Usage

### Set up databases

1. Download and setup CheckM2 database as instructed in [https://github.com/chklovski/CheckM2](https://github.com/chklovski/CheckM2).

1. Download and setup GTDBTk database (release 226, recommended) as instructed in [https://ecogenomics.github.io/GTDBTk/installing/index.html](https://ecogenomics.github.io/GTDBTk/installing/index.html).

### Set up sample data

Sample data file must have three columns: *sample*, *fastq_1*, *fastq_2*
e.g.

```bash
sample,fastq_1,fastq_2
EV25,data/EV25_L1_resampled_1.fq.gz,data/EV25_L1_resampled_2.fq.gz
EV25,data/EV25_L3_resampled_1.fq.gz,data/EV25_L3_resampled_2.fq.gz
```

In the example above the sample EV25 was spread to two lines during sequencing and is merged by the workflow before trimming. Assigning them unique sample ids keeps them separate if desired.

### Running

```bash
nextflow pull tpall/aftekas
```

Locally, this workflow can be run using docker, assuming that java, nextflow and docker are running:

```bash
nextflow run tpall/aftekas --input my_samples.csv --checkm2_db "<path to directory with>/CheckM2_database/uniref100.KO.1.dmnd" --gtdbtk_db "<path to directory with>/release226" --fastp_dedup --fastp_trim_polyg
```

In a slurm cluster: activate required software (java, nextflow, singularity) and run e.g.:

```bash
nextflow run tpall/aftekas -profile cluster --input my_samples.csv --array_size 4 --queue main -resume --checkm2_db "<path to directory with>/CheckM2_database/uniref100.KO.1.dmnd" --gtdbtk_db "<path to directory with>/release226" --fastp_dedup --fastp_trim_polyg 
```

Above, we have four unique sample ids in my_samples.csv, therefore `--array_size 4`.

## Todo

Further analyses are open-ended. The original idea was to integrate DRAM for metabolic characterisation. Still, since its 2nd version is going to be a Nextflow workflow itself, it's probably easier to run it separately using the final bins as input. The same holds for the virus identification workflow and the AMR workflow (nf-core/funcscan).

