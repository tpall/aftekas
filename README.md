# Aftekas (meta) -- a minimal metagenomic workflow

This [Nextflow](https://www.nextflow.io/docs/latest/index.html) workflow was conceived to implement more efficient execution of the [nf-core/mag](https://github.com/nf-core/mag) workflow in slurm hpc cluster.

1. The main feature is the [**use of array jobs**](https://www.nextflow.io/docs/latest/reference/process.html#array) when executed in slurm cluster.
2. Many data munging steps are moved into the customised local modules.
3. The workflow is simple and has one choice point, the assembler:  
      - uses **only double-ended short reads**,
      - uses **fastp** to trim raw reads,
      - always maps reads to host genome and phix
      - assembles with **MEGAHIT** (default) or **metaSPAdes** (`--assembler spades`)
      - uses **MAXBIN2**, **METABAT2**, **VAMB**, **CONCOCT** for binning
      - uses **BINETTE** for bin refinement
      - assigns taxonomy to refined bins using **GTDBTk**.

By default, *aftekas* uses the "*human-t2t-hla masked with 150mers for 985 FDA-ARGOS bacterial, 18,719 RefSeq viral, and 26,928 Millard Lab phage genomes*" reference genome and bowtie2 index from the [bede/hostile](https://github.com/bede/hostile).

A short Methods description (English and Estonian) with software versions and references is available in [METHODS.md](METHODS.md).

## Usage

### Set up databases

1. Download and setup CheckM2 database as instructed in [https://github.com/chklovski/CheckM2](https://github.com/chklovski/CheckM2).

1. Download and setup GTDBTk database (release 226, recommended) as instructed in [https://ecogenomics.github.io/GTDBTk/installing/index.html](https://ecogenomics.github.io/GTDBTk/installing/index.html).
   Release 232 needs GTDB-Tk ≥ 2.7: point `--gtdbtk_db` at it, run the classify step in the `ecogenomic/gtdbtk:2.7.1` container (a `withName: '.*CLASSIFY.*' { container = ... }` line in your config) and set `--gtdbtk_skip_ani_screen false`, since 2.7 has no such flag.

1. Optional, for read-based profiling (`--profile_reads`): download a sylph GTDB database, ideally the same release as GTDB-Tk so read-level and genome-level taxonomy agree, and the sylph-tax metadata.

   ```bash
   wget -P databases/sylph https://storage.googleapis.com/sylph-stuff/gtdb-r226-c200-dbv2.syl2db   # 14 GB; needs sylph >= 1.0
   mkdir -p databases/sylph-tax && sylph-tax download --download-to databases/sylph-tax          # ~50 MB, or let the task fetch it
   ```

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

### Choosing the assembler

| flag | default | effect |
|---|---|---|
| `--assembler` | `megahit` | `megahit` or `spades`; `spades` runs SPAdes with `--meta` and also publishes the assembly graph |
| `--save_bams` | off | publish the sorted BAM and its `.bai` to `assembly/bam/` |

**MEGAHIT is the default because it is enough for plain MAG recovery and much cheaper.**
metaSPAdes costs considerably more memory and wall time, and buys one thing MEGAHIT cannot
give: the assembly graph, published as
`assembly/graph/metaSPAdes_<id>_assembly_graph.gfa.gz`. Use `--assembler spades` only for
runs that need it — graph-based recruitment of unbinned accessory and mobile contigs is the
case it was added for. The graph **cannot be recovered after the run**, so a run that may
need it has to ask for it up front.

Two consequences of the metaSPAdes path worth knowing before you choose it:

- **The binned sequences are SPAdes *scaffolds*, not `contigs.fasta`.** This is deliberate:
  the graph's `P` (path) lines name scaffolds, and the two files carry independent `NODE_*`
  numbering, so binning contigs would leave every contig name unresolvable against the
  graph.
- **BBNorm normalisation runs on both paths**, which is what keeps metaSPAdes' memory
  tractable on a gut metagenome.

`--save_bams` is the same kind of up-front decision. The publish gate is evaluated at
publish time, so a finished run without it cannot be re-published — only re-run. Each arm
records what it was given in `pipeline_info/run_manifest.json`, which is how a later
consumer tells "no BAM was requested" apart from "a BAM was requested and is missing".

### Output layout

The published tree splits at the assembler. Read-level outputs are produced before the
assembler branch and sit at the top of `--outdir`; everything derived from an assembly
sits under `MEGAHIT/` or `metaSPAdes/`:

```
results/
├── processed_reads/              shared: trimmed, host- and phiX-removed reads
├── profile/                      shared: sylph read profiling
├── pipeline_info/
│   └── input_samplesheet.csv     the sheet these results were built from
└── MEGAHIT/                      or metaSPAdes/, per --assembler
    ├── assembly/                 contigs, covstats, graph (metaSPAdes only), bam, prodigal
    ├── binning/                  per-binner bins and their QC
    ├── binrefine/                Binette-refined MAGs, contig_to_bin, QC
    ├── taxonomy/                 GTDB-Tk
    ├── genome_set/               dereplicated reference set
    ├── pipeline_info/
    │   └── run_manifest.json     assembler, --save_bams, revision, run name, start time
    └── multiqc_report.html
```

This means **both assemblers can share one `--outdir`**, and the second arm reuses the
first arm's preprocessing instead of repeating it. The subdirectory is what makes that
safe: filenames under `assembly/` and `binning/` carry an assembler token, but
`binrefine/`, `taxonomy/` and `genome_set/` do not, so without it a second arm would
overwrite the first in place — on exactly the files that record which contig belongs to
which bin, and with no error.

Because the read-level outputs are shared, the samplesheet is kept beside them and a run
that passes a different one into the same `--outdir` is refused rather than silently
overwriting another read set's results.

On a `--skip_assembly` snapshot run no assembler subdirectory is created; the report and
`versions.yml` stay at the top.

### Common reference set for a sample series

Per-sample runs give every sample its own MAGs, so two samples of the same environment cannot be
compared genome by genome: a genome "absent" in one sample may only have been binned differently.
With `--genome_set_id`, the refined MAGs are dereplicated into one versioned genome set
([galah](https://github.com/wwood/galah), using the CheckM2 statistics Binette already produced), and
each sample's host- and phiX-removed reads are mapped to it with
[CoverM](https://github.com/wwood/CoverM). Abundance always comes from that mapping, never from the
sample's own assembly. Design and identity rules:
[docs/design-common-reference-set.md](docs/design-common-reference-set.md).

| flag | default | effect |
|---|---|---|
| `--genome_set_id` | off | id of the series' genome set, e.g. `CGS-SITE-A`; turns the reference set on (needs the assembly branch) |
| `--genome_set` | none | directory of the previous version to extend (`results/<assembler>/genome_set/genome_set/` of the run that built it); without it the run builds v1 |
| `--genome_set_ani` | 95 | clustering threshold; a different value is a different set, never a new version |

The first sample builds version 1, and every later sample extends the latest version:

```bash
nextflow run tpall/aftekas --input sample1.csv --genome_set_id CGS-SITE-A ...
nextflow run tpall/aftekas --input sample2.csv --genome_set_id CGS-SITE-A \
    --genome_set /path/to/first-run/results/MEGAHIT/genome_set/genome_set ...
```

Output under `results/<assembler>/genome_set/` (see [Output layout](#output-layout)):

- `genome_set/representatives/<id>.fa`: one genome per representative, named by its id;
- `genome_set/members.tsv`: every MAG ever binned in the series and the representative it belongs to;
- `genome_set/quality.tsv`: CheckM2 statistics per representative, which the next version needs;
- `genome_set/set.json`: id, version, parent version, ANI, genomes added this version, sha256 of each representative;
- `abundance/<sample>.coverm.tsv`: relative abundance, mean coverage and read count per genome, with the set version.

Rules that keep identity meaningful: an id is never reused; a representative is never replaced by a
better MAG of a later sample (the later MAG is recorded in `members.tsv`); a later MAG that falls
within the ANI threshold of two representatives keeps both and is reported in `set.json` as a bridge.

Only the samples in the run are mapped. When a new version adds genomes, the earlier samples of the
series have to be mapped again against the new representatives (CoverM `genome`, same settings as the
`COVERM_GENOME` module) before their abundances are compared with the new sample's; the pipeline does
not do that yet.

### Read-based profiling and the read-only snapshot

Two independent switches:

| flag | default | effect |
|---|---|---|
| `--profile_reads` | off | sketch the cleaned reads with [sylph](https://sylph-docs.github.io/) and profile them against `--sylph_db` (GTDB); [sylph-tax](https://sylph-docs.github.io/sylph-tax/) adds lineages (`--sylph_tax_id`, default `GTDB_r226`) |
| `--skip_assembly` | off | stop after preprocessing: no assembly, binning, refinement or GTDB-Tk |

`--profile_reads` alone is a MAG run that also reports what binning did not recover. Both together is the read-only taxonomic snapshot, which needs no assembly-grade memory:

```bash
nextflow run tpall/aftekas --input my_samples.csv --profile_reads --skip_assembly \
    --sylph_db databases/sylph/gtdb-r226-c200-dbv2.syl2db --sylph_tax_dir databases/sylph-tax
```

Output under `results/profile/` (shared by both assembler arms): `sylph_profile.tsv` (genome level: ANI, effective coverage, taxonomic and sequence abundance per sample), `taxprof/<sample>.sylphmpa` (MetaPhlAn-style lineage profiles), and `sylph_relative_abundance.tsv` / `sylph_sequence_abundance.tsv` (clade x sample matrices). Profiling reads the host/phiX-removed reads, not the BBNorm-normalised ones, because normalisation flattens the coverage the abundances are estimated from.
