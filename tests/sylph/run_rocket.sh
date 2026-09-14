#!/usr/bin/env bash
# Real (non-stub) test of the read-only snapshot on rocket:
#   --profile_reads --skip_assembly on a 500k-pair Gimenells digester subsample.
# Run from a login node inside tmux after prepare_rocket.sh has finished.
# Expect: Methanothrix / Methanosarcinaceae among the archaea in
# results/profile/taxprof/Gimenells_sub.sylphmpa, since the demo report's
# Gimenells profile is Methanothrix-dominated.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
DB=${DB:-/gpfs/space/projects/preterm/databases}
module load singularity/3.8.5 2>/dev/null || true
command -v nextflow >/dev/null || { echo "!! nextflow not on PATH (conda activate java?)" >&2; exit 1; }
[ -e "$DB/sylph/gtdb-r226-c200-dbv2.syl2db.done" ] || { echo "!! sylph db not ready; run prepare_rocket.sh" >&2; exit 1; }
cd "$HERE"
# Host removal needs some bowtie2 index; the real human+phage index is a multi-GB
# download that this test does not need, so default to phiX as a stand-in host.
if [ -z "${HOST_FASTA:-}" ]; then
  HOST_FASTA="$HERE/reads/host_stub_phix.fa"
  [ -s "$HOST_FASTA" ] || curl -sL "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_001422.1&rettype=fasta" > "$HOST_FASTA"
  HOST_INDEX=""
fi
# host_index_url must be an empty string (falsy, so SETUP builds the index from the fasta): an empty CLI value is read as a flag (true)
# under the strict parser, so the run-specific params go through a YAML file.
cat > params.yml <<YML
input: $HERE/samples.csv
array_size: 2
queue: ${PARTITION:-main}
profile_reads: true
skip_assembly: true
sylph_db: $DB/sylph/gtdb-r226-c200-dbv2.syl2db
sylph_tax_dir: $DB/sylph-tax
gtdbtk_db: $DB/release226
checkm2_db: $DB/CheckM2_database/uniref100.KO.1.dmnd
host_fasta_url: $HOST_FASTA
host_index_url: "${HOST_INDEX:-}"
${PHIX:+phix_reference: $PHIX}
YML
nextflow run "$ROOT" -profile cluster -resume -params-file params.yml "$@"
echo; echo "== profile:"; head -5 results/profile/sylph_profile.tsv | cut -f1-6
echo "== archaea in taxprof:"; grep -E "d__Archaea" results/profile/taxprof/*.sylphmpa | grep -E "s__" | sort -t$'\t' -k2,2gr | head -10
