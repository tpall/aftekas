#!/usr/bin/env bash
#SBATCH --job-name=sylph-prep
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=4:00:00
#SBATCH --output=%x-%j.out
# One-off preparation for testing aftekas --profile_reads on rocket:
#   1. sylph GTDB r226 two-stage database (14 GB) -> $DB/sylph/
#   2. sylph-tax taxonomy metadata (~50 MB)       -> $DB/sylph-tax/
#   3. a 500k-pair subsample of one Gimenells pig-manure digester run
#      (SRR5079467, PRJNA356200, the plant behind the demo report) -> reads/
# Resumable: every step skips work that is already there.
set -euo pipefail
DB=${DB:-/gpfs/space/projects/preterm/databases}
# under sbatch the script runs from the slurm spool, so anchor on the submit dir
HERE="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
module load singularity/3.8.5 2>/dev/null || true

# 1. database: curl -C - resumes a partial download
mkdir -p "$DB/sylph"
f="$DB/sylph/gtdb-r226-c200-dbv2.syl2db"
if [ ! -s "$f.done" ]; then
  curl -L -C - -o "$f" https://storage.googleapis.com/sylph-stuff/gtdb-r226-c200-dbv2.syl2db
  ls -l "$f" && touch "$f.done"
fi

# 2. taxonomy metadata, via the same container the workflow uses
mkdir -p "$DB/sylph-tax"
# the metadata host answers 504 now and then; retry, and do not let it block the reads
if [ ! -s "$DB/sylph-tax/gtdb_r226_metadata.tsv.gz" ]; then
  for i in 1 2 3; do
    singularity exec -B "$DB/sylph-tax":/tax docker://quay.io/biocontainers/sylph-tax:1.9.1--pyhdfd78af_0 \
      sylph-tax --no-config --taxonomy-dir /tax download && break
    echo "!! sylph-tax download attempt $i failed; retrying in 60 s" >&2; sleep 60
  done || true
fi
ls -l "$DB/sylph-tax"

# 3. reads: first 500k pairs, taken from the head of the gzip streams so only
#    the first ~100 MB of each file is transferred; R1 and R2 stay in step
mkdir -p "$HERE/reads"
N=$((500000*4))
for r in 1 2; do
  out="$HERE/reads/SRR5079467_sub_${r}.fq.gz"
  [ -s "$out" ] && continue
  # head closes the pipe early, so curl and gzip exit non-zero by design: the
  # pipeline must not trip pipefail. The read count below is the real check.
  ( curl -sL -r 0-150000000 "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR507/007/SRR5079467/SRR5079467_${r}.fastq.gz" \
    | gzip -dc 2>/dev/null | head -n $N | gzip > "$out" ) || true
  n=$(( $(gzip -dc "$out" | wc -l) / 4 ))
  echo "reads $r: $n reads"
  [ "$n" -eq 500000 ] || { echo "!! expected 500000 reads in $out" >&2; rm -f "$out"; exit 1; }
done
printf 'sample,fastq_1,fastq_2\nGimenells_sub,%s/reads/SRR5079467_sub_1.fq.gz,%s/reads/SRR5079467_sub_2.fq.gz\n' "$HERE" "$HERE" > "$HERE/samples.csv"
cat "$HERE/samples.csv"
echo "done $(date)"
