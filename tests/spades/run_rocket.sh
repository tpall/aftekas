#!/usr/bin/env bash
# Wiring test for --assembler spades on rocket. Run from a login node inside tmux.
#
# What it checks, in order of what would actually hurt:
#   1. SPAdes ran with --meta
#   2. scaffolds and the assembly graph both came out
#   3. FIX_FASTA_HEADER left the NODE_* names alone
#   4. the graph's P (path) lines name the SAME sequences that were binned.
#      This is the assumption the whole graph evidence channel rests on: scaffolds.fasta
#      and contigs.fasta carry independent NODE_* numbering, so if the subworkflow were
#      emitting contigs instead of scaffolds, the names here would not line up.
#   5. the sorted BAM is named after the assembler, not hardcoded MEGAHIT
#
# Reads: 500k-pair subsample of SRR5079467, same one tests/sylph uses.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
module load singularity/3.8.5 2>/dev/null || true
command -v nextflow >/dev/null || { echo "!! nextflow not on PATH (conda activate java?)" >&2; exit 1; }

# --- reads ---------------------------------------------------------------------------
mkdir -p "$HERE/reads"
N=$((500000*4))
for r in 1 2; do
  out="$HERE/reads/SRR5079467_sub_${r}.fq.gz"
  [ -s "$out" ] && continue
  # head closes the pipe early, so curl and gzip exit non-zero by design
  ( curl -sL -r 0-150000000 "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR507/007/SRR5079467/SRR5079467_${r}.fastq.gz" \
    | gzip -dc 2>/dev/null | head -n $N | gzip > "$out" ) || true
  n=$(( $(gzip -dc "$out" | wc -l) / 4 ))
  echo "reads $r: $n"
  [ "$n" -eq 500000 ] || { echo "!! expected 500000 reads in $out" >&2; rm -f "$out"; exit 1; }
done

# --- run -----------------------------------------------------------------------------
# Launch from the repo root: Nextflow auto-loads nextflow.config from the launch
# directory, and that is where the 'cluster' profile lives (via conf/base.config).
# Running from this directory instead gives "Unknown configuration profile: 'cluster'".
# work/ is kept here so the test leaves nothing behind in the repo root.
WORK="$HERE/work"
cd "$ROOT"
nextflow run tests/spades/test_assembly.nf \
  -profile cluster -c tests/spades/test.config -resume -work-dir "$WORK" \
  --assembler spades --save_bams \
  --sample spades_test \
  --r1 "$HERE/reads/SRR5079467_sub_1.fq.gz" \
  --r2 "$HERE/reads/SRR5079467_sub_2.fq.gz" "$@"

# --- checks --------------------------------------------------------------------------
echo
echo "================ checks ================"
fail=0
pass() { echo "  PASS $1"; }
nope() { echo "  FAIL $1"; fail=1; }

tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

spades_cmd=$(grep -rl "spades.py" "$WORK"/*/*/.command.sh 2>/dev/null | head -1 || true)
if [ -z "$spades_cmd" ]; then
  echo "!! no SPAdes task found under $WORK" >&2; exit 1
fi
spades_dir=$(dirname "$spades_cmd")
echo "SPAdes task: $spades_dir"

if grep -q -- "--meta" "$spades_cmd"; then pass "spades.py called with --meta"
else nope "spades.py called with --meta"; fi

scaf=$(ls "$spades_dir"/*.scaffolds.fa.gz 2>/dev/null | head -1 || true)
gfa=$(ls "$spades_dir"/*.assembly.gfa.gz 2>/dev/null | head -1 || true)
[ -n "$scaf" ] && pass "scaffolds emitted" || nope "scaffolds emitted"
[ -n "$gfa" ]  && pass "assembly graph emitted" || nope "assembly graph emitted"
if [ -z "$scaf" ] || [ -z "$gfa" ]; then
  echo "cannot continue without both"; exit 1
fi

fixed=$(ls -t "$WORK"/*/*/*_fixed.fa.gz 2>/dev/null | head -1 || true)
if [ -n "$fixed" ]; then pass "FIX_FASTA_HEADER output found"
else nope "FIX_FASTA_HEADER output found"; echo "cannot continue"; exit 1; fi

if [ "$(gzip -dc "$fixed" | grep -c '^>NODE_' || true)" -gt 0 ]; then
  pass "fixed headers still NODE_*"
else nope "fixed headers still NODE_* (header rewriting would break the graph channel)"; fi

# the load-bearing check: do the graph's path names name the sequences that get binned?
gzip -dc "$gfa"   | awk '$1=="P"{print $2}'        | sort -u > "$tmp/pnames"
gzip -dc "$fixed" | awk '/^>/{print substr($1,2)}' | sort -u > "$tmp/seqnames"
np=$(wc -l < "$tmp/pnames")
ns=$(wc -l < "$tmp/seqnames")
shared=$(comm -12 "$tmp/pnames" "$tmp/seqnames" | wc -l)
echo "  P-line names: $np   sequence names: $ns   shared: $shared"
if [ "$np" -gt 0 ] && [ "$shared" -gt $(( np * 90 / 100 )) ]; then
  pass "graph path names match the binned sequence names (>90%)"
else
  nope "graph path names match the binned sequence names (>90%)"
  echo "    ^ the scaffolds-vs-contigs assumption in subworkflows/assembly.nf is wrong"
  comm -13 "$tmp/pnames" "$tmp/seqnames" | head -3 | sed 's/^/    only in graph: /'
  comm -23 "$tmp/pnames" "$tmp/seqnames" | head -3 | sed 's/^/    only in seqs:  /'
fi

bam=$(ls -t "$WORK"/*/*/metaSPAdes_*_sorted.bam 2>/dev/null | head -1 || true)
if [ -n "$bam" ]; then
  pass "sorted BAM named after the assembler ($(basename "$bam"))"
  # Look for the index itself, not beside $bam: downstream tasks stage the BAM alone,
  # so ls -t can land on a staged symlink whose directory has no index.
  bai=$(ls -t "$WORK"/*/*/metaSPAdes_*_sorted.bam.bai 2>/dev/null | head -1 || true)
  if [ -n "$bai" ]; then pass "BAM index present ($(basename "$bai"))"
  else nope "BAM index present (metaSPAdes_*_sorted.bam.bai)"; fi
else
  nope "sorted BAM named after the assembler (metaSPAdes_*)"
fi

echo "========================================"
if [ "$fail" = 0 ]; then echo "all checks passed"; else echo "SOME CHECKS FAILED"; exit 1; fi
