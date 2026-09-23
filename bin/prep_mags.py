#!/usr/bin/env python3
"""Rename one sample's MAGs so they are unique across the series, and restate Binette's quality
report in CheckM2 column names under those new names.

Bin names collide between samples: every sample has a "bin_1". Clustering staged files from several
samples without renaming would merge two different genomes into one cluster member list, which is
the one error the genome set exists to prevent. Stdlib only.
"""
import argparse, csv, os, shutil, sys

FASTA_EXT = ('.fa', '.fna', '.fasta')
# Binette and CheckM2 both carry these, under different names
COMPLETENESS = ('completeness', 'Completeness')
CONTAMINATION = ('contamination', 'Contamination')
NAME = ('bin_id', 'name', 'Name', 'genome', 'Bin Id')


def pick(header, options, path):
    for o in options:
        if o in header:
            return o
    sys.exit(f'!! {path}: no column among {options} (has {header})')


a = argparse.ArgumentParser()
a.add_argument('--sample', required=True)
a.add_argument('--bins', required=True)
a.add_argument('--quality', required=True)
a.add_argument('--out', required=True)
a.add_argument('--quality-out', required=True)
a = a.parse_args()

os.makedirs(a.out, exist_ok=True)
renamed = {}
for f in sorted(os.listdir(a.bins)):
    stem, ext = os.path.splitext(f)
    if ext not in FASTA_EXT:
        continue
    new = f'{a.sample}__{stem}'
    shutil.copyfile(os.path.join(a.bins, f), os.path.join(a.out, f'{new}.fa'))
    renamed[stem] = new
if not renamed:
    sys.exit(f'!! {a.bins}: no MAG fasta for sample {a.sample}')

with open(a.quality) as fh:
    rows = list(csv.DictReader(fh, delimiter='\t'))
if not rows:
    sys.exit(f'!! {a.quality}: empty quality report')
head = rows[0].keys()
kn, kc, kx = pick(head, NAME, a.quality), pick(head, COMPLETENESS, a.quality), pick(head, CONTAMINATION, a.quality)

def contig_stats(path):
    """Genome size and contig N50, read from the fasta itself rather than trusted from elsewhere."""
    lengths, n = [], 0
    with open(path) as fh:
        for line in fh:
            if line.startswith('>'):
                if n:
                    lengths.append(n)
                n = 0
            else:
                n += len(line.strip())
    if n:
        lengths.append(n)
    total = sum(lengths)
    half, run = total / 2, 0
    for L in sorted(lengths, reverse=True):
        run += L
        if run >= half:
            return total, L
    return total, 0


# galah reads this file with the checkm crate, which wants the full CheckM2 column set: a
# three-column table makes it panic on an unwrap. Completeness and contamination come from
# Binette, size and N50 are measured here, and the rest are the constants CheckM2 writes.
CHECKM2_COLS = ('Name', 'Completeness', 'Contamination', 'Completeness_Model_Used',
                'Translation_Table_Used', 'Coding_Density', 'Contig_N50', 'Average_Gene_Length',
                'Genome_Size', 'GC_Content', 'Total_Coding_Sequences', 'Additional_Notes')

written = 0
with open(a.quality_out, 'w') as out:
    out.write('\t'.join(CHECKM2_COLS) + '\n')
    for r in rows:
        stem = os.path.splitext(os.path.basename(r[kn]))[0]
        if stem not in renamed:
            continue                      # a bin the refinement dropped: no fasta, no quality row
        size, n50 = contig_stats(os.path.join(a.out, f'{renamed[stem]}.fa'))
        out.write('\t'.join((renamed[stem], str(r[kc]), str(r[kx]), 'Neural Network (Specific Model)',
                              '11', 'NA', str(n50), 'NA', str(size), 'NA', 'NA', '')) + '\n')
        written += 1

missing = set(renamed) - {os.path.splitext(os.path.basename(r[kn]))[0] for r in rows}
if missing:
    sys.exit(f'!! {a.sample}: {len(missing)} MAG(s) without a quality row ({sorted(missing)[:3]}). '
             f'galah would pick representatives blind, so this stops here.')
print(f'{a.sample}: {len(renamed)} MAGs prepared, {written} quality rows')
