#!/usr/bin/env python3
"""Behaviour tests for bin/genome_set.py: the identity rules of a client genome set.

    python3 tests/test_genome_set.py

Every rule here exists because breaking it would silently change a number already delivered to a
client. See docs/design-common-reference-set.md.
"""
import json, os, shutil, subprocess, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPT = os.path.join(HERE, '..', 'bin', 'genome_set.py')
ok = fail = 0


def case(label, cond):
    global ok, fail
    if cond:
        ok += 1
        print(f'ok   {label}')
    else:
        fail += 1
        print(f'FAIL {label}')


def run(work, clusters, sample, out, previous=None, set_id='CGS-T', ani='95'):
    cmd = [sys.executable, SCRIPT, '--clusters', os.path.join(work, clusters), '--set-id', set_id,
           '--ani', ani, '--sample', sample, '--fasta-dir', os.path.join(work, 'mags'),
           '--out', os.path.join(work, out)]
    if previous:
        cmd += ['--previous', os.path.join(work, previous)]
    p = subprocess.run(cmd, capture_output=True, text=True)
    return p.returncode, p.stdout + p.stderr


def members(work, version):
    with open(os.path.join(work, version, 'members.tsv')) as fh:
        head = fh.readline().rstrip('\n').split('\t')
        return [dict(zip(head, l.rstrip('\n').split('\t'))) for l in fh if l.strip()]


def meta(work, version):
    with open(os.path.join(work, version, 'set.json')) as fh:
        return json.load(fh)


def clusters(work, name, pairs):
    with open(os.path.join(work, name), 'w') as fh:
        for rep, mem in pairs:
            fh.write(f'mags/{rep}.fa\tmags/{mem}.fa\n')


work = tempfile.mkdtemp(prefix='genome_set_test_')
os.makedirs(os.path.join(work, 'mags'))
for g in ('s1_a', 's1_b', 's2_a_better', 's2_new', 's3_bridge'):
    with open(os.path.join(work, 'mags', f'{g}.fa'), 'w') as fh:
        fh.write(f'>{g}\nACGTACGTACGT\n')

# v1: two genomes from the first sample
clusters(work, 'c1.tsv', [('s1_a', 's1_a'), ('s1_b', 's1_b')])
rc, out = run(work, 'c1.tsv', 'S1', 'v1')
m1 = members(work, 'v1')
case('v1 builds from the first sample', rc == 0 and meta(work, 'v1')['version'] == 1)
case('v1 has both genomes as representatives', sorted(r['genome_id'] for r in m1) == ['CGS-0001', 'CGS-0002'])
case('v1 records the contributing sample', all(r['sample'] == 'S1' for r in m1))

# v2: a better MAG of genome 1 arrives with sample 2, and one genuinely new genome
clusters(work, 'c2.tsv', [('s1_a', 's1_a'), ('s1_a', 's2_a_better'), ('s1_b', 's1_b'),
                          ('s2_new', 's2_new')])
rc, out = run(work, 'c2.tsv', 'S2', 'v2', previous='v1')
m2 = members(work, 'v2')
byf = {r['fasta']: r for r in m2}
case('v2 adds only the genome that matched nothing', meta(work, 'v2')['added_this_version'] == ['CGS-0003'])
case('a better MAG does not replace the representative',
     byf['s1_a.fa']['is_representative'] == 'yes' and byf['s2_a_better.fa']['is_representative'] == 'no'
     and byf['s2_a_better.fa']['genome_id'] == 'CGS-0001')
case('the sequence behind an id does not change between versions',
     meta(work, 'v1')['representatives_sha256']['CGS-0001']
     == meta(work, 'v2')['representatives_sha256']['CGS-0001'])
case('the id of an existing genome keeps its original version_added',
     byf['s1_a.fa']['version_added'] == '1' and byf['s2_new.fa']['version_added'] == '2')
case('who found a MAG is not rewritten by a later round',
     byf['s1_a.fa']['sample'] == 'S1' and byf['s2_new.fa']['sample'] == 'S2')
case('a version that added a genome requires a remap', meta(work, 'v2')['remap_required'] is True)
case('the new version knows its parent', meta(work, 'v2')['parent_version'] == 1)

# v3: a MAG bridges two existing representatives. Both ids stay, and the bridge is reported.
clusters(work, 'c3.tsv', [('s1_a', 's1_a'), ('s1_a', 's1_b'), ('s1_a', 's3_bridge')])
rc, out = run(work, 'c3.tsv', 'S3', 'v3', previous='v2')
m3 = members(work, 'v3')
b3 = meta(work, 'v3')
case('a bridge does not merge two ids',
     rc == 0 and {r['genome_id'] for r in m3 if r['is_representative'] == 'yes'}
     == {'CGS-0001', 'CGS-0002', 'CGS-0003'})
case('a bridge is recorded for review', len(b3['bridged_representatives']) == 1
     and b3['bridged_representatives'][0]['kept'] == 'CGS-0001')
case('a bridge is loud on stderr', 'bridged representatives kept apart' in out)
case('the bridging MAG joins one id only',
     len([r for r in m3 if r['fasta'] == 's3_bridge.fa']) == 1)
case('no fasta is listed under two ids',
     len({r['fasta'] for r in m3}) == len([(r['genome_id'], r['fasta']) for r in m3]))

# guards: a different set id or a different threshold is a different set, not a new version
rc, out = run(work, 'c2.tsv', 'S2', 'vx', previous='v1', set_id='CGS-OTHER')
case('a different set id stops', rc != 0 and 'never mixed' in out)
rc, out = run(work, 'c2.tsv', 'S2', 'vx', previous='v1', ani='97')
case('a different ANI threshold stops', rc != 0 and 'new set id' in out)

# every representative in the set has a fasta on disk, named by its id
reps = {r['genome_id'] for r in m3 if r['is_representative'] == 'yes'}
case('every representative has a fasta named by its id',
     all(os.path.exists(os.path.join(work, 'v3', 'representatives', f'{g}.fa')) for g in reps))

# v2 as REFSET really stages it: the stored representatives arrive named by their id (CGS-0001.fa,
# from <previous>/representatives/), not by the MAG they came from. The first real extend
# (2026-10-13) gave every stored genome a second id this way; the cases above name them by MAG.
clusters(work, 'c2s.tsv', [('CGS-0001', 'CGS-0001'), ('CGS-0001', 's2_a_better'),
                           ('CGS-0002', 'CGS-0002'), ('s2_new', 's2_new')])
rc, out = run(work, 'c2s.tsv', 'S2', 'v2s', previous='v1')
m2s = members(work, 'v2s')
case('staged by id: only the new genome is added', rc == 0 and meta(work, 'v2s')['added_this_version'] == ['CGS-0003'])
case('staged by id: the set counts each genome once', meta(work, 'v2s')['genomes'] == 3)
case('staged by id: no member row for the staged file names',
     not any(r['fasta'].startswith('CGS-') for r in m2s))
case('staged by id: the new MAG joins the stored genome',
     [(r['genome_id'], r['is_representative']) for r in m2s if r['fasta'] == 's2_a_better.fa'] == [('CGS-0001', 'no')])
case('staged by id: the stored representative keeps its history row',
     [(r['genome_id'], r['sample'], r['is_representative']) for r in m2s if r['fasta'] == 's1_a.fa'] == [('CGS-0001', 'S1', 'yes')])

shutil.rmtree(work)
print(f'{ok}/{ok + fail} passed')
sys.exit(1 if fail else 0)
