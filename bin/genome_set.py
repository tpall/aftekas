#!/usr/bin/env python3
"""Build or extend a client's genome set from a galah cluster definition.

The set is the thing that makes a series comparable, so identity is the only rule that matters here:

  - a representative id is stable and is never reused;
  - an existing representative is never replaced by a "better" MAG from a later sample. Replacing the
    sequence would silently change every abundance already delivered. The better MAG is recorded as a
    member and waits for an explicit rebuild, which is a new set id;
  - a cluster that holds no existing representative gets a new id, and the set version increases;
  - two existing representatives that a new genome bridges into one cluster are NOT merged. Both ids
    stay, and the bridge is recorded for review. Identity beats parsimony.

See docs/design-common-reference-set.md. Stdlib only.

Usage:
  genome_set.py --clusters clusters.tsv --set-id CGS-DEMO --ani 95 \
      --new-mags mags/ --sample S2 [--previous genome_set/v3] --out genome_set/v4
"""
import argparse, hashlib, json, os, shutil, sys, time

REP_PREFIX = 'CGS'


def read_clusters(path):
    """galah writes 'representative<TAB>member', one line per member. Returns {rep: [members]}."""
    clusters = {}
    with open(path) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) < 2:
                sys.exit(f'!! {path}: cannot read cluster line {line!r}')
            rep, member = os.path.basename(parts[0]), os.path.basename(parts[1])
            clusters.setdefault(rep, []).append(member)
    if not clusters:
        sys.exit(f'!! {path}: no clusters')
    return clusters


def read_members(path):
    """members.tsv of the previous version: genome_id, mag, sample, version_added, is_representative."""
    rows = []
    if not path or not os.path.exists(path):
        return rows
    with open(path) as fh:
        header = fh.readline().rstrip('\n').split('\t')
        for line in fh:
            if line.strip():
                rows.append(dict(zip(header, line.rstrip('\n').split('\t'))))
    return rows


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def next_id(used):
    n = 1
    while f'{REP_PREFIX}-{n:04d}' in used:
        n += 1
    return f'{REP_PREFIX}-{n:04d}'


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--clusters', required=True, help='galah cluster definition over previous reps + new MAGs')
    p.add_argument('--set-id', required=True)
    p.add_argument('--ani', type=float, required=True)
    p.add_argument('--sample', required=True, help='sample that contributed the new MAGs')
    p.add_argument('--fasta-dir', required=True, help='directory holding every fasta named in the clusters')
    p.add_argument('--previous', help='previous version directory (absent: this builds v1)')
    p.add_argument('--out', required=True)
    p.add_argument('--tool-versions', default='', help='versions.yml content to record in set.json')
    a = p.parse_args()

    clusters = read_clusters(a.clusters)
    prev_members = read_members(os.path.join(a.previous, 'members.tsv') if a.previous else None)
    prev_set = {}
    if a.previous:
        with open(os.path.join(a.previous, 'set.json')) as fh:
            prev_set = json.load(fh)
        if prev_set.get('set_id') != a.set_id:
            sys.exit(f"!! previous set is {prev_set.get('set_id')}, not {a.set_id}: two set ids are "
                     f"never mixed. A rebuild is a new set id, not a new version.")
        if float(prev_set.get('ani', a.ani)) != a.ani:
            sys.exit(f"!! previous set clustered at ANI {prev_set['ani']}, this run at {a.ani}: "
                     f"changing the threshold makes a new set id, not a new version.")

    # fasta file name -> existing genome id, for every representative the set already has
    fasta_to_id = {r['fasta']: r['genome_id'] for r in prev_members if r['is_representative'] == 'yes'}
    # REFSET stages the stored representatives from <previous>/representatives/, where each is named
    # by its id (CGS-0001.fa), not by the MAG it came from. Those names are the same genomes: map them
    # to their id, and never give them a member row of their own (the historical row is carried
    # forward below). Without this every stored genome came back under a second id (2026-10-13).
    staged = {f"{r['genome_id']}.fa": r['genome_id'] for r in prev_members if r['is_representative'] == 'yes'}
    fasta_to_id.update(staged)
    prev_by_fasta = {r['fasta']: r for r in prev_members}
    added_in = {r['genome_id']: r['version_added'] for r in prev_members}
    used_ids = {r['genome_id'] for r in prev_members}
    version = int(prev_set.get('version', 0)) + 1

    members, added, bridged = [], [], []
    for rep, mem in sorted(clusters.items()):
        known = [f for f in mem if f in fasta_to_id]
        known_ids = sorted({fasta_to_id[f] for f in known})
        if known_ids:
            gid = known_ids[0]                      # the oldest id wins; the sequence does not change
            rep_fasta = next(f for f in known if fasta_to_id[f] == gid)
            if len(known_ids) > 1:
                # a new genome bridged two existing representatives: keep both, never merge silently
                bridged.append({'kept': gid, 'also_in_cluster': known_ids[1:],
                                'reason': 'a later MAG falls within ANI of both'})
        else:
            gid = next_id(used_ids)
            used_ids.add(gid)
            rep_fasta = rep
            added.append(gid)
        for f in mem:
            if f in staged:
                continue
            known_row = prev_by_fasta.get(f)
            # a MAG that is already a representative keeps its own id, even when the cluster holds
            # another one: a bridge must not absorb an id that abundances were reported against
            gid_f = fasta_to_id.get(f, gid)
            members.append({'genome_id': gid_f, 'fasta': f,
                            # the sample that contributed this MAG, kept from the history when the
                            # set has seen it before: a later round must not rewrite who found it
                            'sample': known_row['sample'] if known_row else a.sample,
                            # when the genome id entered the set, constant for the life of the id
                            'version_added': added_in.get(gid_f, str(version)),
                            'is_representative': 'yes' if (f == rep_fasta or f in fasta_to_id) else 'no'})

    # carry forward every historical member that this run did not see again
    seen = {(m['genome_id'], m['fasta']) for m in members}
    for m in prev_members:
        if (m['genome_id'], m['fasta']) not in seen:
            members.append(m)

    os.makedirs(os.path.join(a.out, 'representatives'), exist_ok=True)
    reps = sorted({(m['genome_id'], m['fasta']) for m in members if m['is_representative'] == 'yes'})
    for gid, fasta in reps:
        src = os.path.join(a.fasta_dir, fasta)
        if not os.path.exists(src) and a.previous:
            src = os.path.join(a.previous, 'representatives', f'{gid}.fa')
        if not os.path.exists(src):
            sys.exit(f'!! representative {gid}: no fasta at {src}')
        shutil.copyfile(src, os.path.join(a.out, 'representatives', f'{gid}.fa'))

    cols = ['genome_id', 'fasta', 'sample', 'version_added', 'is_representative']
    with open(os.path.join(a.out, 'members.tsv'), 'w') as fh:
        fh.write('\t'.join(cols) + '\n')
        for m in sorted(members, key=lambda x: (x['genome_id'], x['fasta'])):
            fh.write('\t'.join(m.get(c, '') for c in cols) + '\n')

    meta = {
        'set_id': a.set_id,
        'version': version,
        'parent_version': prev_set.get('version'),
        'ani': a.ani,
        'genomes': len(reps),
        'added_this_version': added,
        'sample_added': a.sample,
        'bridged_representatives': bridged,
        'remap_required': bool(added) or bool(bridged),
        'built': time.strftime('%Y-%m-%dT%H:%M:%S%z'),
        'tool_versions': a.tool_versions,
        'representatives_sha256': {gid: sha256(os.path.join(a.out, 'representatives', f'{gid}.fa'))
                                   for gid, _ in reps},
    }
    with open(os.path.join(a.out, 'set.json'), 'w') as fh:
        json.dump(meta, fh, indent=2, sort_keys=True)

    print(f"{a.set_id} v{version}: {len(reps)} genomes, {len(added)} added, "
          f"{len(bridged)} bridge(s) for review, remap {'required' if meta['remap_required'] else 'not required'}")
    if bridged:
        print('!! bridged representatives kept apart (review before the next report): '
              + ', '.join(b['kept'] for b in bridged), file=sys.stderr)


if __name__ == '__main__':
    main()
