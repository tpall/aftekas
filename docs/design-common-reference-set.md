# Design: one common reference set for a sample series

_Decided 2026-09-23. Owner: TP. Status: design agreed, implementation started._

## The problem

`aftekas` today is a per-sample pipeline: each sample is assembled, binned and refined on its own, so
each sample produces its own MAG set. That is correct for a one-off health check and wrong for a
series.

Two samples from the same reactor give two genome sets that cannot be compared genome by genome:

- a genome "absent at T2" can be a genome that T2's assembly split in two, or that binning cut on a
  different boundary;
- the abundance of "the same" MAG is a share of a different sequence in each sample, because the two
  bins are not the same sequence;
- a route that "appeared" can be a marker gene that landed in an unbinned contig in the first sample.

Every quantity in the monitoring report (`magrittr-business/report/MONITORING-PLAN.md`) rests on
genome identity across samples. Without one reference set the report cannot be produced honestly, and
this is its first refusal condition.

## The decision

**One versioned genome set per client series. Built by dereplication, grown incrementally, never
rebuilt from scratch.**

1. Round 1: assemble and bin the first sample as today. Dereplicate its MAGs at 95 % ANI. That is
   genome set **v1**.
2. Every later round: assemble and bin the new sample as today. Cluster the new MAGs against the
   current set. A new MAG that falls inside an existing cluster does **not** enter the set: the
   existing representative stands for it. A new MAG that matches nothing becomes a new representative,
   and the set version increases.
3. After any version increase, **every earlier sample's reads are mapped again** against the new set.
   Mapping is cheap. Assembly and binning are not repeated.
4. Abundances always come from mapping reads to the set, never from the sample's own assembly.

### Why not co-assembly

Co-assembly of the whole series gives one genome set by construction and better bins for low-abundance
organisms. It fails on the property this product needs:

- every new sample forces a re-co-assembly of the whole series, so the cost per round grows with the
  series;
- every number already delivered can change, because the bins change. A client who received a report
  in March would get different March numbers in June, with no error on our side and no way to explain
  it;
- a strain mixture that shifts between samples produces chimeric bins, and the shift is exactly what
  the series is meant to detect.

Dereplication keeps the delivered numbers stable: a genome that was in v1 is the same sequence in v4,
so its abundance series is a series of the same thing. When v4 adds genomes, earlier samples get
abundances for them as well, from remapping — the series extends backwards, and no earlier number
moves except by the small amount that competition for multi-mapping reads changes it. That amount is
reported, not hidden (see "What the report must state").

Co-assembly stays the better choice for a one-off deep characterisation of a single reactor, and the
pipeline keeps supporting it through the existing per-sample path with all samples given as one
sample. It is not the subscription default.

## The genome set

A directory plus one table, versioned, one per client series:

```
<client>/genome_set/
  v4/
    representatives/            # one fasta per representative genome
    members.tsv                 # every MAG ever binned -> its representative
    quality.tsv                 # CheckM2 completeness/contamination per representative
    taxonomy.tsv                # GTDB-Tk per representative
    set.json                    # id, version, ANI threshold, tool versions, parent version, hash
```

Rules that keep identity meaningful:

- a representative id is stable and never reused. `CGS-0007` is the same sequence in every version of
  the set that contains it;
- a representative is never replaced by a "better" MAG from a later sample, even when the later MAG is
  more complete. Replacing the sequence would silently change every earlier abundance. A better MAG is
  recorded in `members.tsv` and the improvement waits for an explicit set rebuild, which is a new set
  id, not a new version;
- `members.tsv` keeps the whole history, so we can always say which sample contributed which MAG to a
  representative, and how many samples support it;
- a rebuild (new ANI threshold, a decision to re-derive representatives) creates a **new set id**. Two
  set ids are never compared in one report.

## Pipeline shape

Two new subworkflows, both after `BINREFINE`:

| subworkflow | tool | what it does |
|---|---|---|
| `REFSET` | galah (dRep-equivalent, reuses the CheckM2 stats Binette already produced) | clusters MAGs at 95 % ANI, picks representatives by CheckM2 score, writes `representatives/`, `members.tsv`, `set.json` |
| `ABUNDANCE` | CoverM genome mode (minimap2/bwa-mem2 + samtools) | maps every sample's host-removed reads to the representatives, writes relative abundance, mean coverage and **read counts** per genome per sample |

Read counts are not decoration: the noise floor in the monitoring plan needs them to decide when a
route is too thinly covered to call a change (Poisson limit).

New parameters:

- `--genome_set <dir>` — an existing set to extend. Absent: build v1 from this run's MAGs.
- `--genome_set_ani 95` — clustering threshold, recorded in `set.json`; changing it makes a new set id.
- `--remap_all` — map every sample in the set's history, not only this run's samples. Default true
  when the set version increased.

GTDB-Tk and the digestome panel then run **once per representative**, not once per sample. Route calls
become identical across samples by construction, which is the whole point: a route difference between
two samples can then only come from abundance, never from a different genome.

## What the report must state

- `abundance_basis`: "reads mapped to client genome set `<id>` v`<n>` (`<k>` genomes)";
- the set id and version in provenance, with the hash of `set.json`;
- how many genomes the version added, and that earlier samples were remapped;
- for a genome added in a later version: that it was not in the set when the earlier report was
  written. Its earlier abundance is real (the reads were always there), but it was not reportable
  then. This is the honest form of "appeared", and it replaces the censoring wording for that case.

The over-time render refuses when two samples carry different set ids, different versions without a
remap, or a different ANI threshold.

## Limits to write down now

1. **Multi-mapping reads.** Two representatives at 95–99 % ANI compete for the same reads. CoverM's
   handling is documented and deterministic, but a genome's share depends on what else is in the set,
   so adding a genome can move an earlier genome's share slightly. We report the magnitude on remap
   (the diff between old and new abundance for every genome) rather than claiming it is zero.
2. **95 % ANI is species level.** Strain succession inside a species is invisible by design. The
   report says so; strain-level tracking is a different product.
3. **A genome the series never assembles well** stays absent, however abundant. The read-level profile
   (sylph, already in the pipeline) is the cross-check, and a large gap between the two is a flag.
4. **The set is client-specific.** Genome sets are not shared between clients: it would leak one
   client's assemblies into another's report.

## Validation before this is used for a client

The two public BP3 samples (T1 and T2, Campanaro et al. 2020, PRJNA311380) are today exactly the
defective case: two independent MAG sets from two studies, and they are on the public demo pages.
Rebuild them on one set and compare:

1. how many representatives the two samples share, versus how many MAGs each carries now;
2. what happens to the finding the demo page reports ("one acetate-route genome at T2, at the
   detection floor"): the same genome below the cutoff at T1, a genuinely new genome, or a binning
   artefact;
3. how far the route fractions within the mcrA carriers move;
4. the remap diff from limit 1 above, measured rather than assumed.

Whatever the answer, the demo pages get corrected and the manuscript gets the comparison. If the
result is that per-sample binning invented a change, that is the strongest possible argument for this
work, and it belongs in the paper.

## Effort

Three to four working days, apart from the renderer's four:

| | days |
|---|---|
| galah + CoverM modules, `REFSET` and `ABUNDANCE` subworkflows, wiring and params | 1.5 |
| set store (`members.tsv`, `set.json`, id stability, incremental add, remap logic) | 1 |
| BP3 T1/T2 validation on HPC and the comparison table | 1 |
| METHODS.md, this document, tests | 0.5 |
