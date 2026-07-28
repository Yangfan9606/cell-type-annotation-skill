---
name: cell-type-annotation
description: Interpret single-cell RNA-seq cluster-annotation outputs and decide each cluster's cell type. Reads the results of the sc_004 annotation pipeline (consensus table + threshold log + cluster_resolve logs + marker_viz summary), applies a margin / vote-split / QC / coexp decision logic, and proposes a per-cluster label with evidence. Also drives the sc_006 marker/DEG/heterogeneity compute + plotting family (consensus one-vs-all markers, FindMarkers A-vs-B + pseudobulk DESeq2, subcluster drill-down, and 11 standalone plot scripts). Generates the exact Rscript commands to run next (does NOT auto-run them). Also reads clustree logs to decide whether to raise clustering resolution. Universal (human/mouse, any tissue); ships a built-in CNS marker library. Use when the user wants to annotate clusters, interpret sc_004 / cluster_resolve / thresholds logs, resolve "Unknown" clusters, compute subtype markers / condition DEGs / pseudobulk / subcluster heterogeneity, generate the matching volcano/dotplot/heatmap/violin/dimplot, or judge clustering resolution.
---

# Cell-type annotation (interpret + command-generation)

This skill **reads sc_004 pipeline outputs and decides how to annotate each cluster.** It is an
**interpretation + command-generation** skill — it never runs the R scripts itself. It tells you
exactly which command to run, then reads the resulting log/table and gives a per-cluster verdict
with the evidence behind it.

> **You are responsible for reviewing the marker genes.** The built-in marker library
> (`references/cell_markers.cns.txt` + `marker_library.md`) is a starting point with cited traps and
> literature references — not a substitute for checking that each panel is correct for your species,
> tissue, and the specific biology you're studying. Always verify a call against the *specific* genes
> before trusting it (see `references/decision_logic.md` §3), and correct/extend the marker library
> when you find something wrong — see "Self-iteration" below.

Two jobs:
1. **Decide cluster identity** from the annotation engine's consensus + threshold log, and from
   targeted `cluster_resolve` / `marker_viz` outputs for ambiguous clusters.
2. **Decide clustering resolution** from clustree logs — whether to split finer for sharper markers.

## Architecture (read this first)

The skill = **reusable decision logic + a marker library + the R engine scripts**.

```
skill (this dir, reusable):
  SKILL.md                       this file — when to use + workflow
  README.md                      install, usage, examples, FAQ
  CHANGELOG.md                   version history of the skill itself
  ITERATION_LOG.md               迭代反馈报告 — every problem + diagnosis + fix, growing over time (the skill's self-iteration mechanism)
  references/
    marker_library.md            the built-in marker library (CNS human + mouse homolog), dated, how to extend
    cell_markers.cns.txt         the actual seed marker file (feed to sc_004 as markers=)
    decision_logic.md            HOW to read the outputs: margins, vote-split, presence_floor, coexp, QC, marker promiscuity
    scripts_reference.md         the R scripts: paths, args, outputs
    report_template.md           shape of the final per-cluster annotation report
  examples/
    example_cns_dataset.md       full worked example (illustrative CNS dataset, 20 clusters, every decision type)
  scripts/
    sc_004.CellType_Annotation.R      annotation engine (MOD+UCell, z-score argmax, margin threshold)
    sc_004.02.cluster_resolve.R       ambiguous-cluster deep-dive (markers/anno/pair/sub/coexp + summary TSV)
    sc_004.03.marker_viz.R            DotPlot/Violin + *_summary.txt numeric tables
    sc_004.04.recode_celltype.R       recode a cluster column into a proper cell_type factor from a map file
    sc_005.SubClass_extract.R         clean cell/cluster extraction + plain re-cluster/UMAP (no integration)
    sc_006.find_marker.DEG.Heter.R    subtype marker / DEG / heterogeneity compute (4 switches; see scripts_reference §4)
    sc_006.plot.<job>.<type>.R        11 standalone plot scripts (find_marker/deg/pseudobulk/hetero × their plots)
    sc_006.plot.dotplot_engine.R      shared 3-style DotPlot engine (sourced by both dotplot wrappers)
```

Adapting this to your own project: point the commands you generate at wherever you deploy
`scripts/` (a fixed absolute path is normal — most users symlink or copy `scripts/` next to their
data pipeline and hardcode that path in the commands they run). `scripts_reference.md` documents the
exact flags each script takes; update its path block if you move the scripts.

## Workflow

**Step 0 — orient.** Confirm species (human/mouse) and tissue. Switch the QC regexes for mouse
(`^mt-`, `^Rp[ls]`, `^Hb`) and use mouse gene symbols. Point at the marker file in use, and — if this
is a new tissue/species — **review the built-in marker panels against your own literature before
trusting them** (the CNS library ships with citations and known traps, but every tissue has its own).

**Step 1 — read the first-pass result.** Read `*_annotation_consensus.csv` + `*_thresholds_log.txt`.
For every cluster note: `MOD_call`, top1/top2 ids + z, `MOD_margin`, reject reason, and whether UCell
agrees. Apply `references/decision_logic.md`:
- A clean call (margin ≥ threshold, both methods agree) → accept.
- `Unknown` from a **near-neighbour vote-split** (top1 is the right lineage, strong z, top2 is a
  related type) → hard-assign top1 if biology + UCell support it; flag the splitter.
- `Unknown` from a **genuine tie** (top1 ≈ top2, both moderate) → needs a `cluster_resolve` run.
- Watch the known **marker-promiscuity traps** (VCAN/LHFPL3 → false OPC on neurons; C1Q → false
  myeloid on ambient; CSPG4 → pericyte vs OPC). See decision_logic.md.

**Step 2 — generate the resolve commands.** For each genuinely-ambiguous cluster, emit a ready-to-run
`nohup Rscript sc_004.02.cluster_resolve.R … > <named>.log 2>&1` command (see scripts_reference.md).
Prefer `coexp='lineageA|lineageB'` to adjudicate ties; `sub=` to split suspected mixtures; `pair=`
for within-lineage splits; `anno=` for pseudobulk similarity. Name each log after its clusters.

**Step 3 — read the resolve outputs and decide.** Read the per-cluster log (top40 pos / top15 neg /
QC / coexp) **and** the new `*_resolve_summary.tsv` (one machine-readable row per cluster). Give each
cluster a verdict + the decisive evidence. For finer numeric checks use `marker_viz`'s `*_summary.txt`
(gene×cluster pct.exp + scaled mean + QC medians).

**Step 4 — resolution decision (optional).** If asked, read the clustree log: per-resolution cluster
sizes + adjacent-resolution split contingency tables. Messy off-diagonal reshuffling = over-splitting
a continuous manifold → do NOT raise resolution globally; fix via markers + targeted `sub=` instead.

**Step 5 — report.** Produce the per-cluster annotation table per `references/report_template.md`:
cluster → label → confidence → evidence, plus a list of marker-file edits the run motivated.
**Before presenting, self-check every row** (see report_template.md "Self-check" — this is mandatory,
not optional): re-query the source CSV keyed by that row's own cluster id and confirm the printed
logFC/pct.1/pct.2 for every cited gene actually belongs to that cluster. Adjacent-cluster row swaps
(writing cluster N's evidence next to cluster M's label) are the single most common failure mode when
hand-building a multi-cluster table from looped output — they look completely plausible and pass a
casual read, so only a keyed re-query catches them. See `ITERATION_LOG.md` (cl8↔cl10 swap entry).

**Step 6 — iterate.** Whenever a verdict is wrong or a trap is found, append a dated entry to
`ITERATION_LOG.md` (problem → diagnosis → fix) and, if it's a marker issue, edit the marker library.
This is how the skill improves with use — see "Self-iteration" below.

## Alternate track — sc_006 marker / DEG / heterogeneity (compute + plot)

When the object is already subset/annotated and the user wants **subtype markers, condition DEGs,
pseudobulk, or subcluster heterogeneity** (not first-pass whole-object annotation), use the sc_006
family instead of / alongside sc_004. Full flags + traps in `references/scripts_reference.md §4`.

- **Pick the switch(es)**: `--do_find_marker` (annotation-grade markers, all groups) · `--do_deg`
  (A-vs-B) · `--do_pseudobulk` (sample-level DESeq2, needs `--sample_column`) · `--do_hetero`
  (drill one group). They are independent and combinable.
- **Set "who vs who" with `-c` + `--subset_value`**: the usual ask "compare condition A vs B within
  one cell type" = `--subset_column cell_type --subset_value c1 -c condition --contrasts 'A,B'`
  (e.g. `--contrasts 'Disease,Control'`).
- **Generate `nohup Rscript … > named.log 2>&1` commands** (per scripts_reference §4a), then read the
  `*_find_marker_strict_final.csv` / `*_DEG_*_vs_*.csv` / `*_resolve`-style outputs and apply the same
  decision logic. For figures, hand over the matching standalone `sc_006.plot.<job>.<type>.R` commands
  (dotplot for annotation, volcano for DEG — §4b) — compute once, plot/edit separately.
- Do **not** re-trip the baked-in traps (SCT PrepSCTFindMarkers, `--outfile` not `--out`, geom_violin
  instead of VlnPlot on ggplot2≥4.0, optparse must be installed) — see §4c.

## Self-iteration

`ITERATION_LOG.md` is a living log of every problem this skill has hit — a mistaken verdict, a
marker-promiscuity trap, a script bug — with **Problem → Diagnosis → Fix → Where it landed** for each.
It is not documentation you write once; it is meant to grow every time you (or the assistant) use this
skill and find something wrong. Two rules make this work:
1. **Read it on any new problem** — it's the first place to check whether a failure mode has already
   been diagnosed.
2. **Append to it whenever a verdict was wrong, a trap was found, or a marker/script was changed** —
   include the deciding evidence (the actual numbers), not just the conclusion, so the next reader can
   verify it rather than take it on faith.

The marker library and decision-logic rules in this repo were themselves built this way — every trap
documented in `references/decision_logic.md` §3 started as a mis-call caught and logged here first.

## Hard rules (learned, do not relearn)
- **Never trust a single marker.** A panel "winning" can be one promiscuous gene. Always check the
  panel's *specific* members (e.g. OPC = real only if PDGFRA **and** CSPG4 fire, not VCAN/LHFPL3 alone).
- **Low absolute z across ALL panels = low quality, not an identity.** Check nFeature/%ribo.
- **`Proliferating` is a state, not a lineage.** Strip cell-cycle genes and read the underlying lineage.
- **High coexp double% is doublet only if BOTH panels are high.** One-sided high = shared markers.
- **SCT objects need `PrepSCTFindMarkers` before any FindMarkers**, including on subclustered subsets.
- **Never hand-transcribe a multi-cluster table without a keyed re-query pass.** Looping over clusters
  to gather evidence, then writing labels into a table by hand/memory, silently swaps adjacent rows
  (real incident: one cluster's fibroblast evidence written next to another cluster's label and vice
  versa — both looked internally consistent). Verify every row against the source keyed by its own
  cluster id before presenting.
