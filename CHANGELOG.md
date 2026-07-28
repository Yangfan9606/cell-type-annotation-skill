# Changelog — cell-type-annotation

## v0.3.3 — 2026-07-28 (new script: check_markers_at_cluster_level.sh)
- New `scripts/check_markers_at_cluster_level.sh`: wraps `sc_006.find_marker.DEG.Heter.R
  --do_find_marker` + the per-cluster strict→sensitive top50 aggregation
  (`references/scripts_reference.md §4a`) into one command, for the common "run once, hand the
  results folder to the skill" workflow.
- Aggregation gains two new columns: `pct_diff` (pct.1 − pct.2) and `specificity_flag` — any positive
  marker whose background expression (pct.2) exceeds a configurable threshold (default 0.5) is tagged
  `"high_bg_check"` instead of being silently trusted, so the skill can weight it more skeptically
  before treating it as decisive evidence.
- `references/scripts_reference.md §4a` updated to document the wrapper script and the new columns.

## v0.3.2 — 2026-07-22 (new script: sc_004.04.recode_celltype.R)
- New `scripts/sc_004.04.recode_celltype.R` + `references/scripts_reference.md §3a`: recodes a
  cluster column into a proper `cell_type` factor from a required external TSV map file
  (cluster id / name / evidence, first 2 cols mandatory), with `Unassigned` fallback for unmapped
  clusters, a hard error on writing over an existing `celltype_col=` (unless `overwrite=true`), and
  typo warnings for map/data cluster-id mismatches.
- Verified end-to-end locally against a synthetic tiny Seurat object before handing off.

## v0.3.1 — 2026-07-22 (Pericyte ATP1A2 fix; SLC4A4 added to the marker seed)
- Literature audit found a likely copy-paste error in a common Pericyte panel: ATP1A2 has no
  pericyte citation and is instead a well-established astrocyte marker. Corrected to use
  **NDUFA4L2** instead — cited to Mesa-Ciller et al. 2023 *J Cereb Blood Flow Metab*
  (pericyte-vs-arterial-SMC discriminator, mouse+human).
- `references/cell_markers.cns.txt`: added SLC4A4 to the `Meningeal_fibro_pia` seed panel, cited to
  Pietilä et al. 2023 *Neuron* (human leptomeningeal fibroblast marker).
- `references/marker_library.md` + `ITERATION_LOG.md` updated accordingly.

## v0.3.0 — 2026-07-22 (literature audit methodology; Melanocyte panel narrowed)
- Ran a real (web-searched, not recalled) per-gene literature audit of a custom multi-type marker
  file. Findings folded into `marker_library.md`'s literature-audit section: several panels fully or
  partially sourced to published atlases (Vanlandewijck 2018 *Nature*, Pietilä 2023 *Neuron*), one
  likely copy-paste error found and fixed (see v0.3.1), and one activation-state panel set flagged as
  having no external literature source — documented as project-internal/data-derived rather than
  literature-backed.
- Melanocyte panel: added a narrower 3-gene consensus-run variant (`TYR,KIT,MITF`) alongside the
  existing 7-gene seed panel; documented that KIT should be checked against raw pct.exp before
  trusting it in any given dataset (it under-performs the other two genes in snRNA-seq).
- `references/cell_markers.cns.txt`: Melanocyte trap comment updated.

## v0.2.9 — 2026-07-22 (pair= vs coexp= distinction)
- `references/decision_logic.md` new §6c: a `pair=A:B` "relatively higher" gene set is not the same
  question as "this cluster's dominant identity" — always follow up with `coexp=` on the actual gene
  panels before promoting a pair= finding to a label. Worked example added showing a pair= result
  that looked like a nerve-associated identity call but was actually a 4%-of-cluster minority signal
  once checked with `coexp=` on the real gene panels.

## v0.2.8 — 2026-07-22 (meningeal fibroblast layer + melanocyte marker panels)
- `references/marker_library.md` + `cell_markers.cns.txt`: new `Meningeal_fibro_pia`
  (`SLC47A1,PTGDS,S100A6`) and `Melanocyte` (`TYR,TYRP1,DCT,PMEL,MLANA,OCA2,PAX3`) panels. SLC47A1
  documented as a near-exclusive single-gene pia-layer call in a worked example; PTGDS/S100A6/DCN/
  CST3/MGP/CLU flagged as a shared fibroblast-wide baseline (confirm, don't call, with them).
  Melanocyte call backed by literature (leptomeningeal melanocytes are a real population concentrated
  in the ventral cervical spinal cord + base of brain, extending into perivascular space) — not
  artifact/doublet.
- New note: antigen-presentation ("APC-like") fibroblast calls can't be settled from a one-vs-rest
  marker CSV alone; absence of HLA-II genes from the significance-filtered table means "not
  differentially high anywhere," not "absent" — needs a targeted dotplot on an explicit panel.

## v0.2.7 — 2026-07-22 (sc_005: stale reduction stripped from extracted object)
- Fixed a real latent bug in `sc_005.SubClass_extract.R`: Seurat's `subset()` preserves every existing
  `@reductions` entry, so a prior integration's reduction (e.g. `cca`) can ride unchanged into a freshly
  extracted subset object even when it's never actually read downstream — a footgun for anyone later
  plotting that reduction and assuming it reflects the subset. Fix: right after `subset()`, all
  inherited reductions are enumerated and stripped before the script's own pca/umap are computed.
- `references/scripts_reference.md §3b`: new trap paragraph documenting this.

## v0.2.6 — 2026-07-21 (sc_005 registered + R `$` partial-matching trap fixed)
- `references/scripts_reference.md` registers `sc_005.SubClass_extract.R` into the skill for the
  first time (clean keep/exclude cell-cluster extraction + plain re-cluster/UMAP, deliberately no
  integration).
- Fixed a real bug, reproduced and verified directly in R: base R's `$` operator partial-matches list
  element names, so `opt$celltype` silently resolved to `opt$celltype_column`'s value whenever
  `--celltype` wasn't passed — making `--exclude`-only invocations falsely trip a mutual-exclusivity
  error. All `opt$x` rewritten to `opt[["x"]]` (exact match) throughout the script. Documented as a
  **distinct** bug class from the `--out`/`--outdir` optparse CLI-abbreviation trap (§4c) — this one is
  list-indexing, not flag-parsing.

## v0.2.5 — 2026-07-21 (ambient-vs-reactive diagnostic via pair=)
- `references/decision_logic.md` new §4a: when a QC-flagged cluster might instead be a real reactive/
  activated state, run `pair=` against its closest clean same-lineage neighbour(s). Ambient
  contamination signature = "higher in" genes are the object's numerically-dominant population's
  markers at minority pct.1 (10-40%); real activation = lineage-appropriate activation genes at
  meaningful pct. Confirmed on a worked example (a QC-flagged EC cluster settled as `EC_lowQC`,
  ambient fibro soup, over `EC_reactive`).

## v0.2.4 — 2026-07-21 (coexp interpretation refinements)
- `references/decision_logic.md §6`: same-lineage coexp double-positive (e.g. SMC vs Pericyte, both
  mural) is a transitional/continuum state, not a doublet — only cross-lineage double-high is
  doublet-suspicious.
- `references/decision_logic.md` new §6b: negative one-vs-rest logFC across every canonical marker
  does **not** mean a cluster is unresolved/markerless — it means "less than crowded neighbours."
  coexp's raw %-positive is the correct tiebreaker and should run before declaring a cluster
  unresolved.

## v0.2.3 — 2026-07-21 (mandatory table self-check)
- `SKILL.md`: Step 5 now requires a self-check pass before presenting any multi-cluster table
  (re-query the source keyed by each row's own cluster id); new hard rule added to the list.
- `references/report_template.md`: new "Self-check" section, using a real row-swap incident as the
  worked example — two confidently-labeled rows had their evidence/label transposed during
  hand-transcription; every query behind them was actually correct.

## v0.2.2 — 2026-07-21 (strict-tier myAUC-NA per-cluster trap)
- `references/scripts_reference.md §4a`: reworded the "sensitive tier if strict is empty" note — it's
  a **per-cluster** fallback, not per-file. Documented a concrete failure where a naive aggregation
  one-liner silently dropped 10/24 clusters (myAUC entirely NA for those clusters) with no error.
  Added a ready-to-run replacement aggregation snippet with a per-cluster strict→sensitive fallback
  and a post-run sanity check.

## v0.2.1 — 2026-07-09 (dotplot 3-style engine)
- New `sc_006.plot.dotplot_engine.R`: one 3-style DotPlot engine, `source()`d by both
  `sc_006.plot.find_marker.dotplot.R` and `sc_006.plot.hetero.dotplot.R` (thin wrappers now).
- `--style A|B|C` (block/grouped/gap; default gap), curated-panel flags (`--marker_list`,
  `--baseline`, `--group_order`, `--group_map`), gene source priority `--features` > `--genes` >
  marker CSV, palette override `--cols`.
- Engine trap logged: style C + `--transpose` rebuilt natively (no `coord_flip`, which breaks facet
  free-scale dropping).
- All 11 plot scripts + engine: `--help`/zero-arg now print usage + per-script examples; only
  `optparse` loads up front (Seurat/ggplot2/… deferred to plot time) → usage returns ~0.2 s.

## v0.2.0 — 2026-07-07 (sc_006 family)
- Registered the `sc_006` marker/DEG/heterogeneity family (compute + 11 plot scripts) into the skill.
- `references/scripts_reference.md`: added §4 (sc_006 driver's 4 switches, the who-vs-who mental
  model, the 11 plot scripts by job, and traps: SCT prep, `--outfile` vs `--out`, geom_violin vs
  VlnPlot on ggplot2≥4.0, hetero UMAP ≥2 dims, optparse dependency).
- `SKILL.md`: description now covers subtype markers / DEG / pseudobulk / heterogeneity + plotting.
- Five real bugs found by running it end-to-end against a synthetic Seurat object, all fixed and
  documented in scripts_reference §4c (optparse prefix-match collision, ggplot2≥4.0/Seurat5.4
  VlnPlot crash, DoHeatmap needing populated scale.data, hetero UMAP needing ≥2 dims, pseudobulk
  normcounts read by the wrong key).

## v0.1.0 — 2026-06-24 (initial)
- Created the skill: interpretation + command-generation for the `sc_004` annotation pipeline.
- `SKILL.md`: 6-step workflow (orient → read first pass → generate resolve commands → read & decide →
  resolution decision → report) + 5 hard rules.
- `references/decision_logic.md`: scoring model, vote-split vs genuine-tie, marker-promiscuity traps
  (VCAN/LHFPL3/C1Q/CSPG4/PDGFRA), quality-before-identity, state-vs-lineage, coexp doublet semantics,
  clustree resolution logic, one-cluster decision flow.
- `references/marker_library.md` + `cell_markers.cns.txt`: CNS (brain+spinal) marker library, human +
  mouse homolog, dated, iterable; OPC = `PDGFRA,CSPG4,VCAN`; Mural split into pericyte + SMC; added
  Neuroblast / NPC_RadialGlia state panels.
- `references/scripts_reference.md`: the R scripts by path (referenced, not copied) + outputs.
- `references/report_template.md`: Chinese table-first per-cluster report.
- `examples/example_cns_dataset.md`: full 20-cluster worked example.
- `ITERATION_LOG.md`: seeded with the 6 problems/fixes that built the logic.

### Companion script changes (in `scripts/sc_004.02.cluster_resolve.R`)
- Fixed `sub=` DE on SCT (PrepSCTFindMarkers on the subset).
- coexp doublet warning made conditional (`dblwarn=`, default 20) + per-line lean tag.
- Added `*_resolve_summary.tsv` machine-readable per-cluster output.
- (marker_viz already emits `*_summary.txt`; OPC marker library tightened to drop LHFPL3.)
