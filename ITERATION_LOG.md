# Iteration & feedback log

**This is the skill's self-iteration mechanism.** Living record of every problem this skill hit, the
diagnosis, and the fix. Append a dated entry whenever a verdict was wrong, a trap was found, or a
script/marker change was made. Newest on top. Each entry: **Problem → Diagnosis → Fix → Where it
landed.**

When you (Claude) finish an annotation round and something was non-obvious, add an entry here before
closing out. When the user reports a future problem, **this is the first file to read** and to extend.
The marker library and decision-logic rules in this repo were themselves built this way — nearly every
trap documented in `references/decision_logic.md` and `references/marker_library.md` started as a
mis-call caught and logged here first.

> This log is scoped to the scripts and decision logic actually shipped in this repo
> (`scripts/sc_004.*`, `sc_005.*`, `sc_006.*` + the reference docs). Entries have been generalized to
> remove identifying details of the specific research projects the skill was developed against —
> the technical lessons (bugs, traps, methodology) are unchanged and real.

---

## Illustrative round — subtype/state decomposition on a real astrocyte compartment

A 6-cluster one-vs-all subtype run separated cleanly into an **origin × activation** 2×2 structure
(two lineage/origin axes, each with its own homeostatic and activated member) rather than 4
independent types — full methodology and findings folded into `references/marker_library.md`
("Astrocyte subtypes/states"). Two process lessons worth keeping separate from the marker content:

**① A PDF "read" via text extraction has no visual data — describing dot size/color from it is
fabrication.** A rendered ggplot dotplot PDF, opened through a text-extraction tool, returns axis
labels/gene names/legend values as plain text; vector dot geometry (size, color, position) simply
isn't present in that output. A first pass nonetheless described specific per-dot visual properties as
if the image had been seen — caught on re-reading the tool output and noticing it was a flat text
block with no image marker (contrast with a real PNG read earlier in the same session, which did come
back with genuine pixel data). **Rule: before describing any per-element visual property from a Read'd
PDF, confirm the tool result actually contains image data, not just extracted text.** When this
happens, discard the visual claim explicitly to the user and re-derive the same conclusion from the
underlying numeric table instead.

**② Before trusting a `coexp` panel gene, check its own pct.2 (rest-of-object baseline).** A candidate
discriminator gene had pct.2 ≈ 0.85 — near-ubiquitous across the whole compartment — and inflated a
`coexp` double-positive% almost identically in an already-confirmed-clean cluster and the cluster
under test. **A coexp result that looks the same in a confirmed-clean cluster and the cluster under
test is itself evidence *against* a real discriminator/doublet signal** — a real signal should differ
between them. New diagnostic added to `references/decision_logic.md`.

**Where it landed:** `references/marker_library.md` (Astrocyte subtypes/states subsection); this entry.

---

## Command-generation mistakes across the sc_004/sc_006 script families

**Mistake 1 — generated commands from memory instead of re-reading `scripts_reference.md`'s flag
table.** Gave optparse-style flags (`-m`, `-c`) to a script that actually takes `key=value` args, and
used a non-existent flag name for a required argument. Both commands failed with usage errors before
the user caught it. **Fix: re-read the actual flag table in `scripts_reference.md` before generating
any command — don't reuse a remembered flag name across script families,** since sc_004 uses
`key=value` and sc_005/sc_006 use optparse-style flags.

**Mistake 2 — passed cluster IDs as `coexp=` gene panels.** Gave `coexp='16|1'` intending "compare
cluster 16 vs cluster 1," but `coexp=` takes two comma-separated **gene panels**, not cluster ids. The
script correctly caught this itself (a warning that panel1 had no genes found in the object) and
reported a meaningless 0%/0%/0% result — harmless because it self-flagged, but the whole coexp block
in that log was wasted. **Fix: `coexp=` panels must always be real gene symbols; use `pair=` (not
`coexp=`) when the comparison is cluster-vs-cluster.**

**Where it landed:** `references/scripts_reference.md`; this entry.

---

## sc_005.SubClass_extract.R: stale inherited reduction rides along into the extracted object

**Problem:** running `--exclude` against a CCA-integrated parent object, the output object printed
`3 dimensional reductions calculated: pca, umap, cca` — a prior integration step's `cca` embedding was
still present even though this run never used CCA.

**Diagnosis (traced through the script, not guessed):** Seurat's `subset()` subsets every existing
`@reductions` entry by cell rows but never drops one — so `cca` from the parent's integration run rides
into the subset object untouched. However `RunPCA()`/`RunUMAP()` both overwrite `pca`/`umap` from
scratch using only the extracted cells' own (unintegrated) data — neither ever reads the old `cca`, so
no computation is actually contaminated. But the old `cca` matrix is never removed — it survives as
inert baggage: correct numbers, wrong scope (fit across the pre-exclusion cell population), and would
look indistinguishable from a real per-subset integration to a later `DimPlot(reduction="cca")` call.
Same bug family as an earlier stale-SCT-assay issue in the same script: a slot computed on the parent
object outlives the subset without being recomputed or stripped.

**Fix:** right after `subset()`, enumerate `Reductions(subset_obj)` and strip every one of them
(`subset_obj[[rd]] <- NULL`) before the script's own pca/umap get computed later, logging what was
removed.

**Where it landed:** `scripts/sc_005.SubClass_extract.R`; `references/scripts_reference.md §3b`;
`CHANGELOG.md` v0.2.7.

---

## sc_005.SubClass_extract.R: R's `$` partial-matching silently swapped --exclude for --celltype

**Problem:** running the `--exclude` mode:
```bash
Rscript sc_005.SubClass_extract.R -i <input.rds> -o ./clean_extract \
  --celltype_column seurat_clusters --exclude 10,17,21,22,23,24,15 -r 0.3,0.5,0.8 -d 30
```
produced a "--celltype and --exclude cannot both be used" error — despite never passing `--celltype`.

**Diagnosis (reproduced directly in R, not guessed):** built the exact `option_list` + args in
isolation and confirmed `parse_args()` itself was correct — `opt` had no `celltype` element, only
`celltype_column`. Yet `is.null(opt$celltype)` returned `FALSE`. Root cause: **R's `$` operator on a
plain `list` performs unambiguous-prefix partial matching** (documented base R behavior, easy to
forget) — since `opt` had no exact `"celltype"` name, `opt$celltype` silently resolved to
`opt$celltype_column`'s value, making the mutual-exclusivity check see two non-NULL values. Verified
the fix mechanism: `opt[["celltype"]]` (double-bracket) uses exact matching and correctly returns
`NULL` in the same scenario. This is a **different bug class** from an existing `--out`/`--outdir`
optparse-CLI-abbreviation trap (that one is optparse matching typed flags at parse time; this one is R
silently matching *list element names* every time `$` is used afterward) — don't conflate the two.

**Fix:** mechanically rewrote every `opt$x` in the script to `opt[["x"]]`. Re-ran the exact failing
command's argument logic in isolation post-fix: mutual-exclusivity check now passes correctly.

**Generalizable rule:** any optparse script with one option name that is a plain prefix of another (e.g.
`celltype`/`celltype_column`, `out`/`outdir`, `dim`/`dims`) is at risk the moment the shorter-named
option is legitimately left unset — `opt$shorter` will silently return the longer option's value with
no warning. Use `opt[["x"]]` everywhere in this script family going forward, not just where a collision
is currently suspected — it's free and permanently closes the bug class.

**Where it landed:** `scripts/sc_005.SubClass_extract.R` (all `opt$` → `opt[["…"]]`);
`references/scripts_reference.md §3b`.

---

## Ambient contamination vs genuine reactive state — `pair=` against the closest clean neighbour

Asked whether a QC-flagged cluster could be a real reactive state instead of just low-quality. Answered
with `pair=` against its closest clean same-lineage neighbours rather than guessing from the one-vs-rest
table alone. Two checks settled it: (1) the flagged cluster's interferon/MHC-I markers sat in the *same
range* as its clean neighbours — not a distinctly activated cluster, just the object's baseline tone;
(2) genes "higher in" the flagged cluster relative to its neighbours were the object's dominant
population's ECM/structural genes, but only at minority pct.1 (10-41%) — ambient soup on a
low-complexity capture, not a coordinated activation program. New reusable technique added to
`references/decision_logic.md §4a`: when a QC-flagged cluster might instead be a real reactive state,
`pair=` it against its closest clean same-lineage neighbour(s) and check whether the "extra" genes are
the object's dominant-population markers at minority pct (ambient) vs lineage-appropriate activation
genes at meaningful pct (real).

**Where it landed:** `references/decision_logic.md` §4a.

---

## Multi-cluster table row-swap when hand-transcribing (self-check rule added)

**Problem:** in a multi-cluster annotation report, the user flagged two rows as apparently wrong.
Re-querying both clusters directly showed the table had their evidence and labels **swapped** — a
fibroblast cluster's row carried a myeloid label and vice versa.

**Diagnosis:** not a data or extraction bug — every query behind the table was correct (verified by
re-running the exact same queries). The error was introduced purely when hand-writing the final table
from looped per-cluster output: the two clusters were processed adjacent to each other and their rows
got transposed during transcription. Both swapped rows looked completely internally consistent (a real
fibroblast marker set, a real myeloid marker set) — nothing about reading the wrong row in isolation
would look wrong. This is the "vote-split looks plausible" trap, but at the transcription layer instead
of the biology layer.

**Fix:** added a mandatory **self-check pass** before presenting any multi-cluster table — re-query the
source CSV keyed by each row's own cluster id and confirm the cited genes' logFC/pct.1/pct.2 actually
belong there, for the whole table in one sweep (not spot-checks — this incident hit two of the most
*confident* rows, not an ambiguous one). New `SKILL.md` hard rule + Step 5 addendum;
`references/report_template.md` new "Self-check" section with this incident as the canonical example.

**Where it landed:** `SKILL.md` (Step 5 addendum + new hard rule), `references/report_template.md`
(new Self-check section), this log entry.

---

## Per-cluster myAUC=NA silently deletes whole clusters from the strict-tier top-N aggregate

**Problem:** an aggregation one-liner filtering out `myAUC=="NA"` rows before taking the top-N per
cluster, run on `sc_006.find_marker.DEG.Heter.R --do_find_marker` output for a 24-cluster dataset,
produced rows for only **14 of 24 clusters** — 10 clusters had zero rows entirely, not few, zero.

**Diagnosis:** `myAUC` comes from a `left_join` of the Wilcox-significant gene list against a
`FindAllMarkers(test.use="roc", features=cand)` result. ROC applies its own default `min.pct`/
`logfc.threshold` **per cluster**; when a cluster's Wilcox-significant genes are all weak-effect
(common for transcriptionally similar/transitional clusters — confirmed here: the affected clusters
were exactly the ones with continuous, closely-related identities), **none of them clear the ROC bar**,
so every row for that cluster gets `myAUC=NA`. A `$8!="NA"` filter applied *before* the per-cluster
top-N cut means a cluster with 0 non-NA rows contributes 0 rows to the aggregate — it doesn't even
appear as an empty group, it's just silently absent. Verified: the *sensitive* tier had non-NA myAUC
for all 24 clusters — confirming this is a strict-tier ROC-power artifact, not absence of real markers.

**Fix (workflow rule, no script change needed):** the existing "sensitive tier if strict is empty" note
was being read as *file-level* (whole strict file has 0 rows → use sensitive file). **It must be read
per-cluster**: before trusting any strict-tier aggregate, diff the cluster set in the aggregate against
the full cluster list from the object; any cluster missing or under a sanity floor (~10-15 rows) falls
back to the sensitive tier for *that cluster only*, not the whole run.

**Where it landed:** `references/scripts_reference.md §4a` (note reworded to say "per-cluster", not
"per-file", plus a ready-to-run replacement aggregation snippet); this log entry.

---

## Added Schwann_cell panel (MPZ,PRX,EGR2) to the CNS marker library

**Problem:** spinal-cord dissections that include nerve roots pick up PNS Schwann cells, which show up
as "Unknown" clusters if the marker library has no Schwann panel. Asked which genes from a 10-marker
literature Schwann-ID panel belong in this library, given it has to compete against the existing
Oligodendrocyte/OPC/Astrocyte panels under argmax.

**Diagnosis (per-gene, against this library's existing panels):** several candidates directly collide
with other panels already in the library — MBP and SOX10 with Oligodendrocyte (shared CNS+PNS myelin
gene, and the pan-oligodendroglial-lineage TF respectively), S100(B) with Astrocyte, SOX2 with
NPC_RadialGlia. GAP43/NCAM1 are broad axon-growth/adhesion genes, not Schwann-specific — confirmed by
finding GAP43 as a top marker of an unrelated Proliferating-OPC cluster in a real dataset.
NGFR(p75NTR)/POU3F1(Oct6) are real but state-restricted (immature/promyelinating only), left out of
the core panel to avoid diluting margin for myelinating Schwann cells.

**Fix:** added `Schwann_cell	MPZ,PRX,EGR2` to `cell_markers.cns.txt` + rationale/exclusions in
`marker_library.md` (PRX wasn't in the literature panel checked but was independently confirmed as a
clean, specific marker in a real worked cluster's own top-marker table, so added as a 3rd panel member
alongside MPZ).

**Where it landed:** `references/cell_markers.cns.txt`, `references/marker_library.md`.

---

## cluster_resolve summary files silently overwritten across runs + coexp= silent-fail trap

**① `_ambig_topmarkers.csv` / `_resolve_summary.tsv` overwritten on every separate `clusters=` run.**
- Problem: running `sc_004.02.cluster_resolve.R` three times on the same rds for different ambiguous
  cluster batches; each run's summary files are named only from the rds basename, not from
  `clusters=`, so the second and third runs silently clobbered the first's summary rows. Per-cluster
  `_cl<k>_markers.csv` is unaffected (named by cluster id).
- Fix: both write blocks (long-table + summary TSV) now read the existing file first, drop any row
  whose `cluster` matches the current run's clusters, `rbind` the fresh rows in, and write back — so
  batched runs on the same rds accumulate into one running table instead of clobbering.
- **Caveat (not fully fixed):** this read-modify-write is not atomic/locked. If multiple invocations on
  the same rds run concurrently, interleaved read-modify-write can lose rows. **Advise running batched
  `cluster_resolve` invocations on the same rds sequentially, not in parallel**, until/unless a lockfile
  is added.

**② `coexp=` silently returns 0%/0% when given cell-type/group names instead of gene symbols.**
- Problem: passing cell-type labels (e.g. `coexp='OPC|Oligodendrocyte'`) instead of comma-separated
  gene symbols. The gene-lookup function's `intersect(genes, rownames(expr))` came back empty, so it
  fell through to reporting every cluster as panel1+ 0%, panel2+ 0%, dbl 0%, with no error. Looked like
  a real (if confusing) biological result but was pure user-input mismatch; confirmed as a bug not
  biology because the same run's own FindMarkers table showed real expression for one of the panel's
  genes in a cluster the coexp result claimed was 0% positive.
- Fix: the gene-lookup function now takes a label and warns explicitly when zero genes match, plus a
  softer note when only some genes of a panel are missing.

**Where it landed:** `scripts/sc_004.02.cluster_resolve.R` (both fixes).

---

## DotPlot unified into a 3-style engine (find_marker + hetero)

**Problem:** wanted three dotplot looks (a red-block style with key/subtype vlines, a spectral style
with colored group bars, and a default-blue style with real gaps between groups) plus a gene-X↔Y
rotation switch, all as the default dotplot downstream of `sc_006.find_marker`.

**Fix / design:**
- New `sc_006.plot.dotplot_engine.R` holds one shared engine; both dotplot wrapper scripts
  `source()` it. One engine → no divergent copies. Engine must deploy alongside the wrappers.
- `--style A|B|C` (block/grouped/gap), default C. New default orientation genes-on-X @45° (was
  genes-on-Y); `--transpose` flips. Gene source priority `--features` > `--genes(file)` > CSV.
- Trap found while building: `coord_flip` defeats facet free-scale dropping (every panel showed all
  genes). **Style C + `--transpose` is rebuilt natively** from `DotPlot$data` with genes on Y and
  `facet_grid(rows, scales/space=free_y)` — no `coord_flip`. Verified all 3 styles × transpose on a
  synthetic 3-cluster object (block-diagonal enrichment recovered).

**Where it landed:** `scripts/sc_006.plot.dotplot_engine.R` (new) + both dotplot wrappers rewritten;
`references/scripts_reference.md §4d`; `CHANGELOG.md` v0.2.1.

---

## Skill created; three script fixes + OPC marker fix

**① sub= subclustering produced no DE on SCT.**
- Problem: subclustering one ambiguous cluster produced "No DE genes identified" + a
  "Run PrepSCTFindMarkers() before FindMarkers()" warning → subcluster markers silently empty.
- Diagnosis: the subset inherits multiple SCT models with unequal library sizes; FindAllMarkers needs
  PrepSCTFindMarkers re-run *on the subset*, which the script skipped.
- Fix: in `sc_004.02.cluster_resolve.R`'s sub block, run PrepSCTFindMarkers on the subset (if SCT)
  before FindAllMarkers.

**② coexp doublet warning misfired.**
- Problem: the doublet warning printed unconditionally — even at 0% double, and one real case (28%
  double) was shared stromal genes, not a doublet.
- Diagnosis: warning had no threshold and didn't require *both* panels high.
- Fix: per-line lean tag (→ leans panel1/2) + doublet warning only when double% exceeds a threshold
  (`dblwarn=`, default 20) **and** both panels exceed it.

**③ No machine-readable cross-cluster summary.**
- Problem: deciding all clusters meant reading long per-cluster log blocks by hand.
- Fix: `cluster_resolve` now writes `*_resolve_summary.tsv` (one row/cluster: QC + top5 + coexp%);
  `marker_viz` already writes `*_summary.txt` (gene×cluster pct.exp + scaled mean + QC). Read these
  first.

**④ OPC marker promiscuity (LHFPL3).**
- Problem: two clearly-GABAergic clusters (GAD1 64–82%) tied/called OPC.
- Diagnosis: the OPC signal came entirely from LHFPL3 (up to 82% in a GABAergic cluster) + VCAN
  (56–82%), while specific PDGFRA/CSPG4 were ≈ 0. LHFPL3 lights up neurons.
- Fix: OPC panel `PDGFRA,CSPG4,LHFPL3,VCAN` → `PDGFRA,CSPG4,VCAN`. Rule added: OPC is real only if
  PDGFRA AND CSPG4 fire.

**⑤ Wrong prior: "proliferating microglia".**
- Problem: predicted a cycling cluster was proliferating microglia.
- Diagnosis: the cluster's microglia markers were ≈0; top non-cycle markers were
  SOX2/PTPRZ1/PANTR1/NFIB → cycling neural progenitor instead.
- Fix: rule — `Proliferating` is a state; strip cycle genes + coexp lineages before naming. Logged as
  the canonical example in `examples/example_cns_dataset.md`.

**⑥ Earlier microglia reversal (pre-skill, kept for the record).**
- Problem: initially argued a cluster group were excitatory neurons, against the user's microglia
  hypothesis.
- Diagnosis: re-read negative-marker pct.2 + the CNS marker set; the clusters were unambiguous
  microglia (P2RY12/CSF1R/C1Q 90%+), a large fraction of the dataset (consistent with a
  neuroinflammatory/microgliosis condition).
- Fix: rule — trust the *specific* myeloid markers over a neuron prior; vote-split (Microglia vs
  Perivasc_Macro) explains the low margin.

---

## sc_006 marker/DEG/heterogeneity family added (compute + 11 plot scripts)

Built and registered `sc_006.find_marker.DEG.Heter.R` (4 combinable switches: find_marker / deg /
pseudobulk / hetero) + 11 standalone `sc_006.plot.<job>.<type>.R`. Validated end-to-end on a synthetic
Seurat object; the find_marker→dotplot chain correctly recovered injected per-group markers. Five real
bugs surfaced only by running it — all fixed, all now documented in `scripts_reference.md §4c`:

**① optparse prefix-match collision.** `--out` is a prefix of `--outdir`; passing `--outdir` populated
the wrong option → a downstream "supply filename with extension" failure.
- Fix: renamed the plot override flag to `--outfile`.

**② ggplot2 ≥4.0 + Seurat 5.4 breaks VlnPlot/RidgePlot** (an internal S7-dispatch error).
DotPlot/DimPlot/DoHeatmap unaffected.
- Fix: both violin scripts rebuilt with `FetchData` + `geom_violin` (env-independent, more editable).
  Env alt: pin ggplot2 to 3.5.2.

**③ DoHeatmap needs populated scale.data.** RNA assay without ScaleData → DoHeatmap errors.
- Fix: heatmap plots run `GetResidual` (SCT) or `ScaleData(features)` (RNA) before DoHeatmap.

**④ hetero UMAP with 1 significant PC** aborted the whole job (UMAP needs ≥2 dims), losing markers+rds.
- Fix: expand dims to ≥2, wrap RunUMAP in tryCatch (non-fatal); dimplot falls back to PCA.

**⑤ pseudobulk normcounts read by the wrong key** (row-name indexing grabbed a sample column) →
corrupted matrix.
- Fix: read by the `gene` column explicitly.

Design principle carried in: **compute once (CSV/rds/manifest) → plot separately** so plot style edits
never trigger recompute; object-based plots auto-detect the grouping column from the CSV labels so a
stale manifest ident can't mis-group.

**Where it landed:** `scripts/sc_006.*.R` (12 files); `references/scripts_reference.md §4`,
`SKILL.md` (description + referenced-scripts + alternate track), this log.

---

## Template for new entries
```
## YYYY-MM-DD — <short title>
**Problem:** what looked wrong / what the user reported.
**Diagnosis:** the actual cause (cite the deciding numbers).
**Fix:** marker edit / script change / logic rule.
**Where it landed:** file(s) changed.
```
