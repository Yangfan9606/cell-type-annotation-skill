# Scripts reference

The R scripts ship inside this repo under `scripts/`. Deploy that folder wherever your pipeline runs
and generate commands against that path — the examples below use `scripts/<name>.R` (relative to the
repo root); substitute your own deployment path when generating real commands.

```
scripts/sc_004.CellType_Annotation.R      # annotation engine (z-score/margin)
scripts/sc_004.02.cluster_resolve.R       # ambiguous-cluster deep-dive
scripts/sc_004.03.marker_viz.R            # DotPlot/Violin + numeric summary
scripts/sc_006.find_marker.DEG.Heter.R    # subtype marker / DEG / heterogeneity (compute)
scripts/sc_006.plot.<job>.<type>.R        # 11 standalone plot scripts (see §4)
scripts/sc_006.plot.dotplot_engine.R      # shared 3-style DotPlot engine (sourced by both dotplot wrappers; see §4d)
```

The skill's job is to **read their outputs and generate their commands**, never to reimplement them.
If you move the scripts to a different deployment path (e.g. a shared server `r_scripts/` directory),
update the path this file's commands use — the skill always regenerates commands from whatever path
you tell it, so keeping this file current is what keeps generated `nohup Rscript …` commands correct.

---

## 1. sc_004.CellType_Annotation.R — annotation engine

Scores marker panels per cell (AddModuleScore = MOD primary, UCell = reference), aggregates per
cluster, per-panel z-scores across clusters, argmax → top1, margin = top1_z − top2_z, adaptive margin
threshold, calls Unknown below it.

Auto-detects the cluster column (order: `seurat_clusters`, `integrated_snn_res.*`, `RNA_snn_res.*`).
Multi-assay/layer fallback; Assay5 `JoinLayers`; SCT `PrepSCTFindMarkers`. Permutation null with
stratified subsampling (cap 2000/cluster).

Key flags: `--assay`, `--ucell_assay`, `--presence_floor 0.4` (drop panels not expressed anywhere —
defeats the z-score "false winner" trap), `--min_markers`, `--hard_floor`.

**Outputs to read:** `*_annotation_consensus.csv` (per-cluster MOD/UCell call, top1/top2 z, margin,
reject reason, agreement, final `cell_type`) and `*_thresholds_log.txt` (margin distribution, the
three candidate thresholds + which won, and the list of Unknown clusters with their top1/top2).

---

## 2. sc_004.02.cluster_resolve.R — ambiguous-cluster deep-dive

`key=value` args; only `rds=` and `clusters=` required.

| arg | meaning |
|---|---|
| `rds=` | Seurat .rds (required) |
| `clusters=` | comma list of ambiguous clusters (required) |
| `clu=` | cluster column (auto: astro_clu / seurat_clusters / harmony_clusters) |
| `assay=` | default SCT else DefaultAssay |
| `anno=` | a celltype column → pseudobulk Spearman similarity of each ambiguous cluster vs named types |
| `pair=` | `14:1` or `14:1,18:6` — within-lineage A-vs-B FindMarkers |
| `sub=` | subcluster one cluster (re-runs PrepSCTFindMarkers on the subset) |
| `subres=` / `subdims=` | sub resolution (0.6) / dims (30) |
| `coexp=` | two comma-separated gene panels, pipe-joined (quote in shell) — %positive each + double%, with per-line lean tag |
| `topN=` | positive markers printed (40) |
| `neg=` | print negative markers (true) |
| `dblwarn=` | coexp double% above which a doublet warning fires (20) |

**Outputs to read:**
- the **log** (stdout): per-cluster QC line, top40 positive (by pct.diff), top15 negative, coexp block.
- `*_cl<k>_markers.csv` — full FindMarkers per cluster (audit; usually not needed, log top40 suffices).
- `*_ambig_topmarkers.csv` — long table of all ambiguous clusters' top markers.
- `*_resolve_summary.tsv` — **machine-readable one row per cluster** (cluster, n, nFeature, nCount,
  pct_mt, pct_ribo, pct_hb, qc_flag, top5_pos, coexp_p1, coexp_p2, coexp_dbl). **Read this first.**

**Command template (always nohup + named log):**
```bash
nohup Rscript scripts/sc_004.02.cluster_resolve.R rds=<rds> clu=<clusterCol> \
    clusters=<k1,k2> coexp=<panelA-pipe-panelB, quoted> \
    > sc_004.02.cluster_resolve.R.c<k1>_<k2>.log 2>&1
```
Name the log after the clusters so it is self-documenting.

---

## 3. sc_004.03.marker_viz.R — DotPlot / Violin + numeric summary

`key=value`; only `rds=` required.

| arg | meaning |
|---|---|
| `rds=` | Seurat .rds (required) |
| `markers=` | marker file → faceted DotPlot **and** the DotPlot numeric tables in the summary |
| `clu=` | cluster column (auto SCT_snn_res.0.5 else seurat_clusters) |
| `assay=` | SCT else DefaultAssay |
| `genes=` | comma list → per-gene violins + gene×cluster mean/%positive in the summary |
| `clusters=`, `outdir=`, `prefix=` | optional |

**Outputs to read:**
- `<prefix>_summary.txt` — **the text file to read** (no PDF needed): `[QC]` per-cluster medians;
  `[DotPlot]` gene×cluster `pct.exp` + `avg.exp.scaled` (only if `markers=`); `[genes]` gene×cluster
  mean + %positive (only if `genes=`). Wide, non-wrapping, clusters as columns.
- `<prefix>_dotplot.pdf` (if `markers=`), `<prefix>_qc_violin.pdf`, `<prefix>_genes_violin.pdf`.

**Global sweep command (one summary covering all clusters):**
```bash
Rscript scripts/sc_004.03.marker_viz.R rds=<rds> markers=<markerfile> clu=<clusterCol> \
    prefix=full outdir=<dir>
# then read <dir>/full_summary.txt
```

> Note: outputs land in `outdir` (default = dir of the rds). Pass `outdir=` explicitly to control it.

---

## 3a. sc_004.04.recode_celltype.R — cluster → cell_type recode from an external map file

`key=value` args (same family convention). Takes the final per-cluster identity table (built during
annotation) and writes it onto the object as a proper factor column, so downstream plot scripts
(`ident=cell_type`) show real labels instead of cluster numbers.

| arg | meaning |
|---|---|
| `rds=` | input Seurat .rds (required) |
| `clu=` | cluster column name (required) |
| `map=` | TSV mapping file, required (see format below) |
| `celltype_col=` | output column name [default: `cell_type`] — **errors if this column already exists**, unless `overwrite=true` |
| `overwrite=` | `true` to allow clobbering an existing `celltype_col=` |
| `out=` | output rds path [default: input filename + `_labeled.rds`, never overwrites the input] |

**Map file**: tab-separated, header row required (skipped), **first two columns mandatory** — 簇
(cluster id) and 命名 (name); an optional 3rd column (依据/evidence) is for human documentation only,
never parsed. Cluster ids present in the object but absent from the map become `Unassigned` (warned,
not an error — lets you label incrementally). Map rows whose cluster id doesn't exist in the object are
warned as likely typos, not fatal. Duplicate cluster ids in the map are a hard error. Factor level order
= first-appearance order in the map file, with `Unassigned` appended last if present.

```bash
Rscript scripts/sc_004.04.recode_celltype.R rds=my_object.rds \
  clu=seurat_clusters map=celltype_map.tsv
# -> my_object_labeled.rds, adds a `cell_type` factor column
```

Verified end-to-end against a synthetic tiny Seurat object: map/data cluster-id mismatch warnings,
`Unassigned` fallback, the `celltype_col=` conflict error (with and without `overwrite=true`),
missing-required-arg usage printout, and duplicate-cluster-id rejection all behave as specified.

---

## 3b. sc_005.SubClass_extract.R — clean cell/cluster extraction + simple UMAP (no integration)

Extracts a keep-list or exclude-list of cell-type/cluster values from `--celltype_column`, does a
**plain** re-normalize→PCA→cluster→UMAP pass (**no Harmony/CCA** — deliberately out of scope; run your
own integration script downstream for real integration), and saves `subclustered.rds` (fixed filename;
log records the run's timestamp for traceability) plus per-resolution UMAP PNGs.

Key flags: `-i/--input`, `-o/--output`, `--celltype_column` (which metadata column holds the type/cluster
labels), `-c/--celltype` (keep-list) **or** `-e/--exclude` (exclude-list, auto-resolved to
`setdiff(all values, exclude list)`) — mutually exclusive, exactly one required. `-r/--resolution`
(comma list), `-d/--dims`, `-k/--neighbors`, `-f/--features`, `--use_sct` (must be explicit; the only
automatic fallback is when the input object has no RNA assay at all).

```bash
nohup Rscript scripts/sc_005.SubClass_extract.R \
  -i <input.rds> -o <outdir> --celltype_column seurat_clusters \
  --exclude 10,17,21,22,23,24,15 -r 0.3,0.5,0.8 -d 30 \
  > sc_005.SubClass_extract.R.exclude.log 2>&1 &
```

**Trap: R's `$` operator on a `list` silently partial-matches**, and this is *independent of optparse's
CLI parsing* — it happens purely in R after `parse_args()` already returned a correct list. `opt$celltype`
will resolve to `opt$celltype_column`'s value whenever `--celltype` is not passed on the CLI (so `opt` has
no exact `"celltype"` element) and `"celltype"` is an unambiguous prefix of another present element name
(`"celltype_column"`). Reproduced directly: `list(celltype_column="x")$celltype` returns `"x"`,
`is.null(...)` is `FALSE` — with no warning, no error, indistinguishable from a real value. This is a
**different mechanism** from the historical `--out`/`--outdir` optparse CLI-matching trap in §4c below
(that one is optparse matching abbreviated flags at parse time; this one is base R silently matching *list
element names* every time `$` is used afterward) — don't conflate the fixes. Concretely: any script with
an option whose name is a plain prefix of another option's name (e.g. `celltype`/`celltype_column`,
`out`/`outdir`, `dim`/`dims`) is at risk **whenever the shorter-named option is legitimately left
NULL/unset** — the moment code checks `opt$shorter_name`, it silently reads the longer option's value
instead. **Fix applied throughout `sc_005.SubClass_extract.R`: every `opt$x` access rewritten to
`opt[["x"]]`** (double-bracket indexing on a list uses exact matching by default — verified
`list(celltype_column="x")[["celltype"]]` correctly returns `NULL`). Apply the same `opt[["x"]]` convention
in any *new* optparse-based script in this family, not just when a collision is suspected — it's zero-cost
and permanently closes this bug class regardless of what options get added later.

**Trap: stale inherited reductions survive the extraction.** Seurat's `subset()` keeps every existing
`@reductions` entry (pca/umap, plus anything from upstream integration like `cca`), just row-filtered to
the retained cells — it does not drop any of them. This script then computes a genuinely fresh
`pca`→`umap` on the extracted cells alone via `RunPCA()`/`RunUMAP()` (no `reduction=` argument given
anywhere, so both default-consume the freshly-computed `pca`, never the old one) — so the new pca/umap are
correct and not contaminated. But any *other* inherited reduction (typically `cca` from a prior
integration run on the full pre-exclusion object) is never touched or recomputed, so it rides along into
`subclustered.rds` unchanged — coordinates aren't wrong, just scoped to the wrong (pre-exclusion) cell
population, and indistinguishable from a real per-subset integration to a later `DimPlot(reduction="cca")`
call. Same class of issue as the stale-SCT-assay fix above (§3b flags/assay logic): a slot computed on the
*parent* object outlives the subset without being recomputed or removed. **Fix applied**: right after
`subset()`, `Reductions(subset_obj)` is enumerated and every one of them is stripped
(`subset_obj[[rd]] <- NULL`) before the fresh pca/umap are computed later — so the saved object's
`@reductions` only ever contains what *this specific run* produced.

## 4. sc_006 family — subtype marker / DEG / heterogeneity (compute + plotting)

For work on an **already-subset/annotated** object: find annotation-grade markers, contrast conditions,
or drill into one type. **Compute and plotting are separate** — compute once → CSV/rds/manifest, then
run standalone plot scripts (edit their CONFIG block, re-run, no recompute). `--help` or no-arg prints
a full usage table with examples.

### 4a. `sc_006.find_marker.DEG.Heter.R` — the compute driver (`getopt`/optparse flags)

Four **independent, combinable** switches (turn on any subset, incl. all four):

| switch | computes | output |
|---|---|---|
| `--do_find_marker` | one-vs-all consensus markers, **all groups** (Wilcox∩ROC→MAST) | `*_find_marker_{strict,sensitive}_final.csv` (+ wilcox) |
| `--do_deg` | `FindMarkers(A vs B)` per contrast | `*_DEG_<A>_vs_<B>.csv` (+ `_sig`) |
| `--do_pseudobulk` | DESeq2 on sample-aggregated counts | `*_pseudobulk_<A>_vs_<B>.csv` (+ `_normcounts`,`_coldata`) |
| `--do_hetero` | subcluster one group (re-PCA/cluster/UMAP) | `*_hetero_<X>_subcluster_markers.csv` + `*_hetero_<X>.rds` |

Every run also writes `<prefix>_manifest.txt` (rds/ident/outdir/prefix/contrasts/hetero_rds…) that the
plot scripts read automatically.

Common flags: `-i` rds, `-o` outdir (default `results`), `--prefix`, `-c/--ident` (grouping column),
`--assay`, `--subset_column`/`--subset_value` (scope to one population **before** contrasting),
`--contrasts 'A,B;A,C'` (`;`=contrasts, `,`=sides, right side `ALL`=one-vs-rest),
`--sample_column` (pseudobulk), `--compare_group` (hetero target), `--no_mast`, `--padj`, `--logfc`.

**"Who vs who" is set by `-c` + `--subset_value`** (this is the key mental model):
- `c1 vs c2` between types → `-c cell_type --contrasts 'c1,c2'`
- **`c1` internal condition-A vs condition-B (the usual one)** → `--subset_column cell_type
  --subset_value c1 -c condition --contrasts 'Disease,Control'`
- global condition-A vs condition-B → `-c condition --contrasts 'Disease,Control'`

Command templates (nohup + named log; SCT handled once at top):
```bash
# annotation markers (all groups)
nohup Rscript scripts/sc_006.find_marker.DEG.Heter.R -i <rds> -c cell_type --do_find_marker > fm.log 2>&1 &
# cell-type-specific condition DEG + pseudobulk (batch contrasts)
nohup Rscript scripts/sc_006.find_marker.DEG.Heter.R -i <rds> --subset_column cell_type --subset_value <T> \
  -c condition --do_deg --do_pseudobulk --sample_column <samp> --contrasts 'Disease,Control' > deg.log 2>&1 &
# heterogeneity drill-down of one type
nohup Rscript scripts/sc_006.find_marker.DEG.Heter.R -i <rds> -c cell_type --do_hetero --compare_group <T> > het.log 2>&1 &
```

### 4b. The 11 plot scripts `sc_006.plot.<job>.<type>.R`

The compute driver's `--plotdir` defaults to `./scripts` (same folder as itself); at the end of every
run it **prints all 11 ready-to-run plot commands** grouped by job, with `--outdir` pre-filled and the
jobs computed this run tagged as already done. So you don't hand-write plot commands — copy the lines
you want from the compute log (or write them into a `make_plots.sh`). If you deploy the plot scripts to
a different directory than the compute driver, pass `--plotdir` explicitly.

Standalone; default-read the newest `results/*_manifest.txt` + the job's CSV/rds; PNG+PDF out; each has
a CONFIG block (TOPN, thresholds, colors) at the top to tweak. Override with `--outdir --prefix
--input --rds --outfile`. **Every plot script prints usage+examples on `--help` or a zero-arg run**, and
loads only `optparse` up front (Seurat/ggplot2/… deferred until it actually plots) so usage returns
instantly (~0.2 s). **A plot exists only where it's the right plot for that job:**

| job | plot scripts | notes |
|---|---|---|
| find_marker | `dotplot` `heatmap` `violin` | **dotplot = annotation workhorse**; no volcano (no 2-sided contrast) |
| deg | `volcano` `heatmap` `violin` | volcano = the DEG plot; loops over every `_DEG_*_vs_*.csv` |
| pseudobulk | `volcano` `heatmap` | heatmap = top DEG × samples; no per-cell violin |
| hetero | `dimplot` `dotplot` `heatmap` | dimplot = subcluster UMAP (falls back to PCA if no UMAP) |

Which plot when: **annotate → dotplot** (+violin to confirm, heatmap for overview); **condition
contrast → volcano** (+heatmap/violin). Object-based plots **auto-detect the grouping column** from the
CSV's labels, so a stale manifest ident can't mis-group them.

### 4d. The dotplot 3-style engine (find_marker.dotplot + hetero.dotplot)

Both `sc_006.plot.find_marker.dotplot.R` and `sc_006.plot.hetero.dotplot.R` are **thin wrappers** that
`source()` a shared `sc_006.plot.dotplot_engine.R` (same folder — it **must be deployed alongside** the
wrappers, or they fail to load). One engine, so edits land in both dotplots at once (no divergent
copies). `--style` picks the look:

| `--style` | palette | group separation | notes |
|---|---|---|---|
| `C` / `gap` **(default)** | grey→blue (Seurat blue) | **facet** — physical whitespace gap + group-name strip per block | the standard annotation dotplot |
| `B` / `grouped` | multi-stop spectral | colored group **label-bars + segment above** + dashed vlines; `--group_colors 'name=#hex,...'` | publication-style |
| `A` / `block` | grey→`#B2182B` (red) | **one** prominent dashed line after the baseline block + dotted between the rest; **no gridlines**, right-side legend; `dot.min=0.1` | baseline/"shared" block set by `--baseline <group>` (see below) |

Shared behaviour: **default genes on X @45°, groups on Y**; `--transpose` flips to genes-on-Y (legend
right). Gene source priority **`--features "A,B,C"` (inline, no groups) > `--marker_list`/`--genes <file>`
(gene,group[,block]) > the marker CSV** (top-N per `cluster`; each gene assigned to its single best cluster
to avoid cross-block dup). Palette overridable with `--cols "low,high"` or a name (`spectral|blue|red|viridis`).
Other knobs: `--topn`, `--dot_min` (**default 0.1 for all styles** — dots below 10% expressed are hidden;
pass `--dot_min 0` to show all), `--drop_zero`/`--keep_zero`, `--keep_neg`, `--dot_scale`, `--width/--height`.
CONFIG block (STYLE, TOPN, DOT_MIN_PCT) at the top of each **wrapper** for quick default edits.

**Curated panels — `--marker_list` / `--baseline` / `--group_order`:**
- **`--marker_list <file>`** (alias `--genes`): a curated `gene,group` table (comma/tab/space-delimited;
  optional 3rd `block` col). When given, the sc_006 CSV is **not read** — but the Seurat **rds is still
  required** (via `--rds`/manifest) because a dotplot needs expression values. Group order in the file is
  the default block/axis order.
- **`--baseline "Shared[,Pan]"`**: names the shared/pan-marker group(s). They are floated to the front,
  get the **single major dashed line** after them (style A), and are **X-column-only — excluded from the Y
  axis and from ident detection**. In styles B/C the baseline is just another bar/facet; the
  major-dashed-line is A-only.
- **`--group_order "c1,c2,c3"`**: Y-axis cell-group order (first = top) **and** X non-baseline block order
  (kept in sync → block-diagonal); baseline stays first. When set (or when `--marker_list` names groups),
  the object is **subset to exactly those Y groups**. Default Y order = the marker_list's group order.
- **`--group_map "Excitatory=6,10,11; Inhibitory=3,7; ..."`** (or a file `lineage,cluster`): use when the
  marker `group` names are **lineages** but you want **clusters on Y grouped by lineage**. Pass `--ident
  <cluster_col>` too (explicit `--ident` now disables auto-detect). Y = the mapped clusters, **faceted by
  lineage** (left strip = lineage name, panel gap = separator), ordered to match the X gene-blocks
  (→ block-diagonal); object is subset to the mapped clusters. In **style C** this becomes a full grid
  (lineage rows × gene-group cols). Lineage-faceting is **non-transpose only** (transpose falls back to
  ordering).

```bash
# a curated panel: Shared baseline block + one dashed line, dotted between subtypes
Rscript scripts/sc_006.plot.find_marker.dotplot.R --marker_list panel.csv --style A --baseline Shared
# force the Y (and X block) order
Rscript scripts/sc_006.plot.find_marker.dotplot.R --marker_list panel.csv --group_order MotorNeuron,Excitatory,Inhibitory
```

```bash
# annotation workhorse (default C, blue+gaps)
Rscript scripts/sc_006.plot.find_marker.dotplot.R
# spectral bars + colored groups
Rscript scripts/sc_006.plot.find_marker.dotplot.R --style B --group_colors c1=#28A828,c2=#1773B3
# red-block style, 8/group, transposed
Rscript scripts/sc_006.plot.find_marker.dotplot.R --style A --topn 8 --transpose
# custom panel from a gene,group[,block] file
Rscript scripts/sc_006.plot.find_marker.dotplot.R --genes panel.csv --style A
```

**Engine trap (verified):** `coord_flip` defeats facet free-scale dropping, so **style C + `--transpose`
is rebuilt natively** from `DotPlot$data` (genes on Y, `facet_grid(rows, scales/space=free_y)`) instead of
`coord_flip` — that path shows all-genes-in-every-panel if done with `facet_wrap+coord_flip`.

### 4c. Traps baked into these scripts (do not re-trip)

- **SCT**: `PrepSCTFindMarkers` is run once at the top of the compute driver; heatmap plots use
  `GetResidual` (SCT) or `ScaleData` (RNA) before `DoHeatmap` — `DoHeatmap` needs a populated
  `scale.data` or it errors.
- **`--outfile` not `--out`** on plot scripts — `--out` collides with `--outdir` under optparse
  prefix-matching (silently saves to the directory → "no file extension" error).
- **Violins avoid Seurat `VlnPlot`** (built from `FetchData` + `geom_violin`) because
  ggplot2 ≥4.0 + Seurat 5.4 crashes `VlnPlot`/`RidgePlot` with `'S4SXP': should not happen`.
  `DotPlot`/`DimPlot`/`DoHeatmap` are unaffected.
- **hetero UMAP** needs ≥2 dims; with 1 significant PC it expands dims and is non-fatal (markers+rds
  still save).
- **`optparse` must be installed** in the R used (`install.packages("optparse")`) — it's the only
  non-Seurat dep and is easy to miss on a fresh R.
- **pseudobulk** needs `--sample_column` and ≥2 samples/side; right side `ALL` is unsupported there.

**Outputs to read for annotation:** `*_find_marker_strict_final.csv` (specific markers; sensitive tier
if strict is empty) → feed the same margin/promiscuity/coexp decision logic; the dotplot PNG for a
visual specificity check.

> **Trap (per-cluster, not per-file):** `myAUC` is NA whenever a cluster's Wilcox-significant genes don't
> clear ROC's own min.pct/logfc bar — this happens **per cluster**, not just when the whole strict file is
> empty. An aggregation one-liner like `awk '$8!="NA"' | sort -k6,6 -k8,8gr | <top-N per cluster>` will
> silently drop any cluster whose myAUC is *entirely* NA — it won't show up as an empty group, it just
> has zero rows. Before trusting a strict-tier top-N aggregate, check the cluster set in the output against
> the object's full cluster list; any cluster missing (or under ~10-15 rows) needs the **sensitive** tier
> for that cluster specifically, not a switch of the whole run. See `ITERATION_LOG.md`.

**Easiest path: `scripts/check_markers_at_cluster_level.sh`** wraps the compute call and the
aggregation below into one script, and adds a **specificity flag**: any positive marker (avg_log2FC>0)
whose background expression (pct.2) exceeds `HIGHBG` (default 0.5) is tagged `"high_bg_check"` in a new
`specificity_flag` column instead of being silently trusted — alongside a `pct_diff` column (pct.1 −
pct.2). Usage:
```bash
bash scripts/check_markers_at_cluster_level.sh <rds> [cluster_col] [thresh] [high_bg_pct2]
# defaults: cluster_col=seurat_clusters  thresh=10  high_bg_pct2=0.5
```
Output is `results/sc006_top50_per_cluster.csv` — hand this straight to the skill for a per-cluster
judgment; when reading it, weight any `high_bg_check` row more skeptically than an `"ok"` row before
treating it as decisive evidence (same principle as the promiscuity traps in `decision_logic.md §3`).

**Recommended replacement aggregation** (the raw awk pipeline the script above wraps, per-cluster
strict→sensitive fallback, tagged, nothing silently dropped — run from the `results/` dir after
`--do_find_marker`):
```bash
THRESH=10   # min ROC-surviving (myAUC!=NA) rows a cluster needs before its strict tier is trusted
awk -F',' -v th="$THRESH" 'NR>1 && $8!="NA"{cnt[$6]++} END{for (c in cnt) if (cnt[c]>=th) print c}' \
  sc006_find_marker_strict_final.csv > .strict_ok_clusters.tmp
{
  head -n1 sc006_find_marker_strict_final.csv | sed 's/$/,"tier"/'
  awk -F',' 'NR==FNR{ok[$1]=1; next} FNR==1{next} $8!="NA" && ($6 in ok){print $0",\"strict\""}' \
    .strict_ok_clusters.tmp sc006_find_marker_strict_final.csv \
    | sort -t',' -k6,6 -k8,8gr | awk -F',' '{c=$6; if(c!=prev){n=0;prev=c}; n++; if(n<=50) print}'
  awk -F',' 'NR==FNR{ok[$1]=1; next} FNR==1{next} $8!="NA" && !($6 in ok){print $0",\"sensitive\""}' \
    .strict_ok_clusters.tmp sc006_find_marker_sensitive_final.csv \
    | sort -t',' -k6,6 -k8,8gr | awk -F',' '{c=$6; if(c!=prev){n=0;prev=c}; n++; if(n<=50) print}'
} > sc006_top50_per_cluster.csv
rm -f .strict_ok_clusters.tmp
```
A trailing `"tier"` column marks which source each cluster's rows came from — read it as a confidence
signal (`strict` = ROC-confirmed, more decisive; `sensitive` = Wilcox-only fallback, still useful but
weight it against QC/coexp more). Sanity-check after running: `tail -n +2 sc006_top50_per_cluster.csv |
awk -F',' '{print $6}' | sort -u | wc -l` must equal the object's real cluster count — if it's lower,
some cluster is still missing (e.g. `THRESH` too high, or a cluster is NA in *both* tiers → genuinely
too weak, flag it, don't force a label).
