<p align="right"><b>English</b> | <a href="README.zh-CN.md">简体中文</a></p>

# cell-type-annotation

A [Claude Code Skill](https://docs.claude.com/en/docs/claude-code/skills) that interprets single-cell
RNA-seq cluster-annotation output and helps decide each cluster's cell type — plus the R engine it
drives (scoring, ambiguous-cluster resolution, subtype/DEG/heterogeneity analysis, and matching plots).
It reads your pipeline's numbers and proposes a per-cluster label with the deciding evidence; it does
**not** run the R scripts on its own — it generates the exact command for you to run, reads the result,
and iterates.

> **⚠️ You must review the marker genes yourself.** The built-in marker library
> (`references/cell_markers.cns.txt` + `references/marker_library.md`) is a well-sourced starting
> point with citations and documented traps — it is not a certified ground truth for your species,
> tissue, or disease context. Always check a call against the *specific* genes before trusting it, and
> correct the library when something is wrong. See [Marker review](#marker-review) below.

---

## Contents
- [What this is](#what-this-is)
- [Install](#install)
- [Usage](#usage)
- [Simple usage — one-shot marker check](#simple-usage--one-shot-cluster-level-marker-check)
- [Worked example](#worked-example)
- [Marker review](#marker-review)
- [Self-iteration](#self-iteration)
- [Repo layout](#repo-layout)

---

## What this is

Single-cell cluster annotation usually goes through an engine that scores marker panels per cluster
and calls a label — but a chunk of clusters always come back `Unknown`, or get a label that's wrong
for a subtle reason (a promiscuous marker gene, an ambient-contamination artifact, a cell-cycle state
masquerading as a lineage). This skill is the decision layer on top of that engine: it reads the
consensus table, the threshold log, and any follow-up diagnostic output, and tells you — with the
specific genes as evidence — what each cluster actually is, or exactly which follow-up command to run
next. It ships:
- an **annotation engine** (`scripts/sc_004.*.R`): per-cluster marker-panel scoring (AddModuleScore +
  UCell), margin/threshold logic, ambiguous-cluster deep-dives, and a numeric-summary/DotPlot tool;
- a **subtype/DEG/heterogeneity family** (`scripts/sc_006.*.R`): one-vs-all subtype markers, A-vs-B
  differential expression, sample-level pseudobulk DESeq2, subcluster drill-down, and 11 matching plot
  scripts;
- a **decision-logic doc** (`references/decision_logic.md`) that turns those numbers into a call
  without falling into the known traps;
- a **CNS marker library** (`references/marker_library.md` + `.txt`) as a starting point, with
  citations and documented failure modes;
- an **iteration log** (`ITERATION_LOG.md`) that is the skill's memory — every mistake it made and how
  it was fixed, so the same mistake doesn't happen twice.

It's built for CNS tissue (brain + spinal cord, human + mouse) but the decision logic and script family
are tissue-agnostic — swap in your own marker panels for a different tissue.

---

## Install

**As a Claude Code skill.** Copy (or clone) this repo into Claude Code's skill discovery directory,
either globally or per-project, so `SKILL.md` sits at `<dir>/cell-type-annotation/SKILL.md`:

```bash
# global (all projects)
git clone https://github.com/Yangfan9606/cell-type-annotation-skill ~/.claude/skills/cell-type-annotation

# or project-local
git clone https://github.com/Yangfan9606/cell-type-annotation-skill <your-project>/.claude/skills/cell-type-annotation
```

Restart/open a new Claude Code conversation in that project and describe an annotation task — Claude
will pick the skill up automatically from `SKILL.md`'s description, or you can name it directly.

**Standalone (no Claude Code, just the R scripts).** Everything under `scripts/` is plain R, runnable
with `Rscript` on its own — the skill (`SKILL.md`, `references/`) is the decision layer, not a
dependency of the scripts. Requires R with `Seurat` (≥ 4; some plot scripts assume Seurat 5 layers),
`optparse`, `UCell`, and — only for the features that use them — `MAST`, `DESeq2`, `ggrepel`,
`pheatmap`. See `references/scripts_reference.md` for exact flags per script.

---

## Usage

### Track 1 — first-pass whole-object annotation
```bash
# 1. run the engine against your marker file (start from the built-in seed panel)
Rscript scripts/sc_004.CellType_Annotation.R rds=<your.rds> \
    markers=references/cell_markers.cns.txt --presence_floor 0.4

# 2. tell Claude (or read yourself, per references/decision_logic.md):
#    "annotate the clusters in Harmony_thresholds_log.txt / annotation_consensus.csv"
```
Claude reads `*_annotation_consensus.csv` + `*_thresholds_log.txt`, classifies every `Unknown` as a
vote-split / genuine tie / low-QC / state issue, and for genuine ties generates a ready-to-run
`cluster_resolve` command:
```bash
nohup Rscript scripts/sc_004.02.cluster_resolve.R rds=<your.rds> clu=<clusterCol> \
    clusters=<k1,k2> coexp='<geneA,geneB,...>|<geneC,geneD,...>' \
    > sc_004.02.cluster_resolve.R.c<k1>_<k2>.log 2>&1
```
Read the resulting `*_resolve_summary.tsv` + log, and Claude gives each cluster a verdict with the
deciding genes — following the self-check pass in `references/report_template.md` before presenting
the final table.

### Track 2 — subtype markers / DEG / pseudobulk / heterogeneity

For an object that's already subset/annotated. Full flags in `references/scripts_reference.md §4`.

Example — condition contrast within one cell type:
```bash
nohup Rscript scripts/sc_006.find_marker.DEG.Heter.R -i <your.rds> \
    --subset_column cell_type --subset_value Fibroblast \
    -c condition --contrasts 'Disease,Control' \
    --do_deg --do_pseudobulk --sample_column orig.ident \
    > deg.log 2>&1 &
```
The compute driver prints all 11 matching plot commands at the end of its run — copy the ones you want
(dotplot for annotation markers, volcano for DEG/pseudobulk, heatmap/violin/dimplot as needed).

### Simple usage — one-shot cluster-level marker check

The most common day-to-day pattern: run the compute once, then just hand Claude the `results/` output
and ask for a judgment. `scripts/check_markers_at_cluster_level.sh` wraps the compute call and a
per-cluster top50 aggregation (strict→sensitive fallback, per `scripts_reference.md §4a`) into one
script, and adds a **specificity flag**: any positive marker whose background expression (pct.2)
exceeds a threshold gets tagged `high_bg_check` instead of being silently trusted.

```bash
bash scripts/check_markers_at_cluster_level.sh <your.rds> [cluster_col] [thresh] [high_bg_pct2]
# defaults: cluster_col=seurat_clusters  thresh=10  high_bg_pct2=0.5
# to background it: nohup bash scripts/check_markers_at_cluster_level.sh <your.rds> > check_markers.log 2>&1 &
```

Then just say: *"here's `results/sc006_top50_per_cluster.csv`, tell me what each cluster is."* Claude
reads the `tier`/`pct_diff`/`specificity_flag` columns alongside the marker table and applies the same
promiscuity-trap logic from `references/decision_logic.md` — any `high_bg_check` row gets scrutinized
before being trusted as decisive evidence.

---

## Worked example

`examples/example_cns_dataset.md` walks a full 20-cluster round end-to-end: reading the first-pass
table, classifying all 7 `Unknown` clusters (vote-splits, a genuine tie, a low-quality cluster), the
marker-promiscuity catch that fixed the OPC panel, four `coexp` resolve runs, and the final annotation
— useful as a template for how to work through your own dataset.

---

## Marker review

The shipped marker library (`references/cell_markers.cns.txt` for the machine-readable panels,
`references/marker_library.md` for rationale + citations + traps) was built and validated against real
CNS datasets — every panel has a stated source or an explicit "unverified" flag, and every documented
trap (e.g. `VCAN`/`LHFPL3` producing a false OPC call, `ATP1A2` being miscited as a pericyte marker
when it's actually astrocytic) was a real mistake caught and fixed. That track record is exactly why
you shouldn't take it on faith for your own data:
- **Different species/tissue/disease context can break a panel that worked elsewhere.** Check every
  gene in a panel against your own literature before trusting a call built on it.
- **Verify against your own data before trusting a panel member**, the way the worked examples do:
  check pct.1/pct.2/log2FC for each gene, not just whether the panel "won."
- **When you find a wrong or missing marker, fix it** — edit `cell_markers.cns.txt` /
  `marker_library.md`, and log the observation in `ITERATION_LOG.md` (see below) so the fix is durable
  and explained for the next reader.

Never treat a cell-type call as final because a script produced it — the scripts scores panels
mechanically; deciding whether the panel itself is trustworthy for your data is on you.

---

## Self-iteration

`ITERATION_LOG.md` is not documentation you write once — it's meant to grow every time you (or Claude)
use this skill and find something wrong. Every trap in `decision_logic.md` and every correction in
`marker_library.md` started as a dated entry here: **Problem → Diagnosis → Fix → Where it landed**,
with the actual deciding numbers, not just the conclusion. Two rules keep this working:
1. **Read it on any new problem** — check whether the failure mode has already been diagnosed before
   re-deriving it from scratch.
2. **Append to it whenever a verdict was wrong, a trap was found, or a marker/script was changed** —
   include the evidence, so the next reader (human or Claude) can verify it rather than take it on
   faith.

This is how the skill was actually built: nearly every rule in `references/decision_logic.md` §3 (the
marker-promiscuity traps) and every correction in `references/marker_library.md` exists because a wrong
call was caught, logged, and fixed first. Keep extending it the same way for your own data.

---

## Repo layout

```
SKILL.md                       entry point — when to use + the workflow + hard rules
README.md                      this file
README.zh-CN.md                Chinese translation
CHANGELOG.md                   version history of the skill itself
ITERATION_LOG.md               self-iteration log — problem → diagnosis → fix, grows with use
references/
  marker_library.md            built-in marker library, rationale + citations + traps
  cell_markers.cns.txt         the machine-readable seed marker file (feed to sc_004 as markers=)
  decision_logic.md            how scores become calls; the promiscuity traps; coexp semantics
  scripts_reference.md         every script's args, outputs, and baked-in traps
  report_template.md           shape of the final per-cluster annotation report + self-check rule
examples/
  example_cns_dataset.md       full worked example — 20 clusters, every decision type
scripts/
  sc_004.*.R                   annotation engine + ambiguous-cluster resolve + marker_viz + recode
  sc_005.SubClass_extract.R    clean cell/cluster extraction + plain re-cluster/UMAP
  sc_006.*.R                   subtype marker / DEG / pseudobulk / heterogeneity + 11 plot scripts
  check_markers_at_cluster_level.sh   one-shot compute + aggregation, see Simple usage above
```

## License

MIT — see [LICENSE](LICENSE).
