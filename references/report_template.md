# Output report template

Default output is a **Chinese** per-cluster annotation report. Lead with the table; expand evidence on
request. Shape:

## 1. One-line status
`N 簇中 M 个干净直接接受, K 个 Unknown(其中 X 个可硬标, Y 个需 resolve)。`

## 2. Per-cluster table
| 簇 | n | 关键证据(%阳性 / top marker) | 判定 | 置信 |
|---|---|---|---|---|
| 0 | 13438 | P2RY12 94, CSF1R 91, C1QB 92 | Microglia(反应性) | 高 |
| … | | | | |

- **判定** = the cell-type label (use the marker-library names; add a state/subtype suffix where
  relevant: `_lowQC`, `_cycling`, `Mural-pericyte`, `Mural-SMC`).
- **置信** = 高 / 中 / 低. 高 = both methods agree or decisive specific markers; 中 = hard-assigned
  vote-split or one decisive panel; 低 = needs more evidence / flagged.
- **关键证据** = the *specific* genes that decided it (not the whole panel). For a vote-split note the
  splitter ("被 X 压低 margin"). For low-QC give nFeature/%ribo.

## 3. Resolved Unknowns (why each is no longer Unknown)
A short line per Unknown: cause (vote-split / tie / low-QC / state) → decisive evidence → label.

## 4. Marker-file edits this run motivated
Bullet list of changes to the marker library + the observation behind each (so the next run is better).
Mirror each into `ITERATION_LOG.md`.

## 5. Commands to run next (if any clusters still open)
Ready-to-run `nohup Rscript … > <named>.log 2>&1` blocks, named per cluster.

## 6. Resolution note (only if clustree was read)
Raise / don't-raise + why (clean split vs messy reshuffle), and any targeted `sub=` instead.

---

## Self-check (mandatory, run before showing the table to the user)

Multi-cluster tables are built by looping over clusters, gathering evidence, then hand-writing labels
into rows — and hand-writing is where **adjacent-cluster row swaps** happen (cl N's evidence ends up
next to cl M's label). A swapped row reads as completely plausible on its own; only a keyed re-query
catches it. Real incident (a 24-cluster CNS vascular dataset): cluster A (fibroblast: LAMA2 pct.1=.995,
DCN=.97, PDGFRA=.87, all myeloid markers negative) and cluster B (myeloid: PTPRC=.73, CSF1R=.51,
CD163=.49, all fibro markers negative) were written to the *opposite* rows -- "Perivascular
macrophage(CD163+)" landed on cluster A, "Vasc_fibro" landed on cluster B -- caught only when the user
asked "why is cluster A CD163+?" and a direct re-query immediately flipped it back.

**Before presenting the table**, for every row, re-query the source CSV **keyed by that row's own
cluster id** and confirm each cited gene's logFC/pct.1/pct.2 sign and magnitude match what's printed.
Do this for the whole table in one pass (e.g. one awk/grep sweep across all cited genes x their
claimed clusters), not spot-checks on a few rows you feel uncertain about -- the swap incident above
was on two of the *most confident* rows in the table, not an ambiguous one.

---

### Confidence calibration (do not over-hedge)
- If the specific markers are decisive, say the label plainly -- no "可能/或许".
- If genuinely unresolved, say what single piece of evidence would settle it.
- Never report a clean call as uncertain, and never report a tie as resolved.
