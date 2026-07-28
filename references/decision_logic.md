# Decision logic — how to read sc_004 outputs and assign cell types

This is the core of the skill. The annotation engine produces *scores*; this file is *how to turn
scores into calls* without being fooled by the known traps.

## 1. The scoring model (what the numbers mean)

`sc_004.CellType_Annotation.R` scores each marker panel per cell with **AddModuleScore (MOD, primary)**
and **UCell (reference)**, aggregates per cluster, **z-scores each panel across clusters**, then per
cluster takes **argmax** (top1) and the **margin = top1_z − top2_z**.

A cluster is called **Unknown** when `margin < threshold`. The threshold is adaptive:
`max(MAD-based [median − 3·MAD], permutation-95%, hard_floor 0.1)` — usually the permutation 95%
(~1.1–1.3). A **negative MAD threshold** means the margin distribution is symmetric / no clearly-bad
cluster (good — high separability), so the permutation value governs.

Read these columns from `*_annotation_consensus.csv`:
`MOD_call, MOD_top1_id, MOD_top1_z, MOD_top2_id, MOD_top2_z, MOD_margin, MOD_reject_reason`,
the UCell mirror, `mod_ucell_agree`, and final `cell_type`.

## 2. Why clusters go Unknown — two very different causes

| Cause | Signature | Action |
|---|---|---|
| **Near-neighbour vote-split** | top1 is the right lineage with **strong z**, top2 is a *related* type (Microglia vs Perivasc_Macro; Vasc_endo vs Mural; Mural vs Astrocyte). UCell often passes it. | **Hard-assign top1** if biology + the marker table support it. The margin is an artifact of overlapping panels, not real ambiguity. |
| **Genuine tie** | top1 ≈ top2, **both moderate** z (e.g. OPC 2.67 vs Inhibitory 2.40). | **Run `cluster_resolve` with a `coexp='A\|B'`** to adjudicate by the *specific* markers. |

Vote-splitting is worsened by **overlapping / hierarchical panels**: a pan-Neuron panel
(SYT1/SNAP25/RBFOX3) competes with Excitatory/Inhibitory/MotorNeuron and compresses every neuronal
margin → false Unknowns. If a real neuron subtype is lost, check the threshold log for
top1=Neuron_pan / top2=<subtype>; consider dropping the pan panel.

## 3. The marker-promiscuity traps (memorize these)

These are real cases where a panel "won" or "tied" on a gene that is **not specific**:

- **VCAN and LHFPL3 → false OPC on neurons.** Measured: LHFPL3 82% positive in a GABAergic cluster,
  VCAN 56–82% across many clusters, while the *specific* OPC genes PDGFRA/CSPG4 were ~0. **OPC is real
  only if PDGFRA AND CSPG4 fire.** Fix applied: OPC panel = `PDGFRA,CSPG4,VCAN` (LHFPL3 removed).
- **C1QA/C1QB → false myeloid on ambient.** C1Q is shared by microglia and perivascular macrophages
  and leaks as ambient. Pair with `--presence_floor 0.4` and verify CSF1R/TYROBP/AIF1 before calling
  microglia on C1Q alone.
- **CSPG4 → pericyte, not just OPC.** CSPG4 (NG2) is ~37–50% positive in pericyte/SMC clusters. Don't
  read CSPG4 as OPC unless PDGFRA co-fires and RGS5/PDGFRB are absent.
- **PDGFRA → fibroblast too.** Fibroblasts/VLMC express PDGFRA (~31%). PDGFRA high + DCN/COL1A1 high =
  fibroblast, not OPC.
- **SLC1A2 → also neuronal.** High SLC1A2 without AQP4/GFAP is not astrocyte.
- **Do NOT put SOX10 / OLIG1 / OLIG2 / CNP in the OPC panel** — pan-oligo-lineage / mature genes that
  steal votes from Oligodendrocyte.

## 4. Quality before identity

A cluster where **every panel's absolute z is low** is almost always **low quality, not a rare type**.
Check QC: low `nFeature` (e.g. ~800), high `%ribo` (e.g. ~27%), elevated `%MT`. Example: a neuroblast
cluster tied Excitatory(0.82) vs Inhibitory(0.58) — the universally-low z *was* the low-quality
fingerprint. Label such clusters `<lineage>_lowQC` and drop them from downstream quantification.

Hard QC flags (script auto-prints): nFeature < 500 → low quality; %MT > 15 → apoptotic/low quality;
%HB > 5 → red-blood-cell contamination.

### 4a. Ambient contamination vs genuine reactive/activated state — use `pair=` against the closest clean neighbour

A QC-flagged cluster whose own one-vs-rest top markers are generic housekeeping/ribosomal genes can look
like it might instead be a real "reactive/activated" biological state (especially if it also carries some
IFN/MHC-I signal, which reads as plausible activation). Don't guess — run `pair='<flagged>:<neighbour1>,
<flagged>:<neighbour2>'` against its closest same-lineage clean cluster(s) and check two things:
1. **Is the flagged cluster's IFN/MHC-I level actually elevated vs its neighbours, or just similar to
   them?** If similar, it's not a distinct activated state — it's the object's general baseline.
2. **What is "higher in the flagged cluster" made of, and at what pct.1?** If it's the object's
   *numerically-dominant* population's secreted/structural markers (e.g. collagens/OGN/THBS in a
   fibroblast-heavy object) at **minority pct.1 (10-40%, not majority)**, that's ambient soup riding on a
   low-complexity capture — low-UMI cells dilute their own transcriptome less, so background ambient
   reads become relatively more visible. A **real** reactive/activated state instead shows lineage-
   appropriate activation genes (e.g. VCAM1/ICAM1/SELE/CXCL9-11 for EC) at meaningful, often majority pct.

Worked example (a 24-cluster CNS vascular dataset): a QC-flagged EC cluster's own markers were
housekeeping/ribosomal; IFI27/HLA-B/HLA-E/IFITM2 were within the *same range* as its clean EC neighbours
(not elevated — ruled out "more activated"); `pair=` against those neighbours showed the flagged cluster
"higher in" PDGFRA/OGN/COL6A3/COL12A1/ABCA8/ABCA10/THBS1/THBS2/F3 — all fibroblast/ECM genes, i.e. the
dominant population in that object — but at **pct.1 = 10-41%**, and these same genes were *negative* in
the flagged cluster's own one-vs-rest table (the dedicated fibro clusters have even higher pct.2). Classic
ambient-fibro-soup signature on a low-quality EC capture, not a doublet (pct well under majority, coexp
double% was already under `dblwarn`) and not a reactive EC state (no VCAM1/ICAM1/SELE/CXCL9-11 anywhere in
either the one-vs-rest or pair tables). Final call: `EC_lowQC`, not `EC_reactive`.

## 5. States vs lineages

`Proliferating` / `Cycling` (MKI67/TOP2A/CENPF/HELLS/DIAPH3) is a **state axis**, not a cell type.
A cycling cluster scores `Proliferating` overwhelmingly because cell-cycle genes mask the lineage.
To name it: strip cell-cycle genes and read the **non-cycle** top markers + a `coexp` of candidate
lineages. Example: a cycling cluster had microglia markers ≈0 (NOT proliferating microglia) and
SOX2/PTPRZ1/PANTR1/NFIB high → **cycling neural progenitor**. Name it `<lineage>_cycling`.

## 6. coexp / doublet adjudication (post-fix semantics)

`coexp='panelA|panelB'` reports, per cluster, `panel1+ %`, `panel2+ %`, and `double %`. A cell is
"panel-positive" if it detects (>0) at least half the panel's genes. New interpretation rules
(script now annotates each line and only warns when justified):

- **double% high requires BOTH panels high → real doublet suspicion** (confirm with scDblFinder) —
  **but only when the two panels are different lineages.** If both panels sit inside the **same broad
  lineage** (e.g. SMC vs Pericyte — both mural), high double-positivity is expected biology for a
  continuous differentiation zone, not a doublet. Worked example (a 19-cluster CNS vascular dataset):
  panel1=ACTA2/MYH11/TAGLN/CNN1/MYL9 (SMC) 47%+, panel2=RGS5/PDGFRB/NOTCH3/ABCC9 (Pericyte) 79%+,
  double 39% — script's own rule-of-thumb flagged "查 doublet", but this is the well-documented
  arteriolar SMC↔pericyte transition continuum (Vanlandewijck 2018-style mural spectrum), confirmed by
  CARMN (mural-lineage lncRNA) + NOTCH3/ABCC9/PDGFRB + ACTA2/LMOD1/DMD/SYNPO2 all genuinely co-expressed
  in the same cells. **Call it `Mural_transitional` (or `SMC-Pericyte transition`), not a doublet.**
  Cross-lineage double-high (e.g. EC panel + immune panel both high) is the real doublet-suspicion case.
- **One panel high, the other low, double low → that lineage** (the "→ 偏 panel1/2" tag). The shared
  genes inflate nothing; this is a clean call.
- **double% can look high (e.g. 28%) yet be shared markers, not doublets** — e.g. fibroblast vs
  perivascular sharing stromal genes. Only the both-high case (and only cross-lineage) is a doublet.

The `dblwarn=` threshold (default 20%) gates the doublet message.

## 6b. Relative-negative logFC ≠ absent — coexp % positive is the tiebreaker for "markerless" clusters

A cluster can show **negative logFC for every canonical marker of every lineage** in a one-vs-rest
FindMarkers/FindAllMarkers table and still have a real, callable identity. Negative logFC there means
"lower **than the other clusters**," not "absent in absolute terms" — a cluster sitting in a crowded
region of highly-expressing neighbours (e.g. six other fibroblast clusters) will read as fibro-depleted
by relative ranking even while it is, in absolute terms, still the most fibroblast-like population in
the object. **`coexp='panelA|panelB'`'s % positive is a raw-presence calculation** (per docstring: "at
least half the panel's genes detected"), not a cross-cluster comparison — it is the correct tool to
resolve a cluster that looks "markerless" by the one-vs-rest logFC alone.

Worked example (the same 19-cluster CNS vascular dataset): two clusters showed negative logFC for
COL1A1/COL1A2/DCN/PDGFRA/RGS5/ACTA2/MYH11/PECAM1/CLDN5/NOTCH3/ABCC9 in the one-vs-rest table — read
(wrongly, first pass) as "no canonical identity, flag as unresolved." A `coexp='COL1A1,COL1A2,COL3A1,
DCN,PDGFRA|RGS5,PDGFRB,NOTCH3,ABCC9'` run showed **panel1(fibro) 34-35% positive vs panel2(pericyte)
7-11%**, and the negative-marker top15 in the same log decisively excluded EC (VWF/PTPRB/ESAM/TIE1/
S1PR1/FLT1/CD93/EGFL7/NOS3 all strongly negative) — i.e. both clusters **are** fibroblast-lineage, just
an atypical/low-collagen substate, not unresolved. **Before calling a cluster "unresolved" from
one-vs-rest logFC alone, run the coexp check** — it is the deciding evidence, not a formality.

### 6c. `pair=A:B` "relatively higher" genes ≠ that cluster's dominant identity — coexp is still the arbiter

`pair=A:B` answers "what's different between these two specific clusters," not "what is cluster A's
main identity." A gene set that is higher in A than in B can still be a small minority signal within A
overall — `pair=` alone cannot tell you that. **Always follow up a `pair=` finding with `coexp=` on the
actual candidate gene panels (real genes, not cluster ids — see ITERATION_LOG.md for the
cluster-ids-passed-as-coexp-panels trap) before promoting the pair='s "relatively elevated" genes to a
cluster's identity label.**

Worked example (the same 19-cluster CNS vascular dataset): `pair=A:B` on two fibroblast-lineage clusters
showed axon-guidance/neural-adhesion genes (GPM6A/CNTN3/NCAM1/SEMA3A/ROBO2/CACNA1A) elevated in cluster A
relative to cluster B — read (wrongly, first pass) as "cluster A is a nerve-associated/perineurial
fibroblast substate." The follow-up `coexp=` on the actual neural-adhesion panel vs cluster B's plain ECM
panel (COL11A1/THBS2/FBLN5/LTBP2) showed the neural-adhesion panel was only **4% positive** in cluster A
versus **65%** for the ECM panel — i.e. the pair= "winner" was a minority signal all along; cluster A's
actual dominant identity was the same ECM lineage as cluster B, just a distinct (quiescent) activation
state within it (independently confirmed by a formal MOD/UCell consensus run — the cleanest call in that
entire run). The lesson generalizes: a `pair=` result tells you a *direction*, never a *magnitude relative
to the cluster's own baseline* — that second number only comes from `coexp=` (or the plain one-vs-rest
pct.1 of the cluster itself).

## 7. Reading the machine-readable summaries (preferred)

- `*_resolve_summary.tsv` (from cluster_resolve): one row per cluster —
  `cluster, n, nFeature, nCount, pct_mt, pct_ribo, pct_hb, qc_flag, top5_pos, coexp_p1, coexp_p2, coexp_dbl`.
  Read this first to sweep all clusters at once.
- `*_summary.txt` (from marker_viz): QC medians + gene×cluster **pct.exp** and **avg.exp.scaled** +
  per-gene means. Use to confirm exactly which panel member fires where (catches promiscuity).

## 8. Resolution decision (clustree)

The clustree log has [1] per-resolution cluster sizes and [2] adjacent-resolution split contingency
tables. **Clean splits** (a parent cluster cleanly divides into children, near-diagonal) = real
substructure → raising resolution is safe. **Messy off-diagonal reshuffling** (cells scatter across
many children) = over-splitting a continuous manifold → do NOT raise globally; instead keep the
coarse resolution and pull substructure with targeted `sub=<cluster>`. Never raise resolution just to
chase one ambiguous cluster — it fragments the neuron manifold and worsens vote-splitting.

## 9. Decision flow (one cluster)
```
clean call + both methods agree                      → accept label
Unknown, top1 strong + top2 related + UCell agrees    → hard-assign top1 (note splitter)
Unknown, top1 ≈ top2 both moderate                    → coexp resolve, decide by specific genes
all panels low z                                      → check QC → <lineage>_lowQC or drop
Proliferating wins                                    → strip cycle genes, coexp lineages → <lineage>_cycling
coexp both panels high + double high                  → doublet (scDblFinder)
```
