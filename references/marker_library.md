# Built-in marker library — CNS (human brain + spinal cord) + mouse homolog

**Version date: 2026-07-17.** This is the skill's seed marker library. The machine-readable file is
[`cell_markers.cns.txt`](cell_markers.cns.txt) — feed it to the engine as `markers=`. This `.md` is
the human-readable rationale + caveats. **Iterate both together** and bump the date.

> **Review before you trust it.** This library was built and validated against specific real datasets
> (see the worked examples below) — it is a well-sourced starting point, not a certified ground truth
> for your tissue. Before using any panel, check its citations, and re-verify against your own data
> (pct.1/pct.2, log2FC) the same way the worked examples below do.

Format of the `.txt` (same as the pipeline expects): `celltype<TAB>gene,gene,gene`, `#` comments,
mixed tab/comma tolerated. One panel per line.

## Species switch
Human symbols below. For **mouse**: title-case the genes (`Syt1, Slc17a7, Chat, Mki67…`) and switch
the QC regexes in the scripts to `^mt-`, `^Rp[ls]`, `^Hb`. Homologs are 1:1 for essentially all genes
here except where noted.

## Panels (with rationale & traps)

### Neurons
- **Neuron_pan** `SYT1,SNAP25,RBFOX3` — pan-neuronal. ⚠ Competes with the subtype panels under
  argmax+margin and can compress a real subtype to Unknown. Include only when you want a neuron/
  non-neuron split; drop it if a subtype (esp. MotorNeuron) is being lost to vote-splitting.
- **Excitatory** `SLC17A6,SLC17A7,SATB2` — SLC17A7 (VGLUT1) cortical, SLC17A6 (VGLUT2) subcortical/
  spinal; SATB2 callosal/excitatory.
- **Inhibitory** `GAD1,GAD2,SLC32A1,SLC6A5` — GAD1/2 + VGAT(SLC32A1) GABAergic; SLC6A5 (GlyT2)
  glycinergic (spinal/brainstem). Immature interneurons add ARX/DLX6-AS1/NXPH1/ERBB4.
- **MotorNeuron** `CHAT,ISL1,MNX1,SLC5A7,SLC18A3` — spinal/brainstem cholinergic MN. ISL1 alone is
  not specific. If CHAT/SLC5A7/SLC18A3 are all 0 across the dataset, MNs are absent in the sample.

### Glia
- **Astrocyte** `SLC1A2,AQP4,GFAP,ATP1A2` — AQP4/GFAP most specific; SLC1A2 also neuronal. Of the
  common "pan-astrocyte" panel genes, only genes that behave *uniformly* across every astrocyte
  cluster and drop sharply outside the astrocyte compartment should be trusted as pan markers — some
  textbook "pan" genes (e.g. SOX9, ALDH1L1) turn out on inspection to be shared with the
  ependymal/radial-glia lineage or to carry their own subcluster gradient, so always spot-check a
  "presumed uniform" gene's actual dotplot/AverageExpression rather than assuming absence from a
  DE-significance table means it's flat.
- **Astrocyte subtypes/states — illustrative case study: origin × activation, not N independent
  types.** A real 6-cluster one-vs-all subtype run on a CNS astrocyte compartment separated into two
  **lineage/origin** axes that each carried their own **homeostatic vs activated** member — confirmed
  by lineage-identity genes (independent of activation state) splitting cleanly the same way the
  clusters merged at lower clustering resolution. This is a good template for reasoning about any
  "clusters that merge at low resolution but stay separate at high resolution" question:
  1. Pull genes that should carry **lineage identity** (transporters, ECM, receptor tyrosine kinases —
     not activation/stress genes) and check whether they split the clusters into two internally-
     consistent groups, cleanly, with no exceptions.
  2. Separately check genes that should carry **activation state** (e.g. CHI3L1, a reliable
     pan-activation gene in this system) — confirm they split *orthogonally* to the lineage genes, not
     redundantly with them.
  3. **Watch for a promiscuous gene inflating a `coexp` panel.** One candidate discriminator gene had
     pct.2 (rest-of-object baseline) near 0.85 — i.e. near-ubiquitous across the whole compartment —
     and including it in a `coexp` panel inflated the double-positive% almost identically in both a
     confirmed-clean homeostatic cluster and the cluster actually under test. **A `coexp` result that
     looks the same in an already-confirmed-clean cluster and the cluster under test is itself evidence
     *against* a real doublet/discriminator signal** — always check a candidate coexp gene's own pct.2
     before trusting it, and compare against a same-lineage-confirmed clean cluster whenever one is
     available.
  4. **Don't let a superficially "reactive-sounding" gene set overturn a settled activation-axis call.**
     A separate literature-style panel surfaced additional genes (heat-shock/ECM-family) that split the
     same way as the lineage axis, tempting a re-read as "lineage A is the reactive one." Re-checking
     against the genes that most directly define the activation axis (a canonical homeostatic connexin
     that should go *down* under activation) showed the opposite of what that re-read would predict —
     the correct model held: those extra genes are a lineage-level *baseline* difference (independent of
     activation state), not evidence that one lineage's quiescent member is secretly reactive. **General
     rule: when a new gene panel seems to contradict a settled call, re-run its actual genes against your
     own data rather than pattern-matching on category names from someone else's atlas** — and always
     re-check the genes that most directly contradict the new reading before accepting or rejecting it.
  5. **A PDF "read" of a dotplot via text extraction has no size/color/position data.** A ggplot PDF
     opened through a text-extraction tool returns axis labels, gene names, and legend values as plain
     text — vector dot geometry (size, color, position) is simply not present in that output. Describing
     specific per-dot visual properties from such a result is fabrication, not observation; before
     describing any visual property from a Read'd PDF, confirm the tool result actually contains image
     data (a genuine `[Image: …]`-style result), not just extracted text tokens. When this happens,
     discard the visual claim explicitly and re-derive the same conclusion from the underlying numeric
     table instead — real image reads (e.g. a PNG) come back with actual pixel data and can be trusted.
- **Oligodendrocyte** `PLP1,MOBP,MBP,MOG,MAG,CNP,CLDN11,OPALIN` — mature/myelin. MOBP/MOG/MAG/OPALIN
  near-100% in true oligo; MBP leaks (ambient myelin) so don't read MBP alone.
- **OPC** `PDGFRA,CSPG4,VCAN` — ⚠ **real only if PDGFRA AND CSPG4 fire.** VCAN is a promiscuous filler
  (lights neurons/fibroblast); LHFPL3 was removed (82% in a GABA cluster). CSPG4 also marks pericytes;
  PDGFRA also marks fibroblasts. Do NOT add SOX10/OLIG1/2/CNP (steal votes from Oligodendrocyte).

### Immune
- **Microglia** `P2RY12,CX3CR1,C1QB,C1QA,CSF1R,APBB1IP` — P2RY12/CX3CR1/CSF1R most specific. ⚠ C1QA/B
  shared with macrophages + ambient; verify CSF1R/TYROBP/AIF1. Reactive/DAM states up CD83/CCL3/TNF/
  SPP1/CLEC7A; homeostatic keeps P2RY12/CX3CR1/GPR34.
- **Perivasc_Macro** `MRC1,CD163,LYVE1,F13A1` — CD163/MRC1 distinguish from microglia.
- **Tcells** `CD3E,CD3D,CD3G`.

### Vascular / stromal
- **Vasc_endo** `CLDN5,PECAM1,FLT1` — endothelial.
- **Mural_pericyte** `RGS5,PDGFRB,NOTCH3,KCNJ8,HIGD1B,NDUFA4L2` — pericyte; KCNJ8/HIGD1B/NDUFA4L2 are
  the pericyte-specific add-ons vs SMC. NDUFA4L2 is sourced from Mesa-Ciller et al. 2023, *J Cereb
  Blood Flow Metab*, "Unique expression of the atypical mitochondrial subunit NDUFA4L2 in cerebral
  pericytes fine tunes HIF activity in response to hypoxia" — reported as specifically distinguishing
  pericytes from arterial SMC along the arteriole-artery axis in both mouse and human brain (building
  on Vanlandewijck et al. 2018), a clean fit alongside HIGD1B (also hypoxia/HIF-linked).
  ⚠ **A commonly-used "pericyte" gene, ATP1A2, has no pericyte citation anywhere in the literature —
  it is instead a well-established astrocyte marker** (this library's own `Astrocyte` panel uses it).
  If you find ATP1A2 in a custom pericyte panel you're auditing, treat it as a likely copy-paste
  mix-up and check whether removing it clears an unresolved pericyte/other-mural-type tie.
- **Mural_SMC** `ACTA2,MYH11,TAGLN,MYL9,MYLK` — contractile vascular smooth muscle. Splits cleanly
  from pericyte by the contractile genes (ACTA2/MYH11/TAGLN).
- **Vasc_fibro** `DCN,COL1A1,COL1A2,LUM` — fibroblast/VLMC; collagen + decorin. Shares PDGFRA with OPC.

### Literature-audit methodology (worked example)
When auditing a custom marker panel, the useful move is a real, per-gene literature search (not
recalled/guessed) against every gene in every panel — grading each gene as confirmed / plausible /
unverified, and checking the panel's own dataset for whether each gene actually clears a
significance/pct.1 threshold anywhere. A worked pass over a custom vascular/perivascular marker file
(EC arterial/capillary/venous, Pericyte, SMC, several fibroblast activation-state panels, meningeal
fibroblast) found:
- **Solid, single-paper-sourced**: EC_arterial (SEMA3G/EFNB2/BMX/GJA5, all 4 from Vanlandewijck et al.
  2018 *Nature*, "A molecular atlas of cell types and zonation in the brain vasculature"); SMC
  (ACTA2/MYH11/MYL9/TAGLN, vascular-biology consensus, present in that paper too but not exclusive to
  it); EC_capillary (MFSD2A/CA4 confirmed from the same paper, SLC2A1 near-certain/textbook BBB gene,
  one panel member unconfirmed as capillary-specific there).
- **Partial**: EC_venous (NR2F2 confirmed same paper, ACKR1 canonical venous-EC marker broadly but not
  paper-specific, remaining members unverified). Meningeal fibroblast (SLC4A4 and SLC47A1 both
  confirmed — Pietilä et al. 2023 *Neuron*, "Molecular anatomy of adult mouse leptomeninges", reports
  SLC4A4 directly as a human leptomeningeal-fibroblast marker and SLC47A1 for mouse dural-border-cell/
  leptomeningeal fibroblast; remaining panel members unverified).
- **A likely copy-paste error, found and fixed**: see the ATP1A2/NDUFA4L2 pericyte note above.
- **No external literature source found at all** for a 5-panel, ~20-gene set of fibroblast
  activation-state labels (inflammatory/myofibroblast/antigen-presenting-like substates). The
  topically-closest published atlas of the relevant tissue's own key markers didn't match any of these
  genes — the most defensible conclusion for a scheme like this, when it doesn't match the closest
  available published atlas, is that it's a **project's own empirically-derived DE-based
  classification**, not a citable external panel. Label it as such (data-derived, not literature-backed)
  in any downstream write-up rather than presenting it as atlas-confirmed.
This is a reusable methodology, not a one-off result: grade every custom panel gene as
confirmed/plausible/unverified with a real citation search, and flag any panel with zero external
sources as project-internal/data-derived so it's described accurately later.

### Meninges (dural/arachnoid/pia fibroblast layers)
- **Meningeal_fibro_pia** `SLC47A1,PTGDS,S100A6,SLC4A4` — pia-mater-layer fibroblast. SLC4A4 is
  sourced from Pietilä et al. 2023 *Neuron* "Molecular anatomy of adult mouse leptomeninges" (Betsholtz
  lab) — reports SLC4A4 directly as a marker of *human* leptomeningeal fibroblasts ("high expression of
  Lama2, Slc4a4 and Slc7a2"); this is currently the best-sourced gene in the whole panel
  (human-confirmed, not just mouse-orthology). (DeSisto et al. 2020 mouse meninges atlas nomenclature;
  human orthologs used here.) **SLC47A1 is the load-bearing gene**: in a worked-example 19-cluster CNS
  vascular subset it was the single strongest marker of exactly one cluster (avg_log2FC +7.5, pct.1
  .81) and **negative in every other cluster** (avg_log2FC -1.7 to -7.5, pct.1 .002-.08) — about as
  clean a single-gene layer call as this library has seen. PTGDS/S100A6 are less exclusive (moderate
  pct.1 across most fibroblast clusters as a pan-stromal baseline) but rise well above that baseline
  (pct.1 .7-.9 vs pct.2 .35-.57) in the true pia cluster — use them to confirm, SLC47A1 to call.
  ⚠ Don't confuse with `Vasc_fibro` (DCN/COL1A1/LUM) — DCN/MGP/CST3/CLU are broadly expressed across
  *all* fibroblast subtypes in CNS vascular subsets (near-baseline pct.1 .7-.9 everywhere), so a cluster
  "winning" on those alone is not decisive for meninges specifically — always check SLC47A1.
- **A cluster with elevated PTGDS/CLU/CST3/MGP/DCN/S100A6 together with elevated mitochondrial/ribosomal
  genes in its own top markers is not automatically "low quality."** These are also the pan-meninges
  panel — before calling `_lowQC`, check whether the "extra" genes are QC genes only (→ lowQC) or
  include this meningeal panel elevated *above* the fibroblast-wide baseline (→ genuine pia-fibroblast,
  possibly with somewhat higher ambient/mito background, which does not by itself invalidate the call).
  Worked example: a cluster in the dataset above was first called `Fibro_lowQC` from MT-/RPL-gene
  loading alone; re-examination showed PTGDS/CLU/CST3/MGP/S100A6 all decisively above the
  fibroblast-wide baseline — reclassified to `Meningeal_fibro_pia`.
- **Antigen-presentation ("APC-like") fibroblast is NOT decidable from a one-vs-rest marker CSV alone.**
  If HLA class II genes (HLA-DRA/DRB1/DPA1/DQA1/DQB1, CD74, CIITA) don't appear anywhere in the
  significance-filtered marker table, that means they never cleared the per-cluster min.pct/logFC bar
  in *any* cluster — it does not prove they're absent, only that no cluster is differentially high for
  them relative to the rest of that same subset. Settle this with a targeted `marker_viz` `genes=` run
  (raw mean/%positive per cluster, from the `[genes]` block of `*_summary.txt` — NOT `markers=`, which
  expects a marker *file*), not by scanning the one-vs-rest output.
  **Worked example (same 19-cluster CNS vascular dataset):** a targeted `genes=` check across all 19
  clusters showed one cluster leading on every executor gene — HLA-DRA mean 0.59/%pos 32.0 (next-highest
  0.26/19.7), **CD74 mean 0.94/%pos 46.1 (next-highest 0.56/37.7)**, HLA-DPA1 mean 0.33/%pos 18.5,
  HLA-DQA1 mean 0.08/%pos 5.25 — all dataset-wide maxima, well above the EC clusters that carry the
  next-highest values. This *is* real, decisive evidence for an APC-like/immunomodulatory fibroblast
  subtype (this cluster was independently called `Meningeal_fibro_pia` from the SLC47A1-family panel —
  the two identities are not mutually exclusive; MHC-II+ fibroblast/CAF populations are a described
  phenomenon in fibrosis and tumor-stroma literature). ⚠ **Caveat**: CIITA — the actual master TF that
  switches on the whole MHC-II program — was *not* elevated in this cluster (was highest in a different
  cluster entirely). Don't let the CIITA mismatch overturn the call: HLA-DRA/CD74/DPA1/DQA1 are the
  direct downstream readout and are unambiguous here; CIITA transcript counts in snRNA-seq are
  sparse/dropout-prone and don't always track the downstream genes 1:1. Report both numbers rather than
  silently picking the one that agrees.

### Neural-crest pigment lineage
- **Melanocyte** `TYR,TYRP1,DCT,PMEL,MLANA,OCA2,PAX3` — leptomeningeal/perivascular melanocyte, a
  **real, normally-occurring CNS population** (not ambient/artifact/doublet). Histology/anatomy
  literature (not a single-cell atlas — see caveat below) places the highest normal concentration at
  the **ventrolateral medulla oblongata and upper/ventral cervical spinal cord**, with documented
  extension into perivascular (Virchow-Robin) spaces — i.e. exactly the compartment a spinal-cord
  *vascular* subset would capture. Neural-crest lineage, developmentally distinct from both fibroblast
  (mesenchyme) and Schwann cell (also neural-crest, but a separate terminal fate) — do not fold into
  either. A narrower 3-gene panel (`TYR,KIT,MITF`) is also usable for a MOD/UCell consensus run;
  ⚠ **worked-example check**: TYR and MITF were both real and decisive in a real 19-cluster dataset
  (TYR avg_log2FC +9.9, pct.1 .72, pct.2 .004 — dataset-wide cleanest hit; MITF avg_log2FC +4.4, pct.1
  .78, higher baseline than TYR but pct.diff still clearly favoring the melanocyte cluster), but **KIT
  never appeared in any of the marker output tables for any cluster** — never cleared even the loosest
  significance/pct filter. KIT (c-Kit) is a well-established general melanocyte marker in the broader
  literature, but is more tied to melanocyte precursors/stem cells than fully pigmented, differentiated
  melanocytes — plausibly just not well-captured in snRNA-seq (low counts/dropout), or genuinely low in
  a mature population. **Lesson: before trusting any marker gene to contribute discriminating power in
  a specific run, check its raw pct.exp directly** — a textbook marker absent from every DE table in
  your own data is a real finding to flag, not something to override with the literature.
  ⚠ **Literature caveat**: checked two relevant single-cell/single-nucleus atlases (adult human spinal
  cord, eLife 2023; human leptomeninges, Nat Commun 2023) — neither reports a discrete melanocyte
  cluster (the spinal cord atlas didn't subcluster its meningeal/vascular population at all; the
  leptomeninges atlas found only endothelial/mural/fibroblast/immune as its 4 major types). This is not
  evidence against the population's existence (both are old/postmortem-tissue snRNA-seq, where sparse,
  pigment-laden melanocytes are known to be under-captured) — the identity call rests on marker
  specificity + histology/anatomy literature, not on a matching transcriptomic reference. Say so plainly
  if writing this up rather than implying an atlas-confirmed cluster.

### Peripheral (nerve-root contamination in spinal cord samples)
- **Schwann_cell** `MPZ,PRX,EGR2` — PNS myelinating glia (common contaminant in spinal-cord dissections
  that include nerve roots). MPZ/PRX are PNS-exclusive structural myelin genes (CNS oligodendrocytes
  don't express either); EGR2 (Krox20) is the myelinating-SC master TF, distinct from the oligo program.
  ⚠ From a 10-marker Schwann-ID literature panel (S100/p75NTR/Sox10/Sox2/GAP43/NCAM/Krox20/Oct6/MBP/MPZ),
  **do not add**: MBP/SOX10 (shared with Oligodendrocyte — SOX10 is the pan-oligodendroglial-lineage TF,
  MBP is CNS+PNS myelin both), S100(B) (this library's Astrocyte marker), SOX2 (shared with
  NPC_RadialGlia), GAP43/NCAM1 (broad axon-growth/neural-adhesion genes — GAP43 was directly observed
  firing in a Proliferating-OPC cluster in a real worked dataset, not Schwann-specific). NGFR(p75NTR)/
  POU3F1(Oct6) mark immature/promyelinating SC substates — useful for later SC-heterogeneity drill-down,
  but leave out of the core identity panel (only a subset of SC states express them, would dilute margin).

### Ependyma & states
- **Ependymal** `FOXJ1,PIFO,TMEM212,CCDC153` — ciliated; central canal / ventricular.
- **Proliferating** `MKI67,TOP2A,CENPF,HELLS,DIAPH3` — ⚠ STATE not lineage; strip + read underlying.
- **Neuroblast** `DCX,STMN2,SOX11,SOX4,CD24,NNAT` — immature/migrating neuron (often low-QC if sparse).
- **NPC_RadialGlia** `SOX2,PTPRZ1,FABP7,VIM,HES1` — neural progenitor / radial glia; pairs with
  Proliferating in a cycling progenitor pool.

## How to extend
1. Add/adjust a line in `cell_markers.cns.txt`; document the rationale + any trap here; bump the date.
2. Prefer **specific** genes; if you must include a filler (e.g. VCAN), note it as a filler so the
   reader knows a lone hit on it is not evidence.
3. Record the change in `../ITERATION_LOG.md` with the observation that motivated it.
4. For a non-CNS tissue, copy this file to `cell_markers.<tissue>.txt` and swap panels; keep the
   states (Proliferating/Neuroblast removed where irrelevant) and the vascular/immune backbone.
