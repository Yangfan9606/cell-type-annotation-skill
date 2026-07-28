# Worked example — illustrative CNS dataset, 20 clusters

> This is an illustrative worked example built from a real CNS single-cell annotation round, with
> project- and disease-identifying details generalized out. The decision types, marker traps, and
> numbers are real — read it as a template for how to work through your own dataset, not as a
> published result.

Input tissue: adult CNS (spinal cord). The data turned out heavily neurogenic (immature/cycling
neural populations + massive microgliosis) — a useful stress-test for the decision logic because so
many clusters landed near the margin threshold. Shows every decision type. Inputs:
`*_annotation_consensus.csv`, `*_thresholds_log.txt`, `full_summary.txt`, and four `cluster_resolve`
logs.

## First pass
Marker file = the CNS seed panel (`cell_markers.cns.txt`) plus a few dataset-specific additions;
`--presence_floor 0.40`; permutation threshold **1.148** (MAD term negative → symmetric margins,
good). 13/20 clean, **7 Unknown**.

## Reading the 7 Unknowns (consensus + thresholds)
| cluster | top1/top2 (z) | margin | cause | verdict |
|---|---|---|---|---|
| 0 | Microglia 1.37 / PvMacro 0.45 | 0.92 | vote-split (UCell passes 1.26) | **Microglia (reactive)** hard-assign |
| 8 | Mural 3.12 / Astro 2.19 | 0.93 | vote-split, strong z | **Mural** hard-assign |
| 19 | Vasc_endo 2.74 / Mural 1.91 | 0.83 | vote-split, strong z | **Vasc_endo** hard-assign |
| 13 | Excit 0.82 / Inhib 0.58 | 0.24 | **all-low z = low quality** | **Neuroblast_lowQC** |
| 2 | OPC 2.67 / Inhib 2.40 | 0.26 | genuine tie | needs resolve |
| 12 | Inhib 1.78 / OPC 1.50 | 0.28 | genuine tie | needs resolve |
| 4 | Mural 1.39 / Vasc_fibro 0.90 | 0.49 | tie (perivascular) | needs resolve |

So only 2/4/12 truly needed a resolve run; 0/8/13/19 decided from the table + biology.

## The promiscuity catch (full_summary.txt)
The `marker_viz` numeric summary exposed why 2 & 12 "tied" OPC:
- cl2: GAD1 **82%**, GAD2 75% (clear inhibitory) — its "OPC" was **LHFPL3 82% + VCAN 56%** while the
  *specific* PDGFRA = 0.26%, CSPG4 = 0%.
- cl12: GAD2 73%, GAD1 64% — OPC again only VCAN 82% / LHFPL3 24%, PDGFRA/CSPG4 ≈ 0.
- **No cluster was PDGFRA+ CSPG4+ and non-vascular/astro → OPC absent as a distinct cluster.**
- → marker-library fix: OPC `PDGFRA,CSPG4,LHFPL3,VCAN` → **`PDGFRA,CSPG4,VCAN`** (drop LHFPL3).

## The resolve runs
`coexp` adjudicated by specific genes:
- **cl2/cl12**: `GAD1,GAD2,SLC32A1 | PDGFRA,CSPG4,OLIG1,OLIG2,SOX10` → panel1 64%/46%, panel2 **0%/0%**
  → both **Inhibitory** (immature GABAergic; cl2 ARX/DLX6-AS1/NXPH1, cl12 MEIS2/DCC/GABRB3).
- **cl4 vs cl8 vs cl14**: `ACTA2,MYH11,TAGLN,MYL9,KCNJ8,RGS5 | DCN,COL1A1,COL1A2,LUM,PDGFRA`
  - cl4: ACTA2 .98/TAGLN .99/MYH11 .94/MYLK/LMOD1 → **Mural-SMC**; coexp 99%/6%.
  - cl8: HIGD1B/NDUFA4L2/KCNJ8/PDGFRB, no contractile → **Mural-pericyte**.
  - cl14: collagens DCN/COL1A1/COL4A1 → **Vasc_fibro/VLMC**. (double 28% = shared stromal, NOT doublet.)
- **cl10** (Proliferating): `GAD1,GAD2,DCX,STMN2,SOX11 | PDGFRA,CSPG4,OLIG1,OLIG2,SOX10` → microglia
  markers ≈0, OPC 1%, neuronal 37%; top non-cycle markers **SOX2/PTPRZ1/PANTR1/NFIB/MEIS2** →
  **Cycling neural progenitor** (NOT proliferating microglia — corrected a wrong prior). The `sub=10`
  DE failed here ("No DE genes" + PrepSCTFindMarkers warning) → motivated the script fix.
- **cl18**: `AQP4,GFAP,SLC1A2 | OPC` → astro 40%/OPC 8%, double 1% → **Astrocyte** (minor OPC mixed in).

## Final annotation (20/20)
0/1/3 Microglia · 16 Microglia_lowQC · 15 Perivasc_Macro · 11 Tcells · 9/17 Oligodendrocyte ·
18 Astrocyte · 7 Vasc_endo · 19 Vasc_endo(small/mixed) · 8 Mural-pericyte · 4 Mural-SMC ·
14 Vasc_fibro · 6 Excitatory · 5/2/12 Inhibitory · 13 Neuroblast_lowQC · 10 Cycling_NPC. **OPC absent.**

## Lessons that became hard rules
1. A "winning" OPC panel was one promiscuous gene (VCAN/LHFPL3); always check PDGFRA+CSPG4.
2. All-low z (cl13) = low quality, not a rare type — read nFeature/%ribo.
3. `Proliferating` (cl10) hid a neural progenitor; strip cycle genes, coexp the lineages.
4. coexp double 28% (cl14) was shared markers, not doublet — only both-high counts.
5. `sub=` on SCT silently produced no DE → must PrepSCTFindMarkers the subset (fixed in the script).
6. MotorNeuron/Ependymal/glycinergic panels were globally 0 → those populations absent; don't invent.
