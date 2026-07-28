# cell-type-annotation

A [Claude Code Skill](https://docs.claude.com/en/docs/claude-code/skills) that interprets single-cell
RNA-seq cluster-annotation output and helps decide each cluster's cell type — plus the R engine it
drives (scoring, ambiguous-cluster resolution, subtype/DEG/heterogeneity analysis, and matching plots).
It reads your pipeline's numbers and proposes a per-cluster label with the deciding evidence; it does
**not** run the R scripts on its own — it generates the exact command for you to run, reads the result,
and iterates.

一个 [Claude Code Skill](https://docs.claude.com/en/docs/claude-code/skills)，用于解读单细胞 RNA-seq
聚类注释的产出并判定每个 cluster 的细胞类型，同时附带驱动它的 R 引擎（打分、疑难 cluster 排查、亚型/
差异表达/异质性分析及配套画图）。它读取你流程产出的数值并给出带证据的逐 cluster 判定；它**不会**擅自
跑 R 脚本——只生成你要跑的确切命令，读结果，再迭代。

> **⚠️ You must review the marker genes yourself.** The built-in marker library
> (`references/cell_markers.cns.txt` + `references/marker_library.md`) is a well-sourced starting
> point with citations and documented traps — it is not a certified ground truth for your species,
> tissue, or disease context. Always check a call against the *specific* genes before trusting it, and
> correct the library when something is wrong. See [Marker review](#marker-review--marker-审核) below.
>
> **⚠️ 内置 marker 需要你自己审核。** 内置 marker 库（`references/cell_markers.cns.txt` +
> `references/marker_library.md`）是有文献引用、记录了已知陷阱的起点，**不是**适用于你的物种/组织/
> 疾病背景的既定真理。任何判定在采信前都要核对*具体基因*，发现问题要修正库文件。见下方
> [Marker review](#marker-review--marker-审核)。

---

## Contents / 目录
- [What this is](#what-this-is--这是什么)
- [Install](#install--安装)
- [Usage](#usage--用法)
- [Worked example](#worked-example--完整案例)
- [Marker review](#marker-review--marker-审核)
- [Self-iteration](#self-iteration--自我迭代)
- [Repo layout](#repo-layout--仓库结构)

---

## What this is / 这是什么

**English.** Single-cell cluster annotation usually goes through an engine that scores marker panels
per cluster and calls a label — but a chunk of clusters always come back `Unknown`, or get a label
that's wrong for a subtle reason (a promiscuous marker gene, an ambient-contamination artifact, a
cell-cycle state masquerading as a lineage). This skill is the decision layer on top of that engine:
it reads the consensus table, the threshold log, and any follow-up diagnostic output, and tells you —
with the specific genes as evidence — what each cluster actually is, or exactly which follow-up command
to run next. It ships:
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

**中文.** 单细胞聚类注释通常有个打分引擎，给每个 marker panel 打分再选出一个标签——但总有一部分 cluster
回来是 `Unknown`，或者标签因为一些不易察觉的原因判错了（某个基因其实不特异、环境污染的假阳性、细胞周期
状态被误读成谱系）。这个 skill 就是打分引擎之上的判定层：读取 consensus 表、阈值日志、以及任何后续诊断
输出，然后**用具体基因当证据**告诉你每个 cluster 到底是什么，或者接下来该跑哪条命令去确认。包含：
- **注释引擎**（`scripts/sc_004.*.R`）：逐 cluster 的 marker panel 打分（AddModuleScore + UCell）、
  margin/阈值判定、疑难 cluster 深挖、以及一个数值汇总/DotPlot 工具；
- **亚型/差异表达/异质性家族**（`scripts/sc_006.*.R`）：one-vs-all 亚型 marker、A-vs-B 差异表达、
  样本级 pseudobulk DESeq2、亚簇下钻，以及配套的 11 个画图脚本；
- **判定逻辑文档**（`references/decision_logic.md`）：把这些数字变成判定，同时避开已知的坑；
- **CNS marker 库**（`references/marker_library.md` + `.txt`）：作为起点，附文献引用与已记录的失败模式；
- **迭代日志**（`ITERATION_LOG.md`）：这个 skill 的记忆——每次犯的错和怎么修的，避免同一个错犯两次。

默认面向 CNS 组织（脑+脊髓，人+鼠），但判定逻辑和脚本家族与组织无关——换成你自己的 marker panel 就能用
在别的组织上。

---

## Install / 安装

**As a Claude Code skill.** Copy (or clone) this repo into Claude Code's skill discovery directory,
either globally or per-project, so `SKILL.md` sits at `<dir>/cell-type-annotation/SKILL.md`:

```bash
# global (all projects)
git clone <this-repo-url> ~/.claude/skills/cell-type-annotation

# or project-local
git clone <this-repo-url> <your-project>/.claude/skills/cell-type-annotation
```

Restart/open a new Claude Code conversation in that project and describe an annotation task — Claude
will pick the skill up automatically from `SKILL.md`'s description, or you can name it directly.

**Standalone (no Claude Code, just the R scripts).** Everything under `scripts/` is plain R, runnable
with `Rscript` on its own — the skill (`SKILL.md`, `references/`) is the decision layer, not a
dependency of the scripts. Requires R with `Seurat` (≥ 4; some plot scripts assume Seurat 5 layers),
`optparse`, `UCell`, and — only for the features that use them — `MAST`, `DESeq2`, `ggrepel`,
`pheatmap`. See `references/scripts_reference.md` for exact flags per script.

**作为 Claude Code skill 安装。** 把本仓库复制（或 clone）进 Claude Code 的 skill 发现目录，全局或按
项目均可，让 `SKILL.md` 落在 `<目录>/cell-type-annotation/SKILL.md`：

```bash
# 全局（所有项目可用）
git clone <this-repo-url> ~/.claude/skills/cell-type-annotation

# 或项目内
git clone <this-repo-url> <你的项目>/.claude/skills/cell-type-annotation
```

在该项目下重新打开/开始一个 Claude Code 对话，描述你的注释任务即可——Claude 会根据 `SKILL.md` 的
description 自动识别，也可以直接点名调用。

**独立使用（不用 Claude Code，只用 R 脚本）。** `scripts/` 下全是普通 R 脚本，`Rscript` 直接可跑——
skill 本体（`SKILL.md`、`references/`）是判定层，不是脚本的运行依赖。需要 R + `Seurat`（≥4；部分画图
脚本假设 Seurat 5 的 layer 结构）、`optparse`、`UCell`，以及仅在用到对应功能时才需要的 `MAST`、
`DESeq2`、`ggrepel`、`pheatmap`。每个脚本的确切参数见 `references/scripts_reference.md`。

---

## Usage / 用法

### Track 1 — first-pass whole-object annotation / 第一轨：整体首轮注释

**English.**
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

**中文.**
```bash
# 1. 用你的 marker 文件跑引擎（先从内置种子 panel 开始）
Rscript scripts/sc_004.CellType_Annotation.R rds=<你的.rds> \
    markers=references/cell_markers.cns.txt --presence_floor 0.4

# 2. 告诉 Claude（或按 references/decision_logic.md 自己读）：
#    "帮我判读 Harmony_thresholds_log.txt / annotation_consensus.csv 里的 cluster"
```
Claude 会读 `*_annotation_consensus.csv` + `*_thresholds_log.txt`，把每个 `Unknown` 归类为
vote-split（近邻票分）/ 真平局 / 低质量 / 状态问题，对真平局生成可直接跑的 `cluster_resolve` 命令
（同上）。读回 `*_resolve_summary.tsv` + 日志后，Claude 会给每个 cluster 一个带决定性基因证据的判定——
出最终表格前会按 `references/report_template.md` 的自查步骤过一遍。

### Track 2 — subtype markers / DEG / pseudobulk / heterogeneity / 第二轨：亚型/差异/异质性

For an object that's already subset/annotated. Full flags in `references/scripts_reference.md §4`.

**English example** — condition contrast within one cell type:
```bash
nohup Rscript scripts/sc_006.find_marker.DEG.Heter.R -i <your.rds> \
    --subset_column cell_type --subset_value Fibroblast \
    -c condition --contrasts 'Disease,Control' \
    --do_deg --do_pseudobulk --sample_column orig.ident \
    > deg.log 2>&1 &
```
The compute driver prints all 11 matching plot commands at the end of its run — copy the ones you want
(dotplot for annotation markers, volcano for DEG/pseudobulk, heatmap/violin/dimplot as needed).

**中文举例** —— 某个细胞类型内部的条件对比：
```bash
nohup Rscript scripts/sc_006.find_marker.DEG.Heter.R -i <你的.rds> \
    --subset_column cell_type --subset_value Fibroblast \
    -c condition --contrasts 'Disease,Control' \
    --do_deg --do_pseudobulk --sample_column orig.ident \
    > deg.log 2>&1 &
```
计算脚本跑完会自动打印 11 个对应画图命令——复制你要的那些即可（注释类 marker 用 dotplot，DEG/pseudobulk
用 volcano，heatmap/violin/dimplot 按需）。

---

## Worked example / 完整案例

`examples/example_cns_dataset.md` walks a full 20-cluster round end-to-end: reading the first-pass
table, classifying all 7 `Unknown` clusters (vote-splits, a genuine tie, a low-quality cluster), the
marker-promiscuity catch that fixed the OPC panel, four `coexp` resolve runs, and the final annotation
— useful as a template for how to work through your own dataset.

`examples/example_cns_dataset.md` 完整走一遍 20-cluster 的判读流程：读首轮表格、给全部 7 个 `Unknown`
分类（vote-split、真平局、低质量各一种）、抓出促成 OPC panel 修正的 marker 促进性陷阱、四次 `coexp`
排查，以及最终注释结果——可以直接当作梳理自己数据集时的模板。

---

## Marker review / marker 审核

**English.** The shipped marker library (`references/cell_markers.cns.txt` for the machine-readable
panels, `references/marker_library.md` for rationale + citations + traps) was built and validated
against real CNS datasets — every panel has a stated source or an explicit "unverified" flag, and every
documented trap (e.g. `VCAN`/`LHFPL3` producing a false OPC call, `ATP1A2` being miscited as a pericyte
marker when it's actually astrocytic) was a real mistake caught and fixed. That track record is exactly
why you shouldn't take it on faith for your own data:
- **Different species/tissue/disease context can break a panel that worked elsewhere.** Check every
  gene in a panel against your own literature before trusting a call built on it.
- **Verify against your own data before trusting a panel member**, the way the worked examples do:
  check pct.1/pct.2/log2FC for each gene, not just whether the panel "won."
- **When you find a wrong or missing marker, fix it** — edit `cell_markers.cns.txt` /
  `marker_library.md`, and log the observation in `ITERATION_LOG.md` (see below) so the fix is durable
  and explained for the next reader.

Never treat a cell-type call as final because a script produced it — the scripts scores panels
mechanically; deciding whether the panel itself is trustworthy for your data is on you.

**中文.** 内置 marker 库（机读部分在 `references/cell_markers.cns.txt`，理由+引用+陷阱在
`references/marker_library.md`）是针对真实 CNS 数据集搭建并验证过的——每个 panel 都标注了来源，或者
明确标"未验证"；每一条记录在案的陷阱（比如 `VCAN`/`LHFPL3` 会造出假的 OPC 判定，`ATP1A2` 被误引成
pericyte marker 其实是星形胶质细胞 marker）都是真实踩过的坑、修过之后留下的记录。正因为有这样的记录，
才更不应该对自己的数据无条件采信：
- **物种/组织/疾病背景不同，别处验证过的 panel 可能失效。** 用某个 panel 判定前，先核对它每个基因在你
  自己领域文献里的依据。
- **采信一个 panel 成员前，用自己的数据验证**——像案例里那样，核对每个基因的 pct.1/pct.2/log2FC，
  而不只是看 panel 有没有"赢"。
- **发现基因错了或缺了就去改**——编辑 `cell_markers.cns.txt` / `marker_library.md`，并把观察记录进
  `ITERATION_LOG.md`（见下）,让这次修正长期有效、给后来者留下依据。

不要因为脚本给出了一个判定就当作定论——脚本只是机械地给 panel 打分；这个 panel 对你的数据是否可信，
要靠你自己判断。

---

## Self-iteration / 自我迭代

**English.** `ITERATION_LOG.md` is not documentation you write once — it's meant to grow every time you
(or Claude) use this skill and find something wrong. Every trap in `decision_logic.md` and every
correction in `marker_library.md` started as a dated entry here: **Problem → Diagnosis → Fix → Where it
landed**, with the actual deciding numbers, not just the conclusion. Two rules keep this working:
1. **Read it on any new problem** — check whether the failure mode has already been diagnosed before
   re-deriving it from scratch.
2. **Append to it whenever a verdict was wrong, a trap was found, or a marker/script was changed** —
   include the evidence, so the next reader (human or Claude) can verify it rather than take it on
   faith.

This is how the skill was actually built: nearly every rule in `references/decision_logic.md` §3 (the
marker-promiscuity traps) and every correction in `references/marker_library.md` exists because a wrong
call was caught, logged, and fixed first. Keep extending it the same way for your own data.

**中文.** `ITERATION_LOG.md` 不是写一次就完事的文档——它应该随着你（或 Claude）每次使用这个 skill、
发现问题而持续增长。`decision_logic.md` 里的每一条陷阱、`marker_library.md` 里的每一次修正，最初都是
这里的一条带日期的记录：**问题 → 诊断 → 修复 → 落地位置**，附带实际的决定性数字，而不只是结论。两条
规则让这套机制真正起作用：
1. **遇到新问题先读它**——先看这个失败模式是不是已经被诊断过，而不是从头再推一遍。
2. **每次判定判错了、发现新陷阱、或改了 marker/脚本，就补一条**——附上证据，让下一个读者（人或 Claude）
   能核实，而不是凭空相信。

这个 skill 本来就是这么建起来的：`references/decision_logic.md` §3（marker 促进性陷阱）里几乎每一条
规则、`references/marker_library.md` 里每一次修正，起点都是一次判错被抓到、记录下来、修好。用同样的
方式，在你自己的数据上继续把它养大。

---

## Repo layout / 仓库结构

```
SKILL.md                       entry point — when to use + the workflow + hard rules
README.md                      this file
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
```

## License

MIT — see [LICENSE](LICENSE).
