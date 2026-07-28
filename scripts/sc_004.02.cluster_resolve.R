#!/usr/bin/env Rscript
# ============================================================================
# sc_004.02.cluster_resolve.R —— 第一遍注释后, 对【分不清的暧昧簇】单独深挖并辅助判定
# ----------------------------------------------------------------------------
# 设计目标: 不是简单 dump marker, 而是给"候选 + 证据":
#   1) QC 体检    —— 先排除"低质量/红细胞污染"被误当稀有类型
#   2) 特异性 marker (正向)  —— 按 pct.diff 排, 比纯 log2FC 更可信; 存全量 CSV
#   3) 负向 marker (缺失证据) —— "有谱系但缺 X" 往往才是定性关键
#   4) 相似度排序 (anno=)     —— 暧昧簇最像"哪个已命名类型"(pseudobulk 相关)
#   5) 成对比较  (pair=)      —— 谱系内两簇到底差在哪 (one-vs-rest 看不出)
#   6) 就地重聚类 (sub=)      —— "暧昧"其实是"两类硬凑"时, 看能否裂开
#                               (子集含多 SCT 模型时会自动重跑 PrepSCTFindMarkers, 否则子簇 DE 静默失败)
#   7) 共表达体检 (coexp=)    —— 同时高表达两套谱系 marker → 区分 doublet vs 共享 marker
#                               双阳须 >dblwarn(默认20%) 且两套都高才提示 doublet; 否则只标"偏 panel1/2"
#   每次还落盘 <base>_resolve_summary.tsv: per-cluster QC+top5+coexp 一行, 供 skill 一次读全簇
#
# 正确性前提:
#   - SCT 多样本: 循环前自动 PrepSCTFindMarkers (否则 fold change 偏倚)
#   - Seurat v5 RNA: 自动 JoinLayers
#   - 装了 presto 会被 Seurat 自动调用, Wilcoxon 快数倍 (建议 install.packages 不到)
#     用 remotes::install_github('immunogenomics/presto')
#
# 物种说明: 默认人类基因符号 (MT-/RPL/RPS/HB[AB])。小鼠改成 ^mt-/^Rp[ls]/^Hb 即可。
# ----------------------------------------------------------------------------
# 用法 (key=value, 顺序无关; 只有 rds 与 clusters 必填):
#   Rscript sc_004.02.cluster_resolve.R rds=<path> clusters=<逗号分隔暧昧簇> [选项...]
#
# 参数:
#   rds=        必填  Seurat 对象 .rds (你的 *_final.rds / *_annotated.rds)
#   clusters=   必填  要深挖的暧昧簇, 逗号分隔, 如 14,18,23
#   assay=      可选  默认: 有 SCT 用 SCT, 否则 DefaultAssay
#   clu=        可选  聚类列名, 默认自动探测 astro_clu/seurat_clusters/harmony_clusters
#   anno=       可选  已命名类型所在的列 (如 celltype)。给了就做【相似度排序】
#   pair=       可选  成对比较, 如 14:1 或多对 14:1,18:6  (左=ident.1, 右=ident.2)
#   sub=        可选  对某一个簇就地重聚类, 如 sub=14
#   subres=     可选  sub 的分辨率, 默认 0.6
#   subdims=    可选  sub 用的维度数, 默认 30
#   coexp=      可选  两套 marker 用 | 分隔做共表达检查 (shell 里务必加引号!):
#                     'GFAP,AQP4|PTPRC,CSF1R'
#   topN=       可选  打印的正向 marker 数, 默认 40
#   neg=        可选  是否打印负向 marker, 默认 true
#   minpct=     可选  FindMarkers min.pct, 默认 0.1
#   logfc=      可选  FindMarkers logfc.threshold, 默认 0.25
#   dblwarn=    可选  coexp 双阳 doublet 提示阈值(%), 默认 20
#
# 示例 (见文件末尾或 README, 以及下方 README 段)
# ============================================================================

suppressPackageStartupMessages({ library(Seurat); library(Matrix) })

# ---------- 参数解析 (key=value) ----------
argv <- commandArgs(TRUE)
kv <- list()
for (a in argv) if (grepl("=", a)) {
  k <- sub("=.*", "", a); v <- sub("^[^=]*=", "", a); kv[[k]] <- v
}
opt <- function(name, default = NA) if (!is.null(kv[[name]]) && nzchar(kv[[name]])) kv[[name]] else default

print_usage <- function() {
  cat("\n用法: Rscript sc_004.02.cluster_resolve.R rds=<path> clusters=<14,18,23> [选项]\n",
      "必填: rds=  clusters=\n",
      "可选: assay= clu= anno= pair=14:1 sub=14 subres=0.6 subdims=30\n",
      "      coexp='GFAP,AQP4|PTPRC,CSF1R' topN=40 neg=true minpct=0.1 logfc=0.25\n",
      "示例:\n",
      "  Rscript sc_004.02.cluster_resolve.R rds=sub_astro/astro_R2_final.rds clusters=14,18 anno=celltype\n",
      "  Rscript sc_004.02.cluster_resolve.R rds=astro_R2_final.rds clusters=14 pair=14:1,14:6\n",
      "  Rscript sc_004.02.cluster_resolve.R rds=astro_R2_final.rds clusters=14 sub=14 subres=0.8\n",
      "  Rscript sc_004.02.cluster_resolve.R rds=astro_R2_final.rds clusters=23 coexp='GFAP,AQP4|PTPRC,CSF1R'\n\n",
      sep = "")
}

RDS <- opt("rds")
if (is.na(RDS) || is.na(opt("clusters"))) { print_usage(); stop("缺少必填参数 rds= 或 clusters=") }
if (!file.exists(RDS)) stop(sprintf("找不到文件: %s", RDS))

CLS    <- trimws(strsplit(opt("clusters"), ",")[[1]])
ANNO   <- opt("anno")
PAIR   <- opt("pair")
SUB    <- opt("sub")
SUBRES <- as.numeric(opt("subres", "0.6"))
SUBDIM <- as.integer(opt("subdims", "30"))
COEXP  <- opt("coexp")
TOPN   <- as.integer(opt("topN", "40"))
NEG    <- tolower(opt("neg", "true")) %in% c("true", "t", "1", "yes")
MINPCT <- as.numeric(opt("minpct", "0.1"))
LOGFC  <- as.numeric(opt("logfc", "0.25"))
DBLWARN<- as.numeric(opt("dblwarn", "20"))   # 共表达双阳 > 此% 才提示 doublet (默认 20)

# ---------- 读对象 + 选列/assay ----------
obj <- readRDS(RDS)

CLU <- opt("clu")
if (is.na(CLU)) {
  CLU <- if ("astro_clu"       %in% colnames(obj@meta.data)) "astro_clu" else
         if ("seurat_clusters" %in% colnames(obj@meta.data)) "seurat_clusters" else
         if ("harmony_clusters"%in% colnames(obj@meta.data)) "harmony_clusters" else NA
}
if (is.na(CLU) || !CLU %in% colnames(obj@meta.data)) {
  cands <- grep("seurat_clusters|_snn_res\\.|clu", colnames(obj@meta.data),
                value = TRUE, ignore.case = TRUE)
  stop(sprintf("聚类列 '%s' 不存在。可用候选: %s", CLU, paste(cands, collapse = ", ")))
}

ASSAY <- opt("assay")
if (is.na(ASSAY)) ASSAY <- if ("SCT" %in% Assays(obj)) "SCT" else DefaultAssay(obj)
DefaultAssay(obj) <- ASSAY

# ---------- 正确性预处理 (SCT / v5 RNA) ----------
if (ASSAY == "SCT") {
  obj <- tryCatch(PrepSCTFindMarkers(obj),
                  error = function(e) { cat("[warn] PrepSCTFindMarkers 失败, 继续:", conditionMessage(e), "\n"); obj })
}
if (ASSAY == "RNA" && inherits(obj[["RNA"]], "Assay5")) {
  obj <- tryCatch(JoinLayers(obj), error = function(e) obj)
}

Idents(obj) <- factor(as.character(obj@meta.data[[CLU]]))
outdir <- dirname(RDS); base <- sub("\\.rds$", "", basename(RDS))

# ---------- 预计算 QC 百分比 (不覆盖已有 percent.mt) ----------
qc_assay <- if ("RNA" %in% Assays(obj)) "RNA" else ASSAY
safe_pct <- function(pat) tryCatch(PercentageFeatureSet(obj, pattern = pat, assay = qc_assay),
                                   error = function(e) rep(NA_real_, ncol(obj)))
obj$.pmt <- if ("percent.mt" %in% colnames(obj@meta.data)) obj$percent.mt else safe_pct("^MT-")
obj$.prb <- safe_pct("^(RPL|RPS)")
obj$.phb <- safe_pct("^HB[AB]")

nf_col <- grep("^nFeature", colnames(obj@meta.data), value = TRUE)[1]
nc_col <- grep("^nCount",   colnames(obj@meta.data), value = TRUE)[1]

cat(sprintf(">>> %s | clu=%s | assay=%s | 暧昧簇: %s\n",
            RDS, CLU, ASSAY, paste(CLS, collapse = ",")))
cat(sprintf(">>> 模式: markers%s%s%s%s\n",
            if (!is.na(ANNO)) " + 相似度(anno)" else "",
            if (!is.na(PAIR)) " + 成对(pair)" else "",
            if (!is.na(SUB))  " + 重聚类(sub)" else "",
            if (!is.na(COEXP)) " + 共表达(coexp)" else ""))

# ============================================================================
# A) 逐簇: QC 体检 + 正向/负向 marker
# ============================================================================
display_filter <- function(g) !grepl("^MT-|^(RPL|RPS)|^HB[AB]", g)   # 显示层去 MT/核糖体/血红蛋白
topmark_long <- list()   # 汇总宽表用
summ <- list()           # 机器可读 per-cluster 汇总 (供 skill 一次读全簇)

for (k in CLS) {
  if (!k %in% levels(Idents(obj))) { cat("\n==", k, "不存在==\n"); next }
  sel <- obj@meta.data[[CLU]] == k
  n   <- sum(sel)
  cat(sprintf("\n######### cluster %s (n=%d) #########\n", k, n))

  # --- QC 体检 ---
  med <- function(x) round(median(x[sel], na.rm = TRUE), 1)
  cat(sprintf("[QC] %s中位=%s | %s中位=%s | %%MT=%.1f | %%核糖体=%.1f | %%HB=%.1f\n",
              nf_col, med(obj@meta.data[[nf_col]]),
              nc_col, med(obj@meta.data[[nc_col]]),
              median(obj$.pmt[sel], na.rm = TRUE),
              median(obj$.prb[sel], na.rm = TRUE),
              median(obj$.phb[sel], na.rm = TRUE)))
  flags <- c()
  if (median(obj@meta.data[[nf_col]][sel], na.rm = TRUE) < 500) flags <- c(flags, "基因数偏低(疑低质量)")
  if (median(obj$.pmt[sel], na.rm = TRUE) > 15) flags <- c(flags, "%MT偏高(疑凋亡/低质量)")
  if (median(obj$.phb[sel], na.rm = TRUE) > 5)  flags <- c(flags, "%HB偏高(疑红细胞污染)")
  if (length(flags)) cat("   !! 警示:", paste(flags, collapse = "; "), "—— 注释前先考虑是否丢弃\n")

  # --- 一次 FindMarkers (双向), 拆正/负 ---
  mk <- tryCatch(FindMarkers(obj, ident.1 = k, only.pos = FALSE,
                             min.pct = MINPCT, logfc.threshold = LOGFC),
                 error = function(e) { cat("[skip]", conditionMessage(e), "\n"); NULL })
  if (is.null(mk)) next
  mk$pct.diff <- round(mk$pct.1 - mk$pct.2, 3)
  csv <- file.path(outdir, sprintf("%s_cl%s_markers.csv", base, k))
  write.csv(mk[order(-mk$avg_log2FC), ], csv, row.names = TRUE)

  sig <- mk[mk$p_val_adj < 0.05, ]
  pos <- sig[sig$avg_log2FC > 0, , drop = FALSE]
  pos <- pos[order(-pos$pct.diff, -pos$avg_log2FC), , drop = FALSE]
  pos_disp <- pos[display_filter(rownames(pos)), , drop = FALSE]
  cat(sprintf("[存] %s (%d genes) | 显著正向 %d\n", csv, nrow(mk), nrow(pos)))
  cat(sprintf("[正向 top%d | 去MT/核糖体/HB | 按 pct.diff 排 (特异性)]\n", TOPN))
  print(head(pos_disp[, c("avg_log2FC","pct.1","pct.2","pct.diff","p_val_adj")], TOPN))

  if (NEG) {
    neg <- sig[sig$avg_log2FC < 0, , drop = FALSE]
    neg <- neg[order(neg$avg_log2FC), , drop = FALSE]
    neg_disp <- neg[display_filter(rownames(neg)), , drop = FALSE]
    cat("[负向 top15 | 该簇相对缺失的基因 → 用于排除候选类型]\n")
    print(head(neg_disp[, c("avg_log2FC","pct.1","pct.2","p_val_adj")], 15))
  }

  # 汇总宽表数据
  tp <- head(pos_disp, TOPN)
  if (nrow(tp)) topmark_long[[k]] <- data.frame(
    cluster = k, rank = seq_len(nrow(tp)), gene = rownames(tp),
    avg_log2FC = round(tp$avg_log2FC, 3), pct.diff = tp$pct.diff,
    p_val_adj = signif(tp$p_val_adj, 3), row.names = NULL)

  # 机器可读汇总行 (coexp 列稍后在 E 段补)
  summ[[k]] <- list(
    cluster = k, n = n,
    nFeature = med(obj@meta.data[[nf_col]]), nCount = med(obj@meta.data[[nc_col]]),
    pct_mt = round(median(obj$.pmt[sel], na.rm = TRUE), 1),
    pct_ribo = round(median(obj$.prb[sel], na.rm = TRUE), 1),
    pct_hb = round(median(obj$.phb[sel], na.rm = TRUE), 1),
    qc_flag = if (length(flags)) paste(flags, collapse = ";") else "",
    top5_pos = paste(head(rownames(pos_disp), 5), collapse = ","),
    coexp_p1 = NA_real_, coexp_p2 = NA_real_, coexp_dbl = NA_real_)
}

# 所有暧昧簇 top marker 汇总成一张长表 (便于 pivot 横向比对)
if (length(topmark_long)) {
  allmk <- do.call(rbind, topmark_long)
  sumcsv <- file.path(outdir, sprintf("%s_ambig_topmarkers.csv", base))
  # 追加 + 按 cluster 去重: 本次跑到的簇覆盖旧行, 未涉及的簇保留, 让多次分批 resolve 累积成一张完整表
  if (file.exists(sumcsv)) {
    old <- tryCatch(read.csv(sumcsv, stringsAsFactors = FALSE), error = function(e) NULL)
    if (!is.null(old) && "cluster" %in% colnames(old)) {
      old <- old[!as.character(old$cluster) %in% as.character(allmk$cluster), , drop = FALSE]
      allmk <- rbind(old, allmk)
    }
  }
  write.csv(allmk, sumcsv, row.names = FALSE)
  cat(sprintf("\n[汇总] 各暧昧簇 top marker 长表已存(追加/去重): %s (累计 %d 簇)\n",
              sumcsv, length(unique(allmk$cluster))))
}

# ============================================================================
# B) 相似度排序: 暧昧簇最像哪个已命名类型 (pseudobulk Spearman 相关)
# ============================================================================
if (!is.na(ANNO)) {
  if (!ANNO %in% colnames(obj@meta.data)) {
    cat(sprintf("\n[sim 跳过] 注释列 '%s' 不存在\n", ANNO))
  } else {
    cat(sprintf("\n========== 相似度: 暧昧簇 vs 已命名类型 (列=%s) ==========\n", ANNO))
    grp <- as.character(obj@meta.data[[ANNO]])
    clu <- as.character(obj@meta.data[[CLU]])
    grp[is.na(grp) | grp == ""] <- "unassigned"
    amb <- clu %in% CLS
    grp[amb] <- paste0("cl", clu[amb])          # 暧昧簇用 cl<id> 独立成组
    obj$.simgrp <- grp

    avg <- AverageExpression(obj, assays = ASSAY, group.by = ".simgrp")[[ASSAY]]
    feat <- intersect(VariableFeatures(obj), rownames(avg))
    if (length(feat) < 50) feat <- rownames(avg)            # 兜底
    cm <- cor(as.matrix(avg[feat, , drop = FALSE]), method = "spearman")

    amb_cols  <- intersect(paste0("cl", CLS), colnames(cm))
    type_cols <- setdiff(colnames(cm), c(paste0("cl", CLS), "unassigned"))
    for (cc in amb_cols) {
      sims <- sort(cm[type_cols, cc], decreasing = TRUE)
      cat(sprintf("\n[%s] 最像的已命名类型 (Spearman, top5):\n", cc))
      for (i in seq_len(min(5, length(sims))))
        cat(sprintf("   %-22s %.3f\n", names(sims)[i], sims[i]))
    }
  }
}

# ============================================================================
# C) 成对比较: 谱系内两簇差在哪 (one-vs-rest 看不出共享 marker 抵消后的差异)
# ============================================================================
if (!is.na(PAIR)) {
  cat("\n========== 成对比较 (ident.1 vs ident.2) ==========\n")
  for (p in strsplit(PAIR, ",")[[1]]) {
    ids <- trimws(strsplit(p, ":")[[1]])
    if (length(ids) != 2) { cat("[pair 跳过] 格式应为 a:b ——", p, "\n"); next }
    i1 <- ids[1]; i2 <- ids[2]
    if (!all(c(i1, i2) %in% levels(Idents(obj)))) { cat("[pair 跳过] 簇不存在:", p, "\n"); next }
    cat(sprintf("\n--- %s vs %s ---\n", i1, i2))
    mk <- tryCatch(FindMarkers(obj, ident.1 = i1, ident.2 = i2,
                               min.pct = MINPCT, logfc.threshold = LOGFC),
                   error = function(e) { cat("[skip]", conditionMessage(e), "\n"); NULL })
    if (is.null(mk)) next
    mk <- mk[order(-mk$avg_log2FC), ]
    pc <- file.path(outdir, sprintf("%s_pair_%sv%s.csv", base, i1, i2))
    write.csv(mk, pc, row.names = TRUE)
    keep <- display_filter(rownames(mk))
    cat(sprintf("[存] %s\n[偏向 %s (log2FC>0) top12]\n", pc, i1))
    print(head(mk[keep, c("avg_log2FC","pct.1","pct.2","p_val_adj")], 12))
    cat(sprintf("[偏向 %s (log2FC<0) top12]\n", i2))
    print(head(mk[keep, ][order(mk[keep, "avg_log2FC"]), c("avg_log2FC","pct.1","pct.2","p_val_adj")], 12))
  }
}

# ============================================================================
# D) 就地重聚类: "暧昧"是否其实是"两类硬凑" → 看能否裂开
#    (轻量做法: 复用已有 reduction; 严格做法应对子集重新 PCA/harmony)
# ============================================================================
if (!is.na(SUB)) {
  cat(sprintf("\n========== 就地重聚类 cluster %s (res=%.2f, dims=%d) ==========\n", SUB, SUBRES, SUBDIM))
  if (!SUB %in% levels(Idents(obj))) {
    cat("[sub 跳过] 簇不存在:", SUB, "\n")
  } else {
    reds <- Reductions(obj)
    RED  <- if ("harmony" %in% reds) "harmony" else if ("pca" %in% reds) "pca" else reds[1]
    sc <- subset(obj, idents = SUB)
    maxd <- min(SUBDIM, ncol(Embeddings(sc, RED)))
    sc <- FindNeighbors(sc, reduction = RED, dims = 1:maxd, verbose = FALSE)
    sc <- FindClusters(sc, resolution = SUBRES, verbose = FALSE)
    cat(sprintf("[reduction=%s, dims=1:%d] 裂成 %d 个子簇:\n", RED, maxd, length(levels(Idents(sc)))))
    print(table(Idents(sc)))
    if (length(levels(Idents(sc))) > 1) {
      # 修复: 子集含多个 SCT 模型, 必须在子集上重跑 PrepSCTFindMarkers, 否则 FindAllMarkers 静默失败
      if (ASSAY == "SCT")
        sc <- tryCatch(PrepSCTFindMarkers(sc),
                       error = function(e) { cat("[warn] 子集 PrepSCTFindMarkers 失败:", conditionMessage(e), "\n"); sc })
      am <- tryCatch(FindAllMarkers(sc, only.pos = TRUE, min.pct = MINPCT,
                                    logfc.threshold = LOGFC, verbose = FALSE),
                     error = function(e) NULL)
      if (!is.null(am) && nrow(am)) {
        am <- am[display_filter(am$gene), ]
        for (sub_k in levels(Idents(sc))) {
          tt <- am[am$cluster == sub_k, ]
          tt <- head(tt[order(-(tt$pct.1 - tt$pct.2)), "gene"], 12)
          cat(sprintf("  子簇 %s top: %s\n", sub_k, paste(tt, collapse = ", ")))
        }
      }
    } else cat("  只有 1 个子簇 —— 该簇内部较均一, 不像是混杂。\n")
  }
}

# ============================================================================
# E) 共表达体检: 同时高表达两套谱系 marker → 优先怀疑 doublet
#    判定: 某细胞对一套 panel "阳性" = 检测到 (>0) 该 panel 一半以上基因
# ============================================================================
if (!is.na(COEXP)) {
  cat("\n========== 共表达 / doublet 体检 ==========\n")
  panels <- strsplit(COEXP, "\\|")[[1]]
  if (length(panels) != 2) {
    cat("[coexp 跳过] 需要两套 panel 用 | 分隔, 如 'GFAP,AQP4|PTPRC,CSF1R'\n")
  } else {
    expr <- tryCatch(GetAssayData(obj, assay = ASSAY, layer = "data"),
                     error = function(e) GetAssayData(obj, assay = ASSAY, slot = "data"))
    pos_of <- function(genes_str, label) {
      genes_in <- trimws(strsplit(genes_str, ",")[[1]])
      genes <- intersect(genes_in, rownames(expr))
      if (!length(genes)) {
        cat(sprintf("[警告] coexp %s: '%s' 里没有一个基因在对象中找到 —— 该 panel 会被判全阴性(恒为0%%)!\n",
                    label, genes_str))
        cat("       多半是把细胞类型/簇名当成基因列表传了 (coexp= 两边都必须是逗号分隔的真实基因symbol)。\n")
        return(rep(FALSE, ncol(expr)))
      }
      if (length(genes) < length(genes_in))
        cat(sprintf("[提示] coexp %s: %d/%d 个基因未在对象中找到, 已忽略: %s\n",
                    label, length(genes_in) - length(genes), length(genes_in),
                    paste(setdiff(genes_in, genes), collapse = ", ")))
      det <- Matrix::colSums(expr[genes, , drop = FALSE] > 0)
      det >= ceiling(length(genes) / 2)
    }
    p1 <- pos_of(panels[1], "panel1"); p2 <- pos_of(panels[2], "panel2")
    cat(sprintf("panel1 = %s | panel2 = %s\n", panels[1], panels[2]))
    any_dbl <- FALSE
    for (k in CLS) {
      sel <- obj@meta.data[[CLU]] == k
      if (!sum(sel)) next
      a <- 100*mean(p1[sel]); b <- 100*mean(p2[sel]); d <- 100*mean(p1[sel] & p2[sel])
      # 标注: 谁占优 / 是否真双阳 (双阳须同时 >DBLWARN 且两套都不低, 否则只是共享 marker)
      note <- if (d > DBLWARN && min(a, b) > DBLWARN) "← 双阳偏高, 查 doublet"
              else if (a >= b) "→ 偏 panel1"
              else "→ 偏 panel2"
      cat(sprintf("[cluster %s] panel1+ %2.0f%% | panel2+ %2.0f%% | 双阳 %2.0f%%  %s\n", k, a, b, d, note))
      if (d > DBLWARN && min(a, b) > DBLWARN) any_dbl <- TRUE
      if (!is.null(summ[[k]])) { summ[[k]]$coexp_p1 <- round(a,0); summ[[k]]$coexp_p2 <- round(b,0); summ[[k]]$coexp_dbl <- round(d,0) }
    }
    if (any_dbl) cat(sprintf("(有簇双阳 >%g%% 且两套 panel 均高 → 真 doublet 嫌疑, 配 scDblFinder 复核)\n", DBLWARN))
    else         cat(sprintf("(无簇双阳超 %g%%; 偏向某 panel = 该谱系, 高双阳但另一套低 = 共享 marker 而非 doublet)\n", DBLWARN))
  }
}

# ============================================================================
# F) 机器可读 per-cluster 汇总 TSV (供 skill 一次读全簇定标签)
# ============================================================================
if (length(summ)) {
  sdf <- do.call(rbind, lapply(summ, function(r) as.data.frame(r, stringsAsFactors = FALSE)))
  stsv <- file.path(outdir, sprintf("%s_resolve_summary.tsv", base))
  # 追加 + 按 cluster 去重, 逻辑同上面的 ambig_topmarkers.csv
  if (file.exists(stsv)) {
    old <- tryCatch(read.delim(stsv, stringsAsFactors = FALSE), error = function(e) NULL)
    if (!is.null(old) && "cluster" %in% colnames(old)) {
      old <- old[!as.character(old$cluster) %in% as.character(sdf$cluster), , drop = FALSE]
      sdf <- rbind(old, sdf)
    }
  }
  ord <- suppressWarnings(as.numeric(sdf$cluster))
  if (!anyNA(ord)) sdf <- sdf[order(ord), , drop = FALSE]   # 簇号是数字时按数值排序, 否则保持原序
  write.table(sdf, stsv, sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("\n[汇总TSV] per-cluster 摘要(QC+top5+coexp)已存(追加/去重): %s (累计 %d 簇)\n",
              stsv, nrow(sdf)))
}

cat("\n[done]\n")
