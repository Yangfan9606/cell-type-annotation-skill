#!/usr/bin/env Rscript
# ==============================================================================
# sc_006.find_marker.DEG.Heter.R
#   单细胞 marker / DEG / 异质性 —— 一体化【计算】脚本 (画图分离到 sc_006.plot.*.R)
# ------------------------------------------------------------------------------
# 四个开关, 相互独立, 可任意组合同时计算:
#   --do_find_marker  one-vs-all 共识 marker(所有组; Wilcox∩ROC→MAST)  → *_find_marker_*_final.csv   [注释]
#   --do_deg          FindMarkers(A vs B) 逐对比(可批量)               → *_DEG_<A>_vs_<B>.csv         [组间生物学]
#   --do_pseudobulk   DESeq2 样本聚合(需 --sample_column)              → *_pseudobulk_<A>_vs_<B>.csv  [稳健组间]
#   --do_hetero       亚群下钻(--compare_group X, 重聚类)              → *_hetero_<X>_*.csv + .rds    [亚结构]
#
# "谁 vs 谁" 由 -c(ident 列) + --subset 决定:
#   c1 vs c2 (细胞类型之间)          -c cell_type            --contrasts "c1,c2"
#   c1 内部 Disease vs Ctrl (最常用⭐)   --subset_column cell_type --subset_value c1 -c condition --contrasts "Disease,Control"
#   全局 Disease vs Ctrl                 -c condition            --contrasts "Disease,Control"
#   (--contrasts 批量: "A,B;A,C"  用 ; 分对比, 用 , 分两侧; 右侧写 ALL = one-vs-rest)
#
# 画图: 计算完出 CSV / rds / manifest, 再独立跑 11 个画图脚本 (它们默认读 manifest):
#   find_marker : dotplot  heatmap  violin
#   deg         : volcano  heatmap  violin
#   pseudobulk  : volcano  heatmap
#   hetero      : dimplot  dotplot  heatmap
#
# SCT / v5 split-layer 在开头统一处理一次。
# ==============================================================================

suppressPackageStartupMessages(library(optparse))

USAGE <- paste(
"=== sc_006.find_marker.DEG.Heter.R — 计算脚本 (画图见 sc_006.plot.*.R) ===",
"",
"开关(可任意组合):",
"  --do_find_marker   one-vs-all 共识 marker(所有组)      → *_find_marker_*_final.csv   [注释]",
"  --do_deg           FindMarkers(A vs B) 逐对比           → *_DEG_<A>_vs_<B>.csv        [组间]",
"  --do_pseudobulk    DESeq2 样本聚合(需 --sample_column)  → *_pseudobulk_<A>_vs_<B>.csv [稳健组间]",
"  --do_hetero        亚群下钻(--compare_group X)          → *_hetero_<X>_*.csv + .rds   [亚结构]",
"",
"谁 vs 谁 (由 -c 列 + --subset 决定):",
"  c1 vs c2 细胞类型间        -c cell_type --contrasts 'c1,c2'",
"  c1 内部 Disease vs Ctrl ⭐     --subset_column cell_type --subset_value c1 -c condition --contrasts 'Disease,Control'",
"  全局 Disease vs Ctrl           -c condition --contrasts 'Disease,Control'",
"",
"画图脚本(11 个, 独立运行, 默认读 results/*_manifest.txt):",
"  find_marker : dotplot heatmap violin | deg : volcano heatmap violin",
"  pseudobulk  : volcano heatmap        | hetero : dimplot dotplot heatmap",
"",
"示例:",
"  # 注释用 marker(所有组一次算)",
"  Rscript sc_006.find_marker.DEG.Heter.R -i obj.rds -c cell_type --do_find_marker",
"",
"  # 某类型内 条件 DEG + pseudobulk (批量两对比)",
"  Rscript sc_006.find_marker.DEG.Heter.R -i obj.rds --subset_column cell_type --subset_value Fibroblast -c condition --do_deg --do_pseudobulk --sample_column orig.ident --contrasts 'Disease,Control;Disease,Other'",
"",
"  # 异质性下钻某群",
"  Rscript sc_006.find_marker.DEG.Heter.R -i obj.rds -c cell_type --do_hetero --compare_group Fibroblast",
"",
"  # 四个一起算",
"  Rscript sc_006.find_marker.DEG.Heter.R -i obj.rds -c cell_type --do_find_marker --do_deg --do_hetero --compare_group Fibroblast --contrasts 'c1,c2'",
sep="\n")

option_list <- list(
  make_option(c("-i","--input"),  type="character", default=NULL, help="输入 Seurat .rds", metavar="FILE"),
  make_option(c("-o","--outdir"), type="character", default="results", help="输出目录 [%default]"),
  make_option("--prefix",         type="character", default="sc006", help="输出前缀 [%default]"),
  make_option("--plotdir",        type="character", default="./scripts",
              help="画图脚本所在目录(结尾提示用; 默认假设与本脚本同目录) [%default]"),
  make_option(c("-c","--ident"),  type="character", default="cell_type", help="用作 Idents 的 metadata 列 [%default]"),
  make_option("--assay",          type="character", default=NULL, help="assay [默认: DefaultAssay]"),
  make_option("--subset_column",  type="character", default=NULL, help="预筛列(把分析限定在某群内, 如 cell_type)"),
  make_option("--subset_value",   type="character", default=NULL, help="预筛保留值(逗号分隔, 如 'Fibroblast')"),
  # 开关
  make_option("--do_find_marker", action="store_true", default=FALSE, help="计算 one-vs-all 共识 marker"),
  make_option("--do_deg",         action="store_true", default=FALSE, help="计算 FindMarkers(A vs B)"),
  make_option("--do_pseudobulk",  action="store_true", default=FALSE, help="计算 pseudobulk DESeq2(需 --sample_column)"),
  make_option("--do_hetero",      action="store_true", default=FALSE, help="亚群下钻(需 --compare_group)"),
  # deg / pseudobulk
  make_option("--contrasts",      type="character", default=NULL, help="对比 'A,B;A,C'(右侧 ALL=one-vs-rest)"),
  make_option("--sample_column",  type="character", default=NULL, help="样本/重复列(pseudobulk 必需)"),
  make_option("--test",           type="character", default="wilcox", help="FindMarkers 检验 [%default]"),
  # hetero
  make_option("--compare_group",  type="character", default=NULL, help="要下钻的组(--do_hetero)"),
  make_option("--dim",            type="integer", default=30, help="hetero PCA 维数 [%default]"),
  make_option("--n_var_features", type="integer", default=2000, help="hetero 高变基因数 [%default]"),
  make_option("--hetero_res",     type="double",  default=0.2, help="hetero 亚聚类分辨率 [%default]"),
  # find_marker 阈值
  make_option("--strict_min_pct", type="double", default=0.5,  help="STRICT min.pct [%default]"),
  make_option("--strict_pct2",    type="double", default=0.05, help="STRICT max pct.2 [%default]"),
  make_option("--sens_min_pct",   type="double", default=0.25, help="SENS min.pct [%default]"),
  make_option("--sens_min_diff_pct", type="double", default=0.1, help="SENS min.diff.pct [%default]"),
  make_option("--sens_logfc",     type="double", default=0.25, help="SENS logfc.threshold [%default]"),
  make_option("--latent",         type="character", default="nFeature_RNA", help="MAST latent.vars [%default]"),
  make_option("--no_mast",        action="store_true", default=FALSE, help="跳过 MAST(final=Wilcox∩ROC)"),
  make_option("--no_roc",         action="store_true", default=FALSE, help="跳过 ROC(候选=Wilcox(padj); ROC 无 presto 极慢时用)"),
  # 通用过滤
  make_option("--padj",  type="double", default=0.05, help="padj 阈值 [%default]"),
  make_option("--logfc", type="double", default=0.25, help="deg logfc 过滤阈值(存 _sig) [%default]")
)

parser <- OptionParser(option_list=option_list,
  description="\n单细胞 marker/DEG/异质性 一体化计算。无参数或无开关时打印 usage。\n",
  epilogue=USAGE)

argv <- commandArgs(trailingOnly=TRUE)
if (length(argv) == 0) { cat(USAGE, "\n"); quit(status=0) }
opt <- parse_args(parser)

if (is.null(opt$input)) { cat("\n[!] 缺少 --input\n\n"); cat(USAGE, "\n"); quit(status=1) }
if (!any(opt$do_find_marker, opt$do_deg, opt$do_pseudobulk, opt$do_hetero)) {
  cat("\n[!] 没有打开任何开关(--do_find_marker/--do_deg/--do_pseudobulk/--do_hetero)\n\n")
  cat(USAGE, "\n"); quit(status=1)
}

# ------------------------------------------------------------------ packages
suppressPackageStartupMessages({ library(Seurat); library(dplyr) })
`%||%` <- function(a,b) if (is.null(a) || (length(a)==1 && is.na(a))) b else a

# ------------------------------------------------------------------ step timer
.T0 <- Sys.time(); .Tlast <- .T0
step <- function(msg) {
  now <- Sys.time()
  cat(sprintf("[%s] +%7.1fs (累计 %7.1fs) %s\n", format(now, "%H:%M:%S"),
              as.numeric(difftime(now, .Tlast, units="secs")),
              as.numeric(difftime(now, .T0,    units="secs")), msg))
  flush.console(); .Tlast <<- now
}

# ------------------------------------------------------------------ load/prep
cat("[load]", opt$input, "\n")
obj <- readRDS(opt$input)
step("readRDS 读入完成")
if (!is.null(opt$assay)) DefaultAssay(obj) <- opt$assay
ASSAY <- DefaultAssay(obj)

# 预筛到某群
if (!is.null(opt$subset_column) && !is.null(opt$subset_value)) {
  if (!opt$subset_column %in% colnames(obj@meta.data)) stop("subset_column 不存在: ", opt$subset_column)
  vals <- trimws(unlist(strsplit(opt$subset_value, ",")))
  cells <- colnames(obj)[obj@meta.data[[opt$subset_column]] %in% vals]
  obj <- obj[, cells]
  cat(sprintf("[subset] %s ∈ {%s} → %d cells\n", opt$subset_column, paste(vals,collapse=","), ncol(obj)))
}
if (!opt$ident %in% colnames(obj@meta.data))
  stop("ident 列 '", opt$ident, "' 不存在。可用: ", paste(colnames(obj@meta.data), collapse=", "))
Idents(obj) <- opt$ident
cat(sprintf("[assay] %s | cells %d | ident '%s' | 组: %s\n",
            ASSAY, ncol(obj), opt$ident, paste(levels(Idents(obj)), collapse=", ")))

# assay 统一处理一次 (find_marker/deg 都用)
if (ASSAY == "SCT") {
  cat("[prep] SCT → PrepSCTFindMarkers\n"); obj <- PrepSCTFindMarkers(obj); step("PrepSCTFindMarkers 完成")
} else {
  ly <- tryCatch(Layers(obj[[ASSAY]]), error=function(e) NULL)
  if (!is.null(ly) && any(grepl("^counts\\.|^data\\.", ly))) { cat("[prep] JoinLayers\n"); obj <- JoinLayers(obj); step("JoinLayers 完成") }
}
dir.create(opt$outdir, recursive=TRUE, showWarnings=FALSE)
op <- function(x) file.path(opt$outdir, paste0(opt$prefix, x))

# ------------------------------------------------------------------ manifest
man <- c(
  sprintf("rds=%s", normalizePath(opt$input, winslash="/", mustWork=FALSE)),
  sprintf("outdir=%s", opt$outdir),
  sprintf("prefix=%s", opt$prefix),
  sprintf("ident=%s", opt$ident),
  sprintf("assay=%s", ASSAY),
  sprintf("subset_column=%s", opt$subset_column %||% ""),
  sprintf("subset_value=%s", opt$subset_value %||% ""),
  sprintf("contrasts=%s", opt$contrasts %||% ""),
  sprintf("sample_column=%s", opt$sample_column %||% ""),
  sprintf("compare_group=%s", opt$compare_group %||% "")
)

# helper: raw counts (prefer RNA assay)
get_counts <- function(o) {
  a <- if ("RNA" %in% Assays(o)) "RNA" else DefaultAssay(o)
  tryCatch(LayerData(o, assay=a, layer="counts"),
           error=function(e) GetAssayData(o, assay=a, slot="counts"))
}
parse_contrasts <- function(s) {
  if (is.null(s)) return(list())
  lapply(strsplit(s, ";")[[1]], function(x) trimws(strsplit(x, ",")[[1]]))
}

# ==============================================================================
# 1) FIND MARKER  (one-vs-all 共识, 所有组)
# ==============================================================================
if (opt$do_find_marker) tryCatch({
  cat("\n########## [do_find_marker] ##########\n")
  has_mast <- (!opt$no_mast) && requireNamespace("MAST", quietly=TRUE)
  if (!opt$no_mast && !has_mast) cat("[warn] MAST 未装 → final=Wilcox∩ROC\n")
  latent <- trimws(unlist(strsplit(opt$latent, ",")))

  run_tier <- function(tier, min.pct, min.diff.pct, logfc, pct2_max) {
    cat(sprintf("--- %s ---\n", tier)); step(sprintf("[%s] 开始", tier))
    wil <- FindAllMarkers(obj, only.pos=FALSE, test.use="wilcox",
                          min.pct=min.pct, min.diff.pct=min.diff.pct,
                          logfc.threshold=logfc, verbose=FALSE)
    if (!is.null(pct2_max)) wil <- wil %>% filter(pct.2 < pct2_max)
    wil <- wil %>% arrange(cluster, desc(avg_log2FC))
    step(sprintf("[%s] Wilcox 完成 (%d 行)", tier, nrow(wil)))
    wsig <- wil %>% filter(p_val_adj < opt$padj)
    cand <- unique(wsig$gene)                       # Wilcox 显著基因 → ROC 只需算这些
    if (opt$no_roc || length(cand) == 0) {
      ov <- wsig %>% select(cluster, gene) %>% mutate(myAUC = NA_real_)
      step(sprintf("[%s] 跳过 ROC → 候选=Wilcox(padj) %d 基因", tier, length(cand)))
    } else {
      # ROC 只在 Wilcox 显著基因上算(features=cand): 交集结果不变, 但从几千基因降到几十 → 秒级
      roc <- FindAllMarkers(obj, only.pos=FALSE, test.use="roc", features=cand,
                            min.pct=min.pct, min.diff.pct=min.diff.pct, verbose=FALSE)
      if (!is.null(pct2_max)) roc <- roc %>% filter(pct.2 < pct2_max)
      step(sprintf("[%s] ROC 完成 (只算 %d 候选基因, %d 行)", tier, length(cand), nrow(roc)))
      ov <- inner_join(wsig %>% select(cluster, gene),
                       roc %>% select(cluster, gene, myAUC), by=c("cluster","gene"))
    }
    genes <- unique(ov$gene)
    cat(sprintf("    候选(Wilcox%s): %d 基因\n", if (opt$no_roc) "" else "∩ROC", length(genes)))
    if (has_mast && length(genes) > 0) {
      mast <- FindAllMarkers(obj, only.pos=FALSE, features=genes, test.use="MAST",
                             latent.vars=latent, min.pct=0, logfc.threshold=0, verbose=FALSE)
      step(sprintf("[%s] MAST 完成 (候选 %d 基因)", tier, length(genes)))
      final <- mast %>% filter(p_val_adj < opt$padj) %>%
        left_join(ov %>% select(cluster,gene,myAUC), by=c("cluster","gene")) %>%
        arrange(cluster, desc(avg_log2FC))
    } else {
      final <- wil %>% inner_join(ov %>% select(cluster,gene,myAUC), by=c("cluster","gene")) %>%
        arrange(cluster, desc(avg_log2FC))
    }
    write.csv(wil,   op(sprintf("_find_marker_%s_wilcox.csv", tier)), row.names=FALSE)
    write.csv(final, op(sprintf("_find_marker_%s_final.csv",  tier)), row.names=FALSE)
    cat(sprintf("    [saved] %s final: %d 基因\n", tier, nrow(final)))
    step(sprintf("[%s] 全部完成", tier))
  }
  run_tier("strict",    opt$strict_min_pct, 0, 0, opt$strict_pct2)
  run_tier("sensitive", opt$sens_min_pct, opt$sens_min_diff_pct, opt$sens_logfc, NULL)
}, error=function(e) cat("[find_marker ERROR]", conditionMessage(e), "\n"))

# ==============================================================================
# 2) DEG  (FindMarkers A vs B, 批量)
# ==============================================================================
if (opt$do_deg) tryCatch({
  cat("\n########## [do_deg] ##########\n")
  ctr <- parse_contrasts(opt$contrasts)
  if (!length(ctr)) stop("--do_deg 需 --contrasts 'A,B;...'")
  grp <- levels(Idents(obj))
  for (cc in ctr) {
    A <- cc[1]; B <- cc[2]
    id2 <- if (toupper(B)=="ALL") NULL else B
    if (!A %in% grp) { cat("[skip]", A, "不在", opt$ident, "\n"); next }
    if (!is.null(id2) && !B %in% grp) { cat("[skip]", B, "不在", opt$ident, "\n"); next }
    cat(sprintf("--- DEG %s vs %s ---\n", A, B))
    res <- FindMarkers(obj, ident.1=A, ident.2=id2, test.use=opt$test,
                       logfc.threshold=0, min.pct=0)
    res$gene <- rownames(res); res$group1 <- A; res$group2 <- B
    res <- res %>% arrange(desc(avg_log2FC))
    write.csv(res, op(sprintf("_DEG_%s_vs_%s.csv", A, B)), row.names=FALSE)
    sig <- res %>% filter(p_val_adj < opt$padj, abs(avg_log2FC) >= opt$logfc)
    write.csv(sig, op(sprintf("_DEG_%s_vs_%s_sig.csv", A, B)), row.names=FALSE)
    cat(sprintf("    全基因 %d | 显著 %d (↑%d ↓%d)\n", nrow(res), nrow(sig),
                sum(sig$avg_log2FC>0), sum(sig$avg_log2FC<0)))
  }
}, error=function(e) cat("[deg ERROR]", conditionMessage(e), "\n"))

# ==============================================================================
# 3) PSEUDOBULK  (DESeq2, 样本聚合)
# ==============================================================================
if (opt$do_pseudobulk) tryCatch({
  cat("\n########## [do_pseudobulk] ##########\n")
  if (is.null(opt$sample_column)) stop("--do_pseudobulk 需 --sample_column")
  if (!opt$sample_column %in% colnames(obj@meta.data)) stop("sample_column 不存在")
  if (!requireNamespace("DESeq2", quietly=TRUE)) stop("需要 DESeq2: BiocManager::install('DESeq2')")
  ctr <- parse_contrasts(opt$contrasts); if (!length(ctr)) stop("--do_pseudobulk 需 --contrasts")
  cnts <- get_counts(obj)
  md <- obj@meta.data
  for (cc in ctr) {
    A <- cc[1]; B <- cc[2]
    if (toupper(B)=="ALL") { cat("[skip] pseudobulk 不支持 ALL:", A, "\n"); next }
    cat(sprintf("--- pseudobulk %s vs %s ---\n", A, B))
    pb <- list(); meta <- list()
    for (g in c(A,B)) {
      cg <- rownames(md)[md[[opt$ident]] == g]
      for (s in unique(md[cg, opt$sample_column])) {
        cs <- cg[md[cg, opt$sample_column] == s]
        if (length(cs) >= 10) {
          key <- paste0(g, "__", s)
          pb[[key]]   <- Matrix::rowSums(cnts[, cs, drop=FALSE])
          meta[[key]] <- data.frame(sample=s, group=g, n=length(cs))
        }
      }
    }
    if (length(pb) < 4) { cat("[skip] 样本太少(每组需≥2)\n"); next }
    mat <- do.call(cbind, pb); coldata <- do.call(rbind, meta); rownames(coldata) <- colnames(mat)
    if (min(table(coldata$group)) < 2) { cat("[skip] 某组样本<2\n"); next }
    suppressPackageStartupMessages(library(DESeq2))
    dds <- DESeqDataSetFromMatrix(round(mat), coldata, design=~group)
    dds <- dds[rowSums(counts(dds) >= 10) >= 2, ]
    dds$group <- relevel(factor(dds$group), ref=B)
    dds <- DESeq(dds, quiet=TRUE)
    res <- as.data.frame(results(dds, contrast=c("group", A, B)))
    res$gene <- rownames(res); res$group1 <- A; res$group2 <- B
    res <- res[!is.na(res$padj), ]; res <- res[order(res$padj), ]
    write.csv(res, op(sprintf("_pseudobulk_%s_vs_%s.csv", A, B)), row.names=FALSE)
    nc <- as.data.frame(counts(dds, normalized=TRUE)); nc$gene <- rownames(nc)
    write.csv(nc, op(sprintf("_pseudobulk_%s_vs_%s_normcounts.csv", A, B)), row.names=FALSE)
    write.csv(coldata, op(sprintf("_pseudobulk_%s_vs_%s_coldata.csv", A, B)), row.names=TRUE)
    cat(sprintf("    样本 %d | 显著(padj<0.05) %d\n", ncol(mat), sum(res$padj < 0.05)))
  }
}, error=function(e) cat("[pseudobulk ERROR]", conditionMessage(e), "\n"))

# ==============================================================================
# 4) HETERO  (亚群下钻)
# ==============================================================================
if (opt$do_hetero) tryCatch({
  cat("\n########## [do_hetero] ##########\n")
  if (is.null(opt$compare_group)) stop("--do_hetero 需 --compare_group")
  X <- opt$compare_group
  if (!X %in% levels(Idents(obj))) stop("compare_group 不在 ident: ", X)
  sub <- subset(obj, idents = X)
  cat(sprintf("--- hetero %s: %d cells ---\n", X, ncol(sub)))
  if (ncol(sub) < 50) cat("[warn] 细胞<50, 结果不稳\n")
  sub <- FindVariableFeatures(sub, nfeatures=opt$n_var_features, verbose=FALSE)
  sub <- ScaleData(sub, verbose=FALSE)
  sub <- RunPCA(sub, npcs=opt$dim, verbose=FALSE)
  pct <- sub@reductions$pca@stdev / sum(sub@reductions$pca@stdev) * 100
  sig_pcs <- which(pct > 5); if (!length(sig_pcs)) sig_pcs <- 1:min(10, opt$dim)
  cat("[hetero] 显著 PC(>5%):", paste(sig_pcs, collapse=","), "\n")
  emb <- sub@reductions$pca@cell.embeddings[, sig_pcs, drop=FALSE]
  sub$dominant_PC <- paste0("PC", sig_pcs[apply(abs(emb), 1, which.max)])
  sub <- FindNeighbors(sub, dims=sig_pcs, verbose=FALSE)
  sub <- FindClusters(sub, resolution=opt$hetero_res, verbose=FALSE)
  sub$subcluster <- Idents(sub)
  udims <- if (length(sig_pcs) >= 2) sig_pcs else seq_len(min(10, ncol(sub@reductions$pca@cell.embeddings)))
  if (ncol(sub) >= 50) sub <- tryCatch(RunUMAP(sub, dims=udims, verbose=FALSE),
                                        error=function(e){cat("[warn] UMAP 失败(跳过):", conditionMessage(e), "\n"); sub})
  if (DefaultAssay(sub) == "SCT")
    sub <- tryCatch(PrepSCTFindMarkers(sub), error=function(e){cat("[warn] 子集 PrepSCT 失败\n"); sub})
  mk <- FindAllMarkers(sub, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25, verbose=FALSE)
  mk <- mk %>% arrange(cluster, desc(avg_log2FC))
  Xs <- gsub("[^A-Za-z0-9]+", "_", X)
  write.csv(mk, op(sprintf("_hetero_%s_subcluster_markers.csv", Xs)), row.names=FALSE)
  saveRDS(sub, op(sprintf("_hetero_%s.rds", Xs)))
  man <<- c(man, sprintf("hetero_group=%s", X), sprintf("hetero_rds=%s", op(sprintf("_hetero_%s.rds", Xs))))
  cat(sprintf("    亚聚类 %d 个 | markers %d | 存 rds\n", length(unique(sub$subcluster)), nrow(mk)))
}, error=function(e) cat("[hetero ERROR]", conditionMessage(e), "\n"))

# ------------------------------------------------------------------ write manifest + done
writeLines(man, op("_manifest.txt"))
cat("\n[done] 输出目录:", opt$outdir, "| manifest:", op("_manifest.txt"), "\n")

# ---- 画图提示: 列出全部 11 个脚本(各自独立), 用户自选要跑哪些 ----
PD <- sub("/+$", "", opt$plotdir); OD <- opt$outdir
DOTS <- c("find_marker.dotplot", "hetero.dotplot")  # 三风格引擎 (需 sc_006.plot.dotplot_engine.R 同目录)
line <- function(s){
  hint <- if (s %in% DOTS) "   # 三风格 --style C(默认 蓝+facet留白)|B(spectral 组条)|A(红 分隔线); --transpose 转置" else ""
  cat(sprintf("Rscript %s/sc_006.plot.%s.R --outdir %s%s\n", PD, s, OD, hint))
}
tag  <- function(on) if (on) "   <<< 本次已算, 可直接画" else ""
cat("\n========== 画图脚本 (11 个, 独立运行, 默认读上面的 manifest; 按需选择) ==========\n")
cat(sprintf("# 目录: %s   | 每条出 PNG+PDF 到 %s ; 改图细节见各脚本顶部 CONFIG 块或 --help\n", PD, OD))
cat(sprintf("\n# --- find_marker  [注释]%s ---\n", tag(opt$do_find_marker)))
for (s in c("find_marker.dotplot","find_marker.heatmap","find_marker.violin")) line(s)
cat(sprintf("\n# --- deg  [组间 A vs B]%s ---\n", tag(opt$do_deg)))
for (s in c("deg.volcano","deg.heatmap","deg.violin")) line(s)
cat(sprintf("\n# --- pseudobulk  [稳健组间]%s ---\n", tag(opt$do_pseudobulk)))
for (s in c("pseudobulk.volcano","pseudobulk.heatmap")) line(s)
cat(sprintf("\n# --- hetero  [亚结构]%s ---\n", tag(opt$do_hetero)))
for (s in c("hetero.dimplot","hetero.dotplot","hetero.heatmap")) line(s)
cat("\n# 打包: 把想要的行写进 make_plots.sh 一次跑齐\n")
cat("# 详细开关/示例: 每个脚本加 --help; dotplot 三风格与自定义(--genes/--features/--cols/--group_colors)见\n")
cat("#              skill cell-type-annotation → references/scripts_reference.md §4 (§4d = dotplot 三风格引擎)\n")
