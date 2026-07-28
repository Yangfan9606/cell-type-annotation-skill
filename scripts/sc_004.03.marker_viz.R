#!/usr/bin/env Rscript
# ============================================================================
# sc_004.03.marker_viz.R — DotPlot + Violin 可视化, 辅助 cluster 注释定标签
# ----------------------------------------------------------------------------
# 用法 (key=value, 顺序无关; 只有 rds 必填):
#   Rscript sc_004.03.marker_viz.R rds=<path> [markers=cell_marker.v6.txt]
#       [clu=SCT_snn_res.0.5] [assay=SCT] [genes=DCX,STMN2,SLC17A7,PDGFRA]
#       [clusters=13,10,6] [outdir=...] [prefix=...]
# 产出:
#   <prefix>_dotplot.pdf       marker panel × cluster 分组 DotPlot (一眼给每簇定型)
#   <prefix>_qc_violin.pdf     nFeature/nCount/%MT/%ribo 跨簇小提琴 (暴露低质量簇, 如 cl13)
#   <prefix>_genes_violin.pdf  若给 genes= , 指定基因跨簇小提琴 (查具体竞争, 如 Glut vs OPC)
# 说明: 可视化用 SCT(或 RNA) 的 data 层做相对比较; markers 文件格式同 sc_004 (tab + 逗号, # 注释)
# ============================================================================
suppressPackageStartupMessages({ library(Seurat); library(ggplot2)
  if (requireNamespace("patchwork", quietly = TRUE)) library(patchwork) })

argv <- commandArgs(TRUE); kv <- list()
for (a in argv) if (grepl("=", a)) { k <- sub("=.*","",a); v <- sub("^[^=]*=","",a); kv[[k]] <- v }
opt <- function(n, d = NA) if (!is.null(kv[[n]]) && nzchar(kv[[n]])) kv[[n]] else d

RDS <- opt("rds"); if (is.na(RDS)) stop("缺 rds=")
if (!file.exists(RDS)) stop(sprintf("找不到文件: %s", RDS))
obj <- readRDS(RDS)

CLU <- opt("clu")
if (is.na(CLU)) CLU <- if ("SCT_snn_res.0.5" %in% colnames(obj@meta.data)) "SCT_snn_res.0.5" else
                       if ("seurat_clusters" %in% colnames(obj@meta.data)) "seurat_clusters" else
                       stop("找不到聚类列, 请用 clu= 指定")
ASSAY <- opt("assay"); if (is.na(ASSAY)) ASSAY <- if ("SCT" %in% Assays(obj)) "SCT" else DefaultAssay(obj)
DefaultAssay(obj) <- ASSAY

# 簇按数值排序
lv <- unique(as.character(obj@meta.data[[CLU]]))
num <- suppressWarnings(as.numeric(lv)); lv <- if (all(!is.na(num))) lv[order(num)] else sort(lv)
Idents(obj) <- factor(as.character(obj@meta.data[[CLU]]), levels = lv)

outdir <- opt("outdir", dirname(RDS)); prefix <- opt("prefix", sub("\\.rds$","",basename(RDS)))
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
op <- function(s) file.path(outdir, paste0(prefix, s))
cat(sprintf(">>> %s | clu=%s | assay=%s | %d clusters\n", RDS, CLU, ASSAY, length(lv)))

# ---- 文本汇总收集器 (供无图直接判读: 读 *_summary.txt) ----
TXT <- c(sprintf("# %s | clu=%s | assay=%s | %d clusters", basename(RDS), CLU, ASSAY, length(lv)))
addtxt <- function(...) TXT <<- c(TXT, ...)
fmt_df <- function(d, dig = 2) {                 # 数值列四舍五入后转可读文本 (簇为列, 不换行)
  d <- as.data.frame(d, check.names = FALSE)
  num <- vapply(d, is.numeric, logical(1)); d[num] <- lapply(d[num], round, dig)
  o <- options(width = 10000); on.exit(options(o)); capture.output(print(d))
}

# ---- QC 列 (不覆盖已有 percent.mt) ----
qc_assay <- if ("RNA" %in% Assays(obj)) "RNA" else ASSAY
if (!"percent.mt" %in% colnames(obj@meta.data))
  obj$percent.mt   <- tryCatch(PercentageFeatureSet(obj, pattern = "^MT-",       assay = qc_assay), error = function(e) NA_real_)
obj$percent.ribo   <- tryCatch(PercentageFeatureSet(obj, pattern = "^(RPL|RPS)", assay = qc_assay), error = function(e) NA_real_)

# QC 各簇中位数 -> TXT
.cl     <- factor(as.character(obj@meta.data[[CLU]]), levels = lv)
.qccols <- intersect(c("nFeature_RNA","nCount_RNA","percent.mt","percent.ribo"), colnames(obj@meta.data))
.qcm    <- sapply(.qccols, function(cc) tapply(obj@meta.data[[cc]], .cl, median, na.rm = TRUE))
.qcm    <- cbind(n = as.integer(table(.cl)), .qcm)
addtxt("", "================ [QC] 各簇中位数 ================", fmt_df(.qcm))

# ---- marker 文件 -> named list ----
read_markers <- function(f) {
  ln <- readLines(f); ln <- ln[!grepl("^#|^\\s*$", ln)]; ml <- list()
  for (l in ln) { p <- strsplit(l, "\t")[[1]]
    if (length(p) >= 2) { g <- trimws(unlist(strsplit(p[2], ","))); g <- g[nchar(g) > 0]; ml[[trimws(p[1])]] <- g } }
  ml
}
MK <- opt("markers")
if (!is.na(MK) && file.exists(MK)) {
  ml <- read_markers(MK)
  ml <- lapply(ml, function(g) intersect(g, rownames(obj)))
  ml <- ml[sapply(ml, length) > 0]
  p <- DotPlot(obj, features = ml, cluster.idents = FALSE) + RotatedAxis() +
       theme(axis.text.x = element_text(size = 7),
             strip.text  = element_text(size = 8, angle = 90))
  ggsave(op("_dotplot.pdf"), p,
         width = max(12, length(unlist(ml)) * 0.22), height = max(5, length(lv) * 0.35),
         limitsize = FALSE)
  cat("[dot]", op("_dotplot.pdf"), "\n")
  # DotPlot 背后数值 -> TXT (gene×cluster 的 %阳性 与 标准化均值)
  dd  <- p$data
  g2p <- setNames(rep(names(ml), lengths(ml)), unlist(ml))
  gord<- unlist(ml, use.names = FALSE)
  to_mat <- function(val) {
    m <- tapply(val, list(factor(dd$features.plot, levels = gord), dd$id), function(x) x[1])
    data.frame(panel = g2p[rownames(m)], m, check.names = FALSE)
  }
  addtxt("", "================ [DotPlot] %阳性细胞 pct.exp (gene×cluster) ================",
         fmt_df(to_mat(dd$pct.exp)),
         "", "================ [DotPlot] 标准化均值 avg.exp.scaled (gene×cluster) ================",
         fmt_df(to_mat(dd$avg.exp.scaled)))
} else cat("[skip dotplot] 未给 markers= 或文件不存在\n")

# ---- QC 小提琴 ----
qc <- intersect(c("nFeature_RNA","nCount_RNA","percent.mt","percent.ribo"), colnames(obj@meta.data))
pq <- VlnPlot(obj, features = qc, pt.size = 0, ncol = 2) & theme(axis.text.x = element_text(size = 7))
ggsave(op("_qc_violin.pdf"), pq, width = 12, height = 8); cat("[qc]", op("_qc_violin.pdf"), "\n")

# ---- 指定基因小提琴 ----
G <- opt("genes")
if (!is.na(G)) {
  gg <- intersect(trimws(unlist(strsplit(G, ","))), rownames(obj))
  miss <- setdiff(trimws(unlist(strsplit(G, ","))), rownames(obj))
  if (length(miss)) cat("[note] 不在对象中的基因:", paste(miss, collapse = ", "), "\n")
  if (length(gg)) {
    nc <- min(3, length(gg))
    pv <- VlnPlot(obj, features = gg, pt.size = 0, ncol = nc) & theme(axis.text.x = element_text(size = 7))
    ggsave(op("_genes_violin.pdf"), pv,
           width = 5 * nc, height = 4 * ceiling(length(gg) / nc), limitsize = FALSE)
    cat("[genes]", op("_genes_violin.pdf"), "\n")
    # 指定基因背后数值 -> TXT (gene×cluster 均值 与 %阳性, 用 data 层)
    em    <- as.matrix(GetAssayData(obj, assay = ASSAY, layer = "data")[gg, , drop = FALSE])
    gmean <- t(apply(em, 1, function(r) tapply(r, .cl, mean)))
    gpct  <- t(apply(em, 1, function(r) tapply(r, .cl, function(x) 100 * mean(x > 0))))
    addtxt("", "================ [genes] 均值 mean(data) (gene×cluster) ================", fmt_df(gmean),
           "", "================ [genes] %阳性细胞 (gene×cluster) ================", fmt_df(gpct))
  }
}

writeLines(TXT, op("_summary.txt")); cat("[txt]", op("_summary.txt"), "\n")
cat("[done]\n")
