#!/usr/bin/env Rscript
# sc_006.plot.pseudobulk.heatmap.R —— pseudobulk top DEG × 样本 热图 (样本级)
#   读 *_pseudobulk_*_vs_*.csv + _normcounts.csv (+ _coldata.csv 做注释)。批量每对比一张
#   Rscript sc_006.plot.pseudobulk.heatmap.R
# ============================== CONFIG ==============================
TOPN <- 40; PADJ <- 0.05
# ===================================================================
suppressPackageStartupMessages(library(optparse))   # 只加载轻库; 重库延后, 加快 --help/usage
EX <- "
示例:
  Rscript sc_006.plot.pseudobulk.heatmap.R --outdir results   # 对每个 *_pseudobulk_*_vs_*.csv 各出一张
  Rscript sc_006.plot.pseudobulk.heatmap.R --input myrun_pseudobulk_Disease_vs_Control.csv
读同名 _normcounts.csv(+_coldata.csv 做样本注释)。细调 TOPN/PADJ: 改脚本顶部 CONFIG 块。
"
parser <- OptionParser(option_list=list(
  make_option("--manifest", default=NULL, help="manifest 路径(默认取 outdir 最新)"),
  make_option("--input", default=NULL, help="单个 pseudobulk CSV(默认自动扫全部)"),
  make_option("--outdir", default="results"), make_option("--prefix", default=NULL)),
  description="sc_006.plot.pseudobulk.heatmap.R — pseudobulk top DEG × 样本 热图(样本级)", epilogue=EX)
if (length(commandArgs(TRUE))==0) { print_help(parser); quit(status=0) }
opt <- parse_args(parser)
suppressPackageStartupMessages({library(pheatmap); library(dplyr)})
`%||%` <- function(a,b) if (is.null(a)||(length(a)==1&&is.na(a))||(is.character(a)&&length(a)==1&&a=="")) b else a
read_manifest <- function(mp, outdir){
  if (is.null(mp)){ cc<-list.files(outdir,"_manifest\\.txt$",full.names=TRUE); if(length(cc)) mp<-cc[which.max(file.mtime(cc))] }
  m<-list(); if(!is.null(mp)&&file.exists(mp)) for(ln in readLines(mp)) if(grepl("=",ln)){k<-sub("=.*","",ln);m[[k]]<-sub("^[^=]*=","",ln)}; m }
man <- read_manifest(opt$manifest, opt$outdir)
outdir<-opt$outdir%||%man$outdir%||%"results"; prefix<-opt$prefix%||%man$prefix%||%"sc006"
files <- if (!is.null(opt$input)) opt$input else {
  ff <- list.files(outdir, pattern=paste0("^",prefix,"_pseudobulk_.*_vs_.*\\.csv$"), full.names=TRUE)
  ff[!grepl("_(normcounts|coldata)\\.csv$", ff)] }
if (!length(files)) stop("找不到 *_pseudobulk_*_vs_*.csv")

for (f in files) {
  base <- sub("\\.csv$", "", f)
  ncf <- paste0(base, "_normcounts.csv"); cdf <- paste0(base, "_coldata.csv")
  if (!file.exists(ncf)) { cat("[skip] 缺 normcounts:", ncf, "\n"); next }
  res <- read.csv(f); A<-if("group1"%in%names(res)) res$group1[1] else "A"; B<-if("group2"%in%names(res)) res$group2[1] else "B"
  nc <- read.csv(ncf, check.names=FALSE); rownames(nc) <- nc$gene
  genes <- res %>% filter(!is.na(padj), padj<PADJ) %>% arrange(padj) %>% head(TOPN) %>% pull(gene)
  genes <- intersect(genes, rownames(nc))
  if (length(genes) < 2) { cat("[skip] 显著基因<2:", f, "\n"); next }
  mat <- as.matrix(nc[genes, setdiff(colnames(nc), "gene"), drop=FALSE])
  mat <- log2(mat + 1)
  ann <- NA
  if (file.exists(cdf)) { cd <- read.csv(cdf, row.names=1); ann <- cd[colnames(mat), "group", drop=FALSE] }
  out <- file.path(outdir, sprintf("%s_plot_pseudobulk_heatmap_%s_vs_%s.png", prefix, A, B))
  pheatmap(mat, scale="row", annotation_col=ann, show_colnames=TRUE,
           fontsize=12, fontsize_row=9, fontsize_col=11,   # 放大(pheatmap 无粗体选项; main 默认加粗)
           main=paste0("Pseudobulk top", length(genes), " DEG (log2 norm): ", A, " vs ", B),
           filename=out, width=max(6, ncol(mat)*0.5+3), height=max(6, length(genes)*0.15+2))
  cat("[saved]", out, "\n")
}
