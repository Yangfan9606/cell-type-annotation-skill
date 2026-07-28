#!/usr/bin/env Rscript
# sc_006.plot.hetero.dimplot.R —— 亚群 UMAP (亚聚类 + dominant_PC 两版)
#   读 manifest 的 hetero_rds (或 outdir 里最新 *_hetero_*.rds)。 Rscript sc_006.plot.hetero.dimplot.R
# ============================== CONFIG ==============================
GROUPBYS <- c("subcluster", "dominant_PC")   # 各画一版; 有哪个画哪个
# ===================================================================
suppressPackageStartupMessages(library(optparse))   # 只加载轻库; 重库延后, 加快 --help/usage
EX <- "
示例:
  Rscript sc_006.plot.hetero.dimplot.R --outdir results      # 读 manifest 的 hetero_rds(或最新 *_hetero_*.rds)
  Rscript sc_006.plot.hetero.dimplot.R --rds myrun_hetero_c1.rds
按 subcluster 和 dominant_PC 各出一版(有哪个画哪个); 无 UMAP 时退回 PCA 前两维。
改要画的分组: 脚本顶部 CONFIG 的 GROUPBYS。
"
parser <- OptionParser(option_list=list(
  make_option("--manifest", default=NULL, help="manifest 路径(默认取 outdir 最新)"),
  make_option("--rds", default=NULL, help="hetero rds(默认取 manifest 的 hetero_rds)"),
  make_option("--outfile", default=NULL, help="输出 PNG(PDF 自动同名)"),
  make_option("--outdir", default="results"), make_option("--prefix", default=NULL)),
  description="sc_006.plot.hetero.dimplot.R — 亚群 UMAP(subcluster + dominant_PC 两版)", epilogue=EX)
if (length(commandArgs(TRUE))==0) { print_help(parser); quit(status=0) }
opt <- parse_args(parser)
suppressPackageStartupMessages({library(Seurat); library(ggplot2); library(patchwork)})
`%||%` <- function(a,b) if (is.null(a)||(length(a)==1&&is.na(a))||(is.character(a)&&length(a)==1&&a=="")) b else a
read_manifest <- function(mp, outdir){
  if (is.null(mp)){ cc<-list.files(outdir,"_manifest\\.txt$",full.names=TRUE); if(length(cc)) mp<-cc[which.max(file.mtime(cc))] }
  m<-list(); if(!is.null(mp)&&file.exists(mp)) for(ln in readLines(mp)) if(grepl("=",ln)){k<-sub("=.*","",ln);m[[k]]<-sub("^[^=]*=","",ln)}; m }
man <- read_manifest(opt$manifest, opt$outdir)
outdir<-opt$outdir%||%man$outdir%||%"results"; prefix<-opt$prefix%||%man$prefix%||%"sc006"
rds <- opt$rds %||% man$hetero_rds %||% {
  ff <- list.files(outdir, pattern="_hetero_.*\\.rds$", full.names=TRUE); if(length(ff)) ff[which.max(file.mtime(ff))] else NULL }
if (is.null(rds) || !file.exists(rds)) stop("找不到 hetero rds (先跑 --do_hetero)")
cat("[read]", rds, "\n"); so <- readRDS(rds)
redu <- if ("umap" %in% names(so@reductions)) "umap" else "pca"
if (redu == "pca") cat("[warn] 无 UMAP, 改用 PCA 前两维\n")
gbs <- GROUPBYS[GROUPBYS %in% colnames(so@meta.data)]
pl <- lapply(gbs, function(g) DimPlot(so, reduction=redu, group.by=g, label=TRUE, label.size=5) +
             labs(title=g) +
             theme(plot.title=element_text(face="bold", hjust=.5, size=15),
                   axis.title=element_text(face="bold", size=13), axis.text=element_text(face="bold", size=10),
                   legend.text=element_text(size=11, face="bold")))
p <- wrap_plots(pl, ncol=length(pl))
out <- opt$outfile %||% file.path(outdir, paste0(prefix, "_plot_hetero_dimplot.png"))
ggsave(out, p, width=6*length(pl), height=5.5, dpi=300); ggsave(sub("\\.png$",".pdf",out), p, width=6*length(pl), height=5.5)
cat("[saved]", out, "\n")
