#!/usr/bin/env Rscript
# sc_006.plot.hetero.heatmap.R —— 亚聚类 marker 热图
#   读 hetero rds + *_hetero_*_subcluster_markers.csv。 Rscript sc_006.plot.hetero.heatmap.R
# ============================== CONFIG ==============================
TOPN <- 10; CELLS_PER_GRP <- 100
# ===================================================================
suppressPackageStartupMessages(library(optparse))   # 只加载轻库; 重库延后, 加快 --help/usage
EX <- "
示例:
  Rscript sc_006.plot.hetero.heatmap.R --outdir results      # 读 hetero rds + *_hetero_*_subcluster_markers.csv
  Rscript sc_006.plot.hetero.heatmap.R --rds myrun_hetero_c1.rds --input myrun_hetero_c1_subcluster_markers.csv
细调 TOPN / 每组抽样细胞 CELLS_PER_GRP: 改脚本顶部 CONFIG 块。
"
parser <- OptionParser(option_list=list(
  make_option("--manifest", default=NULL, help="manifest 路径(默认取 outdir 最新)"),
  make_option("--rds", default=NULL, help="hetero rds(默认取 manifest 的 hetero_rds)"),
  make_option("--input", default=NULL, help="subcluster markers CSV(默认取最新)"),
  make_option("--outfile", default=NULL, help="输出 PNG(PDF 自动同名)"),
  make_option("--outdir", default="results"), make_option("--prefix", default=NULL)),
  description="sc_006.plot.hetero.heatmap.R — 亚聚类 marker 热图", epilogue=EX)
if (length(commandArgs(TRUE))==0) { print_help(parser); quit(status=0) }
opt <- parse_args(parser)
suppressPackageStartupMessages({library(Seurat); library(ggplot2); library(dplyr)})
`%||%` <- function(a,b) if (is.null(a)||(length(a)==1&&is.na(a))||(is.character(a)&&length(a)==1&&a=="")) b else a
read_manifest <- function(mp, outdir){
  if (is.null(mp)){ cc<-list.files(outdir,"_manifest\\.txt$",full.names=TRUE); if(length(cc)) mp<-cc[which.max(file.mtime(cc))] }
  m<-list(); if(!is.null(mp)&&file.exists(mp)) for(ln in readLines(mp)) if(grepl("=",ln)){k<-sub("=.*","",ln);m[[k]]<-sub("^[^=]*=","",ln)}; m }
man <- read_manifest(opt$manifest, opt$outdir)
outdir<-opt$outdir%||%man$outdir%||%"results"; prefix<-opt$prefix%||%man$prefix%||%"sc006"
rds <- opt$rds %||% man$hetero_rds %||% {
  ff <- list.files(outdir, pattern="_hetero_.*\\.rds$", full.names=TRUE); if(length(ff)) ff[which.max(file.mtime(ff))] else NULL }
csv <- opt$input %||% {
  ff <- list.files(outdir, pattern="_hetero_.*_subcluster_markers\\.csv$", full.names=TRUE); if(length(ff)) ff[which.max(file.mtime(ff))] else NULL }
if (is.null(rds) || !file.exists(rds)) stop("找不到 hetero rds")
if (is.null(csv) || !file.exists(csv)) stop("找不到 subcluster markers csv")
so <- readRDS(rds); Idents(so) <- "subcluster"
mk <- read.csv(csv)
genes <- mk %>% group_by(cluster) %>% slice_max(avg_log2FC, n=TOPN, with_ties=FALSE) %>% pull(gene) %>% unique()
genes <- genes[genes %in% rownames(so)]
set.seed(42)
md <- so@meta.data; md$.cell <- rownames(md)
keep <- md %>% group_by(subcluster) %>% slice_sample(n=CELLS_PER_GRP) %>% pull(.cell)
sh <- so[, keep]; Idents(sh) <- "subcluster"
sh <- if (DefaultAssay(sh)=="SCT") tryCatch(GetResidual(sh, features=genes, verbose=FALSE), error=function(e) sh) else ScaleData(sh, features=genes, verbose=FALSE)
p <- DoHeatmap(sh, features=genes, group.by="subcluster", size=4.5, angle=0) +
  labs(title=paste0("Hetero subcluster heatmap (top ", TOPN, ")")) +
  theme(plot.title=element_text(face="bold", hjust=.5, size=15),
        axis.text.y=element_text(size=8, face="bold"),
        legend.title=element_text(face="bold", size=13), legend.text=element_text(size=11))
out <- opt$outfile %||% file.path(outdir, paste0(prefix, "_plot_hetero_heatmap.png"))
h <- max(8, length(genes)*0.16)
ggsave(out, p, width=12, height=h, dpi=300, limitsize=FALSE); ggsave(sub("\\.png$",".pdf",out), p, width=12, height=h, limitsize=FALSE)
cat("[saved]", out, "\n")
