#!/usr/bin/env Rscript
# sc_006.plot.find_marker.heatmap.R —— top marker × 细胞 热图 (总览)
#   默认读 manifest + *_find_marker_strict_final.csv。 Rscript sc_006.plot.find_marker.heatmap.R
# ============================== CONFIG ==============================
TOPN          <- 10     # 每组 top N marker
CELLS_PER_GRP <- 200    # 每组抽样细胞数(控图幅/速度)
# ===================================================================
suppressPackageStartupMessages(library(optparse))   # 只加载轻库; 重库延后, 加快 --help/usage
EX <- "
示例:
  Rscript sc_006.plot.find_marker.heatmap.R --outdir results              # 读最新 manifest + strict CSV
  Rscript sc_006.plot.find_marker.heatmap.R --input markers.csv --rds obj.rds
  Rscript sc_006.plot.find_marker.heatmap.R --outdir results --prefix myrun
细调 TOPN(每组基因数) / CELLS_PER_GRP(每组抽样细胞): 改脚本顶部 CONFIG 块。
"
parser <- OptionParser(option_list=list(
  make_option("--manifest", default=NULL, help="manifest 路径(默认取 outdir 最新)"),
  make_option("--rds", default=NULL, help="Seurat rds(默认取 manifest)"),
  make_option("--input", default=NULL, help="marker CSV(默认自动挑 strict/sensitive)"),
  make_option("--ident", default=NULL, help="分组列(默认从 CSV 标签自动探测)"),
  make_option("--outfile", default=NULL, help="输出 PNG(PDF 自动同名)"),
  make_option("--outdir", default="results"), make_option("--prefix", default=NULL)),
  description="sc_006.plot.find_marker.heatmap.R — top marker × 细胞 热图(注释总览)", epilogue=EX)
if (length(commandArgs(TRUE))==0) { print_help(parser); quit(status=0) }
opt <- parse_args(parser)
suppressPackageStartupMessages({library(Seurat); library(ggplot2); library(dplyr)})
`%||%` <- function(a,b) if (is.null(a)||(length(a)==1&&is.na(a))||(is.character(a)&&length(a)==1&&a=="")) b else a
read_manifest <- function(mp, outdir){
  if (is.null(mp)){ cc<-list.files(outdir,"_manifest\\.txt$",full.names=TRUE); if(length(cc)) mp<-cc[which.max(file.mtime(cc))] }
  m<-list(); if(!is.null(mp)&&file.exists(mp)) for(ln in readLines(mp)) if(grepl("=",ln)){k<-sub("=.*","",ln);m[[k]]<-sub("^[^=]*=","",ln)}; m }
pick_marker_csv <- function(outdir, prefix){
  cand <- file.path(outdir, paste0(prefix, c("_find_marker_strict_final.csv","_find_marker_sensitive_final.csv")))
  for (f in cand) if (file.exists(f) && nrow(read.csv(f))>0) return(f)
  e <- cand[file.exists(cand)]; if (length(e)) e[1] else cand[1] }
man <- read_manifest(opt$manifest, opt$outdir)
rds<-opt$rds%||%man$rds; ident<-opt$ident%||%man$ident%||%"cell_type"
outdir<-opt$outdir%||%man$outdir%||%"results"; prefix<-opt$prefix%||%man$prefix%||%"sc006"
csv <- opt$input %||% pick_marker_csv(outdir, prefix)
out <- opt$outfile %||% file.path(outdir, paste0(prefix,"_plot_find_marker_heatmap.png"))
stopifnot(file.exists(rds), file.exists(csv))

mk<-read.csv(csv); obj<-readRDS(rds)
detect_ident <- function(obj, want, fb){ if(!length(want)) return(fb); for(cn in colnames(obj@meta.data)){v<-obj@meta.data[[cn]]; if((is.character(v)||is.factor(v))&&all(want %in% unique(as.character(v)))) return(cn)}; fb }
ident <- detect_ident(obj, unique(as.character(mk$cluster)), ident); Idents(obj)<-ident
genes <- mk %>% filter(avg_log2FC>0) %>% group_by(cluster) %>%
  slice_max(avg_log2FC, n=TOPN, with_ties=FALSE) %>% pull(gene) %>% unique()
genes <- genes[genes %in% rownames(obj)]
set.seed(42)
cells <- obj@meta.data; cells$.cell <- rownames(cells)
keep <- cells %>% group_by(.data[[ident]]) %>% slice_sample(n=CELLS_PER_GRP) %>% pull(.cell)
so <- obj[, keep]; Idents(so)<-ident
so <- if (DefaultAssay(so)=="SCT") tryCatch(GetResidual(so, features=genes, verbose=FALSE), error=function(e) so) else ScaleData(so, features=genes, verbose=FALSE)
p <- DoHeatmap(so, features=genes, group.by=ident, size=4.5, angle=45) +
  labs(title=paste0("Find-marker heatmap (top ", TOPN, "/group)")) +
  theme(plot.title=element_text(face="bold", hjust=.5, size=15),
        axis.text.y=element_text(size=8, face="bold"),
        legend.title=element_text(face="bold", size=13), legend.text=element_text(size=11))
h <- max(8, length(genes)*0.14)
ggsave(out, p, width=12, height=h, dpi=300, limitsize=FALSE)
ggsave(sub("\\.png$",".pdf",out), p, width=12, height=h, limitsize=FALSE)
cat("[saved]", out, "\n")
