#!/usr/bin/env Rscript
# sc_006.plot.find_marker.violin.R —— 少数关键 marker 跨组小提琴 (确认用)
#   默认读 manifest + *_find_marker_strict_final.csv。 Rscript sc_006.plot.find_marker.violin.R
# ============================== CONFIG ==============================
PER_CLUSTER <- 1     # 每组取几个 top marker 进小提琴
NCOL        <- 3
# ===================================================================
suppressPackageStartupMessages(library(optparse))   # 只加载轻库; 重库延后, 加快 --help/usage
EX <- "
示例:
  Rscript sc_006.plot.find_marker.violin.R --outdir results              # 读最新 manifest + strict CSV
  Rscript sc_006.plot.find_marker.violin.R --input markers.csv --rds obj.rds
细调 PER_CLUSTER(每组取几个 top marker) / NCOL(分面列数): 改脚本顶部 CONFIG 块。
"
parser <- OptionParser(option_list=list(
  make_option("--manifest", default=NULL, help="manifest 路径(默认取 outdir 最新)"),
  make_option("--rds", default=NULL, help="Seurat rds(默认取 manifest)"),
  make_option("--input", default=NULL, help="marker CSV(默认自动挑 strict/sensitive)"),
  make_option("--ident", default=NULL, help="分组列(默认从 CSV 标签自动探测)"),
  make_option("--outfile", default=NULL, help="输出 PNG(PDF 自动同名)"),
  make_option("--outdir", default="results"), make_option("--prefix", default=NULL)),
  description="sc_006.plot.find_marker.violin.R — 少数关键 marker 跨组小提琴(注释确认)", epilogue=EX)
if (length(commandArgs(TRUE))==0) { print_help(parser); quit(status=0) }
opt <- parse_args(parser)
suppressPackageStartupMessages({library(Seurat); library(ggplot2); library(dplyr); library(patchwork)})
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
out <- opt$outfile %||% file.path(outdir, paste0(prefix,"_plot_find_marker_violin.png"))
stopifnot(file.exists(rds), file.exists(csv))

mk<-read.csv(csv); obj<-readRDS(rds)
detect_ident <- function(obj, want, fb){ if(!length(want)) return(fb); for(cn in colnames(obj@meta.data)){v<-obj@meta.data[[cn]]; if((is.character(v)||is.factor(v))&&all(want %in% unique(as.character(v)))) return(cn)}; fb }
ident <- detect_ident(obj, unique(as.character(mk$cluster)), ident); Idents(obj)<-ident
genes <- mk %>% filter(avg_log2FC>0) %>% group_by(cluster) %>%
  slice_max(avg_log2FC, n=PER_CLUSTER, with_ties=FALSE) %>% pull(gene) %>% unique()
genes <- genes[genes %in% rownames(obj)]
cat("[genes]", length(genes), "\n")
# 用 FetchData+geom_violin 手搓, 绕开 Seurat VlnPlot(在 ggplot2>=4.0 上有 S7 崩溃)
df <- FetchData(obj, vars=c(genes, ident), layer="data")
long <- do.call(rbind, lapply(genes, function(g) data.frame(gene=g, expr=df[[g]], grp=as.character(df[[ident]]))))
long$gene <- factor(long$gene, levels=genes)
p <- ggplot(long, aes(grp, expr, fill=grp)) +
  geom_violin(scale="width", trim=TRUE, linewidth=0.2) +
  facet_wrap(~gene, ncol=NCOL, scales="free_y") +
  labs(x=NULL, y="expression") + theme_classic(base_size=13) +
  theme(legend.position="none",
        axis.text.x=element_text(angle=45, hjust=1, size=11, face="bold"),
        axis.text.y=element_text(size=10, face="bold"),
        axis.title.y=element_text(face="bold", size=13),
        strip.text=element_text(face="bold", size=12))
h <- max(6, ceiling(length(genes)/NCOL)*2.5)
ggsave(out, p, width=4*NCOL, height=h, dpi=300, limitsize=FALSE)
ggsave(sub("\\.png$",".pdf",out), p, width=4*NCOL, height=h, limitsize=FALSE)
cat("[saved]", out, "\n")
