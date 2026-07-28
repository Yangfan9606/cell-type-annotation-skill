#!/usr/bin/env Rscript
# sc_006.plot.deg.violin.R —— DEG top 基因跨 A/B 小提琴; 批量每对比一张
#   需 rds + ident。 Rscript sc_006.plot.deg.violin.R
# ============================== CONFIG ==============================
TOPN <- 9; PADJ <- 0.05; LOGFC <- 0.25; NCOL <- 3
# ===================================================================
suppressPackageStartupMessages(library(optparse))   # 只加载轻库; 重库延后, 加快 --help/usage
EX <- "
示例:
  Rscript sc_006.plot.deg.violin.R --outdir results      # 对每个 *_DEG_*_vs_*.csv 各出一张(需 rds)
  Rscript sc_006.plot.deg.violin.R --rds obj.rds --input myrun_DEG_Disease_vs_Control.csv
细调 TOPN/PADJ/LOGFC/NCOL(分面列数): 改脚本顶部 CONFIG 块。
"
parser <- OptionParser(option_list=list(
  make_option("--manifest", default=NULL, help="manifest 路径(默认取 outdir 最新)"),
  make_option("--rds", default=NULL, help="Seurat rds(默认取 manifest)"),
  make_option("--input", default=NULL, help="单个 DEG CSV(默认自动扫全部)"),
  make_option("--ident", default=NULL, help="分组列(默认从 CSV 组名自动探测)"),
  make_option("--outdir", default="results"), make_option("--prefix", default=NULL)),
  description="sc_006.plot.deg.violin.R — DEG top 基因跨 A/B 小提琴, 批量每对比一张", epilogue=EX)
if (length(commandArgs(TRUE))==0) { print_help(parser); quit(status=0) }
opt <- parse_args(parser)
suppressPackageStartupMessages({library(Seurat); library(ggplot2); library(dplyr); library(patchwork)})
`%||%` <- function(a,b) if (is.null(a)||(length(a)==1&&is.na(a))||(is.character(a)&&length(a)==1&&a=="")) b else a
read_manifest <- function(mp, outdir){
  if (is.null(mp)){ cc<-list.files(outdir,"_manifest\\.txt$",full.names=TRUE); if(length(cc)) mp<-cc[which.max(file.mtime(cc))] }
  m<-list(); if(!is.null(mp)&&file.exists(mp)) for(ln in readLines(mp)) if(grepl("=",ln)){k<-sub("=.*","",ln);m[[k]]<-sub("^[^=]*=","",ln)}; m }
man <- read_manifest(opt$manifest, opt$outdir)
rds<-opt$rds%||%man$rds; ident<-opt$ident%||%man$ident%||%"cell_type"
outdir<-opt$outdir%||%man$outdir%||%"results"; prefix<-opt$prefix%||%man$prefix%||%"sc006"
stopifnot(file.exists(rds)); obj<-readRDS(rds)
detect_ident <- function(obj, want, fb){ if(!length(want)) return(fb); for(cn in colnames(obj@meta.data)){v<-obj@meta.data[[cn]]; if((is.character(v)||is.factor(v))&&all(want %in% unique(as.character(v)))) return(cn)}; fb }
files <- if (!is.null(opt$input)) opt$input else {
  ff <- list.files(outdir, pattern=paste0("^",prefix,"_DEG_.*_vs_.*\\.csv$"), full.names=TRUE)
  ff[!grepl("_sig\\.csv$", ff)] }
if (!length(files)) stop("找不到 *_DEG_*_vs_*.csv")

for (f in files) {
  d <- read.csv(f); A<-d$group1[1]; B<-d$group2[1]
  want <- if (toupper(B)=="ALL") A else c(A,B)
  ic <- detect_ident(obj, want, ident); Idents(obj) <- ic
  genes <- d %>% filter(p_val_adj<PADJ, abs(avg_log2FC)>=LOGFC) %>%
    arrange(desc(abs(avg_log2FC))) %>% head(TOPN) %>% pull(gene)
  genes <- genes[genes %in% rownames(obj)]
  if (!length(genes)) { cat("[skip] 无显著:", f, "\n"); next }
  grps <- if (toupper(B)=="ALL") levels(Idents(obj)) else c(A,B)
  so <- obj[, WhichCells(obj, idents=grps)]; Idents(so)<-ic
  # FetchData+geom_violin, 绕开 Seurat VlnPlot(ggplot2>=4.0 上 S7 崩溃)
  df <- FetchData(so, vars=c(genes, ic), layer="data")
  long <- do.call(rbind, lapply(genes, function(g) data.frame(gene=g, expr=df[[g]], grp=as.character(df[[ic]]))))
  long$gene <- factor(long$gene, levels=genes)
  p <- ggplot(long, aes(grp, expr, fill=grp)) +
    geom_violin(scale="width", trim=TRUE, linewidth=0.2) +
    facet_wrap(~gene, ncol=NCOL, scales="free_y") +
    labs(title=paste0("DEG violin: ",A," vs ",B), x=NULL, y="expression") + theme_classic(base_size=13) +
    theme(legend.position="none", strip.text=element_text(face="bold", size=12),
          axis.text.x=element_text(angle=45, hjust=1, size=11, face="bold"),
          axis.text.y=element_text(size=10, face="bold"),
          axis.title.y=element_text(face="bold", size=13),
          plot.title=element_text(face="bold", hjust=.5, size=15))
  out <- file.path(outdir, sprintf("%s_plot_deg_violin_%s_vs_%s.png", prefix, A, B))
  h <- max(6, ceiling(length(genes)/NCOL)*2.5)
  ggsave(out, p, width=4*NCOL, height=h, dpi=300, limitsize=FALSE); ggsave(sub("\\.png$",".pdf",out), p, width=4*NCOL, height=h, limitsize=FALSE)
  cat("[saved]", out, "\n")
}
