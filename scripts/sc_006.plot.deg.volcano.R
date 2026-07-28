#!/usr/bin/env Rscript
# sc_006.plot.deg.volcano.R —— DEG 火山图 (A vs B); 批量: 对每个 *_DEG_*_vs_*.csv 各出一张
#   只需 CSV, 不需 rds。 Rscript sc_006.plot.deg.volcano.R
# ============================== CONFIG ==============================
LOGFC <- 0.25; PADJ <- 0.05; LABEL_N <- 10
UP_COL <- "red3"; DN_COL <- "blue3"; NS_COL <- "grey70"
# ===================================================================
suppressPackageStartupMessages(library(optparse))   # 只加载轻库; 重库延后, 加快 --help/usage
EX <- "
示例:
  Rscript sc_006.plot.deg.volcano.R --outdir results     # 对 outdir 里每个 *_DEG_*_vs_*.csv 各出一张
  Rscript sc_006.plot.deg.volcano.R --input myrun_DEG_Disease_vs_Control.csv
只需 CSV, 不需 rds。细调 LOGFC/PADJ/LABEL_N/配色: 改脚本顶部 CONFIG 块。
"
parser <- OptionParser(option_list=list(
  make_option("--manifest", default=NULL, help="manifest 路径(默认取 outdir 最新, 仅用于定位 prefix/outdir)"),
  make_option("--input", default=NULL, help="单个 DEG CSV(默认自动扫全部 *_DEG_*_vs_*.csv)"),
  make_option("--outdir", default="results"), make_option("--prefix", default=NULL)),
  description="sc_006.plot.deg.volcano.R — DEG 火山图(A vs B), 批量每对比一张", epilogue=EX)
if (length(commandArgs(TRUE))==0) { print_help(parser); quit(status=0) }
opt <- parse_args(parser)
suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(ggrepel)})
`%||%` <- function(a,b) if (is.null(a)||(length(a)==1&&is.na(a))||(is.character(a)&&length(a)==1&&a=="")) b else a
read_manifest <- function(mp, outdir){
  if (is.null(mp)){ cc<-list.files(outdir,"_manifest\\.txt$",full.names=TRUE); if(length(cc)) mp<-cc[which.max(file.mtime(cc))] }
  m<-list(); if(!is.null(mp)&&file.exists(mp)) for(ln in readLines(mp)) if(grepl("=",ln)){k<-sub("=.*","",ln);m[[k]]<-sub("^[^=]*=","",ln)}; m }
man <- read_manifest(opt$manifest, opt$outdir)
outdir<-opt$outdir%||%man$outdir%||%"results"; prefix<-opt$prefix%||%man$prefix%||%"sc006"
files <- if (!is.null(opt$input)) opt$input else {
  ff <- list.files(outdir, pattern=paste0("^",prefix,"_DEG_.*_vs_.*\\.csv$"), full.names=TRUE)
  ff[!grepl("_sig\\.csv$", ff)] }
if (!length(files)) stop("找不到 *_DEG_*_vs_*.csv (先跑 --do_deg)")

for (f in files) {
  d <- read.csv(f)
  A <- if("group1"%in%names(d)) d$group1[1] else "A"; B <- if("group2"%in%names(d)) d$group2[1] else "B"
  d$sig <- "NS"
  d$sig[d$avg_log2FC>LOGFC  & d$p_val_adj<PADJ] <- "Up"
  d$sig[d$avg_log2FC< -LOGFC & d$p_val_adj<PADJ] <- "Down"
  lab <- rbind(d %>% filter(sig=="Up") %>% arrange(desc(avg_log2FC)) %>% head(LABEL_N),
               d %>% filter(sig=="Down") %>% arrange(avg_log2FC) %>% head(LABEL_N))
  p <- ggplot(d, aes(avg_log2FC, -log10(p_val_adj+1e-300), color=sig)) +
    geom_point(alpha=.5, size=1.3) +
    scale_color_manual(values=c(Up=UP_COL, Down=DN_COL, NS=NS_COL)) +
    geom_vline(xintercept=c(-LOGFC,LOGFC), linetype="dashed", alpha=.5) +
    geom_hline(yintercept=-log10(PADJ), linetype="dashed", alpha=.5) +
    ggrepel::geom_text_repel(data=lab, aes(label=gene), size=4, fontface="bold", max.overlaps=15, show.legend=FALSE) +
    labs(title=paste0("DEG volcano: ", A, " vs ", B), x="log2FC", y="-log10(padj)", color="") +
    theme_classic(base_size=13) +
    theme(plot.title=element_text(face="bold", hjust=.5, size=15),
          axis.title=element_text(face="bold", size=13), axis.text=element_text(face="bold", size=11),
          legend.text=element_text(size=11))
  out <- file.path(outdir, sprintf("%s_plot_deg_volcano_%s_vs_%s.png", prefix, A, B))
  ggsave(out, p, width=9, height=7, dpi=300); ggsave(sub("\\.png$",".pdf",out), p, width=9, height=7)
  cat("[saved]", out, "\n")
}
