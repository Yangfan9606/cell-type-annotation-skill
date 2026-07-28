#!/usr/bin/env Rscript
# sc_006.plot.find_marker.dotplot.R —— 注释主力 DotPlot (三风格引擎)
#   默认读 results/*_manifest.txt + *_find_marker_{strict,sensitive}_final.csv,
#   每组取 top-N marker 画气泡图, 用于亚类注释的特异性目检。
#
#   风格 (--style): C 默认(蓝+facet留白) | B(spectral+组名色条) | A(红+分隔线, fig2B 风)
#   基因默认在 X 轴 45°; --transpose 翻转成基因在 Y。
#
#   独立运行 (默认 C 风格):        Rscript sc_006.plot.find_marker.dotplot.R
#   换 spectral 组条风:            Rscript sc_006.plot.find_marker.dotplot.R --style B
#   fig2B 红块风, 每组 8 个:       Rscript sc_006.plot.find_marker.dotplot.R --style A --topn 8
#   竖版(基因在 Y):               Rscript sc_006.plot.find_marker.dotplot.R --transpose
#   自定义基因面板(gene,group):    Rscript sc_006.plot.find_marker.dotplot.R --genes panel.csv --style B
#   给组上色:                      ... --style B --group_colors 'c1=#28A828,c2=#1773B3'
#   覆盖: --rds --input(csv) --ident --outfile --outdir --prefix --cols --dot_min --width --height
# ============================== CONFIG (随意改) ==============================
STYLE       <- "C"    # 默认风格: A/block(红分隔线) | B/grouped(spectral组条) | C/gap(蓝+facet留白)
TOPN        <- 5      # 每组 top-N marker
DOT_MIN_PCT <- 0.1    # 低于此表达比例(0-1)的点不画(=dot.min); 设 0 关闭。--dot_min 可覆盖
# ===========================================================================
suppressPackageStartupMessages(library(optparse))   # 只加载轻库; 重库经 run_dotplot() 延后
.here <- function(){ a<-commandArgs(FALSE); m<-grep("^--file=",a,value=TRUE)
  d<-if(length(m)) dirname(normalizePath(sub("^--file=","",m[1]))) else getwd()
  if(file.exists(file.path(d,"sc_006.plot.dotplot_engine.R"))) d else getwd() }
source(file.path(.here(), "sc_006.plot.dotplot_engine.R"))
parser <- OptionParser(option_list=dotplot_option_list(default_style=STYLE, default_topn=TOPN, default_dot_min=DOT_MIN_PCT),
                       description="sc_006.plot.find_marker.dotplot.R — 注释主力 DotPlot(三风格); 读 manifest + find_marker CSV, 每组 top-N marker。",
                       epilogue=DP_EX)
if (length(commandArgs(TRUE))==0) { print_help(parser); quit(status=0) }   # 空运行 → usage
opt <- parse_args(parser)
run_dotplot(opt, job="find_marker")
