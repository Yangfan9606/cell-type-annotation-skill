#!/usr/bin/env Rscript
# sc_006.plot.hetero.dotplot.R —— 亚聚类 marker DotPlot (三风格引擎, 与 find_marker 共用)
#   默认读 hetero rds + *_hetero_*_subcluster_markers.csv, 每个 subcluster 取 top-N marker。
#
#   风格 (--style): C 默认(蓝+facet留白) | B(spectral+组名色条) | A(红+分隔线)
#   基因默认在 X 轴 45°; --transpose 翻转成基因在 Y。
#
#   独立运行 (默认 C 风格):   Rscript sc_006.plot.hetero.dotplot.R
#   换风格 / 竖版:            Rscript sc_006.plot.hetero.dotplot.R --style B --transpose
#   覆盖: --rds --input(csv) --outfile --outdir --prefix --topn --cols --dot_min --width --height
# ============================== CONFIG (随意改) ==============================
STYLE       <- "C"    # 默认风格: A/block | B/grouped | C/gap
TOPN        <- 5      # 每个 subcluster top-N marker
DOT_MIN_PCT <- 0.1    # 低于此表达比例(0-1)的点不画(=dot.min); 设 0 关闭。--dot_min 可覆盖
# ===========================================================================
suppressPackageStartupMessages(library(optparse))   # 只加载轻库; 重库经 run_dotplot() 延后
.here <- function(){ a<-commandArgs(FALSE); m<-grep("^--file=",a,value=TRUE)
  d<-if(length(m)) dirname(normalizePath(sub("^--file=","",m[1]))) else getwd()
  if(file.exists(file.path(d,"sc_006.plot.dotplot_engine.R"))) d else getwd() }
source(file.path(.here(), "sc_006.plot.dotplot_engine.R"))
parser <- OptionParser(option_list=dotplot_option_list(default_style=STYLE, default_topn=TOPN, default_dot_min=DOT_MIN_PCT),
                       description="sc_006.plot.hetero.dotplot.R — 亚聚类 DotPlot(三风格); 读 hetero rds + subcluster markers CSV, 每个 subcluster top-N。",
                       epilogue=DP_EX)
if (length(commandArgs(TRUE))==0) { print_help(parser); quit(status=0) }   # 空运行 → usage
opt <- parse_args(parser)
run_dotplot(opt, job="hetero")
