#!/usr/bin/env Rscript
# sc_006.plot.dotplot_engine.R —— 三风格 DotPlot 引擎 (find_marker / hetero 两个 wrapper 共用)
# ---------------------------------------------------------------------------
#  不单独运行。由 sc_006.plot.find_marker.dotplot.R / sc_006.plot.hetero.dotplot.R
#  用 source() 载入 (三者部署时必须放在同一目录)。
#
#  三种风格 (--style A|B|C, 也接受 block|grouped|gap):
#    A  block  : grey88→#B2182B(红), 虚线(主/key)+点线(次)分隔, dot.min=0.1 去零点   —— fig2B 面板风
#    B  grouped: 多色 spectral, 组名色条+色带在上方 + 分隔虚线, 支持 --group_colors    —— sc_dotplot 风
#    C  gap    : 默认蓝 lightgrey→blue, facet 物理留白 + 组名 strip (**默认**)          —— 标准蓝+留白
#
#  通用: 默认 **基因在 X 轴 45°, 组在 Y**; --transpose 翻转(基因→Y, 组→X, legend 右)。
#  基因来源优先级: --features "A,B,C"  >  --genes <file(gene,group[,block])>  >  marker CSV 每组 top-N。
# ---------------------------------------------------------------------------
suppressPackageStartupMessages(library(optparse))   # 只加载轻库; 重库(Seurat/ggplot2/dplyr)延后到 run_dotplot 内, 加快 --help/usage 首响应

`%||%` <- function(a,b) if (is.null(a)||(length(a)==1&&is.na(a))||(is.character(a)&&length(a)==1&&a=="")) b else a

read_manifest <- function(mp, outdir){
  if (is.null(mp)){ cc<-list.files(outdir,"_manifest\\.txt$",full.names=TRUE); if(length(cc)) mp<-cc[which.max(file.mtime(cc))] }
  m<-list(); if(!is.null(mp)&&file.exists(mp)) for(ln in readLines(mp)) if(grepl("=",ln)){k<-sub("=.*","",ln);m[[k]]<-sub("^[^=]*=","",ln)}
  m
}

# ---- shared CLI option list (both wrappers reuse) ---------------------------
dotplot_option_list <- function(default_style="C", default_topn=5L, default_dot_min=0.1) list(
  make_option("--manifest",     default=NULL, help="manifest path (default: newest results/*_manifest.txt)"),
  make_option("--rds",          default=NULL, help="Seurat rds (default from manifest)"),
  make_option("--input",        default=NULL, help="marker CSV (default: auto-pick per job)"),
  make_option("--outfile",      default=NULL, help="output PNG path (PDF sibling auto)"),
  make_option("--outdir",       default="results", help="[default %default]"),
  make_option("--prefix",       default=NULL),
  make_option("--style",        default=default_style, help="A|block / B|grouped / C|gap  [default %default]"),
  make_option("--features",     default=NULL, help="inline gene list 'A,B,C' (skips CSV, no grouping)"),
  make_option("--genes",        default=NULL, help="gene file: gene,group[,block]  (skips CSV)"),
  make_option("--ident",        default=NULL, help="grouping/ident column (find_marker; auto-detected from CSV)"),
  make_option("--topn",         default=default_topn, type="integer", help="top-N markers per group [default %default]"),
  make_option("--dot_scale",    default=6, type="double", help="max dot size [default %default]"),
  make_option("--transpose",    action="store_true", default=FALSE, help="genes on Y, groups on X, legend right"),
  make_option("--cols",         default=NULL, help="'low,high' 2-color OR name: spectral|blue|red|viridis"),
  make_option("--dot_min",      default=default_dot_min, type="double", help="低于此表达比例(0-1)的点不画 [默认 %default]; 传 0 关闭"),
  make_option("--drop_zero",    action="store_true", default=FALSE, help="drop zero-pct dots"),
  make_option("--keep_zero",    action="store_true", default=FALSE, help="force-keep zero dots (overrides style A)"),
  make_option("--keep_neg",     action="store_true", default=FALSE, help="include negative markers (default pos only)"),
  make_option("--group_colors", default=NULL, help="'NAME=COLOR,NAME=COLOR' color block labels (style B/axis)"),
  make_option("--marker_list",  default=NULL, help="自定义 marker 表 gene,group[,block]; 开启则不读 sc_006 CSV (=--genes 别名)"),
  make_option("--baseline",     default=NULL, help="共有/baseline 组名(逗号); 置最前 + 一条主虚线, 不作 Y 行 (fig2B; style A 最明显)"),
  make_option("--group_order",  default=NULL, help="Y 轴(及 X 块)组顺序 'c1,c2,c3'(首个在顶); 默认随 marker_list 组序"),
  make_option("--group_map",    default=NULL, help="cluster→谱系映射 'Exc=6,10,11;Inh=3,7'(或文件 lineage,cluster); Y=cluster 按谱系分面(左侧 strip)"),
  make_option("--width",        default=NULL, type="double", help="PDF/PNG width in (auto if omitted)"),
  make_option("--height",       default=NULL, type="double", help="PDF/PNG height in (auto if omitted)")
)

# ---- helpers ----------------------------------------------------------------
.parse_group_colors <- function(s){
  if (is.null(s) || !nzchar(trimws(s))) return(character(0))
  out <- character(0)
  for (p in trimws(strsplit(s, ",")[[1]])) {
    if (!nzchar(p)) next
    kv <- strsplit(p, "=", fixed=TRUE)[[1]]
    if (length(kv)==2 && nzchar(trimws(kv[1])) && nzchar(trimws(kv[2]))) out[trimws(kv[1])] <- trimws(kv[2])
  }
  out
}

.parse_group_map <- function(s){
  # 返回命名 list: 谱系 -> cluster 向量。支持内联 'L=1,2;M=3' 或文件(两列 lineage,cluster)。
  if (is.null(s) || !nzchar(trimws(s))) return(NULL)
  if (file.exists(s)){
    raw <- readLines(s, warn=FALSE); raw <- raw[nzchar(trimws(raw))]
    delim <- if (grepl("\t",raw[1])) "\t" else if (grepl(",",raw[1])) "," else " "
    df <- read.table(s, header=TRUE, sep=delim, strip.white=TRUE, stringsAsFactors=FALSE, check.names=FALSE)
    cn <- tolower(colnames(df))
    li <- which(cn %in% c("lineage","group","celltype","type"))[1]; if(is.na(li)) li<-1L
    ci <- which(cn %in% c("cluster","clusters","id"))[1];            if(is.na(ci)) ci<-setdiff(seq_len(ncol(df)),li)[1]
    split(as.character(df[[ci]]), factor(as.character(df[[li]]), levels=unique(as.character(df[[li]]))))
  } else {
    ml <- list()
    for (p in trimws(strsplit(s, ";")[[1]])) { if(!nzchar(p)) next
      kv <- strsplit(p, "=", fixed=TRUE)[[1]]; if (length(kv)!=2) next
      ml[[trimws(kv[1])]] <- trimws(strsplit(kv[2], ",")[[1]]) }
    if (length(ml)) ml else NULL
  }
}

.read_gene_file <- function(path){
  raw <- readLines(path, warn=FALSE); raw <- raw[nzchar(trimws(raw))]
  if(!length(raw)) stop("gene file empty: ", path)
  h <- raw[1]
  delim <- if (grepl("\t",h)) "\t" else if (grepl(",",h)) "," else if (grepl(" ",trimws(h))) " " else NULL
  if (!is.null(delim)){
    df <- tryCatch(read.table(path, header=TRUE, sep=delim, strip.white=TRUE, stringsAsFactors=FALSE,
                              check.names=FALSE, comment.char=""), error=function(e) NULL)
    if (!is.null(df) && ncol(df)>=2){
      cn <- tolower(colnames(df))
      gi <- which(cn %in% c("gene","genes","feature","features"))[1]; if(is.na(gi)) gi<-1L
      pi <- which(cn %in% c("group","groups","cluster","label","celltype"))[1]
      if(is.na(pi)) pi<-setdiff(seq_len(ncol(df)),gi)[1]
      bi <- which(cn %in% c("block","panel"))[1]
      return(data.frame(gene=as.character(df[[gi]]), group=as.character(df[[pi]]),
                        block=if(!is.na(bi)) as.character(df[[bi]]) else NA_character_, stringsAsFactors=FALSE))
    }
  }
  g <- trimws(raw); if (tolower(g[1]) %in% c("gene","genes","feature","features")) g<-g[-1]
  data.frame(gene=g, group=NA_character_, block=NA_character_, stringsAsFactors=FALSE)
}

.resolve_cols <- function(cols, style){
  SPECTRAL <- c("#E8E8EE","#A8C8E8","#4A90D9","#F5C518","#E86A10","#B22222")
  if (!is.null(cols) && nzchar(cols)){
    if (grepl(",", cols)) return(list(type="grad2", val=trimws(strsplit(cols,",")[[1]])))
    nm <- tolower(cols)
    if (nm=="spectral") return(list(type="gradn",   val=SPECTRAL))
    if (nm=="viridis")  return(list(type="viridis", val=NULL))
    if (nm=="blue")     return(list(type="grad2",   val=c("lightgrey","blue")))
    if (nm=="red")      return(list(type="grad2",   val=c("grey88","#B2182B")))
    return(list(type="grad2", val=c("lightgrey", nm)))            # single color name = grey→color
  }
  if (style=="A") list(type="grad2", val=c("grey88","#B2182B"))
  else if (style=="B") list(type="gradn", val=SPECTRAL)
  else list(type="grad2", val=c("lightgrey","blue"))
}

# ---- shared examples (used as --help/usage epilogue by both wrappers) --------
DP_EX <- "
示例 (find_marker.dotplot / hetero.dotplot 通用):
  # 默认 C 风格(蓝+facet留白), 读最新 manifest
  Rscript sc_006.plot.find_marker.dotplot.R --outdir results
  # B 风格 spectral + 组名色条, 给组上色
  Rscript sc_006.plot.find_marker.dotplot.R --style B --group_colors 'c1=#28A828,c2=#1773B3'
  # A 风格 红块(fig2B), 每组 8 个, 转置(基因在 Y)
  Rscript sc_006.plot.find_marker.dotplot.R --style A --topn 8 --transpose
  # 自定义 marker 表(gene,group), 不读 sc_006 CSV; --baseline 那组=共有(主虚线, 不占 Y 行)
  Rscript sc_006.plot.find_marker.dotplot.R --marker_list panel.csv --style A --baseline Shared
  # 指定 Y 轴(及 X 块)组顺序(首个在顶, X/Y 同排保持对角)
  Rscript sc_006.plot.find_marker.dotplot.R --marker_list panel.csv --group_order 'MotorNeuron,Excitatory,Inhibitory'
  # Y=cluster 但按谱系分组: --group_map(左侧谱系 strip + 间隙); X 基因块=谱系, 与 Y 对齐
  Rscript sc_006.plot.find_marker.dotplot.R --marker_list panel.csv --style A --baseline Shared \
    --ident seurat_clusters --group_map 'Excitatory=6,10,11,13,24,25; Inhibitory=3,7,9; MotorNeuron=15'
  # 内联基因 / 亚聚类 dotplot 同理
  Rscript sc_006.plot.find_marker.dotplot.R --features 'CD3D,MS4A1,LYZ'
  Rscript sc_006.plot.hetero.dotplot.R --style C --transpose

风格: A/block=红+分隔线(dot.min去零) | B/grouped=spectral+组条 | C/gap=蓝+facet留白(默认)
基因来源优先级: --features > --marker_list/--genes(gene,group[,block]) > marker CSV(每组 top-N)
marker 表格式: 逗号/制表/空格分隔; 列 gene + group(细胞组); --baseline 命名的组=共有块(仅 X 列, 一条主虚线)
--group_order 定 Y(首个在顶)且 X 非baseline 块同排; 默认随 marker_list 组序。默认基因在 X 45°, --transpose 转 Y。
"

# =============================================================================
run_dotplot <- function(opt, job=c("find_marker","hetero")){
  job <- match.arg(job)
  suppressPackageStartupMessages({library(Seurat); library(ggplot2); library(dplyr)})  # 重库按需加载
  man    <- read_manifest(opt$manifest, opt$outdir)
  outdir <- opt$outdir %||% man$outdir %||% "results"
  prefix <- opt$prefix %||% man$prefix %||% "sc006"

  st <- toupper(substr(opt$style,1,1))
  if (!st %in% c("A","B","C")) st <- switch(tolower(opt$style), block="A", grouped="B", gap="C", "C")

  # ---- object + default group column --------------------------------------
  if (job=="find_marker"){
    rds <- opt$rds %||% man$rds
    if (is.null(rds) || !file.exists(rds)) stop("找不到 rds (给 --rds 或确保 manifest 有 rds=)")
    obj <- readRDS(rds); default_group <- opt$ident %||% man$ident %||% "cell_type"
  } else {
    rds <- opt$rds %||% man$hetero_rds %||% {
      ff<-list.files(outdir,"_hetero_.*\\.rds$",full.names=TRUE); if(length(ff)) ff[which.max(file.mtime(ff))] else NULL }
    if (is.null(rds) || !file.exists(rds)) stop("找不到 hetero rds")
    obj <- readRDS(rds); default_group <- opt$ident %||% "subcluster"
  }
  cat("[rds]", rds, "\n")
  a <- DefaultAssay(obj)
  if (inherits(obj[[a]], "Assay5")) obj[[a]] <- tryCatch(JoinLayers(obj[[a]]), error=function(e) obj[[a]])

  # ---- gene / group table -------------------------------------------------
  mlist    <- opt$marker_list %||% opt$genes                          # 自定义 marker 表(开启则不读 CSV)
  base_grp <- if (!is.null(opt$baseline)    && nzchar(opt$baseline))    trimws(strsplit(opt$baseline,",")[[1]])    else character(0)
  gorder   <- if (!is.null(opt$group_order) && nzchar(opt$group_order)) trimws(strsplit(opt$group_order,",")[[1]]) else NULL
  gmap     <- .parse_group_map(opt$group_map)                         # cluster→谱系 (Y=cluster 按谱系分面)
  if (!is.null(opt$features) && nzchar(opt$features)){
    g <- trimws(strsplit(opt$features,",")[[1]]); g<-g[nzchar(g)]
    gt <- data.frame(gene=g, group=NA_character_, block=NA_character_, stringsAsFactors=FALSE)
    cat("[genes] inline:", length(g), "\n")
  } else if (!is.null(mlist) && file.exists(mlist)){
    gt <- .read_gene_file(mlist); cat("[genes] marker_list:", mlist, " (跳过 sc_006 CSV)\n")
  } else {
    csv <- opt$input
    if (is.null(csv)){
      if (job=="find_marker"){
        cand <- file.path(outdir, paste0(prefix, c("_find_marker_strict_final.csv","_find_marker_sensitive_final.csv")))
        for (f in cand) if (file.exists(f) && nrow(read.csv(f))>0){ csv<-f; break }
        if (is.null(csv)) { e<-cand[file.exists(cand)]; if(length(e)) csv<-e[1] }
      } else {
        ff <- list.files(outdir, "_hetero_.*_subcluster_markers\\.csv$", full.names=TRUE)
        if (length(ff)) csv <- ff[which.max(file.mtime(ff))]
      }
    }
    if (is.null(csv) || !file.exists(csv)) stop("找不到 marker CSV (给 --input / --genes / --features)")
    cat("[csv]", csv, "\n")
    mk <- read.csv(csv)
    if (!"cluster" %in% names(mk)) stop("CSV 缺 cluster 列")
    if (!opt$keep_neg && "avg_log2FC" %in% names(mk)) mk <- mk %>% filter(avg_log2FC > 0)
    mk <- mk %>% group_by(cluster) %>% slice_max(avg_log2FC, n=opt$topn, with_ties=FALSE) %>% ungroup()
    mk <- mk %>% group_by(gene) %>% slice_max(avg_log2FC, n=1, with_ties=FALSE) %>% ungroup()  # 1 gene → 1 best group
    gt <- data.frame(gene=as.character(mk$gene), group=as.character(mk$cluster),
                     block=NA_character_, stringsAsFactors=FALSE)
  }
  has_groups <- !all(is.na(gt$group))
  all_grp  <- if (has_groups) unique(as.character(gt$group[!is.na(gt$group)])) else character(0)
  base_grp <- intersect(base_grp, all_grp)          # 只保留确实出现的 baseline 组
  nonbase  <- setdiff(all_grp, base_grp)            # 非 baseline 组 = 对应 Y 细胞组

  # ---- resolve ident (探测承载 非baseline 组名 的列; baseline 不参与) --------
  #      显式 --ident 时直接用它, 不再自动探测(支持 Y=cluster 而 X 基因组=谱系 的场景)
  group_col <- default_group
  if (job=="find_marker" && is.null(opt$ident)){
    want <- if (!is.null(gmap)) unique(unlist(gmap, use.names=FALSE)) else nonbase   # gmap 时按 cluster 值探测
    detect <- function(obj,want,fb){ if(!length(want)) return(fb)
      for(cn in colnames(obj@meta.data)){v<-obj@meta.data[[cn]]
        if((is.character(v)||is.factor(v))&&all(want %in% unique(as.character(v)))) return(cn)}; fb }
    group_col <- detect(obj, want, default_group)
  }
  if (!group_col %in% colnames(obj@meta.data)) stop("分组列不存在: ", group_col)
  Idents(obj) <- group_col
  ident_vals <- unique(as.character(obj@meta.data[[group_col]]))

  # ---- Y 轴组顺序: group_map(cluster→谱系) > group_order > marker_list 组序 > 全部
  ylin <- NULL; lin_ord <- NULL                     # cluster→谱系 (group_map 分面用)
  if (!is.null(gmap)){
    miss_lin <- setdiff(names(gmap), nonbase)          # gmap 谱系名在 marker_list 里找不到同名基因组
    if (length(miss_lin)) cat(
      "[warn] --group_map 谱系名在 marker_list 里【无同名基因组】→ 其 Y 行会排到末尾, 与基因块错位! 请核对拼写:\n",
      "       对不上:", paste(miss_lin, collapse=", "), "\n",
      "       marker_list 组名:", paste(nonbase, collapse=", "), "\n")
    lin_ord <- c(intersect(nonbase, names(gmap)), setdiff(names(gmap), nonbase))   # 谱系序跟 panel(X) 对齐
    yflat   <- unlist(lapply(lin_ord, function(L) gmap[[L]]), use.names=FALSE)
    ylin    <- setNames(rep(lin_ord, vapply(lin_ord, function(L) length(gmap[[L]]), 1L)), yflat)
    y_src   <- yflat
    if (opt$transpose) cat("[warn] group_map 的谱系分面暂仅非转置支持; 本次仅按谱系排序/子集\n")
  } else {
    y_src <- if (!is.null(gorder)) gorder else if (!is.null(mlist) && has_groups) nonbase else NULL
  }
  if (!is.null(y_src)){
    yg <- y_src[y_src %in% ident_vals]
    ymiss <- setdiff(y_src, ident_vals)
    if (length(ymiss)) cat("[warn] group_order/group_map/marker_list 中不在对象里的组(跳过):", paste(ymiss, collapse=", "), "\n")
    if (!length(yg)) stop("没有可显示的 Y 组(检查 --group_order / --group_map / marker_list 组名是否匹配 ", group_col, ")")
    obj <- subset(obj, idents = yg)
    lev <- if (opt$transpose) yg else rev(yg)       # 非转置: Y 从下往上, rev 使首个在顶
    obj@meta.data[[group_col]] <- factor(as.character(obj@meta.data[[group_col]]), levels=lev)
    Idents(obj) <- group_col
    if (!is.null(ylin)) ylin <- ylin[yg]            # 只留实际显示的 cluster
  }
  cat("[ident]", group_col, "  [style]", st,
      if(length(base_grp)) paste0("  [baseline:", paste(base_grp,collapse=","), "]") else "",
      if(opt$transpose)"  [transpose]" else "", "\n")

  # ---- X 基因块顺序: baseline 在前; 非baseline 跟 group_order(否则文件序) -----
  if (has_groups){
    nb_order <- if (!is.null(gorder)) c(intersect(gorder, nonbase), setdiff(nonbase, gorder)) else nonbase
    xlev <- c(base_grp, nb_order)
    gt$group <- factor(gt$group, levels=xlev)
    if (!all(is.na(gt$block))) gt$block <- factor(gt$block, levels=unique(gt$block))
    gt <- gt[order(gt$group), , drop=FALSE]
  }

  # ---- filter genes to object, drop duplicated symbols --------------------
  have <- rownames(obj); dup <- unique(have[duplicated(have)])
  miss <- setdiff(unique(gt$gene), have)
  if (length(miss)) cat("[warn] 缺失基因跳过:", paste(head(miss,20),collapse=", "), if(length(miss)>20)"..." else "", "\n")
  gt <- gt[gt$gene %in% have & !(gt$gene %in% dup), , drop=FALSE]
  gt <- gt[!duplicated(gt$gene), , drop=FALSE]
  if (!nrow(gt)) stop("过滤后无可用基因")
  if (has_groups) gt$group <- droplevels(gt$group)   # 去掉被基因过滤清空的组
  genes <- gt$gene
  gene2group <- setNames(as.character(gt$group), gt$gene)
  cat("[genes] usable:", length(genes), "\n")

  # ---- palette / limits / dot_min / drop_zero -----------------------------
  pal     <- .resolve_cols(opt$cols, st)
  col_min <- if (st=="A") -1 else if (st=="B") -1.5 else -2.5
  col_max <- if (st=="A")  2 else if (st=="B")  2.5 else  2.5
  dot_min <- opt$dot_min %||% 0.1                     # 三风格统一默认 DOT_MIN_PCT=0.1; 传 0 关闭
  drop_zero <- if (opt$keep_zero) FALSE else if (opt$drop_zero) TRUE else (st=="A")
  base_cols <- if (pal$type=="grad2") pal$val else c("lightgrey","blue")

  p <- DotPlot(obj, features=genes, group.by=group_col, cols=base_cols,
               dot.scale=opt$dot_scale, dot.min=dot_min, col.min=col_min, col.max=col_max)
  if (pal$type=="gradn")   p <- suppressMessages(p + scale_color_gradientn(colors=pal$val, name="Avg\nexpr"))
  if (pal$type=="viridis") p <- suppressMessages(p + scale_color_viridis_c(name="Avg\nexpr"))
  if (drop_zero) p$data <- p$data[is.na(p$data$pct.exp) | p$data$pct.exp > 0, , drop=FALSE]
  if (has_groups) p$data$.grp <- factor(gene2group[as.character(p$data$features.plot)], levels=levels(gt$group))

  # ---- style C + transpose: rebuild natively (genes on Y) ------------------
  # coord_flip defeats facet free-scale dropping, so for the faceted-gap look in
  # transposed orientation we rebuild from p$data with genes on Y and facet_grid
  # rows (scales/space free_y) — real per-block gaps, only each block's own genes.
  rebuilt <- FALSE
  if (st=="C" && has_groups && opt$transpose){
    d <- p$data
    q <- ggplot(d, aes(x=id, y=features.plot, size=pct.exp, color=avg.exp.scaled)) + geom_point()
    if (pal$type=="grad2")        q <- q + scale_color_gradient(low=pal$val[1], high=pal$val[2],
                                        limits=c(col_min,col_max), oob=scales::squish, name="Avg\nexpr")
    else if (pal$type=="gradn")   q <- q + scale_color_gradientn(colors=pal$val,
                                        limits=c(col_min,col_max), oob=scales::squish, name="Avg\nexpr")
    else                          q <- q + scale_color_viridis_c(limits=c(col_min,col_max),
                                        oob=scales::squish, name="Avg\nexpr")
    q <- q + scale_radius(range=c(0, opt$dot_scale), name="Percent Expressed") +
      facet_grid(rows=vars(.grp), scales="free_y", space="free_y", switch="y") +
      theme_classic(base_size=13) +
      theme(axis.title=element_blank(), plot.title=element_text(face="bold", hjust=.5, size=15),
            legend.title=element_text(face="bold", size=13), legend.text=element_text(size=11),
            panel.grid.major=element_line(color="grey92", linewidth=0.3),
            axis.text.x=element_text(angle=45,hjust=1,vjust=1,face="bold",size=14),
            axis.text.y=element_text(face="bold.italic",size=13),
            strip.text.y.left=element_text(angle=0, face="bold", size=13),
            strip.background=element_rect(fill="grey92", color=NA),
            panel.spacing=unit(0.5,"lines"), legend.position="right")
    p <- q; rebuilt <- TRUE
  }

  if (!rebuilt){
  # ---- base theme + orientation (文字统一加粗放大) -------------------------
  p <- p + theme_classic(base_size=13) +
    theme(axis.title=element_blank(),
          plot.title=element_text(face="bold", hjust=.5, size=15),
          legend.title=element_text(face="bold", size=13), legend.text=element_text(size=11),
          strip.text=element_text(face="bold", size=13),
          panel.grid.major=element_line(color="grey92", linewidth=0.3))
  if (st=="A") p <- p + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())  # fig2B: 去网格线, 只留虚线/点线
  if (opt$transpose){
    p <- p + theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1,face="bold",size=14),
                   axis.text.y=element_text(face="bold.italic",size=13), legend.position="right")
  } else {
    p <- p + theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1,face="bold.italic",size=13),
                   axis.text.y=element_text(face="bold",size=14), legend.position="top")
  }
  if (st=="A") p <- p + theme(legend.position="right")   # fig2B: 宽图右侧图例, 避免顶部裁切

  ug   <- .parse_group_colors(opt$group_colors)
  glev <- if (has_groups) levels(gt$group) else character(0)

  # ---- style C: facet with gaps + strip labels ----------------------------
  if (st=="C" && has_groups && is.null(gmap)){       # gmap 时由下方谱系分面统一处理
    if (opt$transpose){
      # facet_wrap(ncol=1) frees the gene axis PER panel (facet_grid rows can't free x);
      # coord_flip (added later) then stacks the per-group gene blocks vertically with gaps.
      p <- p + facet_wrap(~.grp, ncol=1, scales="free_x", strip.position="right") +
        theme(strip.text.y=element_text(angle=0, face="bold"),
              strip.background=element_rect(fill="grey92", color=NA),
              panel.spacing=unit(0.5,"lines"))
    } else {
      p <- p + facet_grid(cols=vars(.grp), scales="free_x", space="free_x") +
        theme(strip.text.x=element_text(face="bold"),
              strip.background=element_rect(fill="grey92", color=NA),
              panel.spacing=unit(0.5,"lines"))
    }
  }

  # ---- style A: 一条主虚线(baseline/共有 之后) + 点线(各 cluster 之间) --------
  if (st=="A" && has_groups){
    gf <- gt$group
    allb    <- which(gf[-length(gf)] != gf[-1]) + 0.5      # 所有相邻组边界
    is_base <- as.character(gf) %in% base_grp              # baseline 已排在最前
    major <- numeric(0)
    if (any(is_base) && any(!is_base)) {
      major <- sum(is_base) + 0.5                          # baseline→非baseline 处一条主虚线
    } else if (!length(base_grp) && !all(is.na(gt$block))) {
      bf <- gt$block; major <- which(bf[-length(bf)] != bf[-1]) + 0.5   # 无 --baseline 时退回 block 列
    }
    minor <- setdiff(allb, major)                          # 其余组边界 = 点线
    if (length(minor)) p <- p + geom_vline(xintercept=minor, linetype="dotted", color="grey60", linewidth=0.4)
    if (length(major)) p <- p + geom_vline(xintercept=major, linetype="dashed", color="grey25", linewidth=0.6)  # Shared 长虚线
  }

  # ---- style B: colored group bars + labels above (non-transpose) ----------
  if (st=="B" && has_groups){
    gf <- gt$group; gpos <- seq_along(genes)
    sep <- which(gf[-length(gf)] != gf[-1]) + 0.5
    if (length(sep)) p <- p + geom_vline(xintercept=sep, linetype="dashed", color="grey50", linewidth=0.5)
    if (!opt$transpose && is.null(gmap)){            # gmap 分面时省略顶部组条(会逐面重复)
      GP  <- c("#4C72B0","#DD8452","#55A868","#C44E52","#8172B3","#937860","#DA8BC3","#8C8C8C","#CCB974","#64B5CD")
      mid <- tapply(gpos, gf, mean); rng <- tapply(gpos, gf, range)
      gc  <- setNames(rep(NA_character_,length(glev)), glev)
      hit <- intersect(glev, names(ug)); if(length(hit)) gc[hit]<-ug[hit]
      un  <- glev[is.na(gc)]
      if(length(un)){ av<-setdiff(GP, gc[!is.na(gc)]); if(!length(av)) av<-GP; gc[un]<-av[(seq_along(un)-1)%%length(av)+1] }
      nI  <- length(levels(Idents(obj)))
      p <- p +
        annotate("text", x=as.numeric(mid[glev]), y=nI+1.2, label=glev, color=gc[glev],
                 fontface="bold", size=5, hjust=.5) +
        lapply(seq_along(glev), function(i){ g<-glev[i]
          annotate("segment", x=rng[[g]][1]-0.45, xend=rng[[g]][2]+0.45, y=nI+0.6, yend=nI+0.6,
                   color=gc[i], linewidth=3, lineend="round") }) +
        theme(plot.margin=margin(42,16,6,6))
    }
  }

  # ---- group_map: Y 轴按谱系分面(左侧 strip = 谱系名, 面间隙 = 分隔) ----------
  if (!is.null(gmap) && !opt$transpose && !is.null(ylin)){
    p$data$.ylin <- factor(ylin[as.character(p$data$id)], levels=intersect(lin_ord, unique(ylin)))
    if (st=="C" && has_groups)
      p <- p + facet_grid(rows=vars(.ylin), cols=vars(.grp), scales="free", space="free", switch="y")
    else
      p <- p + facet_grid(rows=vars(.ylin), scales="free_y", space="free_y", switch="y")
    p <- p + theme(strip.text.y.left=element_text(angle=0, face="bold"),
                   strip.background=element_rect(fill="grey92", color=NA),
                   panel.spacing=unit(0.4,"lines"))
  }

  # ---- single coord (flip for transpose; clip off for B bars) --------------
  if (opt$transpose) p <- p + coord_flip(clip="off")
  else if (st=="B" && has_groups && !opt$transpose && is.null(gmap)) p <- p + coord_cartesian(clip="off")
  }  # end if(!rebuilt)

  ttl <- if (job=="find_marker") "Find-marker DotPlot" else "Hetero subcluster DotPlot"
  p <- p + labs(title=sprintf("%s  (top %d/group · style %s%s)", ttl, opt$topn, st,
                              if(opt$transpose)" · T" else ""))

  # ---- size + save --------------------------------------------------------
  nG <- length(genes); nI <- length(levels(Idents(obj)))
  if (opt$transpose){
    w <- opt$width  %||% max(5, 2 + nI*0.5 + 1.5)
    h <- opt$height %||% max(5, 1.5 + nG*0.28)
  } else {
    w <- opt$width  %||% max(6, 1.7 + nG*0.32 + 1)
    h <- opt$height %||% max(4, 2.5 + nI*0.4 + (if(st=="B") 1 else 0) +
                            (if(!is.null(ylin)) length(unique(ylin))*0.35 else 0))   # 谱系分面留白
  }
  out <- opt$outfile %||% file.path(outdir, paste0(prefix, "_plot_", job, "_dotplot.png"))
  dir.create(dirname(out), showWarnings=FALSE, recursive=TRUE)
  ggsave(out, p, width=w, height=h, dpi=300, limitsize=FALSE)
  ggsave(sub("\\.png$",".pdf",out), p, width=w, height=h, limitsize=FALSE)
  cat("[saved]", out, sprintf("  (%d genes × %d groups, style %s%s)\n",
                              nG, nI, st, if(opt$transpose)", transposed" else ""))
  invisible(p)
}

# ---- 若被【直接】Rscript 运行(而非被 wrapper source): 打印说明并退出 ----
.dpe_invoked <- function(){ a<-commandArgs(FALSE); m<-grep("^--file=",a,value=TRUE)
  if(length(m)) basename(sub("^--file=","",m[1])) else "" }
if (identical(.dpe_invoked(), "sc_006.plot.dotplot_engine.R")) {
  cat("\n本文件是 DotPlot 三风格【共享引擎】, 由以下两个 wrapper 用 source() 调用, 不单独作图:\n",
      "  sc_006.plot.find_marker.dotplot.R   (注释主力 dotplot)\n",
      "  sc_006.plot.hetero.dotplot.R        (亚聚类 dotplot)\n",
      "两个 wrapper 与本文件必须同目录。共享参数 / 示例:\n", sep="")
  print_help(OptionParser(option_list=dotplot_option_list(),
             usage="(本引擎不直接运行; 请跑上面的 wrapper)", epilogue=DP_EX))
  quit(status=0)
}
