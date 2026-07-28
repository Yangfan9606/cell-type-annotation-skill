#!/usr/bin/env Rscript
#------------------------------------------#
# author:    Yangfan Zhou
# brief:     按外部映射文件把 cluster 列重编码为 cell_type 列，另存新rds。
#            部署位置: ~/Data/Bin/r_scripts/sc_004.Marker_gene_check/sc_004.04.recode_celltype.R
#------------------------------------------#

suppressPackageStartupMessages(library(Seurat))

usage <- function() {
  cat("
用法: Rscript sc_004.04.recode_celltype.R rds=<path> clu=<clusterCol> map=<mapfile> [选项]

必填:
  rds=          Seurat .rds 路径
  clu=          簇号所在的 metadata 列名 (如 seurat_clusters)
  map=          簇号->命名的映射文件路径 (TSV, 见下方格式与示例)

可选:
  celltype_col= 最终写入的cell type列名 [默认: cell_type]；若该列已存在于metadata中,
                不会静默覆盖——会直接报错退出,除非同时给 overwrite=true
  overwrite=    true 时允许覆盖已存在的 celltype_col= 列 [默认: false]
  out=          输出rds路径 [默认: 在输入文件名后加 _labeled，同目录]

映射文件格式 (map=)：
  制表符(\\t)分隔，第一行是表头(会被跳过)，前两列必填：
    簇        命名              依据
    1         pvFB-Inflam-Act   THBS2/FBLN5/COL12A1/CXCL12
    2         pvFB-pMyo-Act     CDH19/CSMD3/TIAM1
    19        Melanocyte        TYR/MITF/PMEL/MLANA
  第三列(依据)仅供人工阅读/存档，脚本不解析、不校验，可以留空或省略整列。
  数据里存在、但没出现在 map= 里的簇号，会统一归到 'Unassigned'（不报错，会打印提醒）。
  map= 里出现、但数据里根本没有的簇号，只警告，不中止（大概率是打错簇号，请核对）。

示例：
  Rscript sc_004.04.recode_celltype.R \\
    rds=subset_vasc_cca_d30_res1.1.rds clu=seurat_clusters map=vasc_celltype_map.tsv

  Rscript sc_004.04.recode_celltype.R \\
    rds=astro_R2_final.rds clu=harmony_clusters map=astro_map.tsv \\
    celltype_col=cell_type_v2 out=astro_R2_final_v2.rds

  Rscript sc_004.04.recode_celltype.R \\
    rds=x.rds clu=seurat_clusters map=x_map.tsv overwrite=true   # 允许覆盖已存在的cell_type列
\n")
}

# ---------------- 参数解析 (key=value，全部 [[ ]] 精确取值，不用 $) ----------------
args_raw <- commandArgs(trailingOnly = TRUE)
if (length(args_raw) == 0) { usage(); quit(save = "no", status = 0) }
if (any(args_raw %in% c("-h", "--help"))) { usage(); quit(save = "no", status = 0) }

opt <- list()
for (a in args_raw) {
  kv <- regmatches(a, regexpr("=", a), invert = TRUE)[[1]]
  if (length(kv) != 2) {
    stop(paste0("错误：参数必须是 key=value 形式，无法解析: ", a), call. = FALSE)
  }
  opt[[kv[1]]] <- kv[2]
}

if (is.null(opt[["rds"]]) || is.null(opt[["clu"]]) || is.null(opt[["map"]])) {
  cat("\n错误：rds=、clu=、map= 均为必填参数\n")
  usage()
  quit(save = "no", status = 1)
}

celltype_col <- if (!is.null(opt[["celltype_col"]])) opt[["celltype_col"]] else "cell_type"
overwrite    <- !is.null(opt[["overwrite"]]) && tolower(opt[["overwrite"]]) == "true"
out_path     <- if (!is.null(opt[["out"]])) opt[["out"]] else sub("\\.rds$", "_labeled.rds", opt[["rds"]])

if (!file.exists(opt[["rds"]])) stop(paste0("错误：rds 文件不存在: ", opt[["rds"]]), call. = FALSE)
if (!file.exists(opt[["map"]])) stop(paste0("错误：map 文件不存在: ", opt[["map"]]), call. = FALSE)

cat("==============================================================================\n")
cat("cluster -> cell_type 重编码 (sc_004.04.recode_celltype.R)\n")
cat("开始时间:", as.character(Sys.time()), "\n")
cat("==============================================================================\n\n")

# ---------------- 加载数据 ----------------
cat("步骤 1/5: 加载数据...\n")
obj <- readRDS(opt[["rds"]])
cat("  细胞数:", ncol(obj), "\n")

if (!opt[["clu"]] %in% colnames(obj@meta.data)) {
  stop(paste0("错误：metadata 中不存在列 clu='", opt[["clu"]], "'"), call. = FALSE)
}

if (celltype_col %in% colnames(obj@meta.data) && !overwrite) {
  stop(paste0(
    "错误：metadata 中已存在列 celltype_col='", celltype_col, "'。\n",
    "      如果确认要覆盖，请加上 overwrite=true 重新运行；否则换一个 celltype_col= 名字。"
  ), call. = FALSE)
}
cat("  簇列:", opt[["clu"]], "| 目标cell type列:", celltype_col,
    if (overwrite && celltype_col %in% colnames(obj@meta.data)) "(将覆盖已有列)" else "", "\n\n")

# ---------------- 读取并校验映射文件 ----------------
cat("步骤 2/5: 读取映射文件", opt[["map"]], "...\n")
map_lines <- readLines(opt[["map"]], warn = FALSE)
map_lines <- map_lines[nzchar(trimws(map_lines))]  # 去空行
if (length(map_lines) < 2) {
  stop("错误：map 文件至少需要一行表头 + 一行数据", call. = FALSE)
}

map_split <- strsplit(map_lines, "\t")
n_cols <- lengths(map_split)
if (any(n_cols < 2)) {
  bad <- which(n_cols < 2)
  stop(paste0("错误：map 文件第 ", paste(bad, collapse = ","),
              " 行少于2列（必须制表符\\t分隔，且至少有 簇/命名 两列）"), call. = FALSE)
}

# 第一行当表头跳过（约定：map= 文件第一行永远是表头，即使没有也建议加上占位表头）
header <- map_split[[1]]
data_rows <- map_split[-1]

cluster_ids <- trimws(vapply(data_rows, `[`, character(1), 1))
cell_names  <- trimws(vapply(data_rows, `[`, character(1), 2))

if (any(cluster_ids == "") || any(cell_names == "")) {
  stop("错误：map 文件的第1列(簇)或第2列(命名)存在空值，请检查", call. = FALSE)
}
if (any(duplicated(cluster_ids))) {
  dup <- unique(cluster_ids[duplicated(cluster_ids)])
  stop(paste0("错误：map 文件里簇号重复出现: ", paste(dup, collapse = ", ")), call. = FALSE)
}

celltype_map <- setNames(cell_names, cluster_ids)
cat("  读取到", length(celltype_map), "条映射 (表头:", paste(header[1:2], collapse = " | "), "...)\n")
for (i in seq_along(celltype_map)) {
  cat("    -", names(celltype_map)[i], "->", celltype_map[i], "\n")
}
cat("\n")

# ---------------- 校验映射与数据的簇号是否对得上 ----------------
cat("步骤 3/5: 核对映射表与数据里的簇号...\n")
data_clusters <- as.character(unique(obj@meta.data[[opt[["clu"]]]]))

map_only <- setdiff(names(celltype_map), data_clusters)
if (length(map_only) > 0) {
  cat("  ⚠ 警告：map 文件里以下簇号在数据中不存在(大概率打错簇号，请核对): ",
      paste(map_only, collapse = ", "), "\n")
}

data_only <- setdiff(data_clusters, names(celltype_map))
if (length(data_only) > 0) {
  cat("  ⚠ 数据里以下簇号未出现在 map 文件中，将统一归为 'Unassigned': ",
      paste(data_only, collapse = ", "), "\n")
} else {
  cat("  ✓ 数据里的全部簇号都在 map 文件里找到了对应命名\n")
}
cat("\n")

# ---------------- 赋值 + 设定factor顺序 ----------------
cat("步骤 4/5: 写入", celltype_col, "列...\n")
clu_chr <- as.character(obj@meta.data[[opt[["clu"]]]])
new_col <- unname(celltype_map[clu_chr])
new_col[is.na(new_col)] <- "Unassigned"

level_order <- unique(celltype_map)  # 保持map文件里出现的顺序（同名簇合并、按首次出现去重）
if ("Unassigned" %in% new_col) level_order <- c(level_order, "Unassigned")

obj@meta.data[[celltype_col]] <- factor(new_col, levels = level_order)

cat("\n[", celltype_col, "] 各类型细胞数:\n", sep = "")
print(table(obj@meta.data[[celltype_col]], useNA = "ifany"))
cat("\n")

# ---------------- 保存 ----------------
cat("步骤 5/5: 保存...\n")
if (file.exists(out_path)) {
  old_mtime <- format(file.info(out_path)$mtime, "%Y-%m-%d %H:%M:%S")
  cat("  ⚠ 输出文件已存在(上次写入时间:", old_mtime, ")，本次运行将覆盖它\n")
}
saveRDS(obj, file = out_path)
cat("  ✓ 已保存:", out_path, "\n")

cat("\n==============================================================================\n")
cat("完成！结束时间:", as.character(Sys.time()), "\n")
cat("==============================================================================\n")
