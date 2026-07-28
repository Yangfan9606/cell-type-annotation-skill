#!/usr/bin/env Rscript
#------------------------------------------#
# author:    Yangfan Zhou
# email:     yangfan.zhou@ki.se
# date:      2026-05-26
# version:   2.4
# license:   MIT
# brief:     单细胞数据细胞类型注释 v2.3
#            - 双方法: AddModuleScore (主决策) + UCell (参考对照)
#            - 自适应 margin 阈值: max(MAD-based, permutation-95%, 0.1)
#            - 输出完整诊断日志: top1_z, top2_z, margin, 阈值来源
#            - 置换检验支持分层子采样, 大数据集加速 6-7x
#
# Changelog v2.4 (2026-06-09):
#   - New: --presence_floor (默认 0 = 关闭, 结果与旧版完全一致). 解决 z-score 陷阱:
#          按列 z-score 会给【全样本都不表达】的 marker panel 也凑出一个 argmax
#          假赢家 (例: 纯髓系 subset 里 T/NK/B panel 全 0, 却被标到某 cluster)。
#          开启后要求每种 cell type 至少一个 marker 在某 cluster 平均表达 >= floor,
#          否则判该谱系在本数据中不存在并丢弃。建议值 0.3-0.5 (基于归一化 data 层)。
#          与 min_markers 区别: min_markers 查"基因在不在", presence_floor 查"表不表达"。
#
# Changelog v2.3 (2026-05-26):
#   - New: 置换检验加入分层子采样 (默认开启) 加速大数据集
#          --permutation_subsample_cap   每 cluster 最多采样数 (默认 2000)
#          --permutation_subsample_floor 小 cluster 警告阈值  (默认 50)
#          预期加速: 36 万细胞 / 29 cluster 从 ~15min 降到 ~2min (50 次置换)
#   - Note: 观测 margin 仍用全部细胞计算, 子采样仅影响零分布估计.
#          策略与 scTRIP/DESeq2/MAGMA 等工具一致 (calibration steps 用子样本).
#
# Changelog v2.2 (2026-05-26):
#   - Fix: 移除 StoreRankings_UCell(verbose=FALSE) 的非法参数 (UCell 不同版本
#          API 不一致, 部分版本不接受 verbose), 改用 suppressMessages 兜底
#   - Fix: ScoreSignatures_UCell 同样包 suppressMessages, 避免冗余输出
#   - Change: --min_markers 默认值 3 -> 1 (允许保留小 marker 列表如 2 基因的
#            Vasc_fibro), 但 < 3 时打运行时警告
#
# Changelog v2.1 (2026-05-26):
#   - Fix: UCell 矩阵提取增加空矩阵检测 (原 tryCatch 只接错误, 不接空返回,
#          导致 SCT assay 的空 counts 槽不会回退到 data, 直接报错)
#   - New: 新增 --ucell_assay 参数 (默认 RNA), UCell 优先用 RNA 原始 counts,
#          失败时自动回退到 --assay 指定的 assay (向后兼容)
#   - Fix: 抽取出 get_expression_matrix() 辅助函数, 统一处理 Assay5 多 layer
#          的 JoinLayers 逻辑, 错误信息更清晰
#------------------------------------------#

# ----------------------------------------------------------------------------
# 第一步: 仅加载 optparse 用于命令行参数解析 (加速 --help 响应)
# ----------------------------------------------------------------------------
suppressPackageStartupMessages(library(optparse))

# ----------------------------------------------------------------------------
# 命令行参数定义
# ----------------------------------------------------------------------------
option_list <- list(
  make_option(c("-i", "--input"),
              type = "character", default = NULL,
              help = "输入的 Seurat RDS 文件路径 [必需]",
              metavar = "FILE"),

  make_option(c("-m", "--markers"),
              type = "character", default = NULL,
              help = "Marker 基因列表文件路径 [必需]
                      文件格式: 每行一种细胞类型
                      第一列: 细胞类型名称
                      后续列: 该类型的 Marker 基因
                      支持分隔符: 逗号(,) 分号(;) 制表符(\\t)",
              metavar = "FILE"),

  make_option(c("-o", "--outdir"),
              type = "character", default = ".",
              help = "输出目录 [默认: 当前目录]",
              metavar = "DIR"),

  make_option(c("--prefix"),
              type = "character", default = "output",
              help = "输出文件前缀 [默认: %default]",
              metavar = "STRING"),

  make_option(c("-c", "--cluster_column"),
              type = "character", default = NULL,
              help = "Seurat 对象 meta.data 中的分群列名 [默认: 自动检测]
                      自动检测顺序: seurat_clusters, integrated_snn_res.*, RNA_snn_res.*",
              metavar = "STRING"),

  make_option(c("-a", "--assay"),
              type = "character", default = "RNA",
              help = "用于 AddModuleScore 的 Assay 名称 [默认: %default]",
              metavar = "STRING"),

  make_option("--ucell_assay",
              type = "character", default = "RNA",
              help = "用于 UCell 的 Assay 名称 (UCell 推荐用 RNA 原始 counts);
                      若该 assay 不可用, 自动回退到 --assay [默认: %default]",
              metavar = "STRING"),

  make_option(c("--reduction"),
              type = "character", default = "umap",
              help = "用于可视化的降维方法 [默认: %default]",
              metavar = "STRING"),

  make_option(c("--subset_types"),
              type = "character", default = NULL,
              help = "需要提取的细胞类型 (逗号分隔) [可选]",
              metavar = "STRING"),

  make_option(c("--sep"),
              type = "character", default = "auto",
              help = "Marker 文件分隔符: auto/comma/semicolon/tab/
                      mixed_tab_comma/mixed_tab_semicolon [默认: %default]",
              metavar = "STRING"),

  make_option(c("--n_permutations"),
              type = "integer", default = 500,
              help = "置换零分布的重复次数 [默认: %default]
                      减小可加速, 增大可提升零分布稳定性",
              metavar = "INTEGER"),

  make_option(c("--skip_permutation"),
              type = "logical", default = FALSE, action = "store_true",
              help = "跳过置换零分布计算 (加快运行, 但只用 MAD 阈值, 回退到 0.2)
                      [默认: FALSE]"),

  make_option(c("--permutation_subsample_cap"),
              type = "integer", default = 2000,
              help = "置换检验时, 每个 cluster 最多采样的细胞数 [默认: %default]
                      策略: 按 cluster 分层随机采样, 大 cluster 截断到此值
                      原因: 置换检验估的是零分布, cluster mean 在 ~2000 cell 已收敛,
                            更多细胞收益递减但耗时线性增加
                      重要: 观测 margin (主结果) 仍用全部细胞, 子采样仅影响阈值估计
                      设为 0 关闭子采样 (用全部细胞, 大数据集会很慢)
                      建议: 中等数据 (<10万 cell) 用 1000-2000; 大数据 (>30万) 用 2000-3000",
              metavar = "INTEGER"),

  make_option(c("--permutation_subsample_floor"),
              type = "integer", default = 50,
              help = "置换检验时, cluster 细胞数低于此值会打 warning [默认: %default]
                      逻辑: < floor 的小 cluster 会全保留 (不采样), 但 cluster mean 估计
                            本身就不稳, 该 cluster 在零分布里的贡献可能噪声较大
                      影响: 对最终 95%% 分位影响很小 (因为是 n_perm x n_cluster 个点的合集)
                            但若有大量 < floor 的小 cluster, 建议人工检查注释结果",
              metavar = "INTEGER"),

  make_option(c("--min_markers"),
              type = "integer", default = 1,
              help = "每种细胞类型至少需要的可用 marker 基因数 [默认: %default]
                      低于此值的细胞类型将跳过
                      注: 设为 1 会保留所有非空 marker 列表, 但单基因得分稳定性差,
                      建议 >= 2 (生产环境推荐 3)",
              metavar = "INTEGER"),

  make_option(c("--presence_floor"),
              type = "double", default = 0.0,
              help = "存在性阈值: 每种细胞类型至少一个 marker 在某 cluster 平均表达
                      >= 此值, 否则判该谱系在本数据中不存在并丢弃 [默认: %default = 关闭]
                      目的: 防止按列 z-score 给全样本零表达的 panel 凑出 argmax 假赢家
                      与 min_markers 区别: min_markers 查基因是否存在, 本参数查是否表达
                      建议: 开启时用 0.3-0.5 (基于归一化 data 层); 0 保持旧行为不变",
              metavar = "FLOAT"),

  make_option(c("--seed"),
              type = "integer", default = 1,
              help = "随机种子, 用于 UCell 和置换检验 [默认: %default]",
              metavar = "INTEGER"),

  make_option(c("--save_percell"),
              type = "logical", default = TRUE, action = "store_true",
              help = "同时输出 per-cell UCell 得分 CSV [默认: %default]"),

  make_option(c("--verbose"),
              type = "logical", default = TRUE, action = "store_true",
              help = "显示详细运行信息 [默认: %default]"),

  make_option(c("--debug"),
              type = "logical", default = FALSE, action = "store_true",
              help = "调试模式: 显示更详细诊断 [默认: %default]")
)

opt_parser <- OptionParser(
  option_list = option_list,
  description = "\n单细胞 RNA-seq 数据细胞类型注释工具 v2.0",
  epilogue = "
示例用法:
  # 基本用法
  Rscript sc_004_CellType_Annotation_v2.R -i seurat.rds -m markers.txt -o results --prefix sample1

  # 跳过置换 (急用)
  Rscript sc_004_CellType_Annotation_v2.R -i seurat.rds -m markers.txt --skip_permutation

  # 减少置换次数加速
  Rscript sc_004_CellType_Annotation_v2.R -i seurat.rds -m markers.txt --n_permutations 100

Marker 文件格式示例:
  T cell,CD3D,CD3E,CD8A,CD4
  B cell,CD19,MS4A1,CD79A
  Monocyte,CD14,FCGR3A,CD68
"
)
opt <- parse_args(opt_parser)

# ----------------------------------------------------------------------------
# 参数验证
# ----------------------------------------------------------------------------
if (is.null(opt$input))   { cat("错误: 缺少 --input\n");   print_help(opt_parser); quit(status = 1) }
if (is.null(opt$markers)) { cat("错误: 缺少 --markers\n"); print_help(opt_parser); quit(status = 1) }
if (!file.exists(opt$input))   { cat(sprintf("错误: 输入文件不存在: %s\n", opt$input));   quit(status = 1) }
if (!file.exists(opt$markers)) { cat(sprintf("错误: Marker 文件不存在: %s\n", opt$markers)); quit(status = 1) }
if (opt$n_permutations < 50 && !opt$skip_permutation) {
  cat(sprintf("警告: --n_permutations=%d 偏少, 零分布估计可能不稳, 建议 >= 200\n", opt$n_permutations))
}
if (opt$min_markers < 1) { cat("错误: --min_markers 必须 >= 1\n"); quit(status = 1) }

# ----------------------------------------------------------------------------
# 加载重型依赖 (放在参数验证之后, 加速 --help)
# ----------------------------------------------------------------------------
if (opt$verbose) {
  cat("==================================================\n")
  cat("单细胞 RNA-seq 细胞类型注释工具 v2.0\n")
  cat("==================================================\n\n")
  cat("正在加载 R 包...\n")
}
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(UCell)
})
if (opt$verbose) cat("R 包加载完成\n\n")

if (!dir.exists(opt$outdir)) {
  dir.create(opt$outdir, recursive = TRUE)
  if (opt$verbose) cat(sprintf("创建输出目录: %s\n\n", opt$outdir))
}

# ============================================================================
# 辅助函数
# ============================================================================

# ---- Marker 文件分隔符检测 (沿用 sc_004 v1 逻辑) ----
detect_separator <- function(file_path, sep_option = "auto") {
  if (sep_option != "auto") {
    sep_map <- list("comma" = ",", "semicolon" = ";", "tab" = "\t",
                    "mixed_tab_comma" = "mixed_tab_comma",
                    "mixed_tab_semicolon" = "mixed_tab_semicolon")
    return(sep_map[[sep_option]])
  }
  first_lines <- readLines(file_path, n = 5)
  first_lines <- first_lines[!grepl("^#|^\\s*$", first_lines)]
  if (length(first_lines) == 0) stop("Marker 文件中没有有效数据行")

  comma_count     <- sum(stringr::str_count(first_lines, ","))
  semicolon_count <- sum(stringr::str_count(first_lines, ";"))
  tab_count       <- sum(stringr::str_count(first_lines, "\t"))

  if (tab_count > 0) {
    lines_with_tab <- sum(grepl("\t", first_lines))
    if (lines_with_tab == length(first_lines)) {
      after_first_tab <- sapply(first_lines, function(line) {
        parts <- strsplit(line, "\t")[[1]]
        if (length(parts) > 1) return(paste(parts[-1], collapse = "\t"))
        return("")
      })
      comma_in_genes <- sum(stringr::str_count(after_first_tab, ","))
      semicolon_in_genes <- sum(stringr::str_count(after_first_tab, ";"))
      if (comma_in_genes > tab_count)     return("mixed_tab_comma")
      if (semicolon_in_genes > tab_count) return("mixed_tab_semicolon")
    }
  }
  if (comma_count >= semicolon_count && comma_count >= tab_count) return(",")
  if (semicolon_count >= comma_count && semicolon_count >= tab_count) return(";")
  "\t"
}

# ---- 读取 Marker 列表 ----
read_marker_list <- function(file_path, sep_info, verbose = TRUE) {
  lines <- readLines(file_path)
  lines <- lines[!grepl("^#|^\\s*$", lines)]
  marker_list <- list()

  for (line in lines) {
    if (sep_info %in% c("mixed_tab_comma", "mixed_tab_semicolon")) {
      inner_sep <- if (sep_info == "mixed_tab_comma") "," else ";"
      parts <- strsplit(line, "\t")[[1]]
      if (length(parts) >= 2) {
        cell_type <- trimws(parts[1])
        genes <- trimws(unlist(strsplit(parts[2], inner_sep)))
        genes <- genes[nchar(genes) > 0]
        marker_list[[cell_type]] <- genes
      }
    } else {
      parts <- trimws(unlist(strsplit(line, sep_info)))
      parts <- parts[nchar(parts) > 0]
      if (length(parts) >= 2) marker_list[[parts[1]]] <- parts[-1]
    }
  }
  if (verbose) {
    cat(sprintf("从 Marker 文件读取到 %d 种细胞类型\n", length(marker_list)))
    for (ct in names(marker_list)) {
      cat(sprintf("  [%s] %d 基因: %s\n", ct, length(marker_list[[ct]]),
                  paste(marker_list[[ct]], collapse = ", ")))
    }
    cat("\n")
  }
  marker_list
}

# ---- 过滤掉数据中不存在或基因数过少的细胞类型 ----
filter_marker_list <- function(marker_list, available_genes, min_markers, verbose = TRUE) {
  filtered <- list()
  skipped  <- character(0)
  for (ct in names(marker_list)) {
    avail <- intersect(marker_list[[ct]], available_genes)
    if (length(avail) >= min_markers) {
      filtered[[ct]] <- avail
      if (verbose) cat(sprintf("  [%s] 保留 %d/%d 基因\n",
                               ct, length(avail), length(marker_list[[ct]])))
    } else {
      skipped <- c(skipped, ct)
      if (verbose) cat(sprintf("  [%s] 跳过 (仅 %d 基因, 低于 min_markers=%d)\n",
                               ct, length(avail), min_markers))
    }
  }
  if (length(filtered) < 2)
    stop(sprintf("可用细胞类型不足 2 种 (仅 %d), 无法进行 argmax 注释", length(filtered)))
  if (verbose) cat(sprintf("\n最终保留 %d 种细胞类型%s\n\n",
                           length(filtered),
                           if (length(skipped) > 0) sprintf(" (跳过: %s)", paste(skipped, collapse = ", ")) else ""))
  filtered
}

# ---- 存在性过滤: 丢弃"全样本几乎不表达"的 marker panel (默认关闭) ----
# 动机: 按列 z-score 会给一个【在本数据里根本不表达】的 panel (如纯髓系 subset
#       里的 T/NK/B) 也凑出 argmax 假赢家 -> 微胶质被误判成淋巴。min_markers 只查
#       基因在不在 (T/NK/B 基因都"在", 只是不表达), 挡不住; 故另设存在性闸。
# 判据: 每种 type 的"最强 marker 在最强 cluster 的平均表达" >= floor 才保留。
#       用归一化 data 层 (无则回退 counts)。floor<=0 时直接返回 (旧行为, 结果不变)。
presence_filter_marker_list <- function(marker_list, seurat_obj, cluster_col,
                                        assay_name, floor, verbose = TRUE) {
  if (is.null(floor) || floor <= 0) return(marker_list)
  expr <- tryCatch(GetAssayData(seurat_obj, assay = assay_name, layer = "data"),
                   error = function(e) NULL)
  if (is.null(expr) || nrow(expr) == 0)
    expr <- GetAssayData(seurat_obj, assay = assay_name, layer = "counts")
  cl <- factor(as.character(seurat_obj@meta.data[[cluster_col]]))
  presence <- sapply(marker_list, function(genes) {
    g <- intersect(genes, rownames(expr))
    if (length(g) == 0) return(0)
    gm <- sapply(levels(cl), function(k) {
      sub <- expr[g, cl == k, drop = FALSE]
      if (ncol(sub) == 0) rep(0, length(g)) else Matrix::rowMeans(sub)
    })
    max(gm)
  })
  kept    <- names(marker_list)[presence >= floor]
  dropped <- names(marker_list)[presence <  floor]
  if (verbose) {
    cat(sprintf("存在性过滤 (--presence_floor=%.2f; max marker x cluster 平均表达):\n", floor))
    for (ct in names(sort(presence, decreasing = TRUE)))
      cat(sprintf("  [%s] presence=%.2f%s\n", ct, presence[ct],
                  if (ct %in% dropped) "  -> 丢弃 (该谱系在本数据中不存在)" else ""))
    cat("\n")
  }
  if (length(kept) < 2)
    stop(sprintf("存在性过滤后可用类型不足 2 种 (仅 %d), 请降低 --presence_floor", length(kept)))
  marker_list[kept]
}

# ---- 按列做 z-score; 零方差列置 0 ----
col_zscore <- function(mat) {
  z <- scale(mat)
  z[is.nan(z)] <- 0
  attr(z, "scaled:center") <- attr(z, "scaled:scale") <- NULL
  z
}

# ---- 按 cluster 聚合得分 (得到 cluster x cell_type 矩阵) ----
aggregate_by_cluster <- function(score_mat_cells, cluster_vec) {
  # score_mat_cells: 行=cell, 列=cell_type
  # cluster_vec:     length = ncell, 与 score_mat_cells 行对应
  df <- as.data.frame(score_mat_cells, check.names = FALSE)
  df$.cluster <- as.character(cluster_vec)
  out <- df %>%
    group_by(.cluster) %>%
    summarise(across(everything(), mean), .groups = "drop") %>%
    column_to_rownames(".cluster") %>%
    as.matrix()
  out
}

# ---- 对一个 z-score 矩阵跑 argmax + margin 过滤, 返回 data.frame ----
argmax_call <- function(zmat, margin_threshold, label) {
  out <- do.call(rbind, lapply(seq_len(nrow(zmat)), function(i) {
    v <- zmat[i, ]
    o <- order(v, decreasing = TRUE)
    top1 <- v[o[1]]
    top2 <- if (length(v) >= 2) v[o[2]] else NA_real_
    marg <- top1 - ifelse(is.na(top2), -Inf, top2)
    if (top1 < 0) {
      call_label <- "Unknown"
      reason     <- "top1_z<0"
    } else if (marg < margin_threshold) {
      call_label <- "Unknown"
      reason     <- sprintf("margin<%.3f", margin_threshold)
    } else {
      call_label <- names(v)[o[1]]
      reason     <- ""
    }
    data.frame(
      cluster       = rownames(zmat)[i],
      call          = call_label,
      top1_id       = names(v)[o[1]],
      top1_z        = round(top1, 3),
      top2_id       = if (length(v) >= 2) names(v)[o[2]] else NA_character_,
      top2_z        = round(top2, 3),
      margin        = round(marg, 3),
      reject_reason = reason,
      stringsAsFactors = FALSE
    )
  }))
  names(out)[-1] <- paste0(label, "_", names(out)[-1])
  out
}

# ---- 计算观测 margins (跨所有 cluster, 不做阈值过滤) ----
extract_margins <- function(zmat) {
  apply(zmat, 1, function(v) {
    o <- order(v, decreasing = TRUE)
    if (length(v) >= 2) v[o[1]] - v[o[2]] else NA_real_
  })
}

# ---- 从 Seurat 对象中尝试获取表达矩阵, 多层回退 ----
# 优先级: 1) preferred_assay 2) fallback_assay 3) 其他所有 assay
# 每个 assay 内部: 先试 counts, 再试 data
# Assay5 多 layer 自动 JoinLayers
# 返回: list(matrix, assay, layer)
get_expression_matrix <- function(seurat_obj, preferred_assay,
                                  fallback_assay = NULL, verbose = TRUE) {
  available <- Assays(seurat_obj)
  if (length(available) == 0)
    stop("Seurat 对象中没有任何 Assay")

  # 构建尝试顺序
  to_try <- character(0)
  if (preferred_assay %in% available) {
    to_try <- preferred_assay
  } else if (verbose) {
    cat(sprintf("  注意: 首选 assay '%s' 不存在 (可用: %s)\n",
                preferred_assay, paste(available, collapse = ", ")))
  }
  if (!is.null(fallback_assay) && fallback_assay %in% available &&
      !fallback_assay %in% to_try) {
    to_try <- c(to_try, fallback_assay)
  }
  # 最后保底: 所有剩余 assay
  for (a in available) if (!a %in% to_try) to_try <- c(to_try, a)

  # 逐个尝试
  for (a in to_try) {
    obj_tmp <- seurat_obj
    if (inherits(obj_tmp[[a]], "Assay5") &&
        length(SeuratObject::Layers(obj_tmp[[a]], search = "counts")) > 1) {
      obj_tmp <- tryCatch(JoinLayers(obj_tmp, assay = a),
                          error = function(e) obj_tmp)
    }
    for (layer in c("counts", "data")) {
      m <- tryCatch(GetAssayData(obj_tmp, assay = a, layer = layer),
                    error = function(e) NULL)
      if (!is.null(m) && nrow(m) > 0 && ncol(m) > 0) {
        if (verbose) {
          cat(sprintf("  -> 使用 assay='%s', layer='%s' (%d 基因 x %d 细胞)\n",
                      a, layer, nrow(m), ncol(m)))
        }
        return(list(matrix = m, assay = a, layer = layer))
      } else if (verbose) {
        cat(sprintf("     skip: assay='%s', layer='%s' (空或不可用)\n", a, layer))
      }
    }
  }

  stop(sprintf("所有 assay 都无法提供有效表达矩阵 (尝试过: %s)",
               paste(to_try, collapse = ", ")))
}

# ---- AddModuleScore: 得到 cluster x cell_type 的原始得分矩阵 + per-cell ----
# 优化: 一次性传入所有 cell type 的 features list, 共用基因 bin, 显著加速
score_mod <- function(seurat_obj, marker_list, cluster_col, assay_name, verbose = TRUE) {
  if (verbose) cat("计算 AddModuleScore...\n")
  DefaultAssay(seurat_obj) <- assay_name
  original_cells <- colnames(seurat_obj)

  # 单次调用: features 为命名 list, AddModuleScore 会为每个元素加 "1", "2", ... 后缀
  tmp_prefix <- "MODtmp"
  seurat_obj <- suppressWarnings(AddModuleScore(
    seurat_obj, features = marker_list,
    name = tmp_prefix, assay = assay_name
  ))
  added_cols <- paste0(tmp_prefix, seq_along(marker_list))
  if (!all(added_cols %in% colnames(seurat_obj@meta.data)))
    stop("AddModuleScore 未返回预期列, 请检查 marker_list")

  percell <- as.matrix(seurat_obj@meta.data[original_cells, added_cols, drop = FALSE])
  colnames(percell) <- names(marker_list)
  # 清理临时列
  for (cl in added_cols) seurat_obj@meta.data[[cl]] <- NULL

  # 聚合到 cluster
  cluster_vec <- as.character(seurat_obj@meta.data[[cluster_col]])
  names(cluster_vec) <- colnames(seurat_obj)
  cluster_mat <- aggregate_by_cluster(percell, cluster_vec)

  list(percell = percell, cluster_raw = cluster_mat, seurat_obj = seurat_obj)
}

# ---- UCell: 得到 cluster x cell_type 的原始得分矩阵 + per-cell + 预存 ranking ----
# 注意: UCell 基于排序, 推荐用 RNA 原始 counts;
#       SCT 的 counts 是校正后值, scale.data 是 Pearson 残差 (有负值, 不适合 UCell)
score_ucell <- function(seurat_obj, marker_list, cluster_col,
                        preferred_assay = "RNA", fallback_assay = NULL,
                        seed = 1, verbose = TRUE) {
  if (verbose) cat(sprintf("计算 UCell (首选 assay='%s', 回退 assay='%s')...\n",
                           preferred_assay,
                           if (is.null(fallback_assay)) "无" else fallback_assay))

  # 通过统一辅助函数取表达矩阵
  expr_info <- get_expression_matrix(seurat_obj, preferred_assay,
                                     fallback_assay = fallback_assay,
                                     verbose = verbose)
  expr <- expr_info$matrix
  used_assay <- expr_info$assay
  used_layer <- expr_info$layer

  # 预存 ranking: 置换检验复用, 大幅提速
  # 注: StoreRankings_UCell 不同版本 API 不一致, 不传 verbose 用 suppressMessages 兜底
  set.seed(seed)
  if (verbose) cat("  预计算 UCell ranking (这一步耗时, 但只算一次)...\n")
  ranks <- suppressMessages(StoreRankings_UCell(expr))

  # 计算 per-cell 得分
  if (verbose) cat("  计算 UCell 原始得分...\n")
  set.seed(seed)
  ucell_mat <- suppressMessages(ScoreSignatures_UCell(
    matrix = NULL, features = marker_list,
    precalc.ranks = ranks, name = "_UCell"
  ))
  # UCell 列名是 "{name}_UCell", 去掉后缀对齐
  colnames(ucell_mat) <- sub("_UCell$", "", colnames(ucell_mat))
  ucell_mat <- ucell_mat[colnames(seurat_obj), names(marker_list), drop = FALSE]

  # 聚合到 cluster
  cluster_vec <- as.character(seurat_obj@meta.data[[cluster_col]])
  names(cluster_vec) <- colnames(seurat_obj)
  cluster_mat <- aggregate_by_cluster(ucell_mat, cluster_vec)

  list(percell = ucell_mat, cluster_raw = cluster_mat, ranks = ranks,
       used_assay = used_assay, used_layer = used_layer)
}

# ---- 分层子采样: 用于置换检验加速 ----
# 策略:
#   - 大 cluster (n > cap): 随机采样 cap 个细胞
#   - 中 cluster (floor <= n <= cap): 全部保留
#   - 小 cluster (n < floor): 全部保留, 但打 warning
#   - 返回采样的细胞名 + 统计信息
# 科学性:
#   - cluster mean 的 SE = sigma/sqrt(n), 在 n=2000 时已经 sigma/45,
#     边际收益极小. 不引入偏差 (无偏估计), 只略增方差.
#   - 主结果 (observed margin) 不受影响, 仅置换零分布用此采样.
stratified_subsample <- function(cluster_vec, cap, floor, seed = 1, verbose = TRUE) {
  set.seed(seed)
  cluster_sizes <- table(cluster_vec)
  sampled_cells   <- character(0)
  small_clusters  <- character(0)  # n < floor
  capped_clusters <- character(0)  # n > cap (被截断)

  for (cl in names(cluster_sizes)) {
    cells_in_cl <- names(cluster_vec)[cluster_vec == cl]
    n_cl <- length(cells_in_cl)
    if (n_cl < floor) {
      small_clusters <- c(small_clusters, sprintf("%s(n=%d)", cl, n_cl))
      sampled_cells <- c(sampled_cells, cells_in_cl)
    } else if (n_cl > cap) {
      capped_clusters <- c(capped_clusters, sprintf("%s(%d->%d)", cl, n_cl, cap))
      sampled_cells <- c(sampled_cells, sample(cells_in_cl, cap))
    } else {
      sampled_cells <- c(sampled_cells, cells_in_cl)
    }
  }

  if (verbose) {
    pct <- 100 * length(sampled_cells) / length(cluster_vec)
    cat(sprintf("  分层子采样: %d / %d 细胞 (%.1f%%), cap=%d floor=%d\n",
                length(sampled_cells), length(cluster_vec), pct, cap, floor))
    if (length(capped_clusters) > 0)
      cat(sprintf("    %d 个 cluster 被 cap: %s\n",
                  length(capped_clusters),
                  paste(capped_clusters, collapse = ", ")))
    if (length(small_clusters) > 0)
      cat(sprintf("    ⚠ %d 个 cluster 低于 floor=%d: %s\n",
                  length(small_clusters), floor,
                  paste(small_clusters, collapse = ", ")))
  }

  list(
    sampled_cells   = sampled_cells,
    n_sampled       = length(sampled_cells),
    n_total         = length(cluster_vec),
    small_clusters  = small_clusters,
    capped_clusters = capped_clusters
  )
}

# ---- 置换检验: 打乱基因到细胞类型的映射, 收集 margin 零分布 ----
# 设计:
#   pool = 所有 cell type 的 marker 基因合集 (去重)
#   每次置换: 每种 cell type 从 pool 独立无放回抽样 (size = 原始基因数)
#   v2.3+: 支持分层子采样, 大数据集时仅在子采样上做置换, 6-7x 加速
#     - 观测 margin 仍用全部细胞 (主流程已计算, 不在此函数)
#     - 此函数返回的是子采样上的零分布 margin
#   MOD 路径优化: 一次 AddModuleScore 调用处理所有 type, 共用 bin
#   UCell 路径: 若子采样, 重新计算 ranks (基于子采样矩阵); 否则复用传入 ranks
run_permutation <- function(method,                # "MOD" 或 "UCell"
                            seurat_obj, marker_list, cluster_col, assay_name,
                            ucell_ranks = NULL,    # 不子采样时复用
                            ucell_preferred_assay = NULL,
                            ucell_fallback_assay  = NULL,
                            n_perm = 500, seed = 1, verbose = TRUE,
                            subsample_cap = 2000, subsample_floor = 50) {

  pool  <- unique(unlist(marker_list, use.names = FALSE))
  sizes <- sapply(marker_list, length)

  if (length(pool) < max(sizes))
    stop(sprintf("基因池过小 (%d) 小于最大 type 基因数 (%d), 无法置换",
                 length(pool), max(sizes)))
  if (verbose) cat(sprintf("  基因池大小: %d (各 type 原始大小: %s)\n",
                           length(pool),
                           paste(sizes, collapse = ", ")))

  # ---- 子采样 (默认开启) ----
  cluster_vec_full <- as.character(seurat_obj@meta.data[[cluster_col]])
  names(cluster_vec_full) <- colnames(seurat_obj)

  subsample_info <- list(enabled = FALSE,
                         n_sampled = length(cluster_vec_full),
                         n_total   = length(cluster_vec_full),
                         small_clusters = character(0),
                         capped_clusters = character(0))

  if (subsample_cap > 0) {
    ss <- stratified_subsample(cluster_vec_full,
                               cap = subsample_cap, floor = subsample_floor,
                               seed = seed, verbose = verbose)
    if (ss$n_sampled < ss$n_total) {
      subsample_info <- modifyList(ss, list(enabled = TRUE))
      # 子采样 seurat 对象 (drop 多余细胞)
      seurat_obj  <- seurat_obj[, ss$sampled_cells]
      cluster_vec <- cluster_vec_full[ss$sampled_cells]

      # UCell: 在子采样矩阵上重算 ranks
      if (method == "UCell") {
        if (is.null(ucell_preferred_assay))
          stop("子采样模式下 UCell 置换需要 ucell_preferred_assay 参数")
        if (verbose) cat("  为置换重算 UCell ranking (基于子采样矩阵)...\n")
        expr_info <- get_expression_matrix(seurat_obj,
                                           preferred_assay = ucell_preferred_assay,
                                           fallback_assay  = ucell_fallback_assay,
                                           verbose = FALSE)
        set.seed(seed)
        ucell_ranks <- suppressMessages(StoreRankings_UCell(expr_info$matrix))
      }
    } else {
      cluster_vec <- cluster_vec_full
      if (verbose) cat("  (子采样上限大于实际细胞数, 不做截断, 直接用全部细胞)\n")
    }
  } else {
    cluster_vec <- cluster_vec_full
    if (verbose) cat("  子采样已关闭 (--permutation_subsample_cap=0), 用全部细胞\n")
  }

  unique_clusters <- sort(unique(cluster_vec))

  # 输出: n_perm x n_cluster 的 margin 矩阵
  null_margins <- matrix(NA_real_, nrow = n_perm, ncol = length(unique_clusters),
                         dimnames = list(NULL, unique_clusters))

  if (verbose) cat(sprintf("  开始置换 (%s, n=%d, 基于 %d 细胞):\n",
                           method, n_perm, length(cluster_vec)))
  pb_every <- max(1, n_perm %/% 20)
  tmp_prefix <- "permMODtmp"

  for (i in seq_len(n_perm)) {
    set.seed(seed + i)
    # 每种 type 独立从池中无放回抽样
    perm_list <- lapply(sizes, function(s) sample(pool, s, replace = FALSE))
    names(perm_list) <- names(sizes)

    if (method == "UCell") {
      set.seed(seed + i)
      pmat <- suppressMessages(ScoreSignatures_UCell(
        matrix = NULL, features = perm_list,
        precalc.ranks = ucell_ranks, name = "_UCell"
      ))
      colnames(pmat) <- sub("_UCell$", "", colnames(pmat))
      pmat <- pmat[colnames(seurat_obj), names(perm_list), drop = FALSE]

    } else if (method == "MOD") {
      # 单次 AddModuleScore 调用, 共用 bin -> 快很多
      obj_tmp <- suppressWarnings(AddModuleScore(
        seurat_obj, features = perm_list,
        name = tmp_prefix, assay = assay_name
      ))
      added_cols <- paste0(tmp_prefix, seq_along(perm_list))
      pmat <- as.matrix(obj_tmp@meta.data[colnames(seurat_obj), added_cols, drop = FALSE])
      colnames(pmat) <- names(perm_list)
      rm(obj_tmp)
    } else stop("未知方法: ", method)

    # 聚合 -> z-score -> margin
    cluster_mat <- aggregate_by_cluster(pmat, cluster_vec)
    z <- col_zscore(cluster_mat)
    null_margins[i, rownames(z)] <- extract_margins(z)

    if (verbose && (i %% pb_every == 0))
      cat(sprintf("    %d / %d (%.0f%%)\n", i, n_perm, 100 * i / n_perm))
  }
  if (verbose) cat(sprintf("  %s 置换完成\n", method))

  # 附加子采样元信息到返回值
  attr(null_margins, "subsample_info") <- subsample_info
  null_margins
}

# ---- 计算最终阈值 ----
compute_threshold <- function(margins_observed, null_margins = NULL,
                              skip_perm = FALSE,
                              hard_floor = 0.1, fallback_skip = 0.2) {
  med <- median(margins_observed, na.rm = TRUE)
  ma  <- mad(margins_observed, na.rm = TRUE)
  mad_th <- med - 3 * ma  # SingleR 风格 outlier 下界

  if (skip_perm) {
    candidates <- c(MAD = mad_th, fallback = fallback_skip, hard_floor = hard_floor)
    final <- max(candidates)
    src   <- names(candidates)[which.max(candidates)]
    perm_th <- NA_real_
  } else {
    perm_th <- as.numeric(quantile(as.vector(null_margins), 0.95, na.rm = TRUE))
    candidates <- c(MAD = mad_th, permutation = perm_th, hard_floor = hard_floor)
    final <- max(candidates)
    src   <- names(candidates)[which.max(candidates)]
  }

  list(
    final  = final,
    source = src,
    mad    = mad_th,
    perm   = perm_th,
    stats  = list(median = med, mad_value = ma,
                  min = min(margins_observed, na.rm = TRUE),
                  max = max(margins_observed, na.rm = TRUE),
                  n   = length(margins_observed))
  )
}

# ---- 小工具: 从日志标题里提取方法标签 ----
title_to_label <- function(title) {
  if (grepl("^MOD", title)) "MOD" else "UCell"
}

# ---- 写阈值日志 ----
write_threshold_log <- function(path, opt, info_mod, info_ucell, calls,
                                ucell_assay_info = NULL,
                                subsample_mod = NULL, subsample_ucell = NULL) {
  con <- file(path, "w")
  on.exit(close(con))

  w <- function(...) writeLines(sprintf(...), con)

  w("==============================================================")
  w("细胞类型注释 - Margin 阈值诊断日志")
  w("==============================================================")
  w("生成时间:        %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  w("输入文件:        %s", opt$input)
  w("Marker 文件:     %s", opt$markers)
  w("输出前缀:        %s", opt$prefix)
  w("MOD assay:       %s (--assay)", opt$assay)
  if (!is.null(ucell_assay_info)) {
    w("UCell assay:     %s, layer=%s (--ucell_assay=%s%s)",
      ucell_assay_info$used_assay,
      ucell_assay_info$used_layer,
      opt$ucell_assay,
      if (ucell_assay_info$used_assay != opt$ucell_assay) " [已回退]" else "")
  }
  w("置换次数:        %s", if (opt$skip_permutation) "跳过" else as.character(opt$n_permutations))
  if (!opt$skip_permutation) {
    if (opt$permutation_subsample_cap > 0) {
      w("置换子采样:      启用 (cap=%d, floor=%d)",
        opt$permutation_subsample_cap, opt$permutation_subsample_floor)
      if (!is.null(subsample_mod) && subsample_mod$enabled) {
        w("  MOD   置换基于:  %d / %d 细胞 (%.1f%%)",
          subsample_mod$n_sampled, subsample_mod$n_total,
          100 * subsample_mod$n_sampled / subsample_mod$n_total)
        if (length(subsample_mod$capped_clusters) > 0)
          w("  MOD   被 cap 的 cluster: %s",
            paste(subsample_mod$capped_clusters, collapse = ", "))
        if (length(subsample_mod$small_clusters) > 0)
          w("  ⚠ MOD 小 cluster (< floor): %s",
            paste(subsample_mod$small_clusters, collapse = ", "))
      }
      if (!is.null(subsample_ucell) && subsample_ucell$enabled) {
        w("  UCell 置换基于:  %d / %d 细胞 (%.1f%%)",
          subsample_ucell$n_sampled, subsample_ucell$n_total,
          100 * subsample_ucell$n_sampled / subsample_ucell$n_total)
      }
    } else {
      w("置换子采样:      已关闭 (--permutation_subsample_cap=0)")
    }
  }
  w("min_markers:     %d", opt$min_markers)
  w("presence_floor:  %s", if (!is.null(opt$presence_floor) && opt$presence_floor > 0)
                             sprintf("%.2f", opt$presence_floor) else "关闭")
  w("hard_floor:      0.1")
  w("")

  for (info_pair in list(list("MOD (AddModuleScore - 主决策)", info_mod, calls$mod),
                         list("UCell (参考)",                   info_ucell, calls$ucell))) {
    title  <- info_pair[[1]]
    info   <- info_pair[[2]]
    callsf <- info_pair[[3]]

    w("==============================================================")
    w("[%s]", title)
    w("==============================================================")
    w("观测 margin 分布 (跨 %d 个 cluster):", info$stats$n)
    w("  median     = %.4f", info$stats$median)
    w("  MAD        = %.4f", info$stats$mad_value)
    w("  min / max  = %.4f / %.4f", info$stats$min, info$stats$max)
    w("")
    w("候选阈值:")
    w("  MAD 自适应 (median - 3*MAD):  %.4f", info$mad)
    if (!is.na(info$perm))
      w("  置换 95%% 分位:                %.4f", info$perm)
    else
      w("  置换 95%% 分位:                跳过")
    w("  hard_floor:                   %.4f", 0.1)
    if (opt$skip_permutation)
      w("  skip_permutation 回退:        %.4f", 0.2)
    w("")
    w("最终阈值:    %.4f", info$final)
    w("来源:        %s", info$source)
    w("")

    # 诊断提示
    if (info$mad < 0) {
      w("⚠ MAD 自适应阈值为负:")
      w("  说明 margin 分布偏对称 / 无明显坏 cluster, 通常是好事 (cluster 区分度高);")
      w("  也可能是 cluster 数太少导致 MAD 估计不稳, 请检查 stats 部分的 n.")
    } else if (info$mad > 1.0) {
      w("⚠ MAD 自适应阈值较高 (>1.0):")
      w("  说明 margin 分布很离散 (有少数 cluster 信号特别强);")
      w("  此时阈值会偏严, 可能导致较多 cluster 被标为 Unknown.")
    }
    if (info$stats$n < 8) {
      w("⚠ Cluster 数过少 (n=%d): MAD 估计本身方差大, 建议参考置换阈值.", info$stats$n)
    }
    w("")

    # 被剔除的 cluster
    rejected <- callsf[callsf[[paste0(title_to_label(title), "_call")]] == "Unknown", ]
    if (nrow(rejected) > 0) {
      w("被标为 Unknown 的 cluster (%d 个):", nrow(rejected))
      for (j in seq_len(nrow(rejected))) {
        row <- rejected[j, ]
        prefix <- title_to_label(title)
        w("  cluster %s -> top1=%s (z=%.3f), top2=%s (z=%.3f), margin=%.3f [%s]",
          row$cluster,
          row[[paste0(prefix, "_top1_id")]], row[[paste0(prefix, "_top1_z")]],
          row[[paste0(prefix, "_top2_id")]], row[[paste0(prefix, "_top2_z")]],
          row[[paste0(prefix, "_margin")]],
          row[[paste0(prefix, "_reject_reason")]])
      }
    } else {
      w("无 cluster 被标为 Unknown")
    }
    w("")
  }
  invisible(NULL)
}

# ---- 诊断图: margin 分布直方图 + 阈值线 ----
plot_diagnostics <- function(margins_mod, margins_ucell, info_mod, info_ucell, save_path) {
  df_mod   <- data.frame(margin = margins_mod,   method = "MOD")
  df_ucell <- data.frame(margin = margins_ucell, method = "UCell")
  df <- rbind(df_mod, df_ucell)

  thr_df <- data.frame(
    method = c("MOD", "UCell"),
    final  = c(info_mod$final, info_ucell$final),
    mad    = c(info_mod$mad,   info_ucell$mad),
    perm   = c(info_mod$perm,  info_ucell$perm)
  )

  perm_layer <- if (any(!is.na(thr_df$perm)))
    geom_vline(data = thr_df, aes(xintercept = perm, color = "perm-95%"),
               linetype = "dotted", linewidth = 0.6, na.rm = TRUE)
  else NULL

  p <- ggplot(df, aes(x = margin, fill = method)) +
    geom_histogram(bins = 15, alpha = 0.6, position = "identity") +
    geom_vline(data = thr_df, aes(xintercept = final, color = "final (used)"),
               linetype = "solid", linewidth = 0.8) +
    geom_vline(data = thr_df, aes(xintercept = mad, color = "MAD"),
               linetype = "dashed", linewidth = 0.6) +
    perm_layer +
    facet_wrap(~ method, ncol = 1, scales = "free_y") +
    scale_color_manual(name = "阈值",
                       values = c("final (used)" = "red",
                                  "MAD" = "darkgreen",
                                  "perm-95%" = "blue")) +
    labs(title = "Cluster margin 分布与阈值",
         x = "margin = top1_z - top2_z",
         y = "Cluster 数") +
    theme_bw() +
    theme(legend.position = "bottom")

  ggsave(save_path, p, width = 8, height = 7, dpi = 200)
}

# ---- UMAP 可视化: cluster vs cell_type 并排 ----
save_umap_plots <- function(seurat_obj, output_path, reduction, cluster_col,
                            celltype_col, title_suffix = "") {
  p1 <- DimPlot(seurat_obj, reduction = reduction, group.by = cluster_col,
                label = TRUE, repel = TRUE) +
    ggtitle(paste0("Cluster", title_suffix)) +
    theme_minimal() + theme(legend.position = "right")
  p2 <- DimPlot(seurat_obj, reduction = reduction, group.by = celltype_col,
                label = TRUE, repel = TRUE) +
    ggtitle(paste0("Cell type (MOD-based)", title_suffix)) +
    theme_minimal() + theme(legend.position = "right")
  ggsave(output_path, p1 + p2, width = 16, height = 7, dpi = 200)
}

# ============================================================================
# 主流程
# ============================================================================
tryCatch({

  # ---- 1. 读 Seurat ----
  if (opt$verbose) cat(sprintf("读取 Seurat 对象: %s\n", opt$input))
  seurat_obj <- readRDS(opt$input)
  if (opt$verbose) {
    cat(sprintf("  细胞数: %d  基因数: %d\n", ncol(seurat_obj), nrow(seurat_obj)))
    cat(sprintf("  Assays:  %s\n", paste(Assays(seurat_obj), collapse = ", ")))
    cat("\n")
  }

  # ---- 2. cluster 列 ----
  if (is.null(opt$cluster_column)) {
    candidates <- c("seurat_clusters",
                    grep("integrated_snn_res", colnames(seurat_obj@meta.data), value = TRUE),
                    grep("RNA_snn_res",        colnames(seurat_obj@meta.data), value = TRUE),
                    grep("liger_clusters",     colnames(seurat_obj@meta.data), value = TRUE))
    cluster_col <- NULL
    for (cc in candidates)
      if (cc %in% colnames(seurat_obj@meta.data)) { cluster_col <- cc; break }
    if (is.null(cluster_col))
      stop("自动检测失败, 请用 --cluster_column 指定")
    if (opt$verbose) cat(sprintf("自动检测到 cluster 列: %s\n\n", cluster_col))
  } else {
    cluster_col <- opt$cluster_column
    if (!cluster_col %in% colnames(seurat_obj@meta.data))
      stop(sprintf("meta.data 中找不到列 '%s'", cluster_col))
  }

  if (!opt$assay %in% Assays(seurat_obj))
    stop(sprintf("未找到 Assay '%s'", opt$assay))

  plot_available <- opt$reduction %in% names(seurat_obj@reductions)
  if (!plot_available)
    warning(sprintf("未找到降维 '%s', 跳过 UMAP 可视化", opt$reduction))

  # ---- 3. 读 marker, 过滤 min_markers ----
  if (opt$verbose) cat("读取 Marker 列表...\n")
  sep_info    <- detect_separator(opt$markers, opt$sep)
  marker_list <- read_marker_list(opt$markers, sep_info, opt$verbose)
  if (opt$verbose) cat(sprintf("基因可用性检查 (min_markers = %d):\n", opt$min_markers))
  if (opt$min_markers < 3 && opt$verbose) {
    cat(sprintf("  ⚠ min_markers=%d 较低: 该值下的细胞类型得分稳定性会下降,\n",
                opt$min_markers))
    cat("    margin 阈值的解释可能受影响, 注意人工核查 thresholds_log.txt\n")
  }
  marker_list <- filter_marker_list(marker_list, rownames(seurat_obj),
                                    opt$min_markers, opt$verbose)

  # 存在性过滤 (默认关闭; --presence_floor>0 时丢弃全样本不表达的 panel)
  marker_list <- presence_filter_marker_list(marker_list, seurat_obj, cluster_col,
                                             opt$assay, opt$presence_floor, opt$verbose)

  # ---- 4. 评分: MOD ----
  if (opt$verbose) cat("==================== 计算原始得分 ====================\n")
  mod_res <- score_mod(seurat_obj, marker_list, cluster_col, opt$assay, opt$verbose)
  seurat_obj <- mod_res$seurat_obj
  mod_cluster_raw <- mod_res$cluster_raw
  mod_z <- col_zscore(mod_cluster_raw)
  if (opt$verbose) cat(sprintf("  MOD z-score 矩阵: %d cluster x %d cell_type\n",
                               nrow(mod_z), ncol(mod_z)))

  # ---- 5. 评分: UCell ----
  # UCell 用独立 assay 选择: 默认 RNA, 失败回退到 --assay
  ucell_res <- score_ucell(seurat_obj, marker_list, cluster_col,
                           preferred_assay = opt$ucell_assay,
                           fallback_assay  = opt$assay,
                           seed = opt$seed, verbose = opt$verbose)
  ucell_cluster_raw <- ucell_res$cluster_raw
  ucell_z <- col_zscore(ucell_cluster_raw)
  if (opt$verbose) cat(sprintf("  UCell z-score 矩阵: %d cluster x %d cell_type (实际使用 assay='%s', layer='%s')\n\n",
                               nrow(ucell_z), ncol(ucell_z),
                               ucell_res$used_assay, ucell_res$used_layer))

  # ---- 6. 置换检验 ----
  if (opt$skip_permutation) {
    if (opt$verbose) cat("==================== 置换检验已跳过 ====================\n\n")
    null_mod   <- NULL
    null_ucell <- NULL
  } else {
    if (opt$verbose) cat("==================== 置换检验 ====================\n")
    null_mod <- run_permutation("MOD", seurat_obj, marker_list, cluster_col, opt$assay,
                                n_perm = opt$n_permutations, seed = opt$seed,
                                verbose = opt$verbose,
                                subsample_cap   = opt$permutation_subsample_cap,
                                subsample_floor = opt$permutation_subsample_floor)
    null_ucell <- run_permutation("UCell", seurat_obj, marker_list, cluster_col, opt$assay,
                                  ucell_ranks = ucell_res$ranks,
                                  ucell_preferred_assay = opt$ucell_assay,
                                  ucell_fallback_assay  = opt$assay,
                                  n_perm = opt$n_permutations, seed = opt$seed,
                                  verbose = opt$verbose,
                                  subsample_cap   = opt$permutation_subsample_cap,
                                  subsample_floor = opt$permutation_subsample_floor)
    if (opt$verbose) cat("\n")
  }

  # ---- 7. 计算阈值 ----
  margins_mod   <- extract_margins(mod_z)
  margins_ucell <- extract_margins(ucell_z)
  info_mod   <- compute_threshold(margins_mod,   null_mod,   opt$skip_permutation)
  info_ucell <- compute_threshold(margins_ucell, null_ucell, opt$skip_permutation)

  if (opt$verbose) {
    cat("==================== 阈值计算结果 ====================\n")
    cat(sprintf("  MOD   最终阈值: %.4f (来源: %s; MAD=%.4f, perm=%s)\n",
                info_mod$final, info_mod$source, info_mod$mad,
                if (is.na(info_mod$perm)) "NA" else sprintf("%.4f", info_mod$perm)))
    cat(sprintf("  UCell 最终阈值: %.4f (来源: %s; MAD=%.4f, perm=%s)\n\n",
                info_ucell$final, info_ucell$source, info_ucell$mad,
                if (is.na(info_ucell$perm)) "NA" else sprintf("%.4f", info_ucell$perm)))
  }

  # ---- 8. 调用 argmax + margin 过滤 ----
  mod_calls   <- argmax_call(mod_z,   info_mod$final,   "MOD")
  ucell_calls <- argmax_call(ucell_z, info_ucell$final, "UCell")

  # ---- 9. 合并: MOD 主, UCell 标记冲突 ----
  res <- mod_calls %>%
    left_join(ucell_calls, by = "cluster") %>%
    mutate(
      mod_ucell_agree = MOD_call == UCell_call,
      cell_type = ifelse(
        mod_ucell_agree | MOD_call == "Unknown",
        MOD_call,
        paste0(MOD_call, "(?)")
      )
    )

  if (opt$verbose) {
    cat("==================== 注释结果 ====================\n")
    print(res[, c("cluster", "MOD_call", "MOD_top1_z", "MOD_margin",
                  "UCell_call", "UCell_top1_z", "UCell_margin",
                  "mod_ucell_agree", "cell_type")])
    cat("\n")

    n_agree    <- sum(res$mod_ucell_agree)
    n_disagree <- sum(!res$mod_ucell_agree)
    n_unknown  <- sum(res$MOD_call == "Unknown")
    cat(sprintf("两法一致:    %d / %d\n", n_agree, nrow(res)))
    cat(sprintf("两法不一致:  %d (cell_type 加 '(?)' 后缀)\n", n_disagree))
    cat(sprintf("MOD Unknown: %d\n\n", n_unknown))
  }

  # ---- 10. 写回 meta.data ----
  # 注意: ct_map[cell_clusters] 用字符向量索引带名向量, 返回值会带上 cluster-id
  #       作为 names. 直接 seurat_obj$col <- 该向量会触发 AddMetaData 把 names
  #       当成 cell barcode 去匹配 Cells(obj), 二者零交集 -> "No cell overlap".
  #       必须把 names 设回细胞 barcode (或 unname 走位置赋值).
  cells <- colnames(seurat_obj)
  ct_map    <- setNames(res$cell_type,  res$cluster)
  mod_map   <- setNames(res$MOD_call,   res$cluster)
  ucell_map <- setNames(res$UCell_call, res$cluster)
  cell_clusters <- as.character(seurat_obj@meta.data[[cluster_col]])

  seurat_obj$cell_type       <- setNames(ct_map[cell_clusters],    cells)
  seurat_obj$cell_type_MOD   <- setNames(mod_map[cell_clusters],   cells)
  seurat_obj$cell_type_UCell <- setNames(ucell_map[cell_clusters], cells)

  # per-cell UCell 得分写回
  for (ct in colnames(ucell_res$percell)) {
    col_nm <- paste0("UCell_", make.names(ct))
    seurat_obj@meta.data[[col_nm]] <- ucell_res$percell[colnames(seurat_obj), ct]
  }
  # per-cell MOD 得分写回
  for (ct in colnames(mod_res$percell)) {
    col_nm <- paste0("MOD_", make.names(ct))
    seurat_obj@meta.data[[col_nm]] <- mod_res$percell[colnames(seurat_obj), ct]
  }

  # ---- 11. 保存所有输出文件 ----
  if (opt$verbose) cat("==================== 写入输出文件 ====================\n")

  f <- function(suffix) file.path(opt$outdir, paste0(opt$prefix, "_", suffix))

  # 主结果
  saveRDS(seurat_obj, f("annotated.rds"))
  write.csv(res,             f("annotation_consensus.csv"),   row.names = FALSE, quote = FALSE)
  write.csv(round(mod_z, 3), f("MOD_zscore.csv"),             row.names = TRUE)
  write.csv(round(ucell_z, 3), f("UCell_zscore.csv"),         row.names = TRUE)
  write.csv(round(mod_cluster_raw, 4),   f("MOD_raw.csv"),    row.names = TRUE)
  write.csv(round(ucell_cluster_raw, 4), f("UCell_raw.csv"),  row.names = TRUE)

  # per-cell UCell
  if (opt$save_percell) {
    percell_df <- as.data.frame(ucell_res$percell)
    percell_df <- tibble::rownames_to_column(percell_df, var = "cell")
    write.csv(percell_df, f("ucell_percell_scores.csv"), row.names = FALSE, quote = FALSE)
  }

  # 阈值日志
  write_threshold_log(f("thresholds_log.txt"), opt, info_mod, info_ucell,
                      list(mod = mod_calls, ucell = ucell_calls),
                      ucell_assay_info = list(used_assay = ucell_res$used_assay,
                                              used_layer = ucell_res$used_layer),
                      subsample_mod   = if (!opt$skip_permutation) attr(null_mod,   "subsample_info") else NULL,
                      subsample_ucell = if (!opt$skip_permutation) attr(null_ucell, "subsample_info") else NULL)

  # 置换零分布
  if (!opt$skip_permutation) {
    null_df <- rbind(
      data.frame(method = "MOD",   as.data.frame(null_mod),   check.names = FALSE),
      data.frame(method = "UCell", as.data.frame(null_ucell), check.names = FALSE)
    )
    write.csv(null_df, f("permutation_null.csv"), row.names = FALSE, quote = FALSE)
  }

  # 诊断图
  plot_diagnostics(margins_mod, margins_ucell, info_mod, info_ucell,
                   f("margin_diagnostics.pdf"))

  # UMAP
  if (plot_available)
    save_umap_plots(seurat_obj, f("annotation.pdf"), opt$reduction,
                    cluster_col, "cell_type")

  if (opt$verbose) cat("基本输出已写入\n\n")

  # ---- 12. 子集提取 ----
  if (!is.null(opt$subset_types)) {
    if (opt$verbose) cat("==================== 子集提取 ====================\n")
    subset_types <- trimws(unlist(strsplit(opt$subset_types, ",")))
    # 子集匹配: 允许命中 cell_type, MOD 主标签, 或 (?) 后缀变体
    base_celltype <- sub("\\(\\?\\)$", "", seurat_obj$cell_type)
    keep_cells <- base_celltype %in% subset_types

    missing_types <- setdiff(subset_types, unique(base_celltype))
    if (length(missing_types) > 0)
      cat(sprintf("警告: 未找到细胞类型: %s\n", paste(missing_types, collapse = ", ")))

    if (sum(keep_cells) == 0) {
      cat("错误: 所有请求的细胞类型都未找到, 跳过子集提取\n")
    } else {
      seurat_subset <- subset(seurat_obj, cells = colnames(seurat_obj)[keep_cells])
      if (opt$verbose)
        cat(sprintf("提取 %d / %d 细胞\n", ncol(seurat_subset), ncol(seurat_obj)))
      saveRDS(seurat_subset, f("subset.rds"))
      if (plot_available)
        save_umap_plots(seurat_subset, f("subset.pdf"), opt$reduction,
                        cluster_col, "cell_type", " (Subset)")
      stats <- as.data.frame(table(seurat_subset$cell_type))
      names(stats) <- c("cell_type", "count")
      stats$percentage <- round(stats$count / sum(stats$count) * 100, 2)
      write.csv(stats, f("subset_stats.csv"), row.names = FALSE, quote = FALSE)
      if (opt$verbose) cat("子集已保存\n\n")
    }
  }

  # ---- 13. 最终总结 ----
  if (opt$verbose) {
    cat("==================================================\n")
    cat("分析完成!\n")
    cat("==================================================\n")
    cat("输出文件:\n")
    cat(sprintf("  注释 RDS:               %s\n", f("annotated.rds")))
    cat(sprintf("  注释主表:               %s\n", f("annotation_consensus.csv")))
    cat(sprintf("  MOD z-score / raw:      %s , %s\n", f("MOD_zscore.csv"), f("MOD_raw.csv")))
    cat(sprintf("  UCell z-score / raw:    %s , %s\n", f("UCell_zscore.csv"), f("UCell_raw.csv")))
    if (opt$save_percell)
      cat(sprintf("  UCell per-cell:         %s\n", f("ucell_percell_scores.csv")))
    cat(sprintf("  阈值日志:               %s\n", f("thresholds_log.txt")))
    if (!opt$skip_permutation)
      cat(sprintf("  置换零分布:             %s\n", f("permutation_null.csv")))
    cat(sprintf("  Margin 诊断图:          %s\n", f("margin_diagnostics.pdf")))
    if (plot_available)
      cat(sprintf("  UMAP 注释图:            %s\n", f("annotation.pdf")))
    cat("\n")
  }

}, error = function(e) {
  cat("\n错误发生:\n")
  cat(sprintf("  %s\n", conditionMessage(e)))
  if (opt$debug) { cat("\n详细信息:\n"); print(e) }
  quit(status = 1)
})

if (opt$verbose) cat("程序成功完成!\n")
quit(status = 0)
