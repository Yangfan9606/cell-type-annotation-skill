#!/usr/bin/env Rscript
#------------------------------------------#
# author:    Yangfan Zhou
# email:     yangfan.zhou@ki.se
# date:      2025-12-05
# version:   1.0
# license:   MIT
# brief:     single-cell data subsets Cell extract.
#------------------------------------------#

# ==============================================================================
library(optparse)
# ==============================================================================
# 1. 参数解析
# ==============================================================================

option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="输入RDS文件路径（Seurat对象） [必需]", metavar="PATH"),
  make_option(c("-o", "--output"), type="character", default=NULL,
              help="输出目录路径 [必需]", metavar="PATH"),
  make_option(c("-c", "--celltype"), type="character", default=NULL,
              help="要提取(保留)的细胞类型/簇号，多个用逗号分隔 [与--exclude二选一]", metavar="STRING"),
  make_option(c("-e", "--exclude"), type="character", default=NULL,
              help="要排除的细胞类型/簇号，多个用逗号分隔 [与--celltype二选一；保留名单=全部可用值减去此列表]", metavar="STRING"),
  make_option(c("--celltype_column"), type="character", default="celltype",
              help="存储细胞类型的列名 [默认: %default]", metavar="STRING"),
  make_option(c("-r", "--resolution"), type="character", default="0.3,0.5,0.8",
              help="聚类分辨率，用逗号分隔 [默认: %default]", metavar="STRING"),
  make_option(c("-d", "--dims"), type="integer", default=30,
              help="PCA维度数量 [默认: %default]", metavar="INTEGER"),
  make_option(c("-k", "--neighbors"), type="integer", default=20,
              help="KNN邻居数 [默认: %default]", metavar="INTEGER"),
  make_option(c("-f", "--features"), type="integer", default=2000,
              help="高变基因数量 [默认: %default]", metavar="INTEGER"),
  make_option(c("-t", "--top_markers"), type="integer", default=10,
              help="每个cluster保存的top标记基因数 [默认: %default]", metavar="INTEGER"),
  make_option(c("--use_sct"), action="store_true", default=FALSE,
              help="使用SCT标准化进行重新分析 [默认: FALSE，使用RNA]"),
  make_option(c("--debug"), action="store_true", default=FALSE,
              help="启用调试模式，输出详细信息")
)

opt_parser <- OptionParser(
  option_list=option_list,
  usage = "usage: %prog [options]",
  description = "\n单细胞亚群重注释分析脚本 (支持SCT)\n提取特定细胞类型进行重新聚类和差异分析\n"
)

# 解析参数
opt <- parse_args(opt_parser)

# 检查必需参数
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  print_help(opt_parser)
  quit(save="no", status=0)
}

if (is.null(opt[["input"]]) || is.null(opt[["output"]])) {
  cat("\n错误：必须提供 --input 和 --output 参数\n\n")
  quit(save="no", status=1)
}

if (is.null(opt[["celltype"]]) && is.null(opt[["exclude"]])) {
  cat("\n错误：必须提供 --celltype (保留名单) 或 --exclude (排除名单) 二者之一\n\n")
  quit(save="no", status=1)
}

if (!is.null(opt[["celltype"]]) && !is.null(opt[["exclude"]])) {
  cat("\n错误：--celltype 和 --exclude 不能同时使用，请二选一\n\n")
  quit(save="no", status=1)
}

if (!file.exists(opt[["input"]])) {
  stop(paste0("错误：输入文件不存在: ", opt[["input"]]), call.=FALSE)
}

# 创建输出目录
if (!dir.exists(opt[["output"]])) {
  dir.create(opt[["output"]], recursive=TRUE)
}

# 解析参数
resolutions <- as.numeric(strsplit(opt[["resolution"]], ",")[[1]])
celltypes_requested <- if (!is.null(opt[["celltype"]])) trimws(strsplit(opt[["celltype"]], ",")[[1]]) else NULL
exclude_requested   <- if (!is.null(opt[["exclude"]]))  trimws(strsplit(opt[["exclude"]], ",")[[1]])  else NULL
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

# 加载包
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

# 设置日志
log_file <- file.path(opt[["output"]], paste0("analysis_log_", timestamp, ".txt"))
log_conn <- file(log_file, open="wt")
sink(log_conn, type="output", split=TRUE)
sink(log_conn, type="message")

cat("==============================================================================\n")
cat("单细胞亚群提取 (sc_005.SubClass_extract.R)\n")
cat("开始时间:", as.character(Sys.time()), "\n")
cat("==============================================================================\n\n")

# ==============================================================================
# 2. 数据加载和验证
# ==============================================================================

cat("步骤 1/8: 加载数据...\n")
seurat_obj <- readRDS(opt[["input"]])
cat("  成功加载Seurat对象\n")
cat("  总细胞数:", ncol(seurat_obj), "\n")
cat("  总基因数:", nrow(seurat_obj), "\n")
cat("  可用Assays:", paste(Assays(seurat_obj), collapse=", "), "\n")
cat("  当前Active Assay:", DefaultAssay(seurat_obj), "\n\n")

# 检测是否有SCT assay
has_sct <- "SCT" %in% Assays(seurat_obj)
has_rna <- "RNA" %in% Assays(seurat_obj)
original_assay <- DefaultAssay(seurat_obj)

if (has_sct) {
  cat("  检测到SCT assay\n")
  if (original_assay == "SCT") {
    cat("  ⚠️  当前active assay是SCT，将使用特殊处理流程\n\n")
  }
}

# 验证细胞类型列
if (!opt[["celltype_column"]] %in% colnames(seurat_obj@meta.data)) {
  stop(paste0("错误：元数据中不存在列: ", opt[["celltype_column"]]), call.=FALSE)
}

available_celltypes <- unique(seurat_obj@meta.data[[opt[["celltype_column"]]]])
available_celltypes_chr <- as.character(available_celltypes)
cat("  可用细胞类型/簇号:", paste(available_celltypes_chr, collapse=", "), "\n\n")

# 解析保留名单：--celltype 直接用；--exclude 用"全部可用值减去排除值"反推
if (!is.null(exclude_requested)) {
  missing_exclude <- exclude_requested[!exclude_requested %in% available_celltypes_chr]
  if (length(missing_exclude) > 0) {
    stop(paste0("错误：--exclude 中以下值不存在于列 '", opt[["celltype_column"]], "': ",
                paste(missing_exclude, collapse=", ")), call.=FALSE)
  }
  celltypes <- setdiff(available_celltypes_chr, exclude_requested)
  # 若全部值都可解析为数字（如簇号），按数值排序；否则按字符排序，方便阅读日志
  if (!anyNA(suppressWarnings(as.numeric(celltypes)))) {
    celltypes <- as.character(sort(as.numeric(celltypes)))
  } else {
    celltypes <- sort(celltypes)
  }
  cat("  模式: --exclude 排除模式\n")
  cat("  排除:", paste(exclude_requested, collapse=", "), "\n")
  cat("  保留(自动反推):", paste(celltypes, collapse=", "), "\n\n")
} else {
  celltypes <- celltypes_requested
  missing_celltypes <- celltypes[!celltypes %in% available_celltypes_chr]
  if (length(missing_celltypes) > 0) {
    stop(paste0("错误：以下细胞类型不存在: ",
                paste(missing_celltypes, collapse=", ")), call.=FALSE)
  }
  cat("  模式: --celltype 保留模式\n")
  cat("  保留:", paste(celltypes, collapse=", "), "\n\n")
}

# ==============================================================================
# 3. 细胞提取
# ==============================================================================

cat("步骤 2/8: 提取目标细胞类型...\n")
subset_obj <- subset(seurat_obj, 
                     subset = !!sym(opt[["celltype_column"]]) %in% celltypes)

cat("  提取的细胞数:", ncol(subset_obj), "\n")
celltype_counts <- table(subset_obj@meta.data[[opt[["celltype_column"]]]])
for (ct in celltypes) {
  cat("    -", ct, ":", celltype_counts[ct], "个细胞\n")
}

# 清理陈旧reductions：subset()会把父对象的每个reduction(pca/umap/以及此前
# 整合产生的cca等)原样按行筛选保留下来——坐标数值没错，但语义已过期
# (基于排除前、包含被排除细胞在内的那次拟合)。本脚本随后会在提取出的这批
# 细胞上重新计算全新的pca/umap(步骤5、7)，不会读取这些旧reduction；若不清
# 理，保存的对象会同时留着"本次新算的pca/umap"和"过期的cca"，容易被后续误
# 用(如DimPlot(reduction="cca"))误读成这批细胞自己的整合结果。予以移除。
old_reductions <- Reductions(subset_obj)
if (length(old_reductions) > 0) {
  cat("\n  清理陈旧reductions(继承自提取前的父对象:", paste(old_reductions, collapse=", "), ")...\n")
  for (rd in old_reductions) {
    subset_obj[[rd]] <- NULL
  }
  cat("    已清除，pca/umap将在本次运行中重新计算\n")
}

# 清理factor水平：Seurat的subset()不会自动丢弃不再出现的factor level，
# 被排除的类型/簇号会以0个细胞的形式继续留在level里，导致后续DimPlot图例、
# table()等出现"幽灵分类"。对所有factor列做droplevels()清理。
cat("\n  清理factor水平(droplevels)...\n")
factor_cols <- names(subset_obj@meta.data)[sapply(subset_obj@meta.data, is.factor)]
for (fc in factor_cols) {
  subset_obj@meta.data[[fc]] <- droplevels(subset_obj@meta.data[[fc]])
}
if (length(factor_cols) > 0) {
  cat("    已清理", length(factor_cols), "个factor列:", paste(factor_cols, collapse=", "), "\n")
} else {
  cat("    未发现factor列，跳过\n")
}

# 清理metadata
cat("\n  清理metadata...\n")
score_cols <- grep("^Score_", colnames(subset_obj@meta.data), value = TRUE)
if (length(score_cols) > 0) {
  cat("    删除", length(score_cols), "个'Score_'列\n")
  subset_obj@meta.data <- subset_obj@meta.data[, 
    !colnames(subset_obj@meta.data) %in% score_cols, drop = FALSE]
}

# 重命名细胞类型列
cat("  重命名细胞类型列为PrimaryCellType\n")
if (opt[["celltype_column"]] != "PrimaryCellType") {
  subset_obj@meta.data$original_celltype <- subset_obj@meta.data[[opt[["celltype_column"]]]]
  colnames(subset_obj@meta.data)[colnames(subset_obj@meta.data) == opt[["celltype_column"]]] <- "PrimaryCellType"
} else {
  subset_obj@meta.data$original_celltype <- subset_obj@meta.data$PrimaryCellType
}

rm(seurat_obj)
gc()
cat("\n")

# ==============================================================================
# 4. 决定使用哪个assay进行重新分析
# ==============================================================================

cat("步骤 3/8: 决定分析策略...\n")

# 使用SCT须显式指定 --use_sct；唯一的自动例外是原对象根本没有RNA assay
# (此时没有第二选择，只能退回SCT)。绝不会在有RNA可用时静默选SCT。
if (opt[["use_sct"]] || (!has_rna && has_sct)) {
  analysis_assay <- "SCT"
  if (opt[["use_sct"]]) {
    cat("  ✓ 按 --use_sct 显式指定，使用SCT assay进行重新分析\n\n")
  } else {
    cat("  ✓ 原对象无RNA assay，回退使用SCT assay进行重新分析\n\n")
  }
} else if (has_rna) {
  analysis_assay <- "RNA"
  cat("  ✓ 使用RNA assay进行重新分析（推荐）\n")
  cat("  💡 如果要使用SCT，请添加 --use_sct 参数\n")
  if (has_sct) {
    cat("  🧹 清理残留的SCT assay(该模型基于原始全量对象拟合，对当前子集已过期，予以移除)\n")
    subset_obj[["SCT"]] <- NULL
  }
  cat("\n")
} else {
  stop("错误：没有可用的RNA或SCT assay", call.=FALSE)
}

DefaultAssay(subset_obj) <- analysis_assay

# Seurat v5 Assay5多层检查：若原对象整合流程中各样本的counts/data未Join
# (如 counts.sample1, counts.sample2, ...)，NormalizeData/FindVariableFeatures/
# FindAllMarkers 在多层对象上可能按层各算一遍或直接报错。仅在真正存在
# "同类型多层"(如多个counts.*)时才JoinLayers；已是单层(或v4 Assay)则跳过。
cat("  检查", analysis_assay, "assay是否需要JoinLayers()...\n")
if (!analysis_assay %in% Assays(subset_obj)) {
  # --use_sct 在原对象尚无SCT assay时会到这里：SCT将由下一步SCTransform()新建，
  # 此刻还不存在，无需(也无法)检查
  cat("    尚无", analysis_assay, "assay(将由后续SCTransform新建)，跳过\n\n")
} else {
  assay_obj <- subset_obj[[analysis_assay]]
  if (inherits(assay_obj, "Assay5")) {
    layer_names <- tryCatch(Layers(assay_obj), error=function(e) character(0))
    base_names <- sub("\\..*$", "", layer_names)
    needs_join <- length(layer_names) > 0 && any(table(base_names) > 1)
    if (needs_join) {
      cat("    检测到多层:", paste(layer_names, collapse=", "), "→ 执行 JoinLayers()\n")
      subset_obj[[analysis_assay]] <- JoinLayers(subset_obj[[analysis_assay]])
      cat("    ✓ JoinLayers完成，合并后层:", paste(Layers(subset_obj[[analysis_assay]]), collapse=", "), "\n\n")
    } else {
      cat("    未发现按样本拆分的多层(", paste(layer_names, collapse=", "), ")，跳过\n\n")
    }
  } else {
    cat("    非Assay5(v4 Assay或其他)，无JoinLayers概念，跳过\n\n")
  }
}

# ==============================================================================
# 5. 重新分析流程
# ==============================================================================

cat("步骤 4/8: 数据标准化和高变基因识别...\n")

if (analysis_assay == "SCT") {
  # 使用SCT进行标准化
  cat("  使用SCTransform标准化...\n")
  subset_obj <- SCTransform(subset_obj, variable.features.n = opt[["features"]],
                           verbose = FALSE)
  cat("  SCTransform完成\n")
} else {
  # 使用标准流程
  subset_obj <- NormalizeData(subset_obj, verbose=FALSE)
  subset_obj <- FindVariableFeatures(subset_obj, selection.method="vst",
                                     nfeatures=opt[["features"]], verbose=FALSE)
}

cat("  识别高变基因数:", length(VariableFeatures(subset_obj)), "\n\n")

cat("步骤 5/8: 数据缩放和PCA降维...\n")
subset_obj <- ScaleData(subset_obj, verbose=FALSE)
subset_obj <- RunPCA(subset_obj, features=VariableFeatures(subset_obj),
                    npcs=opt[["dims"]], verbose=FALSE)
cat("  PCA降维完成\n\n")

cat("步骤 6/8: 构建邻居图和聚类分析...\n")
subset_obj <- FindNeighbors(subset_obj, 
                           dims=1:opt[["dims"]], k.param=opt[["neighbors"]],
                           verbose=FALSE)

# 聚类（修复版）
for (res in resolutions) {
  cat("  分辨率", res, "聚类中...\n")
  subset_obj <- FindClusters(subset_obj, resolution=res, verbose=FALSE)
  # 智能检测聚类列
  possible_patterns <- c(
    paste0(analysis_assay, "_snn_res.", res),
    paste0(".*_snn_res.", res), "seurat_clusters" )
  
  cluster_col_found <- NULL
  for (pattern in possible_patterns) {
    matching_cols <- grep(pattern, colnames(subset_obj@meta.data), value = TRUE)
    if (length(matching_cols) > 0) {
      cluster_col_found <- matching_cols[length(matching_cols)]
      break
    }
  }
  
  new_col <- paste0("subcluster_res", res)
  
  if (!is.null(cluster_col_found)) {
    subset_obj@meta.data[[new_col]] <- subset_obj@meta.data[[cluster_col_found]]
    n_clusters <- length(unique(subset_obj@meta.data[[new_col]]))
    
    if (opt[["debug"]]) {
      cat("    [调试] 找到聚类列:", cluster_col_found, "\n")
      cat("    [调试] 聚类分布:", paste(table(subset_obj@meta.data[[new_col]]), collapse=", "), "\n")
    }
  } else {
    n_clusters <- 0
    warning(paste0("警告: 未找到分辨率 ", res, " 的聚类结果列"))
  }
  
  cat("    识别到", n_clusters, "个亚群\n")
}
cat("\n")

cat("步骤 7/8: UMAP可视化...\n")
subset_obj <- RunUMAP(subset_obj, 
                     dims=1:opt[["dims"]],
					 n.neighbors = opt[["neighbors"]],
                     verbose=FALSE)
cat("  UMAP降维完成\n\n")

# ==============================================================================
# 6. 差异基因分析（SCT特殊处理）
# ==============================================================================

cat("步骤 8/8: 差异表达基因分析...\n")

default_res <- paste0("subcluster_res", resolutions[1])

if (!default_res %in% colnames(subset_obj@meta.data)) {
  cat("  警告：未找到默认分辨率的聚类结果，跳过marker分析\n")
  all_markers <- data.frame()
  top_markers <- data.frame()
} else {
  Idents(subset_obj) <- default_res
  n_clusters <- length(unique(Idents(subset_obj)))
  
  if (n_clusters < 2) {
    cat("  警告：只有", n_clusters, "个cluster，无法进行差异分析\n")
    all_markers <- data.frame()
    top_markers <- data.frame()
  } else {
    cat("  准备差异分析...\n")
    
    # 关键修复：如果是SCT assay，需要PrepSCTFindMarkers
    if (analysis_assay == "SCT") {
      cat("  ⚠️  检测到SCT assay，运行PrepSCTFindMarkers()...\n")
      tryCatch({
        subset_obj <- PrepSCTFindMarkers(subset_obj, verbose = FALSE)
        cat("  ✓ PrepSCTFindMarkers完成\n")
      }, error = function(e) {
        cat("  ⚠️  PrepSCTFindMarkers失败，将切换到RNA assay\n")
        if ("RNA" %in% Assays(subset_obj)) {
          DefaultAssay(subset_obj) <<- "RNA"
          cat("  ✓ 已切换到RNA assay\n")
        } else {
          stop("无法切换到RNA assay", call.=FALSE)
        }
      })
    }
    
    # FindAllMarkers
    cat("  寻找差异表达基因...\n")
    tryCatch({
      all_markers <- FindAllMarkers(subset_obj,
                                   only.pos=TRUE,
                                   min.pct=0.25,
                                   logfc.threshold=0.25,
                                   verbose=FALSE)
      
      if (nrow(all_markers) > 0) {
        cat("  ✓ 共识别", nrow(all_markers), "个差异基因\n")
        
        # 获取top markers
        top_markers <- all_markers %>%
          group_by(cluster) %>%
          slice_max(order_by = avg_log2FC, n = opt[["top_markers"]])
        
        cat("  ✓ 每个cluster保留top", opt[["top_markers"]], "个标记基因\n")
      } else {
        cat("  ⚠️  未识别到差异基因\n")
        top_markers <- data.frame()
      }
    }, error = function(e) {
      cat("  ✗ 差异分析失败:", conditionMessage(e), "\n")
      all_markers <<- data.frame()
      top_markers <<- data.frame()
    })
  }
}
cat("\n")

# ==============================================================================
# 7. 结果保存
# ==============================================================================

cat("保存结果...\n")

output_rds <- file.path(opt[["output"]], "subclustered.rds")
if (file.exists(output_rds)) {
  old_mtime <- format(file.info(output_rds)$mtime, "%Y-%m-%d %H:%M:%S")
  cat("  ⚠️  输出文件已存在(上次写入时间:", old_mtime, ")，本次运行将覆盖它\n")
}
saveRDS(subset_obj, file=output_rds)
cat("  ✓ Seurat对象已保存:", output_rds, "\n")
cat("  ℹ️  本次运行时间戳:", timestamp, "— 文件名固定为subclustered.rds(不含时间戳)，",
    "如需追溯是哪次运行生成的，请查看本日志文件名/内容中的时间戳\n")

# 保存差异基因
if (nrow(all_markers) > 0) {
  markers_file <- file.path(opt[["output"]], paste0("markers_all_", timestamp, ".csv"))
  write.csv(all_markers, file=markers_file, row.names=FALSE)
  cat("  ✓ 所有差异基因已保存:", markers_file, "\n")
  
  if (nrow(top_markers) > 0) {
    top_markers_file <- file.path(opt[["output"]], paste0("markers_top", opt[["top_markers"]], "_", timestamp, ".csv"))
    write.csv(top_markers, file=top_markers_file, row.names=FALSE)
    cat("  ✓ Top标记基因已保存:", top_markers_file, "\n")
  }
} else {
  cat("  ⚠️  未生成差异基因文件\n")
}

# 保存可视化
cat("\n生成可视化图...\n")

celltype_label <- if(length(celltypes) == 1) celltypes[1] else paste0(length(celltypes), " types")

# 原始注释
p1 <- DimPlot(subset_obj, group.by="original_celltype", reduction="umap") +
  ggtitle(paste("原始注释:", celltype_label)) +
  theme_minimal()
ggsave(file.path(opt[["output"]], paste0("umap_original_", timestamp, ".png")),
       plot=p1, width=8, height=6, dpi=300)
cat("  ✓ 原始注释UMAP已保存\n")

# 各分辨率
for (res in resolutions) {
  cluster_col <- paste0("subcluster_res", res)
  if (cluster_col %in% colnames(subset_obj@meta.data)) {
    p <- DimPlot(subset_obj, group.by=cluster_col, reduction="umap",
                label=TRUE, label.size=4) +
      ggtitle(paste0(celltype_label, " 亚群 (分辨率 ", res, ")")) +
      theme_minimal()
    ggsave(file.path(opt[["output"]], paste0("umap_res", res, "_", timestamp, ".png")),
           plot=p, width=8, height=6, dpi=300)
    cat("  ✓ 分辨率", res, "UMAP已保存\n")
  }
}

# 对比图
if (default_res %in% colnames(subset_obj@meta.data)) {
  p_combined <- DimPlot(subset_obj,
                       group.by=c("original_celltype", default_res),
                       reduction="umap", label=TRUE, label.size=3, ncol=2)
  ggsave(file.path(opt[["output"]], paste0("umap_comparison_", timestamp, ".png")),
         plot=p_combined, width=16, height=6, dpi=300)
  cat("  ✓ 对比图已保存\n")
}

# ==============================================================================
# 8. 完成
# ==============================================================================

cat("\n==============================================================================\n")
cat("分析完成！\n")
cat("结束时间:", as.character(Sys.time()), "\n")
cat("==============================================================================\n")

cat("\n输出文件:\n")
cat("  1. Seurat对象:", output_rds, "\n")
if (nrow(all_markers) > 0) {
  cat("  2. 差异基因表格:", list.files(opt[["output"]], pattern="markers_.*\\.csv", full.names=TRUE), "\n")
}
cat("  3. UMAP可视化图:", opt[["output"]], "\n")
cat("  4. 分析日志:", log_file, "\n")

cat("\n分析策略:\n")
cat("  • 使用的Assay:", analysis_assay, "\n")
if (analysis_assay == "SCT") {
  cat("  • SCT特殊处理: PrepSCTFindMarkers已应用\n")
}
cat("  • 已删除所有'Score_'列\n")
cat("  • 细胞类型列已重命名为'PrimaryCellType'\n")

sink(type="message")
sink(type="output")
close(log_conn)

cat("\n✓ 分析日志已保存至:", log_file, "\n")
