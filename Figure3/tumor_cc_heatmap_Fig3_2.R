#library(monocle3)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(clusterProfiler)
library(org.Hs.eg.db)
library(grid)

# ==============================================================================
# CC pseudotime dynamic heatmap workflow
# ==============================================================================

mat_cc_scaled <- readRDS("F:/Chordoma/Result/Tumor_cc/mat_cc_scaled.rds")
cell_order_cc <- readRDS("F:/Chordoma/Result/Tumor_cc/cell_order_cc.rds")
cluster_ordered <- readRDS("F:/Chordoma/Result/Tumor_cc/cluster_ordered.rds")
pt_cc <- readRDS("F:/Chordoma/Result/Tumor_cc/pt_cc.rds")
gene_clusters <- read.csv("F:/Chordoma/Result/Tumor_cc/gene_clusters.csv", header = TRUE)

rownames(gene_clusters) <- gene_clusters[, 1]
gene_clusters <- gene_clusters[, 2, drop = FALSE]
colnames(gene_clusters) <- "module"

mat_cc_scaled <- mat_cc_scaled[, cell_order_cc]
pt_ordered <- pt_cc[cell_order_cc]
cc_order <- c(12,11,5,19,9,16,13,  15,17,0,6,8,3,2,  18,1,10,7,  4,14 )
cluster_ordered <- factor(cluster_ordered[cell_order_cc], levels = cc_order)
#cluster_ordered <- factor(cluster_ordered[cell_order_cc], levels = 0:20)
cluster_levels <- levels(cluster_ordered)


zzm60colors2 <- c(
  "#76a2be", "#4b6aa8", "#2d3462", "#e0cfda", "#e6b884",
  "#d69a55", "#64a776", "#cc7f73", "#927c9a", "#efd2c9",
  "#c6adb0", "#df5734", "#6c408e", "#ac6894", "#b7deea",
  "#83ab8e", "#d4c2db", "#ece399", "#cbdaa9", "#b95055",
  "#bc9a7f", "#da6f6d", "#ebb1a4", "#a44e89", "#a9c2cb",
  "#b85292", "#6d6fa0", "#8d689d", "#c8c7e1", "#d25774",
  "#c49abc", "#b05545", "#405993", "#9f8d89", "#72567a",
  "#63a3b8", "#c4daec", "#3674a2", "#537eb7", "#e29eaf",
  "#4490c4", "#e6e2a3", "#de8b36", "#c4612f", "#9a70a8",
  "#408444", "#9d3b62", "#d5bb72", "#d8a0c0", "#61bada"
)

#cluster_colors <- zzm60colors2[1:length(cluster_levels)]
#names(cluster_colors) <- cluster_levels

cluster_colors_all <- zzm60colors2[1:20]
names(cluster_colors_all) <- 0:19
# 按 cc_order 重排
cluster_colors <- cluster_colors_all[as.character(cc_order)]
names(cluster_colors) <- cc_order

col_anno <- HeatmapAnnotation(
  Cluster = cluster_ordered,
  Pseudotime = pt_ordered,
  col = list(
    Cluster = cluster_colors,
    Pseudotime = colorRamp2(
      c(min(pt_ordered), max(pt_ordered)),
      c("#440154FF", "#FDE725FF")
    )
  ),
  annotation_name_side = "left"
)

modules <- gene_clusters[rownames(mat_cc_scaled), 1]
modules <- factor(modules)

module_colors <- c("#B7B7EB", "#EAB883", "#9BBBE1")
module_labels <- c(
  "Module 1\nregulation of T-cell activation",
  "Module 2\ncytoplasmic translation",
  "Module 3\nresponse to oxygen levels"
)

right_anno <- rowAnnotation(
  Module_Info = anno_block(
    gp = gpar(fill = module_colors),
    labels = module_labels,
    labels_gp = gpar(col = "white", fontsize = 9, fontface = "bold"),
    width = unit(5, "cm")
  )
)

pdf(
  "F:/Chordoma/Plots/Heatmap_CC_Pseudotime_Dynamic_Genes_original_backup1.pdf",
  width = 12,
  height = 10
)

Heatmap(
  mat_cc_scaled,
  name = "Z-score",
  col = colorRamp2(c(-2, 0, 2), c("navy", "white", "firebrick")),
  row_split = modules,
  cluster_row_slices = FALSE,
  column_split = cluster_ordered,
  cluster_columns = FALSE,
  cluster_rows = TRUE,
  top_annotation = col_anno,
  right_annotation = right_anno,
  show_column_names = FALSE,
  show_row_names = FALSE,
  row_names_gp = gpar(fontsize = 6),
  use_raster = TRUE,
  raster_quality = 3,
  column_gap = unit(0.5, "mm"),
  column_title = "Pseudotime Evolution within CC Tumor Cells (Early -> Late)",
  row_title = "Dynamic Gene Modules"
)

dev.off()
