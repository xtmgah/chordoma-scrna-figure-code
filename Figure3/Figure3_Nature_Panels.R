script_dir <- if (exists("SCRIPT_DIR", inherits = FALSE)) {
  SCRIPT_DIR
} else {
  dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), mustWork = TRUE))
}
source(file.path(dirname(script_dir), "_figure_style.R"))

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
  library(ggrepel)
  library(patchwork)
  library(ggraph)
  library(igraph)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})
if (requireNamespace("scCustomize", quietly = TRUE)) {
  suppressPackageStartupMessages(library(scCustomize))
}

OUT_DIR <- file.path(script_dir, "nature_panels")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

message("Figure 3: loading tumor Seurat object")
sc_tumor <- readRDS(file.path(script_dir, "sc_tumor.rds"))
meta <- sc_tumor@meta.data
cluster_col <- if ("RNA_snn_res.0.5" %in% colnames(meta)) "RNA_snn_res.0.5" else grep("snn_res", colnames(meta), value = TRUE)[1]
stopifnot(!is.na(cluster_col))
cluster_levels <- sort_numeric_labels(unique(as.character(meta[[cluster_col]])))
sc_tumor@meta.data[[cluster_col]] <- factor(as.character(sc_tumor@meta.data[[cluster_col]]), levels = cluster_levels)
cluster_cols <- setNames(nature_discrete(length(cluster_levels)), cluster_levels)

Idents(sc_tumor) <- cluster_col
p_umap <- if (requireNamespace("scCustomize", quietly = TRUE)) {
  DimPlot_scCustom(
    sc_tumor,
    reduction = "umap",
    group.by = cluster_col,
    label = TRUE,
    repel = TRUE,
    label.size = 3.2,
    pt.size = 0.15,
    colors_use = cluster_cols,
    raster = FALSE
  )
} else {
  DimPlot(
    sc_tumor,
    reduction = "umap",
    group.by = cluster_col,
    label = TRUE,
    repel = TRUE,
    label.size = 3.2,
    pt.size = 0.15,
    cols = cluster_cols,
    raster = FALSE
  )
}
p_umap <- p_umap +
  labs(title = NULL, color = "Tumor cluster") +
  theme_umap_nature() +
  guides(color = guide_legend(ncol = 2, override.aes = list(size = 3.0))) +
  theme(legend.position = "right", legend.text = element_text(size = 7.4))
save_panel_pdf(p_umap, file.path(OUT_DIR, "fig3a_tumor_subcluster_umap.pdf"), 4.8, 4.0)

if ("clinical_group" %in% colnames(meta)) {
  p_group_umap <- if (requireNamespace("scCustomize", quietly = TRUE)) {
    DimPlot_scCustom(
      sc_tumor,
      reduction = "umap",
      group.by = "clinical_group",
      colors_use = group_palette[c("CC", "PC")],
      pt.size = 0.13,
      shuffle = TRUE,
      raster = FALSE
    )
  } else {
    DimPlot(
      sc_tumor,
      reduction = "umap",
      group.by = "clinical_group",
      cols = group_palette[c("CC", "PC")],
      pt.size = 0.13,
      shuffle = TRUE,
      raster = FALSE
    )
  }
  p_group_umap <- p_group_umap +
    labs(title = NULL, color = "Group") +
    theme_umap_nature() +
    theme(legend.position = "right")
  save_panel_pdf(p_group_umap, file.path(OUT_DIR, "fig3b_tumor_umap_by_pathology.pdf"), 4.5, 3.8)

  props <- meta %>%
    count(.data[[cluster_col]], clinical_group, name = "n") %>%
    group_by(.data[[cluster_col]]) %>%
    mutate(freq = n / sum(n)) %>%
    ungroup()
  colnames(props)[1] <- "cluster"
  props$cluster <- factor(props$cluster, levels = cluster_levels)
  props$clinical_group <- factor(props$clinical_group, levels = c("CC", "PC"))
  p_check <- ggplot(props, aes(x = cluster, y = freq, fill = clinical_group)) +
    geom_col(width = 0.72, color = "white", linewidth = 0.2) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey45", linewidth = 0.35) +
    scale_fill_manual(values = group_palette[c("CC", "PC")], na.value = "grey80") +
    percent_axis() +
    labs(x = "Tumor cluster", y = "Cell fraction", fill = "Group") +
    theme_nature() +
    theme(panel.grid.major.x = element_blank())
  save_panel_pdf(p_check, file.path(OUT_DIR, "fig3c_pathology_composition_by_tumor_cluster.pdf"), 5.0, 3.0)
}

if ("Cluster" %in% colnames(meta)) {
  props_patient <- meta %>%
    count(Cluster, .data[[cluster_col]], name = "n") %>%
    group_by(Cluster) %>%
    mutate(freq = n / sum(n)) %>%
    ungroup()
  colnames(props_patient)[2] <- "tumor_cluster"
  props_patient$tumor_cluster <- factor(props_patient$tumor_cluster, levels = cluster_levels)
  p_patient <- ggplot(props_patient, aes(x = factor(Cluster), y = freq, fill = tumor_cluster)) +
    geom_col(width = 0.68, color = "white", linewidth = 0.18) +
    scale_fill_manual(values = cluster_cols, drop = FALSE) +
    percent_axis() +
    labs(x = "Patient cluster", y = "Cell fraction", fill = "Tumor cluster") +
    theme_nature() +
    theme(panel.grid.major.x = element_blank())
  save_panel_pdf(p_patient, file.path(OUT_DIR, "fig3d_tumor_subtype_composition_by_patient_cluster.pdf"), 2.65, 3.0)
}

message("Figure 3: pseudotime trajectory panel")
pt_file <- file.path(script_dir, "pt_cc.rds")
if (file.exists(pt_file) && "umap" %in% names(sc_tumor@reductions)) {
  pseudotime <- readRDS(pt_file)
  umap_embed <- Embeddings(sc_tumor, "umap")
  common_cells <- intersect(names(pseudotime), rownames(umap_embed))
  if (length(common_cells) > 100) {
    trajectory_df <- data.frame(
      cell = common_cells,
      UMAP_1 = umap_embed[common_cells, 1],
      UMAP_2 = umap_embed[common_cells, 2],
      pseudotime = as.numeric(pseudotime[common_cells]),
      tumor_cluster = as.character(sc_tumor@meta.data[common_cells, cluster_col]),
      stringsAsFactors = FALSE
    ) %>%
      filter(is.finite(pseudotime), !is.na(tumor_cluster))
    trajectory_df$tumor_cluster <- factor(trajectory_df$tumor_cluster, levels = cluster_levels)

    path_bins <- min(90, max(30, floor(nrow(trajectory_df) / 1500)))
    trajectory_path <- trajectory_df %>%
      arrange(pseudotime) %>%
      mutate(path_bin = ntile(pseudotime, path_bins)) %>%
      group_by(path_bin) %>%
      summarise(
        UMAP_1 = stats::median(UMAP_1, na.rm = TRUE),
        UMAP_2 = stats::median(UMAP_2, na.rm = TRUE),
        pseudotime = stats::median(pseudotime, na.rm = TRUE),
        n = dplyr::n(),
        .groups = "drop"
      ) %>%
      filter(n >= 10) %>%
      arrange(pseudotime)
    if (nrow(trajectory_path) >= 6) {
      path_index <- seq_len(nrow(trajectory_path))
      smooth_df <- min(14, max(5, floor(nrow(trajectory_path) / 6)))
      trajectory_path$UMAP_1 <- stats::predict(stats::smooth.spline(path_index, trajectory_path$UMAP_1, df = smooth_df), path_index)$y
      trajectory_path$UMAP_2 <- stats::predict(stats::smooth.spline(path_index, trajectory_path$UMAP_2, df = smooth_df), path_index)$y
    }

    cluster_label_df <- trajectory_df %>%
      group_by(tumor_cluster) %>%
      summarise(
        UMAP_1 = stats::median(UMAP_1, na.rm = TRUE),
        UMAP_2 = stats::median(UMAP_2, na.rm = TRUE),
        .groups = "drop"
      )

    trajectory_theme <- theme_umap_nature(base_size = 9.2) +
      theme(
        plot.title = element_text(face = "plain", size = 10.5, hjust = 0.5),
        legend.position = "right",
        legend.title = element_text(face = "plain", size = 8.6),
        legend.text = element_text(size = 7.6),
        plot.margin = margin(3, 4, 3, 4)
      )

    p_cluster_traj <- ggplot(trajectory_df, aes(UMAP_1, UMAP_2)) +
      geom_point(aes(color = tumor_cluster), size = 0.018, alpha = 0.80, stroke = 0) +
      geom_path(
        data = trajectory_path,
        aes(UMAP_1, UMAP_2),
        inherit.aes = FALSE,
        linewidth = 0.30,
        alpha = 0.82,
        lineend = "round",
        linejoin = "round",
        color = "black"
      ) +
      ggrepel::geom_text_repel(
        data = cluster_label_df,
        aes(x = UMAP_1, y = UMAP_2, label = tumor_cluster),
        inherit.aes = FALSE,
        size = 3.3,
        family = FONT_FAMILY,
        fontface = "bold",
        color = "black",
        min.segment.length = 0,
        seed = 7,
        box.padding = 0.25,
        point.padding = 0.15,
        segment.size = 0.15
      ) +
      scale_color_manual(values = cluster_cols, drop = FALSE) +
      coord_equal() +
      labs(title = "Tumor clusters", color = "Tumor cluster") +
      guides(color = guide_legend(ncol = 2, override.aes = list(size = 2.6, alpha = 1))) +
      trajectory_theme

    p_time_traj <- ggplot(trajectory_df, aes(UMAP_1, UMAP_2)) +
      geom_point(aes(color = pseudotime), size = 0.018, alpha = 0.86, stroke = 0) +
      geom_path(
        data = trajectory_path,
        aes(UMAP_1, UMAP_2),
        inherit.aes = FALSE,
        linewidth = 0.30,
        alpha = 0.82,
        lineend = "round",
        linejoin = "round",
        color = "black"
      ) +
      scale_color_viridis_c(option = "plasma", name = "Pseudotime") +
      coord_equal() +
      labs(title = "Pseudotime") +
      guides(color = guide_colorbar(
        title.position = "top",
        title.hjust = 0.5,
        barwidth = unit(24, "mm"),
        barheight = unit(3.0, "mm")
      )) +
      trajectory_theme +
      theme(legend.position = "bottom")

    save_panel_pdf(
      p_cluster_traj + p_time_traj + plot_layout(widths = c(1.22, 1)),
      file.path(OUT_DIR, "fig3e_pseudotime_trajectory.pdf"),
      7.2,
      3.45
    )
  }
}

message("Figure 3: TWIST1 and PTGES expression panels")
genes <- intersect(c("TWIST1", "PTGES"), rownames(sc_tumor))
if (length(genes) > 0) {
  feature_colors <- viridis::magma(60, begin = 0.08, end = 0.95)
  fp <- if (requireNamespace("scCustomize", quietly = TRUE)) {
    FeaturePlot_scCustom(
      sc_tumor,
      features = genes,
      colors_use = feature_colors,
      na_color = "grey92",
      order = TRUE,
      min.cutoff = "q05",
      max.cutoff = "q95",
      pt.size = 0.027,
      raster = FALSE,
      combine = FALSE
    )
  } else {
    FeaturePlot(
      sc_tumor,
      features = genes,
      cols = c("grey92", "#B2182B"),
      order = TRUE,
      min.cutoff = "q05",
      max.cutoff = "q95",
      pt.size = 0.027,
      combine = FALSE,
      raster = FALSE
    )
  }
  fp <- lapply(fp, function(p) p + theme_umap_nature() + theme(plot.title = element_text(face = "bold", hjust = 0.5)))
  for (i in seq_along(fp)) {
    save_panel_pdf(fp[[i]], file.path(OUT_DIR, paste0("fig3e_feature_", genes[[i]], ".pdf")), 3.1, 3.0)
  }
  save_panel_pdf(wrap_plots(fp, nrow = 1), file.path(OUT_DIR, "fig3e_twist1_ptges_featureplots.pdf"), 6.2, 3.0)
}

if (all(c("TWIST1", "PTGES") %in% rownames(sc_tumor))) {
  expr_data <- FetchData(sc_tumor, vars = c("TWIST1", "PTGES"))
  if (sum(expr_data$TWIST1 > 0, na.rm = TRUE) > 10) {
    breaks <- quantile(expr_data$TWIST1[expr_data$TWIST1 > 0], probs = c(0.33, 0.66), na.rm = TRUE)
    expr_data$Group <- case_when(
      expr_data$TWIST1 == 0 ~ "Negative",
      expr_data$TWIST1 <= breaks[1] ~ "Low",
      expr_data$TWIST1 <= breaks[2] ~ "Middle",
      TRUE ~ "High"
    )
    expr_data$Group <- factor(expr_data$Group, levels = c("Negative", "Low", "Middle", "High"))
    p_label_y <- max(expr_data$PTGES, na.rm = TRUE) * 1.015
    p_box <- ggplot(expr_data, aes(x = Group, y = PTGES, fill = Group)) +
      geom_violin(width = 0.80, linewidth = 0, alpha = 0.82, trim = FALSE, scale = "width", color = NA) +
      geom_boxplot(width = 0.18, fill = NA, color = "grey20", outlier.shape = NA, linewidth = 0.22) +
      stat_compare_means(method = "kruskal.test", label.y = p_label_y, size = 2.8, family = FONT_FAMILY) +
      scale_fill_manual(values = c("Negative" = "#D9D9D9", "Low" = "#FCAE91", "Middle" = "#FB6A4A", "High" = "#A50F15")) +
      scale_y_continuous(expand = expansion(mult = c(0.02, 0.06))) +
      labs(x = "TWIST1 expression group", y = "PTGES expression") +
      theme_nature() +
      theme(legend.position = "none", panel.grid.major.x = element_blank())
    save_panel_pdf(p_box, file.path(OUT_DIR, "fig3f_ptges_by_twist1_expression.pdf"), 2.35, 3.1)
  }
}

message("Figure 3: co-expression network panels")
plot_network <- function(graph, out_file, core_color, neighbor_color, core_label_size = 3.0, edge_width_range = c(0.12, 0.75)) {
  vertex_types <- igraph::vertex_attr(graph, "type")
  if (is.null(vertex_types)) vertex_types <- rep("Gene", igraph::vcount(graph))
  vertex_cols <- c("Core (Top 10)" = core_color, "Neighbor" = neighbor_color, "Gene" = core_color)
  layout_weights <- abs(igraph::edge_attr(graph, "weight"))
  positive_weights <- layout_weights[!is.na(layout_weights) & layout_weights > 0]
  fallback_weight <- if (length(positive_weights) > 0) min(positive_weights, na.rm = TRUE) / 10 else 1
  layout_weights[is.na(layout_weights) | layout_weights == 0] <- fallback_weight
  p <- ggraph(graph, layout = "fr", weights = layout_weights, niter = 2500) +
    geom_edge_link(aes(edge_alpha = abs(weight), edge_width = abs(weight), edge_colour = weight), show.legend = FALSE) +
    scale_edge_colour_gradient2(low = "#2166AC", mid = "grey92", high = "#B2182B", midpoint = 0) +
    scale_edge_width(range = edge_width_range) +
    geom_node_point(aes(color = type, size = degree), alpha = 0.94) +
    geom_node_text(
      aes(label = ifelse(type == "Core (Top 10)", name, "")),
      repel = TRUE,
      size = core_label_size,
      fontface = "bold",
      family = FONT_FAMILY,
      point.padding = 0.15
    ) +
    geom_node_text(
      aes(label = ifelse(type == "Core (Top 10)", "", name)),
      repel = TRUE,
      size = 1.65,
      family = FONT_FAMILY,
      point.padding = 0.10
    ) +
    scale_color_manual(values = vertex_cols, drop = FALSE, guide = guide_legend(override.aes = list(size = 4.0))) +
    scale_size(range = c(1.2, 4.2), guide = "none") +
    labs(color = "Gene class") +
    theme_void(base_size = 8.5, base_family = FONT_FAMILY) +
    theme(legend.position = "right", plot.margin = margin(2, 2, 2, 2))
  save_panel_pdf(p, out_file, 4.4, 3.8)
}
graph_cc <- readRDS(file.path(script_dir, "graph_cc.rds"))
graph_pc <- readRDS(file.path(script_dir, "graph_pc.rds"))
plot_network(graph_cc, file.path(OUT_DIR, "fig3g_cc_coexpression_network.pdf"), "#B2182B", "#9ECAE1")
plot_network(graph_pc, file.path(OUT_DIR, "fig3h_pc_coexpression_network.pdf"), "#7E5A9B", "#FDD0A2", core_label_size = 2.65, edge_width_range = c(0.08, 0.55))

message("Figure 3: CC pseudotime dynamic heatmap")
mat_cc_scaled <- readRDS(file.path(script_dir, "mat_cc_scaled.rds"))
cell_order_cc <- readRDS(file.path(script_dir, "cell_order_cc.rds"))
cluster_ordered <- readRDS(file.path(script_dir, "cluster_ordered.rds"))
pt_cc <- readRDS(file.path(script_dir, "pt_cc.rds"))
gene_clusters <- read.csv(file.path(script_dir, "gene_clusters.csv"), header = TRUE, check.names = FALSE)
rownames(gene_clusters) <- gene_clusters[[1]]

cell_order_cc <- intersect(cell_order_cc, colnames(mat_cc_scaled))
mat_cc_scaled <- mat_cc_scaled[, cell_order_cc, drop = FALSE]
pt_ordered <- pt_cc[cell_order_cc]
cluster_ordered <- cluster_ordered[cell_order_cc]
cc_order <- as.character(c(12, 11, 5, 19, 9, 16, 13, 15, 17, 0, 6, 8, 3, 2, 18, 1, 10, 7, 4, 14))
cluster_values <- as.character(cluster_ordered)
cluster_levels_original <- c(intersect(cc_order, unique(cluster_values)), setdiff(unique(cluster_values), cc_order))
cluster_ordered <- factor(cluster_values, levels = cluster_levels_original)
cluster_color_levels <- sort_numeric_labels(unique(cluster_values))
cluster_heat_cols_all <- setNames(nature_discrete(length(cluster_color_levels)), cluster_color_levels)
cluster_heat_cols <- cluster_heat_cols_all[levels(cluster_ordered)]
modules <- factor(gene_clusters[rownames(mat_cc_scaled), 2])
pseudotime_col_fun <- colorRamp2(range(pt_ordered, na.rm = TRUE), c("#2166AC", "#FDD142"))
module_levels <- levels(modules)
module_cols <- setNames(c("#B7B7EB", "#EAB883", "#9BBBE1", "#6AAE75")[seq_along(module_levels)], module_levels)
module_label_map <- setNames(
  c(
    "Module 1\nT-cell activation",
    "Module 2\nCytoplasmic translation",
    "Module 3\nOxygen response",
    paste("Module", module_levels[-seq_len(min(3, length(module_levels)))])
  )[seq_along(module_levels)],
  module_levels
)
row_anno <- rowAnnotation(
  Module_Info = anno_block(
    gp = gpar(fill = module_cols),
    labels = unname(module_label_map[module_levels]),
    labels_rot = 90,
    labels_gp = gpar(col = "white", fontsize = 6.2, fontface = "plain", fontfamily = FONT_FAMILY),
    width = unit(8.0, "mm")
  )
)
heatmap_pdf <- file.path(OUT_DIR, "fig3i_cc_pseudotime_dynamic_gene_heatmap.pdf")
open_panel_pdf(heatmap_pdf, 6.3, 5.8)
bar_columns <- colnames(mat_cc_scaled)
cluster_bar <- Heatmap(
  matrix(as.character(cluster_ordered), nrow = 1, dimnames = list("Cluster", bar_columns)),
  name = "Cluster",
  col = cluster_heat_cols,
  column_split = cluster_ordered,
  cluster_column_slices = FALSE,
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  show_column_names = FALSE,
  column_title = NULL,
  show_row_names = TRUE,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 7.2, fontfamily = FONT_FAMILY),
  row_names_max_width = unit(18, "mm"),
  height = unit(3.8, "mm"),
  column_gap = unit(0.35, "mm"),
  rect_gp = gpar(col = NA),
  border = TRUE,
  border_gp = gpar(col = "black", lwd = 0.02),
  heatmap_legend_param = list(
    title = "Cluster",
    title_gp = gpar(fontsize = 7.2, fontfamily = FONT_FAMILY, fontface = "plain"),
    labels_gp = gpar(fontsize = 7.0, fontfamily = FONT_FAMILY),
    grid_width = unit(3.2, "mm"),
    grid_height = unit(3.2, "mm")
  )
)
pseudotime_bar <- Heatmap(
  matrix(as.numeric(pt_ordered), nrow = 1, dimnames = list("Pseudotime", bar_columns)),
  name = "Pseudotime",
  col = pseudotime_col_fun,
  column_split = cluster_ordered,
  cluster_column_slices = FALSE,
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  show_column_names = FALSE,
  column_title = NULL,
  show_row_names = TRUE,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 7.2, fontfamily = FONT_FAMILY),
  row_names_max_width = unit(18, "mm"),
  height = unit(3.8, "mm"),
  column_gap = unit(0.35, "mm"),
  rect_gp = gpar(col = NA),
  border = TRUE,
  border_gp = gpar(col = "black", lwd = 0.02),
  heatmap_legend_param = list(
    title = "Pseudotime",
    title_gp = gpar(fontsize = 7.2, fontfamily = FONT_FAMILY, fontface = "plain"),
    labels_gp = gpar(fontsize = 7.0, fontfamily = FONT_FAMILY),
    legend_height = unit(22, "mm")
  )
)
ht <- Heatmap(
  mat_cc_scaled,
  name = "Z-score",
  col = colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B")),
  row_split = modules,
  cluster_row_slices = FALSE,
  column_split = cluster_ordered,
  cluster_columns = FALSE,
  cluster_rows = TRUE,
  show_row_dend = FALSE,
  right_annotation = row_anno,
  show_column_names = FALSE,
  show_row_names = FALSE,
  use_raster = TRUE,
  raster_quality = 3,
  column_gap = unit(0.35, "mm"),
  row_title = NULL,
  column_title = NULL,
  heatmap_legend_param = list(
    title = "Z-score",
    title_gp = gpar(fontsize = 7.2, fontfamily = FONT_FAMILY, fontface = "plain"),
    labels_gp = gpar(fontsize = 7.0, fontfamily = FONT_FAMILY),
    legend_height = unit(22, "mm")
  )
)
draw(
  cluster_bar %v% pseudotime_bar %v% ht,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  ht_gap = unit(c(0.15, 0.75), "mm")
)
close_panel_pdf()

export_large_panel_pngs(OUT_DIR)
message("Figure 3 nature panels written to: ", OUT_DIR)
