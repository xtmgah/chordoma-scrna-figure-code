script_dir <- if (exists("SCRIPT_DIR", inherits = FALSE)) {
  SCRIPT_DIR
} else {
  dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), mustWork = TRUE))
}
source(file.path(dirname(script_dir), "_figure_style.R"))

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggpubr)
  library(rstatix)
  library(patchwork)
})
if (requireNamespace("scCustomize", quietly = TRUE)) {
  suppressPackageStartupMessages(library(scCustomize))
}

OUT_DIR <- file.path(script_dir, "nature_panels")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

message("Figure 1: loading atlas Seurat object")
sc <- readRDS(file.path(script_dir, "sc_combined_clean.rds"))
meta <- sc@meta.data
celltype_col <- if ("cell_type" %in% colnames(meta)) "cell_type" else stop("No `cell_type` column found in Figure1 Seurat metadata.")
sample_col <- if ("orig.ident" %in% colnames(meta)) "orig.ident" else stop("No `orig.ident` column found in Figure1 Seurat metadata.")
cluster_col <- if ("RNA_snn_res.0.5" %in% colnames(meta)) "RNA_snn_res.0.5" else grep("snn_res", colnames(meta), value = TRUE)[1]
group_col <- intersect(c("clinical_group", "Subtype", "subtype"), colnames(meta))[1]

cell_types <- sort(unique(as.character(meta[[celltype_col]])))
cell_cols <- palette_for(cell_types, cell_type_palette)
sc@meta.data[[celltype_col]] <- factor(as.character(sc@meta.data[[celltype_col]]), levels = cell_types)
meta[[celltype_col]] <- factor(as.character(meta[[celltype_col]]), levels = cell_types)
celltype_plot_col <- "Cell type"
sc@meta.data[[celltype_plot_col]] <- sc@meta.data[[celltype_col]]
sample_ids <- sort(unique(as.character(meta[[sample_col]])))
sample_cols <- setNames(nature_discrete(length(sample_ids)), sample_ids)

if (!is.na(cluster_col)) {
  cluster_levels <- sort_numeric_labels(unique(as.character(meta[[cluster_col]])))
  sc@meta.data[[cluster_col]] <- factor(as.character(sc@meta.data[[cluster_col]]), levels = cluster_levels)
  cluster_cols <- setNames(nature_discrete(length(cluster_levels)), cluster_levels)
  p_cluster <- if (requireNamespace("scCustomize", quietly = TRUE)) {
    DimPlot_scCustom(
      sc,
      reduction = "umap",
      group.by = cluster_col,
      label = TRUE,
      repel = TRUE,
      label.size = 3.2,
      pt.size = 0.10,
      colors_use = cluster_cols,
      raster = FALSE
    )
  } else {
    DimPlot(
      sc,
      reduction = "umap",
      group.by = cluster_col,
      label = TRUE,
      repel = TRUE,
      label.size = 3.2,
      pt.size = 0.10,
      cols = cluster_cols,
      raster = FALSE
    )
  }
  p_cluster <- p_cluster +
    labs(title = NULL, color = "Cluster") +
    theme_umap_nature() +
    guides(color = guide_legend(ncol = 2, override.aes = list(size = 3.0))) +
    theme(legend.position = "right", legend.text = element_text(size = 7.4))
  save_panel_pdf(p_cluster, file.path(OUT_DIR, "fig1a_umap_by_unsupervised_cluster.pdf"), 4.7, 4.0)

  cluster_sample <- meta %>%
    count(.data[[cluster_col]], .data[[sample_col]], name = "n") %>%
    group_by(.data[[cluster_col]]) %>%
    mutate(freq = n / sum(n)) %>%
    ungroup()
  colnames(cluster_sample)[1:2] <- c("cluster", "sample")
  cluster_sample$cluster <- factor(cluster_sample$cluster, levels = cluster_levels)
  cluster_sample$sample <- factor(cluster_sample$sample, levels = sample_ids)
  p_cluster_sample <- ggplot(cluster_sample, aes(x = cluster, y = freq, fill = sample)) +
    geom_col(width = 0.72, color = "white", linewidth = 0.08) +
    scale_fill_manual(values = sample_cols, drop = FALSE) +
    percent_axis() +
    labs(x = "Cluster", y = "Cell fraction", fill = "Sample") +
    theme_nature() +
    theme(panel.grid.major.x = element_blank(), legend.key.size = unit(2.7, "mm"))
  save_panel_pdf(p_cluster_sample, file.path(OUT_DIR, "fig1b_sample_fraction_by_cluster.pdf"), 5.7, 3.2)
}

message("Figure 1: annotated atlas panels")
p_anno <- if (requireNamespace("scCustomize", quietly = TRUE)) {
  DimPlot_scCustom(
    sc,
    reduction = "umap",
    group.by = celltype_plot_col,
    label = TRUE,
    repel = TRUE,
    label.size = 3.4,
    pt.size = 0.10,
    colors_use = cell_cols,
    raster = FALSE
  )
} else {
  DimPlot(
    sc,
    reduction = "umap",
    group.by = celltype_plot_col,
    label = TRUE,
    repel = TRUE,
    label.size = 3.4,
    pt.size = 0.10,
    cols = cell_cols,
    raster = FALSE
  )
}
p_anno <- p_anno +
  labs(title = NULL, color = "Cell type") +
  theme_umap_nature() +
  guides(color = guide_legend(title = "Cell type", override.aes = list(size = 3.1))) +
  theme(legend.position = "right")
save_panel_pdf(p_anno, file.path(OUT_DIR, "fig1c_annotated_celltype_umap.pdf"), 5.2, 4.3)

p_sample <- if (requireNamespace("scCustomize", quietly = TRUE)) {
  DimPlot_scCustom(
    sc,
    reduction = "umap",
    group.by = sample_col,
    colors_use = sample_cols,
    pt.size = 0.05,
    shuffle = TRUE,
    raster = FALSE
  )
} else {
  DimPlot(
    sc,
    reduction = "umap",
    group.by = sample_col,
    cols = sample_cols,
    pt.size = 0.05,
    shuffle = TRUE,
    raster = FALSE
  )
}
p_sample <- p_sample +
  labs(title = NULL, color = "Sample") +
  theme_umap_nature() +
  theme(legend.position = "right", legend.text = element_text(size = 6.2))
save_panel_pdf(p_sample, file.path(OUT_DIR, "fig1d_umap_by_sample.pdf"), 5.2, 4.3)

sample_cell_props <- meta %>%
  count(.data[[sample_col]], .data[[celltype_col]], name = "n") %>%
  group_by(.data[[sample_col]]) %>%
  mutate(freq = n / sum(n)) %>%
  ungroup()
colnames(sample_cell_props)[1:2] <- c("sample", "cell_type")
sample_cell_props$sample <- factor(sample_cell_props$sample, levels = sample_ids)
sample_cell_props$cell_type <- factor(sample_cell_props$cell_type, levels = cell_types)
p_sample_comp <- ggplot(sample_cell_props, aes(x = sample, y = freq, fill = cell_type)) +
  geom_col(width = 0.75, color = "white", linewidth = 0.08) +
  scale_fill_manual(values = cell_cols, drop = FALSE) +
  percent_axis() +
  labs(x = NULL, y = "Cell fraction", fill = "Cell type") +
  theme_nature(base_size = 9.0) +
  theme(
    axis.title.y = element_text(size = 9.6),
    axis.text.x = element_text(angle = 55, hjust = 1, vjust = 1, size = 7.2),
    axis.text.y = element_text(size = 7.8),
    legend.title = element_text(size = 8.5),
    legend.text = element_text(size = 7.4),
    legend.key.height = unit(4.0, "mm"),
    legend.spacing.y = unit(1.3, "mm"),
    panel.grid.major.x = element_blank()
  )
save_panel_pdf(p_sample_comp, file.path(OUT_DIR, "fig1e_celltype_fraction_by_sample.pdf"), 7.0, 3.7)

markers_to_plot <- c(
  "TBXT", "KRT19", "CD24",
  "MKI67", "TOP2A",
  "CD68", "CD163", "C1QA", "CD86",
  "CTSK", "ACP5",
  "CLEC9A",
  "S100A8",
  "TPSAB1",
  "CD3D", "CD3E",
  "MS4A1", "CD79A",
  "JCHAIN",
  "COL1A1", "DCN",
  "PECAM1", "VWF",
  "ACTA2"
)
markers_to_plot <- intersect(markers_to_plot, rownames(sc))
if (length(markers_to_plot) > 0) {
  Idents(sc) <- celltype_col
  p_dot <- DotPlot(sc, features = markers_to_plot, group.by = celltype_col, cols = c("grey94", "#9D3B62"), dot.scale = 4.3) +
    RotatedAxis() +
    labs(x = NULL, y = "Cell type", color = "Average expression", size = "Expressing cells") +
    theme_nature() +
    theme(
      axis.text.x = element_text(size = 6.6, face = "italic", angle = 90, hjust = 1, vjust = 0.5),
      axis.text.y = element_text(size = 7.2),
      panel.grid = element_blank(),
      legend.position = "right"
    )
  save_panel_pdf(p_dot, file.path(OUT_DIR, "fig1f_celltype_marker_dotplot.pdf"), 5.4, 3.6)
}

feature_genes <- intersect(c("TBXT", "MKI67", "COL1A1", "CD3D", "MS4A1", "JCHAIN", "CD163", "CTSK", "S100A8", "TPSAB1", "VWF", "ACTA2"), rownames(sc))
if (length(feature_genes) > 0) {
  feature_colors <- viridis::magma(60, begin = 0.08, end = 0.95)
  feature_plots <- if (requireNamespace("scCustomize", quietly = TRUE)) {
    FeaturePlot_scCustom(
      sc,
      features = feature_genes,
      colors_use = feature_colors,
      na_color = "grey92",
      order = TRUE,
      min.cutoff = "q05",
      max.cutoff = "q95",
      pt.size = 0.013,
      raster = FALSE,
      combine = FALSE
    )
  } else {
    FeaturePlot(
      sc,
      features = feature_genes,
      cols = c("grey92", "#B2182B"),
      order = TRUE,
      min.cutoff = "q05",
      max.cutoff = "q95",
      pt.size = 0.013,
      combine = FALSE,
      raster = FALSE
    )
  }
  feature_plots <- Map(function(p, gene) {
    p + labs(title = gene) + theme_umap_nature() + theme(plot.title = element_text(face = "bold", hjust = 0.5))
  }, feature_plots, feature_genes)
  for (i in seq_along(feature_plots)) {
    save_panel_pdf(feature_plots[[i]], file.path(OUT_DIR, paste0("fig1g_feature_", feature_genes[[i]], ".pdf")), 2.55, 2.45)
  }
  save_panel_pdf(wrap_plots(feature_plots, ncol = 4), file.path(OUT_DIR, "fig1g_marker_featureplots_combined.pdf"), 8.2, 6.4)
}

if (!is.na(group_col)) {
  group_values <- unique(as.character(meta[[group_col]]))
  group_cols <- group_palette[intersect(names(group_palette), group_values)]
  if (length(group_cols) < length(group_values)) {
    missing <- setdiff(group_values, names(group_cols))
    group_cols <- c(group_cols, nature_named_palette(missing))
  }

  grouped_comp <- meta %>%
    filter(!is.na(.data[[group_col]])) %>%
    count(.data[[group_col]], .data[[celltype_col]], name = "n") %>%
    group_by(.data[[group_col]]) %>%
    mutate(freq = n / sum(n)) %>%
    ungroup()
  colnames(grouped_comp)[1:2] <- c("pathology", "cell_type")
  grouped_comp$cell_type <- factor(grouped_comp$cell_type, levels = cell_types)
  grouped_comp$pathology <- factor(grouped_comp$pathology, levels = intersect(c("CC", "PC", "DC"), unique(grouped_comp$pathology)))
  x_positions <- setNames(c(1.00, 1.36, 1.72)[seq_along(levels(grouped_comp$pathology))], levels(grouped_comp$pathology))
  grouped_comp$xpos <- unname(x_positions[as.character(grouped_comp$pathology)])
  p_grouped <- ggplot(grouped_comp, aes(x = xpos, y = freq, fill = cell_type)) +
    geom_col(width = 0.24, color = "white", linewidth = 0.18) +
    scale_fill_manual(values = cell_cols, drop = FALSE) +
    scale_x_continuous(
      breaks = unname(x_positions),
      labels = names(x_positions),
      limits = range(unname(x_positions)) + c(-0.18, 0.18),
      expand = expansion(mult = 0)
    ) +
    percent_axis() +
    labs(x = NULL, y = "Cell fraction", fill = "Cell type") +
    theme_nature(base_size = 7.4) +
    theme(
      axis.title.y = element_text(size = 8.0),
      axis.text = element_text(size = 7.2),
      legend.title = element_text(size = 7.6),
      legend.text = element_text(size = 6.7),
      legend.key.size = unit(2.8, "mm"),
      legend.spacing.y = unit(0.7, "mm"),
      panel.grid.major.x = element_blank()
    )
  save_panel_pdf(p_grouped, file.path(OUT_DIR, "fig1h_celltype_fraction_by_pathology.pdf"), 2.15, 3.0)

  p_split <- if (requireNamespace("scCustomize", quietly = TRUE)) {
    DimPlot_scCustom(
      sc,
      reduction = "umap",
      group.by = celltype_plot_col,
      split.by = group_col,
      colors_use = cell_cols,
      label = TRUE,
      repel = TRUE,
      label.size = 2.7,
      pt.size = 0.08,
      raster = FALSE
    )
  } else {
    DimPlot(
      sc,
      reduction = "umap",
      group.by = celltype_plot_col,
      split.by = group_col,
      cols = cell_cols,
      label = TRUE,
      repel = TRUE,
      label.size = 2.7,
      pt.size = 0.08,
      raster = FALSE
    )
  }
  p_split <- p_split +
    labs(title = NULL, color = "Cell type", fill = "Cell type") +
    theme_umap_nature() +
    guides(color = guide_legend(title = "Cell type", override.aes = list(size = 3.1))) +
    theme(legend.position = "right", legend.title = element_text(size = 8.8), legend.text = element_text(size = 7.6))
  save_panel_pdf(p_split, file.path(OUT_DIR, "fig1i_celltype_umap_split_by_pathology.pdf"), 7.4, 3.9)
}

sample_props_path <- file.path(script_dir, "sample_props_complete.rds")
if (file.exists(sample_props_path)) {
  sample_props <- readRDS(sample_props_path)
  if (all(c("clinical_group", "freq", "cell_type") %in% colnames(sample_props))) {
    sample_props$clinical_group <- factor(sample_props$clinical_group, levels = c("CC", "PC"))
    stat_test <- sample_props %>%
      filter(clinical_group %in% c("CC", "PC")) %>%
      group_by(cell_type) %>%
      wilcox_test(freq ~ clinical_group) %>%
      adjust_pvalue(method = "BH") %>%
      add_significance("p.adj") %>%
      left_join(
        sample_props %>%
          group_by(cell_type) %>%
          summarise(y.position = max(freq, na.rm = TRUE) * 1.07, .groups = "drop"),
        by = "cell_type"
      ) %>%
      mutate(
        p.adj.label = ifelse(p.adj < 0.001, "FDR < 0.001", paste0("FDR = ", signif(p.adj, 2))),
        p.color = ifelse(p.adj < 0.05, "#BB0E3D", "black")
      )

    point_group_palette <- c("CC" = "#007C92", "PC" = "#4F2A6E")
    p_box <- ggplot(sample_props, aes(x = clinical_group, y = freq, fill = clinical_group)) +
      geom_boxplot(outlier.shape = NA, width = 0.44, linewidth = 0.26, alpha = 0.34) +
      geom_point(
        aes(fill = clinical_group),
        position = position_jitter(width = 0.07, height = 0),
        shape = 21,
        stroke = 0.08,
        color = "grey15",
        size = 1.15,
        alpha = 0.95
      ) +
      facet_wrap(~ cell_type, scales = "free_y", ncol = 4) +
      stat_pvalue_manual(stat_test, "p.adj.label", y.position = "y.position", color = "p.color", tip.length = 0.005, bracket.size = 0.22, size = 2.45, family = FONT_FAMILY) +
      scale_color_identity() +
      scale_fill_manual(values = point_group_palette) +
      scale_y_continuous(labels = label_percent(accuracy = 1), expand = expansion(mult = c(0.04, 0.10))) +
      labs(x = NULL, y = "Cell fraction") +
      theme_nature(base_size = 7.5) +
      theme(
        legend.position = "none",
        axis.title.y = element_text(size = 9.6),
        axis.text.y = element_text(size = 8.4),
        panel.grid.major.x = element_blank(),
        strip.text = element_text(size = 7.2, margin = margin(2.4, 2, 2.4, 2)),
        strip.background = element_rect(fill = "grey92", color = NA)
      )
    save_panel_pdf(p_box, file.path(OUT_DIR, "fig1j_celltype_fraction_cc_vs_pc_boxplots.pdf"), 6.2, 6.15)
  }
}

cell_type_props_path <- file.path(script_dir, "cell_type_props.rds")
if (file.exists(cell_type_props_path)) {
  cell_type_props <- readRDS(cell_type_props_path)
  if (all(c("clinical_group", "freq", "cell_type") %in% colnames(cell_type_props))) {
    cell_type_props$clinical_group <- factor(cell_type_props$clinical_group, levels = c("CC", "PC"))
    p_sorted <- ggplot(cell_type_props, aes(y = cell_type, x = freq, fill = clinical_group)) +
      geom_col(position = "fill", width = 0.72, color = "white", linewidth = 0.18) +
      scale_fill_manual(values = group_palette[c("CC", "PC")]) +
      scale_x_continuous(labels = label_percent(accuracy = 1), expand = expansion(mult = c(0, 0.02))) +
      labs(x = "Cell fraction", y = NULL, fill = "Group") +
      theme_nature() +
      theme(panel.grid.major.y = element_blank(), axis.text.y = element_text(size = 7.1))
    save_panel_pdf(p_sorted, file.path(OUT_DIR, "fig1k_sorted_celltype_fraction_cc_vs_pc.pdf"), 3.8, 2.8)
  }
}

message("Figure 1 nature panels written to: ", OUT_DIR)
