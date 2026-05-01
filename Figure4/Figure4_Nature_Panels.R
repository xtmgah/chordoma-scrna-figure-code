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
  library(ggrepel)
  library(scales)
})
if (requireNamespace("scCustomize", quietly = TRUE)) {
  suppressPackageStartupMessages(library(scCustomize))
}

OUT_DIR <- file.path(script_dir, "nature_panels")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

plot_auc_subset <- function(asset_rds, prefix, group_label, composition_label, marker_features = NULL, color_override = NULL) {
  assets <- readRDS(asset_rds)
  obj <- assets$object
  markers <- assets$markers
  annotation_map <- assets$annotation_map
  group_col <- assets$parameters$group_col
  label_col <- assets$parameters$new_label_col %||% group_col
  fc_col <- assets$parameters$fc_col
  color_map <- color_override %||% assets$color_map
  if (is.null(names(color_map))) {
    color_map <- nature_named_palette(sort(unique(obj@meta.data[[label_col]])))
  }
  color_map <- color_map[!is.na(names(color_map))]
  obj@meta.data[[label_col]] <- factor(obj@meta.data[[label_col]], levels = names(color_map))
  markers$cluster <- factor(markers$cluster, levels = unique(markers$cluster))
  markers$plot_auc <- ifelse(markers$myAUC >= 0.5, markers$myAUC, 1 - markers$myAUC)

  label_df <- markers %>%
    filter(myAUC >= 0.55) %>%
    group_by(cluster) %>%
    arrange(desc(plot_auc), desc(.data[[fc_col]]), .by_group = TRUE) %>%
    slice_head(n = 4) %>%
    ungroup()

  p_umap <- if (requireNamespace("scCustomize", quietly = TRUE)) {
    DimPlot_scCustom(
      obj,
      reduction = "umap",
      group.by = label_col,
      colors_use = color_map,
      label = TRUE,
      repel = TRUE,
      label.size = 2.7,
      pt.size = 0.18,
      raster = FALSE
    )
  } else {
    DimPlot(
      obj,
      reduction = "umap",
      group.by = label_col,
      cols = color_map,
      label = TRUE,
      repel = TRUE,
      label.size = 2.7,
      pt.size = 0.18,
      raster = FALSE
    )
  }
  p_umap <- p_umap +
    labs(title = NULL, color = group_label) +
    theme_umap_nature() +
    theme(legend.position = "right")
  save_panel_pdf(p_umap, file.path(OUT_DIR, paste0(prefix, "_annotated_umap.pdf")), 4.7, 3.7)

  p_butterfly <- ggplot(markers, aes(x = .data[[fc_col]], y = plot_auc)) +
    geom_point(aes(fill = plot_auc, size = plot_auc), shape = 21, color = "black", stroke = 0.05, alpha = 0.82) +
    geom_hline(yintercept = 0.70, linetype = "dashed", color = "grey55", linewidth = 0.35) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.35) +
    geom_text_repel(
      data = label_df,
      aes(label = gene),
      size = 2.9,
      family = FONT_FAMILY,
      box.padding = 0.25,
      segment.alpha = 0.45,
      max.overlaps = Inf
    ) +
    viridis::scale_fill_viridis(option = "D", direction = -1, name = "AUC", guide = guide_colorbar(barheight = unit(28, "mm"), barwidth = unit(4.2, "mm"))) +
    scale_size(range = c(0.65, 3.0), guide = "none") +
    facet_wrap(~ cluster, scales = "free_x", ncol = 3) +
    labs(x = "Log2 fold change", y = "AUC") +
    theme_nature(base_size = 9.5) +
    theme(
      axis.text = element_text(size = 8.7),
      axis.title.x = element_text(size = 10.8),
      axis.title.y = element_text(size = 10.8),
      legend.title = element_text(size = 9.0),
      legend.text = element_text(size = 8.4),
      legend.key.height = unit(5.0, "mm"),
      strip.text = element_text(size = 8.8, margin = margin(2, 3, 2, 3)),
      panel.grid.major.x = element_blank(),
      legend.position = "right"
    )
  save_panel_pdf(p_butterfly, file.path(OUT_DIR, paste0(prefix, "_auc_butterfly.pdf")), 7.2, 5.1)

  if ("clinical_group" %in% colnames(obj@meta.data)) {
    comp_df <- obj@meta.data %>%
      filter(clinical_group %in% c("CC", "PC")) %>%
      count(clinical_group, .data[[label_col]], name = "n") %>%
      group_by(clinical_group) %>%
      mutate(freq = n / sum(n), label = ifelse(freq >= 0.035, paste0(round(freq * 100), "%"), "")) %>%
      ungroup()
    colnames(comp_df)[2] <- "annotation_label"
    comp_df$annotation_label <- factor(comp_df$annotation_label, levels = names(color_map))
    comp_df$clinical_group <- factor(comp_df$clinical_group, levels = c("CC", "PC"))

    p_comp <- ggplot(comp_df, aes(x = clinical_group, y = freq, fill = annotation_label)) +
      geom_col(width = 0.46, color = "white", linewidth = 0.22) +
      geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 2.6, color = "white", family = FONT_FAMILY) +
      scale_fill_manual(values = color_map, drop = FALSE) +
      scale_y_continuous(labels = label_percent(accuracy = 1), expand = expansion(mult = c(0, 0.02))) +
      labs(x = NULL, y = "Cell fraction", fill = composition_label) +
      theme_nature() +
      theme(panel.grid.major.x = element_blank())
    save_panel_pdf(p_comp, file.path(OUT_DIR, paste0(prefix, "_composition_cc_vs_pc.pdf")), 2.8, 3.0)
  }

  if (is.null(marker_features)) {
    marker_features <- lapply(seq_len(nrow(annotation_map)), function(i) {
      cluster_id <- annotation_map$source_cluster[i] %||% annotation_map$subtype[i]
      selected <- annotation_map$selected_gene[i]
      genes <- markers %>%
        filter(cluster == cluster_id) %>%
        arrange(desc(plot_auc), desc(.data[[fc_col]])) %>%
        pull(gene)
      unique(c(selected, genes))[seq_len(min(3, length(unique(c(selected, genes)))))]
    })
    names(marker_features) <- annotation_map$recommended_annotation %||% annotation_map$subtype
  }
  marker_features <- lapply(marker_features, function(x) intersect(x, rownames(obj)))
  marker_features <- marker_features[lengths(marker_features) > 0]
  Idents(obj) <- label_col
  p_dot <- DotPlot(obj, features = marker_features, dot.scale = 4.2) +
    RotatedAxis() +
    scale_color_gradientn(colours = c("grey96", "#FCAEA1", "#B2182B"), name = "Expression") +
    guides(color = guide_colorbar(order = 1), size = guide_legend(order = 2)) +
    labs(x = NULL, y = group_label) +
    theme_nature() +
    theme(
      axis.text.x = element_text(size = 7.2, face = "italic", angle = 90, hjust = 1, vjust = 0.5),
      axis.text.y = element_text(size = 7.2),
      strip.text.x = element_text(size = 6.8, margin = margin(2, 4, 2, 4)),
      strip.clip = "off",
      legend.spacing.y = unit(3.0, "mm"),
      legend.box.spacing = unit(2.0, "mm"),
      panel.grid = element_blank(),
      plot.margin = margin(4, 10, 6, 10)
    )
  if (identical(prefix, "fig4_tcell") && "feature.groups" %in% colnames(p_dot$data)) {
    p_dot <- p_dot + facet_wrap(~ feature.groups, scales = "free_x", nrow = 2)
  }
  save_panel_pdf(p_dot, file.path(OUT_DIR, paste0(prefix, "_marker_dotplot.pdf")), 5.7, if (identical(prefix, "fig4_tcell")) 3.6 else 3.0)
}

`%||%` <- function(x, y) if (is.null(x)) y else x

message("Figure 4: macrophage panels")
plot_auc_subset(
  asset_rds = file.path(script_dir, "Mac", "plot_ready_rds", "macrophage_auc_plot_ready.rds"),
  prefix = "fig4_macrophage",
  group_label = "Macrophage subtype",
  composition_label = "Macrophage subtype",
  color_override = macrophage_palette
)

message("Figure 4: T-cell panels")
tcell_colors <- tcell_palette
tcell_markers <- list(
  "CD4 memory/helper T" = c("IL7R", "CCR7", "CD40LG"),
  "Effector-memory T" = c("PRKCQ-AS1", "PCED1B-AS1", "PLAAT4"),
  "Activated cytotoxic T" = c("LAG3", "GZMK", "CCL4"),
  "NK-like cytotoxic T" = c("FGFBP2", "GNLY", "KLRD1"),
  "Treg" = c("FOXP3", "IL2RA", "CCR8"),
  "MAIT-like T" = c("TRAV1-2", "SLC4A10", "KLRB1")
)
plot_auc_subset(
  asset_rds = file.path(script_dir, "T_cell", "plot_ready_rds", "tcell_auc_plot_ready.rds"),
  prefix = "fig4_tcell",
  group_label = "T-cell subtype",
  composition_label = "T-cell subtype",
  marker_features = tcell_markers,
  color_override = tcell_colors
)

message("Figure 4: CC4 CellChat panels")
if (requireNamespace("CellChat", quietly = TRUE)) {
  suppressPackageStartupMessages(library(CellChat))
  cc4_dir <- file.path(script_dir, "subcell_commu_cc4")
  assets <- readRDS(file.path(cc4_dir, "plot_ready_rds", "subcell_cc4_cellchat_input_assets.rds"))
  color_map <- assets$color_map
  cellchat_overall <- readRDS(file.path(cc4_dir, "subcell_cc4_cellchat_overall.rds"))
  local_colors <- palette_for(levels(cellchat_overall@idents), c(cell_type_palette, macrophage_palette, tcell_palette, color_map))
  names(local_colors) <- levels(cellchat_overall@idents)

  p_scatter <- netAnalysis_signalingRole_scatter(
    cellchat_overall,
    color.use = local_colors,
    do.label = TRUE,
    label.size = 2.7
  ) +
    labs(x = "Outgoing signaling strength", y = "Incoming signaling strength", color = "Cells") +
    theme_nature() +
    theme(legend.position = "right")
  save_panel_pdf(p_scatter, file.path(OUT_DIR, "fig4_cellchat_cc4_signaling_role_scatter.pdf"), 4.9, 3.8)

  group_size <- as.numeric(table(cellchat_overall@idents))
  names(group_size) <- levels(cellchat_overall@idents)
  open_panel_pdf(file.path(OUT_DIR, "fig4_cellchat_cc4_network_count.pdf"), 4.0, 4.0)
  par(family = FONT_FAMILY, mar = c(0.5, 0.5, 1.2, 0.5), xpd = TRUE)
  netVisual_circle(cellchat_overall@net$count, vertex.weight = group_size, weight.scale = TRUE, label.edge = FALSE, title.name = "Number of interactions", color.use = local_colors, vertex.label.cex = 0.58, edge.width.max = 2.6, alpha.edge = 0.32, arrow.size = 0.12, arrow.width = 0.6)
  close_panel_pdf()

  open_panel_pdf(file.path(OUT_DIR, "fig4_cellchat_cc4_network_strength.pdf"), 4.0, 4.0)
  par(family = FONT_FAMILY, mar = c(0.5, 0.5, 1.2, 0.5), xpd = TRUE)
  netVisual_circle(cellchat_overall@net$weight, vertex.weight = group_size, weight.scale = TRUE, label.edge = FALSE, title.name = "Interaction strength", color.use = local_colors, vertex.label.cex = 0.58, edge.width.max = 2.6, alpha.edge = 0.32, arrow.size = 0.12, arrow.width = 0.6)
  close_panel_pdf()

  merged_path <- file.path(cc4_dir, "subcell_cc4_cellchat_CC_vs_PC_merged.rds")
  if (file.exists(merged_path)) {
    merged_cellchat <- readRDS(merged_path)
    open_panel_pdf(file.path(OUT_DIR, "fig4_cellchat_cc4_cc_vs_pc_interactions.pdf"), 3.7, 3.0)
    axis_thin <- theme(
      axis.line = element_line(linewidth = 0.4, color = "black"),
      axis.ticks = element_line(linewidth = 0.4, color = "black"),
      axis.title = element_text(face = "plain"),
      legend.title = element_text(face = "plain"),
      panel.grid = element_blank()
    )
    interaction_fill <- setNames(group_palette[c("CC", "PC")], c(1, 2))
    p1 <- compareInteractions(merged_cellchat, show.legend = FALSE, group = c(1, 2), color.use = group_palette[c("CC", "PC")], width = 0.36, size.text = 8, xlabel = NULL) +
      scale_fill_manual(values = interaction_fill, guide = "none") +
      theme_nature(base_size = 8.2) + axis_thin
    p2 <- compareInteractions(merged_cellchat, show.legend = FALSE, group = c(1, 2), measure = "weight", color.use = group_palette[c("CC", "PC")], width = 0.36, size.text = 8, xlabel = NULL) +
      scale_fill_manual(values = interaction_fill, guide = "none") +
      theme_nature(base_size = 8.2) + axis_thin
    print(p1 + p2)
    close_panel_pdf()
  }
} else {
  message("  CellChat is not installed; skipped CC4 communication panels.")
}

message("Figure 4 nature panels written to: ", OUT_DIR)
