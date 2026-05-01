options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(scales)
})

OUT_DIR <- "F:/Chordoma/Result/T_cell"
PLOT_READY_RDS <- file.path(OUT_DIR, "plot_ready_rds", "tcell_auc_plot_ready.rds")

UMAP_PDF <- file.path(OUT_DIR, "Purified_Tcell_AUC_Annotated_UMAP.pdf")
BUTTERFLY_PDF <- file.path(OUT_DIR, "Purified_Tcell_AUC_Butterfly.pdf")
COMPOSITION_PDF <- file.path(OUT_DIR, "Purified_Tcell_AUC_CC_vs_PC_Proportion.pdf")
DOTPLOT_PDF <- file.path(OUT_DIR, "Purified_Tcell_AUC_DotPlot.pdf")

TOP_LABELS_PER_SUBTYPE <- 4
MIN_AUC_TO_SHOW <- 0.55

ZZM60COLORS2 <- c(
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

TCELL_ZZM_INDEX <- c(36, 48, 20, 46, 47, 35)
TCELL_COLORS <- c(
  "CD4 memory/helper T" = ZZM60COLORS2[TCELL_ZZM_INDEX[1]],
  "Effector-memory T" = ZZM60COLORS2[TCELL_ZZM_INDEX[2]],
  "Activated cytotoxic T" = ZZM60COLORS2[TCELL_ZZM_INDEX[3]],
  "NK-like cytotoxic T" = ZZM60COLORS2[TCELL_ZZM_INDEX[4]],
  "Treg" = ZZM60COLORS2[TCELL_ZZM_INDEX[5]],
  "MAIT-like T" = ZZM60COLORS2[TCELL_ZZM_INDEX[6]]
)

message("Reading T-cell plot-ready assets...")
assets <- readRDS(PLOT_READY_RDS)
sc_tcell <- assets$object
markers <- assets$markers
annotation_map <- assets$annotation_map
color_map <- TCELL_COLORS
composition_df <- assets$composition_df
group_col <- assets$parameters$group_col
fc_col <- assets$parameters$fc_col

sc_tcell@meta.data[[group_col]] <- factor(sc_tcell@meta.data[[group_col]], levels = names(color_map))

markers$cluster <- factor(markers$cluster, levels = names(color_map))
markers$plot_auc <- ifelse(markers$myAUC >= 0.5, markers$myAUC, 1 - markers$myAUC)

label_df <- markers %>%
  filter(myAUC >= MIN_AUC_TO_SHOW) %>%
  group_by(cluster) %>%
  arrange(desc(plot_auc), desc(.data[[fc_col]]), .by_group = TRUE) %>%
  slice_head(n = TOP_LABELS_PER_SUBTYPE) %>%
  ungroup() %>%
  bind_rows(
    annotation_map %>%
      transmute(cluster = factor(subtype, levels = names(color_map)), gene = selected_gene) %>%
      left_join(markers, by = c("cluster", "gene"))
  ) %>%
  distinct(cluster, gene, .keep_all = TRUE) %>%
  filter(!is.na(plot_auc), !is.na(.data[[fc_col]]))

message("Plotting T-cell UMAP...")
p_umap <- DimPlot(
  sc_tcell,
  reduction = "umap",
  group.by = group_col,
  cols = color_map,
  label = TRUE,
  repel = TRUE,
  label.size = 4.5,
  pt.size = 0.18,
  raster = FALSE
) +
  ggtitle("Purified T-cell annotation") +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right"
  )
ggsave(UMAP_PDF, p_umap, width = 10.2, height = 7.1)

message("Plotting T-cell butterfly figure...")
p_butterfly <- ggplot(markers, aes(x = .data[[fc_col]], y = plot_auc)) +
  geom_point(aes(color = plot_auc, size = plot_auc), alpha = 0.75) +
  geom_hline(yintercept = 0.70, linetype = "dashed", color = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_text_repel(
    data = label_df,
    aes(label = gene),
    size = 3.1,
    box.padding = 0.35,
    segment.alpha = 0.5,
    max.overlaps = Inf
  ) +
  scale_color_gradientn(colors = c("grey85", "#C6DBEF", "#6BAED6", "#2171B5")) +
  scale_size(range = c(0.5, 3.5), guide = "none") +
  facet_wrap(~ cluster, scales = "free_x", ncol = 3) +
  theme_classic() +
  labs(
    title = "Purified T-cell subtype AUC butterfly plot",
    subtitle = "One-vs-rest ROC markers across refined T-cell subtypes",
    x = "Log2 fold change",
    y = "AUC"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text = element_text(face = "bold")
  )
ggsave(BUTTERFLY_PDF, p_butterfly, width = 12, height = 8.5)

if (!is.null(composition_df)) {
  message("Plotting T-cell CC vs PC composition...")
  comp_df <- composition_df %>%
    mutate(label = ifelse(freq >= 0.025, paste0(round(freq * 100, 1), "%"), ""))

  p_comp <- ggplot(comp_df, aes(x = clinical_group, y = freq, fill = annotation_label)) +
    geom_col(width = 0.72, color = "white", linewidth = 0.2) +
    geom_text(
      aes(label = label),
      position = position_stack(vjust = 0.5),
      size = 3,
      color = "white"
    ) +
    scale_fill_manual(values = color_map, drop = FALSE) +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    labs(
      title = "Purified T-cell composition in CC and PC",
      x = NULL,
      y = "Proportion",
      fill = "Purified T-cell type"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.x = element_text(face = "bold")
    )
  ggsave(COMPOSITION_PDF, p_comp, width = 8.5, height = 6.5)
}

message("Plotting T-cell dot plot...")
marker_list <- list(
  "CD4 memory/helper T" = c("IL7R", "CCR7", "CD40LG"),
  "Effector-memory T" = c("PRKCQ-AS1", "PCED1B-AS1", "PLAAT4"),
  "Activated cytotoxic T" = c("LAG3", "GZMK", "CCL4"),
  "NK-like cytotoxic T" = c("FGFBP2", "GNLY", "KLRD1"),
  "Treg" = c("FOXP3", "IL2RA", "CCR8"),
  "MAIT-like T" = c("TRAV1-2", "SLC4A10", "KLRB1")
)
marker_list <- lapply(marker_list, function(x) intersect(x, rownames(sc_tcell)))
marker_list <- marker_list[names(color_map)]
marker_list <- marker_list[lengths(marker_list) > 0]

Idents(sc_tcell) <- group_col
p_dot <- DotPlot(sc_tcell, features = marker_list) +
  RotatedAxis() +
  scale_color_gradientn(colours = c("white", "red", "darkred")) +
  labs(title = "Purified T-cell marker dot plot", x = NULL, y = "Subtype") +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(size = 9),
    axis.text.y = element_text(size = 10, face = "bold")
  )
ggsave(DOTPLOT_PDF, p_dot, width = 11, height = 6.5)

message("T-cell figures saved to:")
message("  ", UMAP_PDF)
message("  ", BUTTERFLY_PDF)
message("  ", COMPOSITION_PDF)
message("  ", DOTPLOT_PDF)
