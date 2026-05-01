options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(scales)
})

OUT_DIR <- "F:/Chordoma/Result/Mac"
PLOT_READY_RDS <- file.path(OUT_DIR, "plot_ready_rds", "macrophage_auc_plot_ready.rds")

UMAP_PDF <- file.path(OUT_DIR, "Macrophage_AUC_Annotated_UMAP.pdf")
BUTTERFLY_PDF <- file.path(OUT_DIR, "Macrophage_AUC_Butterfly.pdf")
COMPOSITION_PDF <- file.path(OUT_DIR, "Macrophage_AUC_Composition_CC_vs_PC.pdf")
DOTPLOT_PDF <- file.path(OUT_DIR, "Macrophage_AUC_DotPlot.pdf")

TOP_LABELS_PER_CLUSTER <- 4
MIN_AUC_TO_SHOW <- 0.55

message("Reading macrophage plot-ready assets...")
assets <- readRDS(PLOT_READY_RDS)
sc_mac <- assets$object
markers <- assets$markers
annotation_map <- assets$annotation_map
color_map <- assets$color_map
group_col <- assets$parameters$group_col
label_col <- assets$parameters$new_label_col
fc_col <- assets$parameters$fc_col

markers$cluster <- factor(markers$cluster, levels = annotation_map$source_cluster)
markers$plot_auc <- ifelse(markers$myAUC >= 0.5, markers$myAUC, 1 - markers$myAUC)

label_df <- markers %>%
  filter(myAUC >= MIN_AUC_TO_SHOW) %>%
  group_by(cluster) %>%
  arrange(desc(plot_auc), desc(.data[[fc_col]]), .by_group = TRUE) %>%
  slice_head(n = TOP_LABELS_PER_CLUSTER) %>%
  ungroup() %>%
  bind_rows(
    annotation_map %>%
      transmute(cluster = factor(source_cluster, levels = annotation_map$source_cluster), gene = selected_gene) %>%
      left_join(markers, by = c("cluster", "gene"))
  ) %>%
  distinct(cluster, gene, .keep_all = TRUE) %>%
  filter(!is.na(plot_auc), !is.na(.data[[fc_col]]))

message("Plotting macrophage UMAP...")
p_umap <- DimPlot(
  sc_mac,
  reduction = "umap",
  group.by = label_col,
  cols = color_map,
  label = TRUE,
  repel = TRUE,
  label.size = 4.5,
  pt.size = 0.18,
  raster = FALSE
) +
  ggtitle("Macrophage strict AUC-based annotation") +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right"
  )
ggsave(UMAP_PDF, p_umap, width = 10, height = 7)

message("Plotting macrophage butterfly figure...")
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
  scale_color_gradientn(colors = c("grey85", "#FFEDA0", "#FEB24C", "#F03B20")) +
  scale_size(range = c(0.5, 3.5), guide = "none") +
  facet_wrap(~ cluster, scales = "free_x", ncol = 3) +
  theme_classic() +
  labs(
    title = "Macrophage subtype AUC butterfly plot",
    subtitle = "One-vs-rest ROC markers across macrophage subclusters with strict AUC-based labels",
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

if ("clinical_group" %in% colnames(sc_mac@meta.data)) {
  message("Plotting macrophage CC vs PC composition...")
  comp_df <- sc_mac@meta.data %>%
    filter(clinical_group %in% c("CC", "PC")) %>%
    count(clinical_group, .data[[label_col]], name = "n") %>%
    group_by(clinical_group) %>%
    mutate(freq = n / sum(n), label = ifelse(freq >= 0.03, paste0(round(freq * 100, 1), "%"), "")) %>%
    ungroup()
  colnames(comp_df)[2] <- "annotation_label"
  comp_df$annotation_label <- factor(comp_df$annotation_label, levels = names(color_map))

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
      title = "Macrophage composition in CC and PC",
      x = NULL,
      y = "Proportion",
      fill = "Macrophage subtype"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.x = element_text(face = "bold")
    )
  ggsave(COMPOSITION_PDF, p_comp, width = 8, height = 6.2)
}

message("Plotting macrophage dot plot...")
dot_features <- lapply(seq_len(nrow(annotation_map)), function(i) {
  cluster_id <- annotation_map$source_cluster[i]
  label_name <- annotation_map$recommended_annotation[i]
  genes <- markers %>%
    filter(cluster == cluster_id) %>%
    arrange(desc(plot_auc), desc(.data[[fc_col]])) %>%
    pull(gene)
  genes <- unique(c(annotation_map$selected_gene[i], genes))
  genes <- genes[seq_len(min(3, length(genes)))]
  names(genes) <- NULL
  genes
})
names(dot_features) <- annotation_map$recommended_annotation
dot_features <- dot_features[lengths(dot_features) > 0]

Idents(sc_mac) <- label_col
p_dot <- DotPlot(sc_mac, features = dot_features) +
  RotatedAxis() +
  scale_color_gradientn(colours = c("white", "red", "darkred")) +
  labs(title = "Macrophage subtype marker dot plot", x = NULL, y = "Subtype") +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(size = 9),
    axis.text.y = element_text(size = 10, face = "bold")
  )
ggsave(DOTPLOT_PDF, p_dot, width = 11, height = 6.5)


