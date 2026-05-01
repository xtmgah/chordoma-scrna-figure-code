options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(reshape2)
  library(ggrepel)
})

dc=readRDS("F:/Chordoma/Result/AUC_Reference_Matrix/RNAseq/deconv_17subcells/dc_profile_plot_inputs.rds")
OUT_DIR <- "F:/Chordoma/Result/AUC_Reference_Matrix/RNAseq/deconv_17subcells"
OUT_PDF <- file.path(OUT_DIR, "DC_vs_CC_PC_17subtype_profile_summary.pdf")

resolve_pdf_output <- function(path) {
  if (!file.exists(path)) return(path)
  con <- try(file(path, open = "ab"), silent = TRUE)
  if (inherits(con, "try-error")) {
    return(sub("\\.pdf$", "_latest.pdf", path))
  }
  close(con)
  path
}

celltype_order <- c(
  "Tumor_PC_enriched", "CC_tumor1", "CC_tumor2", "CC_tumor3", "CC_tumor4",
  "FTL+ Mac", "CCL3+ Mac", "FN1+ Mac", "CD1C+ Mac", "CD3E+ Mac", "NAMPT+ Mac",
  "CD4 memory/helper T", "Effector-memory T", "Activated cytotoxic T",
  "NK-like cytotoxic T", "Treg", "MAIT-like T"
)


group_colors <- c(
  "CC" = "#61bada",
  "DC" = "#F28E2B",
  "PC" = "#72567a"
)

celltype_colors <- c(
  "Tumor_PC_enriched" = "#7A3E65",
  "CC_tumor1" = "#D95F5F",
  "CC_tumor2" = "#F29E4C",
  "CC_tumor3" = "#F6D55C",
  "CC_tumor4" = "#8AB17D",
  "FTL+ Mac" = "#8E5EA2",
  "CCL3+ Mac" = "#B565A7",
  "FN1+ Mac" = "#C97B84",
  "CD1C+ Mac" = "#D8A48F",
  "CD3E+ Mac" = "#A06CD5",
  "NAMPT+ Mac" = "#7B2CBF",
  "CD4 memory/helper T" = "#2A9D8F",
  "Effector-memory T" = "#00A6A6",
  "Activated cytotoxic T" = "#118AB2",
  "NK-like cytotoxic T" = "#073B4C",
  "Treg" = "#EF476F",
  "MAIT-like T" = "#5C677D"
)

prediction_colors <- c(
  "CC" = "#4DBBD5",
  "PC" = "#E64B35"
)



theme_pub <- theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, color = "grey30"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  )

prop_fill_scale <- scale_fill_gradientn(
  colours = c("#FFFFFF", "#FFF2E2", "#FDD49E", "#FCAE91", "#FB6A4A", "#CB181D"),
  values = scales::rescale(c(0, 0.01, 0.05, 0.10, 0.20, dc$prop_limit)),
  limits = c(0, dc$prop_limit),
  breaks = dc$prop_breaks,
  labels = scales::label_percent(accuracy = 1),
  name = "Proportion",
  guide = guide_colorbar(
    title.position = "top",
    title.hjust = 0.5,
    barheight = grid::unit(4.0, "cm"),
    barwidth = grid::unit(0.45, "cm"),
    ticks = TRUE,
    label = TRUE
  )
)

p_heat_group <- ggplot(dc$group_mean_long, aes(x = Histo, y = cell_type, fill = mean_prop)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = label), size = 3.0, color = "black") +
  prop_fill_scale +
  labs(
    title = "Group-Level Deconvolution Pattern",
    subtitle = "Absolute mean proportion of each deconvolved subtype",
    x = NULL,
    y = NULL
  ) +
  theme_pub +
  theme(
    axis.text.x = element_text(face = "bold", color = unname(group_colors[c("CC", "DC", "PC")])),
    axis.text.y = element_text(size = 9),
    legend.position = "right"
  )

p_heat_dc <- ggplot(dc$dc_long, aes(x = sample_label, y = cell_type, fill = prop)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = label), size = 2.35, color = "black") +
  prop_fill_scale +
  labs(
    title = "Single-Sample DC Profiles",
    subtitle = "Deconvolved subtype proportions across the 5 DC samples",
    x = NULL,
    y = NULL
  ) +
  theme_pub +
  theme(
    axis.text.x = element_text(size = 9, face = "bold"),
    axis.text.y = element_text(size = 9),
    legend.position = "right"
  )

p_scatter_mac <- ggplot(dc$scatter_df, aes(x = `Tumor_PC_enriched`, y = `CCL3+ Mac`, color = Histo_group)) +
  geom_point(aes(size = point_size, alpha = point_alpha)) +
  ggrepel::geom_text_repel(
    data = dc$scatter_df[dc$scatter_df$Histo_group == "DC", , drop = FALSE],
    aes(label = sample_id),
    size = 3.4,
    box.padding = 0.35,
    point.padding = 0.3,
    segment.color = "grey60",
    max.overlaps = 20,
    show.legend = FALSE
  ) +
  scale_color_manual(values = group_colors) +
  scale_size_identity() +
  scale_alpha_identity() +
  labs(
    title = "PC Signal vs CCL3+ Macrophages",
    subtitle = "DC samples are highlighted against all CC and PC cases",
    x = "Tumor_PC_enriched proportion",
    y = "CCL3+ Mac proportion",
    color = "Pathology"
  ) +
  theme_pub +
  theme(legend.position = "right")

p_scatter_t <- ggplot(dc$scatter_df, aes(x = `Tumor_PC_enriched`, y = `Effector-memory T`, color = Histo_group)) +
  geom_point(aes(size = point_size, alpha = point_alpha)) +
  ggrepel::geom_text_repel(
    data = dc$scatter_df[dc$scatter_df$Histo_group == "DC", , drop = FALSE],
    aes(label = sample_id),
    size = 3.4,
    box.padding = 0.35,
    point.padding = 0.3,
    segment.color = "grey60",
    max.overlaps = 20,
    show.legend = FALSE
  ) +
  scale_color_manual(values = group_colors) +
  scale_size_identity() +
  scale_alpha_identity() +
  labs(
    title = "PC Signal vs Effector-memory T Cells",
    subtitle = "A second immune axis of the DC intermediate state",
    x = "Tumor_PC_enriched proportion",
    y = "Effector-memory T proportion",
    color = "Pathology"
  ) +
  theme_pub +
  theme(legend.position = "right")

final_plot <- (p_heat_group + p_heat_dc) / (p_scatter_mac + p_scatter_t) +
  plot_layout(widths = c(1.05, 1.35), heights = c(1.45, 1), guides = "collect") &
  theme(legend.position = "right")

out_pdf_use <- resolve_pdf_output(OUT_PDF)
ggsave(out_pdf_use, plot = final_plot, width = 16, height = 12, dpi = 320, useDingbats = FALSE)

message("DC vs CC vs PC profile figure written to: ", out_pdf_use)
#message("DC plot-ready inputs written to: ", OUT_RDS)
