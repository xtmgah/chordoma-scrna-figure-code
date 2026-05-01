script_dir <- if (exists("SCRIPT_DIR", inherits = FALSE)) {
  SCRIPT_DIR
} else {
  dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), mustWork = TRUE))
}
source(file.path(dirname(script_dir), "_figure_style.R"))

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(survminer)
  library(fmsb)
})

OUT_DIR <- file.path(script_dir, "nature_panels")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

message("Figure 5: marker count reference matrix")
plot_ready <- readRDS(file.path(script_dir, "plot_ready_rds", "cc4_auc_reference_plot_ready.rds"))
marker_counts <- plot_ready$marker_counts
lineage_levels <- c(intersect(c("Tumor", "Macrophage", "T cell", "T cells"), unique(marker_counts$lineage)), setdiff(unique(marker_counts$lineage), c("Tumor", "Macrophage", "T cell", "T cells")))
lineage_palette <- palette_for(lineage_levels, lineage_palette_shared)
write.csv(marker_counts, file.path(OUT_DIR, "fig5_marker_counts_17_subtypes.csv"), row.names = FALSE, quote = TRUE)

plot_df <- marker_counts[marker_counts$n_markers > 0, , drop = FALSE]
plot_df$pretty_label <- factor(plot_df$pretty_label, levels = plot_df$pretty_label)
plot_df$lineage <- factor(plot_df$lineage, levels = lineage_levels)

rose_df <- plot_df %>%
  arrange(lineage, desc(n_markers)) %>%
  mutate(pretty_label = factor(pretty_label, levels = pretty_label))

p_rose <- ggplot(rose_df, aes(x = pretty_label, y = n_markers, fill = lineage)) +
  geom_segment(aes(xend = pretty_label, y = 0, yend = n_markers, color = lineage), linewidth = 0.35, alpha = 0.8) +
  geom_point(shape = 21, size = 3.0, color = "black", stroke = 0.08) +
  geom_text(aes(y = n_markers + ifelse(n_markers < 10, 0.85, 0.45), label = n_markers), hjust = 0, size = 2.6, family = FONT_FAMILY) +
  scale_fill_manual(values = lineage_palette, drop = FALSE) +
  scale_color_manual(values = lineage_palette, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.16))) +
  coord_flip() +
  labs(x = NULL, y = "Marker genes", fill = "Lineage") +
  theme_nature(base_size = 8.3) +
  theme(
    legend.position = "right",
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(size = 7.1)
  )
save_panel_pdf(p_rose, file.path(OUT_DIR, "fig5a_marker_count_rose.pdf"), 4.6, 3.4)

p_bar <- plot_df %>%
  arrange(lineage, desc(n_markers)) %>%
  mutate(pretty_label = factor(pretty_label, levels = rev(pretty_label))) %>%
  ggplot(aes(x = n_markers, y = pretty_label, fill = lineage)) +
  geom_col(width = 0.72, color = "white", linewidth = 0.18) +
  geom_text(aes(label = n_markers), hjust = -0.15, size = 2.8, family = FONT_FAMILY) +
  scale_fill_manual(values = lineage_palette, drop = FALSE) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.09))) +
  labs(x = "Marker genes", y = NULL, fill = "Lineage") +
  theme_nature() +
  theme(panel.grid.major.y = element_blank())
save_panel_pdf(p_bar, file.path(OUT_DIR, "fig5a_marker_count_bar.pdf"), 4.3, 3.8)

message("Figure 5: model performance")
res <- readRDS(file.path(script_dir, "binary_modeling_results.rds"))
MODEL_ORDER <- c("LogisticRegression", "RandomForest", "SVM_Radial", "DecisionTree", "NaiveBayes")
MODEL_LABELS <- c(
  "LogisticRegression" = "Logistic regression",
  "RandomForest" = "Random forest",
  "SVM_Radial" = "SVM radial",
  "DecisionTree" = "Decision tree",
  "NaiveBayes" = "Naive Bayes"
)
MODEL_COLORS <- c(
  "LogisticRegression" = "#1B9E77",
  "RandomForest" = "#D95F02",
  "SVM_Radial" = "#7570B3",
  "DecisionTree" = "#E7298A",
  "NaiveBayes" = "#66A61E"
)
METRIC_ORDER <- c("auc", "accuracy", "recall", "f1", "composite_score")
METRIC_LABELS <- c("auc" = "AUC", "accuracy" = "Accuracy", "recall" = "Recall", "f1" = "F1 score", "composite_score" = "Composite")

perf <- res$performance
perf <- perf[perf$model %in% MODEL_ORDER, , drop = FALSE]
perf <- perf[match(MODEL_ORDER, perf$model), , drop = FALSE]
perf <- perf[!is.na(perf$model), , drop = FALSE]
metric_long <- perf %>%
  select(model, all_of(METRIC_ORDER)) %>%
  pivot_longer(-model, names_to = "metric", values_to = "value") %>%
  mutate(
    model = factor(model, levels = MODEL_ORDER),
    model_label = factor(unname(MODEL_LABELS[as.character(model)]), levels = unname(MODEL_LABELS[MODEL_ORDER])),
    metric = factor(metric, levels = METRIC_ORDER),
    metric_label = factor(unname(METRIC_LABELS[as.character(metric)]), levels = unname(METRIC_LABELS[METRIC_ORDER]))
  )

p_tile <- ggplot(metric_long, aes(x = metric_label, y = model_label, fill = value)) +
  geom_tile(color = "white", linewidth = 0.35) +
  geom_text(aes(label = sprintf("%.3f", value)), size = 3.0, family = FONT_FAMILY) +
  scale_fill_gradientn(colours = c("#F7FBFF", "#9ECAE1", "#08519C"), limits = c(0, 1), name = "Score") +
  labs(x = NULL, y = NULL) +
  theme_nature() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1), panel.grid = element_blank())
save_panel_pdf(p_tile, file.path(OUT_DIR, "fig5b_model_metric_tileplot.pdf"), 4.7, 3.1)

metric_mat <- perf[, METRIC_ORDER, drop = FALSE]
metric_mat[] <- lapply(metric_mat, as.numeric)
rownames(metric_mat) <- perf$model
radar_df <- rbind(rep(1, ncol(metric_mat)), rep(0, ncol(metric_mat)), metric_mat)
rownames(radar_df) <- c("Max", "Min", rownames(metric_mat))
colnames(radar_df) <- unname(METRIC_LABELS[colnames(radar_df)])
radar_pdf <- file.path(OUT_DIR, "fig5b_model_metric_radar.pdf")
open_panel_pdf(radar_pdf, 4.4, 4.4)
op <- par(family = FONT_FAMILY, mar = c(1.8, 1.8, 2.2, 1.8))
radarchart(
  radar_df,
  axistype = 1,
  pcol = MODEL_COLORS[rownames(metric_mat)],
  plwd = 1.7,
  plty = 1,
  cglcol = "grey82",
  cglty = 1,
  cglwd = 0.55,
  axislabcol = "grey35",
  vlcex = 0.85,
  caxislabels = c("0", "0.25", "0.50", "0.75", "1")
)
legend(
  "bottom",
  legend = unname(MODEL_LABELS[rownames(metric_mat)]),
  col = MODEL_COLORS[rownames(metric_mat)],
  lty = 1,
  lwd = 1.7,
  bty = "n",
  cex = 0.72,
  ncol = 2,
  inset = -0.08,
  xpd = TRUE
)
par(op)
close_panel_pdf()

message("Figure 5: adjusted survival")
surv_res <- readRDS(file.path(script_dir, "pc_signal_prognostic_analysis.rds"))
format_p_label <- function(p, prefix) {
  if (is.na(p)) return(paste0(prefix, " P = NA"))
  if (p < 0.001) return(paste0(prefix, " P < 0.001"))
  paste0(prefix, " P = ", formatC(p, format = "f", digits = 3))
}
format_hr_label <- function(cox_table) {
  rows <- cox_table[grepl("group", cox_table$variable), , drop = FALSE]
  if (nrow(rows) == 0) return(NULL)
  paste0(
    gsub(".*group", "", rows$variable),
    ": HR ",
    sprintf("%.2f", rows$HR),
    " (95% CI ",
    sprintf("%.2f", rows$lower95),
    "-",
    sprintf("%.2f", rows$upper95),
    ")",
    collapse = "\n"
  )
}
plot_adjusted_survival <- function(fit, data, palette, title, subtitle, out_pdf, legend_title, hr_label = NULL, hr_y = 0.10) {
  surv_df <- surv_summary(fit, data = data)
  surv_df$strata <- factor(names(palette)[as.integer(as.character(surv_df$strata))], levels = names(palette))
  p <- ggsurvplot_df(
    surv_df,
    conf.int = FALSE,
    censor = FALSE,
    palette = unname(palette),
    legend.title = legend_title,
    legend.labs = names(palette),
    size = 0.85,
    break.time.by = 24,
    xlab = "Time (months)",
    ylab = "Adjusted overall survival probability",
    ggtheme = theme_nature(base_size = 8.2)
  ) +
    labs(title = title, subtitle = subtitle) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 8)) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 7.2),
      legend.position = c(0.76, 0.76),
      legend.background = element_blank(),
      panel.grid = element_blank()
    )
  if (!is.null(hr_label)) {
    p <- p + annotate("text", x = Inf, y = hr_y, label = hr_label, hjust = 1.02, vjust = 0, size = 2.15, family = FONT_FAMILY)
  }
  save_panel_pdf(p, out_pdf, 3.9, 3.1)
}
survival_palette <- group_palette[c("PC-enriched", "CC-PC_enriched", "CC-enriched")]
cc_only_palette <- group_palette[c("CC-PC_enriched", "CC-enriched")]
plot_adjusted_survival(
  surv_res$fit_age,
  surv_res$newdata_age,
  survival_palette,
  "Tumor PC-enriched signal",
  paste0("Best CC cutoff = ", formatC(surv_res$best_cutoff, format = "f", digits = 4), " | ", format_p_label(surv_res$cox_age_p, "Age-adjusted Cox")),
  file.path(OUT_DIR, "fig5d_survival_pc_signal_age_adjusted.pdf"),
  "Group",
  format_hr_label(surv_res$cox_age_table),
  hr_y = 0.045
)
plot_adjusted_survival(
  surv_res$fit_cc_only_age,
  surv_res$newdata_cc_only_age,
  cc_only_palette,
  "CC-only stratification",
  paste0("Best CC cutoff = ", formatC(surv_res$best_cutoff, format = "f", digits = 4), " | ", format_p_label(surv_res$cox_cc_only_age_p, "Age-adjusted Cox")),
  file.path(OUT_DIR, "fig5c_survival_cc_only_age_adjusted.pdf"),
  "Group",
  format_hr_label(surv_res$cox_cc_only_age_table)
)

message("Figure 5: DC, CC and PC profile panels")
dc <- readRDS(file.path(script_dir, "dc_profile_plot_inputs.rds"))
group_colors <- group_palette[c("CC", "PC")]
prop_fill_scale <- function() ggsci::scale_fill_material(
  "deep-orange",
  limits = c(0, dc$prop_limit),
  breaks = scales::pretty_breaks(n = 5),
  labels = label_percent(accuracy = 1),
  name = "Proportion"
)
p_heat_group <- ggplot(dc$group_mean_long %>% filter(Histo %in% c("CC", "PC")), aes(x = Histo, y = cell_type, fill = mean_prop)) +
  geom_tile(color = "white", linewidth = 0.35) +
  geom_text(aes(label = label), size = 2.3, family = FONT_FAMILY) +
  prop_fill_scale() +
  labs(x = NULL, y = NULL) +
  theme_nature() +
  guides(fill = guide_colorbar(barheight = unit(32, "mm"), barwidth = unit(3.0, "mm"))) +
  theme(axis.text.x = element_text(color = unname(group_colors[c("CC", "PC")])), axis.text.y = element_text(size = 7.2), panel.grid = element_blank())
save_panel_pdf(p_heat_group, file.path(OUT_DIR, "fig5e_group_level_deconvolution_heatmap.pdf"), 2.8, 4.4)

p_heat_dc <- ggplot(dc$dc_long, aes(x = sample_label, y = cell_type, fill = prop)) +
  geom_tile(color = "white", linewidth = 0.35) +
  geom_text(aes(label = label), size = 2.25, family = FONT_FAMILY) +
  prop_fill_scale() +
  labs(x = NULL, y = NULL) +
  theme_nature() +
  guides(fill = guide_colorbar(barheight = unit(32, "mm"), barwidth = unit(3.0, "mm"))) +
  theme(axis.text.x = element_text(size = 7, angle = 30, hjust = 1), axis.text.y = element_text(size = 7.2), panel.grid = element_blank())
save_panel_pdf(p_heat_dc, file.path(OUT_DIR, "fig5f_single_sample_dc_profile_heatmap.pdf"), 4.6, 4.4)

scatter_theme <- theme_nature() +
  theme(
    legend.position = c(0.82, 0.78),
    legend.background = element_blank(),
    panel.grid = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
scatter_df <- dc$scatter_df %>% filter(Histo_group %in% c("CC", "PC"))
p_scatter_mac <- ggplot(scatter_df, aes(x = `Tumor_PC_enriched`, y = `CCL3+ Mac`, fill = Histo_group)) +
  geom_point(aes(size = point_size), shape = 21, color = "black", stroke = 0.1, alpha = 1) +
  scale_fill_manual(values = group_colors) +
  scale_size_identity() +
  scale_x_continuous(labels = label_percent(accuracy = 1), breaks = scales::pretty_breaks(n = 7)) +
  scale_y_continuous(labels = label_percent(accuracy = 1), breaks = scales::pretty_breaks(n = 7)) +
  guides(fill = guide_legend(override.aes = list(size = 2.6, shape = 21, color = "black", stroke = 0.1, alpha = 1))) +
  labs(x = "Tumor PC-enriched proportion", y = "CCL3+ macrophage proportion", fill = "Group") +
  scatter_theme
save_panel_pdf(p_scatter_mac, file.path(OUT_DIR, "fig5g_pc_signal_vs_ccl3_macrophages.pdf"), 3.8, 3.2)

p_scatter_t <- ggplot(scatter_df, aes(x = `Tumor_PC_enriched`, y = `Effector-memory T`, fill = Histo_group)) +
  geom_point(aes(size = point_size), shape = 21, color = "black", stroke = 0.1, alpha = 1) +
  scale_fill_manual(values = group_colors) +
  scale_size_identity() +
  scale_x_continuous(labels = label_percent(accuracy = 1), breaks = scales::pretty_breaks(n = 7)) +
  scale_y_continuous(labels = label_percent(accuracy = 1), breaks = scales::pretty_breaks(n = 7)) +
  guides(fill = guide_legend(override.aes = list(size = 2.6, shape = 21, color = "black", stroke = 0.1, alpha = 1))) +
  labs(x = "Tumor PC-enriched proportion", y = "Effector-memory T-cell proportion", fill = "Group") +
  scatter_theme
save_panel_pdf(p_scatter_t, file.path(OUT_DIR, "fig5h_pc_signal_vs_effector_memory_t.pdf"), 3.8, 3.2)

combined <- (p_heat_group + p_heat_dc) / (p_scatter_mac + p_scatter_t) +
  plot_layout(widths = c(1, 1.25), heights = c(1.35, 1), guides = "collect") &
  theme(legend.position = "right")
save_panel_pdf(combined, file.path(OUT_DIR, "fig5_profile_summary_combined.pdf"), 7.6, 6.5)

message("Figure 5 nature panels written to: ", OUT_DIR)
