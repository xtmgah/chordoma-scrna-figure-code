options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(fmsb)
  library(ggplot2)
})

OUT_DIR <- "F:/Chordoma/Result/AUC_Reference_Matrix/RNAseq/deconv_17subcells"
IN_RDS <- file.path(OUT_DIR, "binary_modeling_results.rds")
OUT_RADAR_PDF <- file.path(OUT_DIR, "combined_binary_model_radarplot_R.pdf")
OUT_TILE_PDF <- file.path(OUT_DIR, "combined_binary_model_metrics_tileplot_R.pdf")
OUT_PLOT_RDS <- file.path(OUT_DIR, "model_comparison_plot_inputs.rds")

MODEL_ORDER <- c("LogisticRegression", "RandomForest", "SVM_Radial", "DecisionTree", "NaiveBayes")
MODEL_LABELS <- c(
  "LogisticRegression" = "Logistic Regression",
  "RandomForest" = "Random Forest",
  "SVM_Radial" = "SVM (Radial)",
  "DecisionTree" = "Decision Tree",
  "NaiveBayes" = "Naive Bayes"
)
MODEL_COLORS <- c(
  "LogisticRegression" = "#1b9e77",
  "RandomForest" = "#d95f02",
  "SVM_Radial" = "#7570b3",
  "DecisionTree" = "#e7298a",
  "NaiveBayes" = "#66a61e"
)
METRIC_ORDER <- c("auc", "accuracy", "recall", "f1", "composite_score")
METRIC_LABELS <- c(
  "auc" = "AUC",
  "accuracy" = "Accuracy",
  "recall" = "Recall",
  "f1" = "F1 score",
  "composite_score" = "Composite"
)

resolve_pdf_output <- function(path) {
  if (!file.exists(path)) return(path)
  con <- try(file(path, open = "ab"), silent = TRUE)
  if (inherits(con, "try-error")) {
    return(sub("\\.pdf$", "_latest.pdf", path))
  }
  close(con)
  path
}

res <- readRDS(IN_RDS)
perf <- res$performance
perf <- perf[perf$model %in% MODEL_ORDER, , drop = FALSE]
perf <- perf[match(MODEL_ORDER, perf$model), , drop = FALSE]
perf <- perf[!is.na(perf$model), , drop = FALSE]

metric_mat <- as.data.frame(perf[, METRIC_ORDER, drop = FALSE])
metric_mat[] <- lapply(metric_mat, as.numeric)
rownames(metric_mat) <- perf$model

radar_df <- rbind(
  rep(1, ncol(metric_mat)),
  rep(0, ncol(metric_mat)),
  metric_mat
)
rownames(radar_df)[1:2] <- c("Max", "Min")
rownames(radar_df)[3:nrow(radar_df)] <- rownames(metric_mat)

metric_long <- do.call(
  rbind,
  lapply(rownames(metric_mat), function(model_name) {
    data.frame(
      model = model_name,
      model_label = unname(MODEL_LABELS[model_name]),
      metric = colnames(metric_mat),
      metric_label = unname(METRIC_LABELS[colnames(metric_mat)]),
      value = as.numeric(metric_mat[model_name, ]),
      stringsAsFactors = FALSE
    )
  })
)
metric_long$model <- factor(metric_long$model, levels = MODEL_ORDER)
metric_long$model_label <- factor(metric_long$model_label, levels = unname(MODEL_LABELS[MODEL_ORDER]))
metric_long$metric <- factor(metric_long$metric, levels = METRIC_ORDER)
metric_long$metric_label <- factor(metric_long$metric_label, levels = unname(METRIC_LABELS[METRIC_ORDER]))

saveRDS(
  list(
    performance = perf,
    metric_matrix = metric_mat,
    radar_df = radar_df,
    metric_long = metric_long,
    model_order = MODEL_ORDER,
    model_labels = MODEL_LABELS,
    model_colors = MODEL_COLORS,
    metric_order = METRIC_ORDER,
    metric_labels = METRIC_LABELS,
    best_model_name = res$best_model_name,
    best_threshold = res$best_threshold,
    n_total = res$n_total,
    n_labeled = res$n_labeled,
    n_pc = res$n_pc,
    n_cc = res$n_cc
  ),
  OUT_PLOT_RDS
)

radar_pdf <- resolve_pdf_output(OUT_RADAR_PDF)
pdf(radar_pdf, width = 8.2, height = 8.2)
radarchart(
  radar_df,
  axistype = 1,
  pcol = MODEL_COLORS[rownames(metric_mat)],
  plwd = 2.2,
  plty = 1,
  cglcol = "grey80",
  cglty = 1,
  cglwd = 0.8,
  axislabcol = "grey40",
  vlcex = 1.05,
  title = paste0(
    "Model Performance for CC vs PC Classification\n",
    "(17 subtype proportions; ", res$n_total, " total samples; ",
    res$n_labeled, " labeled CC/PC samples)"
  )
)
legend(
  "topright",
  legend = unname(MODEL_LABELS[rownames(metric_mat)]),
  col = MODEL_COLORS[rownames(metric_mat)],
  lty = 1,
  lwd = 2,
  bty = "n",
  cex = 0.9
)
dev.off()

tile_pdf <- resolve_pdf_output(OUT_TILE_PDF)
tile_plot <- ggplot(metric_long, aes(x = metric_label, y = model_label, fill = value)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.3f", value)), size = 4.2) +
  scale_fill_gradientn(
    colours = c("#f7fbff", "#6baed6", "#08519c"),
    limits = c(0, 1),
    name = "Score"
  ) +
  labs(
    title = "Performance Comparison Across Binary Classification Models",
    subtitle = paste0(
      "Prediction target: pathological CC vs PC labels using 17 deconvolved subtype proportions",
      "\nBest model: ", unname(MODEL_LABELS[res$best_model_name]), " (threshold = ",
      formatC(res$best_threshold, format = "f", digits = 3), ")"
    ),
    x = NULL,
    y = NULL
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, color = "grey30"),
    axis.text.x = element_text(angle = 25, hjust = 1),
    panel.grid = element_blank(),
    legend.position = "right"
  )

ggsave(tile_pdf, plot = tile_plot, width = 8.8, height = 5.8, useDingbats = FALSE)

message("Model comparison radar plot written to: ", radar_pdf)
message("Model comparison tile plot written to: ", tile_pdf)
message("Model comparison plot inputs written to: ", OUT_PLOT_RDS)
