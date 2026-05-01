options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(survminer)
  library(ggplot2)
})

OUT_DIR <- "F:/Chordoma/Result/AUC_Reference_Matrix/RNAseq/deconv_17subcells"
IN_RDS <- file.path(OUT_DIR, "pc_signal_prognostic_analysis.rds")

OUT_AGE_PDF <- file.path(OUT_DIR, "survival_curve_pc_signal_3groups_age_adjusted.pdf")
OUT_AGE_GENDER_PDF <- file.path(OUT_DIR, "survival_curve_pc_signal_3groups_age_gender_adjusted.pdf")
OUT_CC_ONLY_AGE_PDF <- file.path(OUT_DIR, "survival_curve_cc_only_2groups_age_adjusted.pdf")
OUT_CC_ONLY_AGE_GENDER_PDF <- file.path(OUT_DIR, "survival_curve_cc_only_2groups_age_gender_adjusted.pdf")

GROUP_COLORS <- c(
  "PC-enriched" = "#72567a",
  "CC-PC_enriched" = "#E64B35",
  "CC-enriched" = "#4DBBD5"
)
CC_ONLY_GROUP_COLORS <- c(
  "CC-PC_enriched" = "#E64B35",
  "CC-enriched" = "#4DBBD5"
)

format_p_label <- function(p, prefix) {
  if (is.na(p)) {
    paste0(prefix, " P = NA")
  } else if (p < 0.001) {
    paste0(prefix, " P < 0.001")
  } else {
    paste0(prefix, " P = ", formatC(p, format = "f", digits = 3))
  }
}

plot_adjusted_survival <- function(fit, data, palette, title, subtitle, out_pdf, legend_title) {
  surv_df <- surv_summary(fit, data = data)
  surv_df$strata <- factor(
    names(palette)[as.integer(as.character(surv_df$strata))],
    levels = names(palette)
  )

  p <- ggsurvplot_df(
    surv_df,
    conf.int = FALSE,
    censor = FALSE,
    palette = unname(palette),
    legend.title = legend_title,
    legend.labs = names(palette),
    size = 1.2,
    break.time.by = 24,
    xlab = "Time (months)",
    ylab = "Adjusted overall survival probability",
    ggtheme = theme_bw(base_size = 14)
  )

  p <- p +
    labs(title = title, subtitle = subtitle) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, color = "grey30"),
      legend.position = "right",
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.3)
    )

  ggsave(out_pdf, plot = p, width = 8.6, height = 6.8, useDingbats = FALSE)
}

res <- readRDS(IN_RDS)

subtitle_age <- paste0(
  "Best CC cutoff = ",
  formatC(res$best_cutoff, format = "f", digits = 4),
  " | ",
  format_p_label(res$cox_age_p, "Age-adjusted Cox group effect")
)
subtitle_age_gender <- paste0(
  "Best CC cutoff = ",
  formatC(res$best_cutoff, format = "f", digits = 4),
  " | ", format_p_label(res$cox_age_gender_p, "Age+gender-adjusted Cox group effect")
)
subtitle_cc_only_age <- paste0(
  "Best CC cutoff = ",
  formatC(res$best_cutoff, format = "f", digits = 4),
  " | ", format_p_label(res$cox_cc_only_age_p, "Age-adjusted Cox group effect")
)
subtitle_cc_only_age_gender <- paste0(
  "Best CC cutoff = ",
  formatC(res$best_cutoff, format = "f", digits = 4),
  " | ", format_p_label(res$cox_cc_only_age_gender_p, "Age+gender-adjusted Cox group effect")
)

plot_adjusted_survival(
  fit = res$fit_age,
  data = res$newdata_age,
  palette = GROUP_COLORS,
  title = "Adjusted Survival Stratification by Tumor_PC_enriched Signal",
  subtitle = subtitle_age,
  out_pdf = OUT_AGE_PDF,
  legend_title = "Risk group"
)
plot_adjusted_survival(
  fit = res$fit_age_gender,
  data = res$newdata_age_gender,
  palette = GROUP_COLORS,
  title = "Adjusted Survival Stratification by Tumor_PC_enriched Signal",
  subtitle = subtitle_age_gender,
  out_pdf = OUT_AGE_GENDER_PDF,
  legend_title = "Risk group"
)
plot_adjusted_survival(
  fit = res$fit_cc_only_age,
  data = res$newdata_cc_only_age,
  palette = CC_ONLY_GROUP_COLORS,
  title = "Adjusted Survival Stratification Within CC Samples",
  subtitle = subtitle_cc_only_age,
  out_pdf = OUT_CC_ONLY_AGE_PDF,
  legend_title = "CC group"
)
plot_adjusted_survival(
  fit = res$fit_cc_only_age_gender,
  data = res$newdata_cc_only_age_gender,
  palette = CC_ONLY_GROUP_COLORS,
  title = "Adjusted Survival Stratification Within CC Samples",
  subtitle = subtitle_cc_only_age_gender,
  out_pdf = OUT_CC_ONLY_AGE_GENDER_PDF,
  legend_title = "CC group"
)

message("Final survival PDFs written to: ", OUT_DIR)
