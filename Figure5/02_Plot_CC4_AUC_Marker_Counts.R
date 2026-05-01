options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(ggplot2)
})

OUT_DIR <- "F:/Chordoma/Result/AUC_Reference_Matrix/CC4_Extended"
PLOT_READY_RDS <- file.path(OUT_DIR, "plot_ready_rds", "cc4_auc_reference_plot_ready.rds")
COUNT_CSV <- file.path(OUT_DIR, "AUC_Marker_Counts_17Subtypes.csv")
#OUT_BAR_PDF <- file.path(OUT_DIR, "AUC_Marker_Count_Barplot.pdf")
OUT_ROSE_PDF <- file.path(OUT_DIR, "AUC_Marker_Count_Rose.pdf")

plot_ready <- readRDS(PLOT_READY_RDS)
marker_counts <- plot_ready$marker_counts
lineage_palette <- plot_ready$lineage_palette

write.csv(marker_counts, COUNT_CSV, row.names = FALSE, quote = TRUE)

plot_df <- marker_counts[marker_counts$n_markers > 0, , drop = FALSE]
stopifnot(nrow(plot_df) > 0)

plot_df$pretty_label <- factor(plot_df$pretty_label, levels = plot_df$pretty_label)
plot_df$subtype <- factor(plot_df$subtype, levels = plot_df$subtype)

bar_df <- plot_df[order(plot_df$n_markers, decreasing = TRUE), , drop = FALSE]
bar_df$pretty_label <- factor(bar_df$pretty_label, levels = rev(as.character(bar_df$pretty_label)))

#p_bar <- ggplot(bar_df, aes(x = pretty_label, y = n_markers, fill = lineage)) +
#  geom_col(width = 0.72, color = "white", linewidth = 0.3) +
#  geom_text(aes(label = n_markers), hjust = -0.15, size = 4) +
#  coord_flip() +
#  scale_fill_manual(values = lineage_palette) +
#  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
#  labs(
#    title = "Marker Gene Counts Across Reference Cell States",
#    x = NULL,
#    y = "Number of marker genes",
#    fill = "Lineage"
#  ) +
#  theme_bw(base_size = 14) +
#  theme(
#    plot.title = element_text(hjust = 0.5, face = "bold"),
#    legend.position = "right",
#    panel.grid.minor = element_blank(),
#    panel.grid.major.y = element_blank()
#  )

#ggsave(OUT_BAR_PDF, plot = p_bar, width = 9.5, height = 7.2, useDingbats = FALSE)

rose_df <- plot_df
rose_df$id <- seq_len(nrow(rose_df))
rose_df$label_angle <- 90 - 360 * (rose_df$id - 0.5) / nrow(rose_df)
rose_df$hjust <- ifelse(rose_df$label_angle < -90, 1, 0)
rose_df$label_angle_adj <- ifelse(rose_df$label_angle < -90, rose_df$label_angle + 180, rose_df$label_angle)

p_rose <- ggplot(rose_df, aes(x = factor(id), y = n_markers, fill = lineage)) +
  geom_col(width = 0.95, color = "white", linewidth = 0.25) +
  geom_text(aes(y = n_markers + 0.7, label = n_markers), size = 3.4, color = "black") +
  geom_text(
    aes(
      y = max(n_markers) + 2.2,
      label = pretty_label,
      angle = label_angle_adj,
      hjust = hjust
    ),
    size = 3.3
  ) +
  scale_fill_manual(values = lineage_palette) +
  scale_y_continuous(limits = c(0, max(rose_df$n_markers) + 4)) +
  coord_polar(start = 0) +
  labs(
    title = "Marker Gene Counts Across Reference Cell States",
    fill = "Lineage"
  ) +
  theme_void(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right"
  )

ggsave(OUT_ROSE_PDF, plot = p_rose, width = 10, height = 10, useDingbats = FALSE)

#message("Marker count bar plot written to: ", OUT_BAR_PDF)
message("Marker count rose plot written to: ", OUT_ROSE_PDF)
