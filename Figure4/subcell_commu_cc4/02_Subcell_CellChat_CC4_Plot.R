options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

if (!requireNamespace("CellChat", quietly = TRUE)) {
  stop(
    "The `CellChat` package is not installed. Please install CellChat first.",
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(CellChat)
})

OUT_DIR <- "F:/Chordoma/Result/cellchat/Analysis_CellChat/subcell_commu_cc4"
PLOT_READY_DIR <- file.path(OUT_DIR, "plot_ready_rds")

INPUT_ASSET_RDS <- file.path(PLOT_READY_DIR, "subcell_cc4_cellchat_input_assets.rds")
OVERALL_RDS <- file.path(OUT_DIR, "subcell_cc4_cellchat_overall.rds")
CC_RDS <- file.path(OUT_DIR, "subcell_cc4_cellchat_CC.rds")
PC_RDS <- file.path(OUT_DIR, "subcell_cc4_cellchat_PC.rds")
MERGED_RDS <- file.path(OUT_DIR, "subcell_cc4_cellchat_CC_vs_PC_merged.rds")

TOP_PATHWAYS_TO_PLOT <- 10
FIG_WIDTH <- 14
FIG_HEIGHT <- 10
HEATMAP_WIDTH <- 22
HEATMAP_HEIGHT <- 13
BASE_SIZE <- 16
SCATTER_LABEL_SIZE <- 4.2
LEGEND_POSITION <- c(0.83, 0.18)
COMPARISON_COLORS <- c("CC" = "#61bada", "PC" = "#72567a")

plot_single_cellchat <- function(cellchat_obj, prefix, color_map) {
  plot_obj <- cellchat_obj
  group_size <- as.numeric(table(plot_obj@idents))
  names(group_size) <- levels(plot_obj@idents)

  local_colors <- color_map[levels(plot_obj@idents)]
  names(local_colors) <- levels(plot_obj@idents)

  pdf(file.path(OUT_DIR, paste0(prefix, "_network_circle.pdf")), width = FIG_WIDTH, height = FIG_HEIGHT)
  par(mfrow = c(1, 2), xpd = TRUE, mar = c(1.5, 1.5, 4, 4))
  netVisual_circle(
    plot_obj@net$count,
    vertex.weight = group_size,
    weight.scale = TRUE,
    label.edge = FALSE,
    title.name = "Number of interactions",
    color.use = local_colors,
    vertex.label.cex = 0.85
  )
  netVisual_circle(
    plot_obj@net$weight,
    vertex.weight = group_size,
    weight.scale = TRUE,
    label.edge = FALSE,
    title.name = "Interaction weights/strength",
    color.use = local_colors,
    vertex.label.cex = 0.85
  )
  dev.off()

  pdf(file.path(OUT_DIR, paste0(prefix, "_signaling_role_heatmap.pdf")), width = HEATMAP_WIDTH, height = HEATMAP_HEIGHT)
  p1 <- netAnalysis_signalingRole_heatmap(
    plot_obj,
    pattern = "outgoing",
    color.use = local_colors,
    font.size = 6,
    width = 11,
    height = 12
  )
  p2 <- netAnalysis_signalingRole_heatmap(
    plot_obj,
    pattern = "incoming",
    color.use = local_colors,
    font.size = 6,
    width = 11,
    height = 12
  )
  print(p1 + p2)
  dev.off()

  scatter_title <- switch(
    prefix,
    "subcell_cc4_overall" = "Global Intercellular Signaling (CellChat)",
    "subcell_cc4_CC" = "CC Intercellular Signaling (CellChat)",
    "subcell_cc4_PC" = "PC Intercellular Signaling (CellChat)",
    paste0(prefix, " Intercellular Signaling (CellChat)")
  )

  p_scatter <- netAnalysis_signalingRole_scatter(
    plot_obj,
    color.use = local_colors,
    do.label = TRUE,
    label.size = SCATTER_LABEL_SIZE
  ) +
    labs(title = scatter_title) +
    theme_bw(base_size = BASE_SIZE) +
    theme(
      legend.position = LEGEND_POSITION,
      plot.title = element_text(hjust = 0.5),
      plot.margin = margin(12, 20, 12, 12)
    )

  ggsave(
    filename = file.path(OUT_DIR, paste0(prefix, "_signaling_role_scatter.pdf")),
    plot = p_scatter,
    width = FIG_WIDTH,
    height = FIG_HEIGHT
  )

  pathways <- head(plot_obj@netP$pathways, TOP_PATHWAYS_TO_PLOT)
  pathway_dir <- file.path(OUT_DIR, paste0(prefix, "_pathways"))
  dir.create(pathway_dir, recursive = TRUE, showWarnings = FALSE)

  if (length(pathways) > 0) {
    try({
      pdf(file.path(OUT_DIR, paste0(prefix, "_bubble_top_pathways.pdf")), width = 17, height = 11)
      print(
        netVisual_bubble(
          plot_obj,
          signaling = pathways,
          remove.isolate = FALSE,
          angle.x = 45
        )
      )
      dev.off()
    }, silent = TRUE)
    while (dev.cur() > 1) dev.off()
  }

  for (pathway in pathways) {
    try({
      pdf(file.path(pathway_dir, paste0(pathway, "_circle.pdf")), width = FIG_WIDTH, height = FIG_HEIGHT)
      netVisual_aggregate(
        plot_obj,
        signaling = pathway,
        layout = "circle",
        color.use = local_colors,
        vertex.label.cex = 0.85
      )
      dev.off()

      gg <- netAnalysis_contribution(plot_obj, signaling = pathway)
      ggsave(
        filename = file.path(pathway_dir, paste0(pathway, "_contribution.pdf")),
        plot = gg,
        width = FIG_WIDTH,
        height = FIG_HEIGHT
      )
    }, silent = TRUE)
    while (dev.cur() > 1) dev.off()
  }
}

message("Reading CC4 CellChat plotting assets ...")
assets <- readRDS(INPUT_ASSET_RDS)
color_map <- assets$color_map

message("Plotting overall CC4 CellChat figures ...")
cellchat_overall <- readRDS(OVERALL_RDS)
plot_single_cellchat(cellchat_overall, "subcell_cc4_overall", color_map)

if (file.exists(CC_RDS)) {
  message("Plotting CC-only CC4 CellChat figures ...")
  cellchat_cc <- readRDS(CC_RDS)
  plot_single_cellchat(cellchat_cc, "subcell_cc4_CC", color_map)
}

if (file.exists(PC_RDS)) {
  message("Plotting PC-only CC4 CellChat figures ...")
  cellchat_pc <- readRDS(PC_RDS)
  plot_single_cellchat(cellchat_pc, "subcell_cc4_PC", color_map)
}

if (file.exists(MERGED_RDS)) {
  message("Plotting CC vs PC comparison figures ...")
  merged_cellchat <- readRDS(MERGED_RDS)

  pdf(file.path(OUT_DIR, "subcell_cc4_CC_vs_PC_compareInteractions.pdf"), width = FIG_WIDTH, height = FIG_HEIGHT)
  p1 <- compareInteractions(
    merged_cellchat,
    show.legend = FALSE,
    group = c(1, 2),
    color.use = COMPARISON_COLORS
  )
  p2 <- compareInteractions(
    merged_cellchat,
    show.legend = FALSE,
    group = c(1, 2),
    measure = "weight",
    color.use = COMPARISON_COLORS
  )
  print(p1 + p2)
  dev.off()

  rank_data <- rankNet(
    merged_cellchat,
    mode = "comparison",
    stacked = TRUE,
    do.stat = FALSE,
    return.data = TRUE
  )
  rank_df <- rank_data$signaling.contribution
  keep_stat <- split(rank_df$contribution, rank_df$name)
  keep_stat <- vapply(keep_stat, function(x) max(x, na.rm = TRUE), numeric(1))
  keep_signaling <- names(keep_stat)[keep_stat > 0]

  pdf(file.path(OUT_DIR, "subcell_cc4_CC_vs_PC_rankNet.pdf"), width = FIG_WIDTH, height = FIG_HEIGHT)
  print(
    rankNet(
      merged_cellchat,
      mode = "comparison",
      stacked = TRUE,
      do.stat = FALSE,
      signaling = keep_signaling,
      color.use = COMPARISON_COLORS,
      font.size = 9
    )
  )
  dev.off()
}

message("CC4 CellChat plotting finished.")
