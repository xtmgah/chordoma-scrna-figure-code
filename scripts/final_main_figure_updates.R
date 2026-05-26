script_dir <- if (exists("SCRIPT_DIR", inherits = FALSE)) {
  SCRIPT_DIR
} else {
  dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), mustWork = TRUE))
}
repo_dir <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
source(file.path(repo_dir, "_figure_style.R"))

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(scales)
  library(grid)
  library(survival)
  library(survminer)
  library(Seurat)
})

fig2_dir <- file.path(repo_dir, "Figure2")
fig3_dir <- file.path(repo_dir, "Figure3")
fig2_out <- file.path(fig2_dir, "nature_panels")
fig3_out <- file.path(fig3_dir, "nature_panels")
dir.create(fig2_out, recursive = TRUE, showWarnings = FALSE)
dir.create(fig3_out, recursive = TRUE, showWarnings = FALSE)

message("Regenerating requested Figure 2 and Figure 3 individual panels")

cluster_col_from <- function(meta) {
  if ("RNA_snn_res.0.5" %in% colnames(meta)) {
    "RNA_snn_res.0.5"
  } else {
    grep("snn_res", colnames(meta), value = TRUE)[1]
  }
}

get_tumor_cluster_palette <- function(sc) {
  cluster_col <- cluster_col_from(sc@meta.data)
  stopifnot(!is.na(cluster_col))
  cluster_levels <- sort_numeric_labels(unique(as.character(sc@meta.data[[cluster_col]])))
  setNames(nature_discrete(length(cluster_levels)), cluster_levels)
}

save_requested_pdf <- function(plot, path, width, height, dpi = 600) {
  save_panel_pdf(plot, path, width = width, height = height, dpi = dpi)
  invisible(path)
}

make_fig3d_cytotrace <- function() {
  message("  Figure 3d: CytoTRACE2 ranked developmental-potential violin/box plot")
  sc <- readRDS(file.path(fig3_dir, "sc_tumor_with_CytoTRACE2.rds"))
  cluster_col <- cluster_col_from(sc@meta.data)
  if (!"CytoTRACE2_Score" %in% colnames(sc@meta.data)) {
    stop("CytoTRACE2_Score was not found in sc_tumor_with_CytoTRACE2.rds")
  }
  cluster_cols <- get_tumor_cluster_palette(sc)

  plot_df <- data.frame(
    cluster = as.character(sc@meta.data[[cluster_col]]),
    CytoTRACE2_Score = as.numeric(sc@meta.data$CytoTRACE2_Score),
    stringsAsFactors = FALSE
  ) %>%
    filter(is.finite(CytoTRACE2_Score), !is.na(cluster))

  cluster_order <- plot_df %>%
    group_by(cluster) %>%
    summarise(mean_score = mean(CytoTRACE2_Score, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(mean_score)) %>%
    pull(cluster)

  plot_df$cluster <- factor(plot_df$cluster, levels = cluster_order)
  cluster_cols <- cluster_cols[cluster_order]

  p <- ggplot(plot_df, aes(x = cluster, y = CytoTRACE2_Score, fill = cluster)) +
    geom_violin(width = 0.86, linewidth = 0.18, color = "grey25", scale = "width", trim = TRUE) +
    geom_boxplot(width = 0.095, fill = "white", color = "grey20", linewidth = 0.22, outlier.shape = NA) +
    scale_fill_manual(values = cluster_cols, guide = "none") +
    scale_y_continuous(breaks = pretty_breaks(n = 5), expand = expansion(mult = c(0.02, 0.06))) +
    labs(x = "Tumor cluster", y = "CytoTRACE2 score") +
    theme_nature(base_size = 9.2) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8.7),
      axis.title = element_text(face = "plain", size = 10.0),
      plot.margin = margin(4, 5, 4, 5)
    )

  save_requested_pdf(p, file.path(fig3_out, "fig3d_developmental_potential_by_tumor_cluster.pdf"), 5.3, 3.05)
  save_requested_pdf(p, file.path(fig3_dir, "CytoTRACE2_Boxplot_Ranked.pdf"), 5.3, 3.05)
}

make_fig3e_monocle <- function() {
  message("  Figure 3e: Monocle3 trajectory with shared tumor-cluster and pseudotime colors")
  if (!requireNamespace("monocle3", quietly = TRUE)) {
    stop("The monocle3 package is required for fig3e_pseudotime_trajectory.pdf")
  }
  sc <- readRDS(file.path(fig3_dir, "sc_tumor.rds"))
  cds <- readRDS(file.path(fig3_dir, "cds.rds"))
  cluster_col <- cluster_col_from(sc@meta.data)
  cluster_cols <- get_tumor_cluster_palette(sc)
  cluster_levels <- names(cluster_cols)

  cds_meta <- as.data.frame(SummarizedExperiment::colData(cds))
  if (!cluster_col %in% colnames(cds_meta)) {
    stop("Cluster column ", cluster_col, " was not found in cds.rds")
  }
  SummarizedExperiment::colData(cds)[[cluster_col]] <- factor(as.character(cds_meta[[cluster_col]]), levels = cluster_levels)

  base_traj_theme <- theme_umap_nature(base_size = 9.2) +
    theme(
      plot.title = element_text(face = "plain", size = 10.5, hjust = 0.5),
      legend.title = element_text(face = "plain", size = 8.6),
      legend.text = element_text(size = 7.6),
      plot.margin = margin(3, 4, 3, 4)
    )

  p_cluster <- monocle3::plot_cells(
    cds,
    color_cells_by = cluster_col,
    label_cell_groups = TRUE,
    label_leaves = FALSE,
    label_roots = FALSE,
    label_branch_points = FALSE,
    graph_label_size = 3.0,
    cell_size = 0.45,
    trajectory_graph_segment_size = 0.95
  ) +
    scale_color_manual(values = cluster_cols, drop = FALSE) +
    labs(title = "Tumor clusters", color = "Tumor cluster") +
    guides(color = guide_legend(ncol = 2, override.aes = list(size = 2.7, alpha = 1))) +
    base_traj_theme +
    theme(legend.position = "right")

  p_time <- monocle3::plot_cells(
    cds,
    color_cells_by = "pseudotime",
    label_cell_groups = FALSE,
    label_leaves = FALSE,
    label_roots = FALSE,
    label_branch_points = FALSE,
    cell_size = 0.45,
    trajectory_graph_segment_size = 0.95
  ) +
    scale_color_viridis_c(option = "plasma", name = "Pseudotime") +
    labs(title = "Pseudotime") +
    guides(color = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      barwidth = unit(24, "mm"),
      barheight = unit(3.0, "mm")
    )) +
    base_traj_theme +
    theme(legend.position = "bottom")

  save_requested_pdf(p_cluster + p_time + plot_layout(widths = c(1.18, 1)), file.path(fig3_out, "fig3e_pseudotime_trajectory.pdf"), 7.2, 3.45)
}

make_fig3i_go_pc <- function() {
  message("  Figure 3i: GO enrichment for updated PC-enriched AUC genes")
  if (!requireNamespace("clusterProfiler", quietly = TRUE) || !requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    stop("clusterProfiler and org.Hs.eg.db are required for fig3i_pc_gained_go_enrichment.pdf")
  }

  deg <- read.csv(file.path(fig3_dir, "DEGs_PCvsCC_AUC_Result.csv"), header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  blank_cols <- which(!nzchar(colnames(deg)) | is.na(colnames(deg)))
  if (length(blank_cols) > 0) {
    if ("gene" %in% colnames(deg)) {
      deg <- deg[, -blank_cols, drop = FALSE]
    } else {
      colnames(deg)[blank_cols[1]] <- "gene"
    }
  }
  if (!"gene" %in% colnames(deg)) {
    first_col <- colnames(deg)[1]
    if (!first_col %in% c("myAUC", "avg_diff", "power", "avg_log2FC", "final_AUC", "Group")) {
      colnames(deg)[1] <- "gene"
    }
  }
  pc_genes <- deg %>%
    mutate(myAUC = as.numeric(myAUC), avg_log2FC = as.numeric(avg_log2FC)) %>%
    filter(myAUC > 0.7, avg_log2FC > 0.25) %>%
    pull(gene) %>%
    unique()
  if (length(pc_genes) < 5) {
    stop("Too few PC-enriched genes were found for GO enrichment: ", length(pc_genes))
  }

  gene_convert <- clusterProfiler::bitr(
    pc_genes,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db::org.Hs.eg.db
  )
  ego_pc <- clusterProfiler::enrichGO(
    gene = unique(gene_convert$ENTREZID),
    OrgDb = org.Hs.eg.db::org.Hs.eg.db,
    ont = "BP",
    readable = TRUE
  )
  if (is.null(ego_pc) || nrow(as.data.frame(ego_pc)) == 0) {
    stop("No enriched GO biological-process terms were returned for PC genes.")
  }

  wrap_axis_label <- function(x, width = 34) {
    vapply(x, function(label) paste(strwrap(label, width = width), collapse = "\n"), character(1))
  }

  p <- clusterProfiler::dotplot(ego_pc, showCategory = 20, title = "Functions gained in PC group (122 genes)") +
    scale_color_gradientn(
      colors = c("#B2182B", "#FDDDBC", "#2166AC"),
      trans = "reverse",
      name = "FDR",
      guide = guide_colorbar(barheight = unit(18, "mm"), barwidth = unit(3.2, "mm"))
    ) +
    scale_size(range = c(1.6, 4.4), name = "Count", breaks = pretty_breaks(n = 4)) +
    scale_x_continuous(breaks = pretty_breaks(n = 4), expand = expansion(mult = c(0.03, 0.08))) +
    scale_y_discrete(labels = wrap_axis_label) +
    labs(x = "Gene ratio", y = NULL) +
    theme_nature(base_size = 8.8) +
    theme(
      plot.title = element_text(face = "plain", size = 10.0, hjust = 0.5),
      panel.grid.major.x = element_line(color = "grey92", linewidth = 0.22),
      panel.grid.major.y = element_blank(),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y = element_text(size = 8.4, lineheight = 0.90),
      axis.text.x = element_text(size = 8.6),
      axis.title.x = element_text(face = "plain", size = 9.7),
      legend.position = "right",
      legend.title = element_text(face = "plain", size = 8.3),
      legend.text = element_text(size = 7.7),
      plot.margin = margin(3, 3, 3, 3)
    )

  save_requested_pdf(p, file.path(fig3_out, "fig3i_pc_gained_go_enrichment.pdf"), 4.0, 5.65)
  save_requested_pdf(p, file.path(fig3_dir, "GO_PC_122_Genes.pdf"), 4.0, 5.65)
}

fig2d_to_numeric <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "/", "NA", "NaN", "null")] <- NA_character_
  suppressWarnings(as.numeric(x))
}

fig2d_pick_col <- function(df, pattern) {
  hit <- grep(pattern, names(df), value = TRUE)
  if (length(hit) == 0) stop("Cannot find column matching pattern: ", pattern)
  hit[1]
}

fig2d_group_p <- function(null_fit, full_fit) {
  cmp <- anova(null_fit, full_fit, test = "Chisq")
  p_col <- grep("P", colnames(cmp), value = TRUE)
  if (length(p_col) == 0 || nrow(cmp) < 2) return(NA_real_)
  suppressWarnings(as.numeric(cmp[2, p_col[1]]))
}

fig2d_format_p <- function(p, prefix = "Age-adjusted Cox P") {
  if (is.na(p)) {
    paste0(prefix, " = NA")
  } else if (p < 0.001) {
    paste0(prefix, " < 0.001")
  } else {
    paste0(prefix, " = ", formatC(p, format = "f", digits = 3))
  }
}

make_fig2d_survival <- function() {
  message("  Figure 2d: age-adjusted OS/PFS by composition cluster")
  sample_props <- readRDS(file.path(fig2_dir, "sample_props_complete.rds"))
  clinical_df <- read.csv(file.path(fig2_dir, "clinical.csv"), header = TRUE, stringsAsFactors = FALSE, check.names = TRUE)
  key_cells <- c("Chordoma", "Macrophage", "T cells")

  prop_key_df <- sample_props %>%
    mutate(cell_type = as.character(cell_type)) %>%
    filter(cell_type %in% key_cells) %>%
    select(orig.ident, cell_type, freq) %>%
    pivot_wider(names_from = cell_type, values_from = freq, values_fill = 0) %>%
    arrange(orig.ident)

  prop_key <- as.matrix(prop_key_df[, key_cells, drop = FALSE])
  rownames(prop_key) <- prop_key_df$orig.ident
  hclust_res <- hclust(dist(prop_key, method = "euclidean"), method = "ward.D2")
  sample_clusters <- data.frame(
    Sample = rownames(prop_key),
    Cluster = factor(cutree(hclust_res, k = 3), levels = as.character(1:3)),
    stringsAsFactors = FALSE
  )

  os_time_col <- fig2d_pick_col(clinical_df, "^OS\\.Month")
  pfs_time_col <- fig2d_pick_col(clinical_df, "^PFS\\.Month")
  surv_df <- clinical_df %>%
    left_join(sample_clusters, by = "Sample") %>%
    mutate(
      OS_time = fig2d_to_numeric(.data[[os_time_col]]),
      OS_status = fig2d_to_numeric(OS_event),
      PFS_time = fig2d_to_numeric(.data[[pfs_time_col]]),
      PFS_status = fig2d_to_numeric(PFS_event),
      age = fig2d_to_numeric(age),
      Cluster = factor(Cluster, levels = as.character(1:3))
    ) %>%
    filter(!is.na(Cluster))

  cluster_pal <- setNames(c("#D55E00", "#0072B2", "#009E73"), c("1", "2", "3"))
  age_median <- median(surv_df$age, na.rm = TRUE)

  plot_one <- function(endpoint, time_col, status_col, out_file, ylab, root_alias) {
    plot_df <- surv_df %>%
      filter(!is.na(.data[[time_col]]), !is.na(.data[[status_col]]), !is.na(age), !is.na(Cluster)) %>%
      mutate(Cluster = droplevels(Cluster))
    if (nrow(plot_df) == 0 || length(unique(plot_df$Cluster)) < 2) {
      stop(endpoint, " survival plot needs at least two non-empty clusters.")
    }

    cox_fit <- coxph(as.formula(paste0("Surv(", time_col, ", ", status_col, ") ~ Cluster + age")), data = plot_df)
    cox_null <- coxph(as.formula(paste0("Surv(", time_col, ", ", status_col, ") ~ age")), data = plot_df)
    p_label <- fig2d_format_p(fig2d_group_p(cox_null, cox_fit))

    newdata <- data.frame(
      Cluster = factor(levels(plot_df$Cluster), levels = levels(plot_df$Cluster)),
      age = rep(age_median, length(levels(plot_df$Cluster)))
    )
    fit_adj <- survfit(cox_fit, newdata = newdata)
    pal <- cluster_pal[levels(plot_df$Cluster)]

    surv_plot <- ggsurvplot(
      fit_adj,
      data = newdata,
      conf.int = FALSE,
      risk.table = FALSE,
      censor = FALSE,
      break.time.by = 12,
      palette = unname(pal),
      xlab = "Time (months)",
      ylab = ylab,
      legend.title = "Cluster",
      legend.labs = levels(plot_df$Cluster),
      size = 0.85,
      ggtheme = theme_nature(base_size = 9.4)
    )
    max_time <- max(plot_df[[time_col]], na.rm = TRUE)
    surv_plot$plot <- surv_plot$plot +
      annotate("text", x = max_time * 0.58, y = 0.14, label = p_label, size = 2.7, family = FONT_FAMILY) +
      scale_x_continuous(breaks = pretty_breaks(n = 8), expand = expansion(mult = c(0.02, 0.03))) +
      scale_y_continuous(breaks = pretty_breaks(n = 6), limits = c(0, 1), expand = expansion(mult = c(0.01, 0.03))) +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.80, 0.80),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.title = element_text(face = "plain", size = 8.5),
        legend.text = element_text(size = 8.0),
        axis.title = element_text(face = "plain", size = 9.5),
        axis.text = element_text(size = 8.5),
        plot.margin = margin(4, 5, 4, 5)
      )

    save_requested_pdf(surv_plot$plot, file.path(fig2_out, out_file), 3.75, 3.15)
    save_requested_pdf(surv_plot$plot, file.path(fig2_dir, root_alias), 3.75, 3.15)
  }

  plot_one("OS", "OS_time", "OS_status", "fig2d_overall_survival_by_composition_cluster.pdf", "Overall survival probability", "clusterOS_3cluster_age_adjusted.pdf")
  plot_one("PFS", "PFS_time", "PFS_status", "fig2d_progression_free_survival_by_composition_cluster.pdf", "Progression-free survival probability", "PFS_3cluster_age_adjusted.pdf")
}

make_fig3d_cytotrace()
make_fig3e_monocle()
make_fig3i_go_pc()
make_fig2d_survival()

message("Requested individual panels regenerated.")
