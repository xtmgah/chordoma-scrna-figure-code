script_dir <- if (exists("SCRIPT_DIR", inherits = FALSE)) {
  SCRIPT_DIR
} else {
  dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), mustWork = TRUE))
}
source(file.path(dirname(script_dir), "_figure_style.R"))

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(ComplexHeatmap)
  library(pheatmap)
  library(survival)
  library(survminer)
  library(cluster)
  library(RColorBrewer)
  library(circlize)
})

OUT_DIR <- file.path(script_dir, "nature_panels")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

plot_cellchat_networks <- function(cellchat_obj, prefix, source_name = "Chordoma") {
  group_size <- as.numeric(table(cellchat_obj@idents))
  names(group_size) <- levels(cellchat_obj@idents)
  idents <- levels(cellchat_obj@idents)
  colors <- palette_for(idents, c(cell_type_palette, macrophage_palette, tcell_palette))

  circle_one <- function(mat, out, title, edge_max = NULL) {
    open_panel_pdf(file.path(OUT_DIR, out), 3.9, 3.9)
    par(family = FONT_FAMILY, mar = c(0.4, 0.4, 1.2, 0.4), xpd = TRUE)
    args <- list(
      net = mat,
      vertex.weight = group_size,
      weight.scale = TRUE,
      label.edge = FALSE,
      title.name = title,
      color.use = colors,
      vertex.label.cex = 0.58,
      edge.width.max = 2.4,
      alpha.edge = 0.32,
      arrow.size = 0.12,
      arrow.width = 0.6
    )
    if (!is.null(edge_max)) args$edge.weight.max <- edge_max
    do.call(CellChat::netVisual_circle, args)
    close_panel_pdf()
  }

  circle_one(cellchat_obj@net$count, paste0(prefix, "_network_count.pdf"), "Number of interactions")
  circle_one(cellchat_obj@net$weight, paste0(prefix, "_network_strength.pdf"), "Interaction strength")

  mat <- cellchat_obj@net$weight
  source_name <- intersect(source_name, rownames(mat))[1]
  if (!is.na(source_name)) {
    mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    mat2[source_name, ] <- mat[source_name, ]
    circle_one(mat2, paste0(prefix, "_", make.names(source_name), "_outgoing.pdf"), paste0(source_name, " to others"), max(mat, na.rm = TRUE))
  }

  p_scatter <- CellChat::netAnalysis_signalingRole_scatter(cellchat_obj, color.use = colors, do.label = TRUE, label.size = 2.5) +
    labs(x = "Outgoing signaling strength", y = "Incoming signaling strength", color = "Cell type") +
    theme_nature() +
    theme(legend.position = "right")
  save_panel_pdf(p_scatter, file.path(OUT_DIR, paste0(prefix, "_signaling_role_scatter.pdf")), 4.8, 3.6)

  open_panel_pdf(file.path(OUT_DIR, paste0(prefix, "_signaling_role_heatmap.pdf")), 11.8, 7.4)
  ht_out <- make_cellchat_role_heatmap_original(cellchat_obj, pattern = "outgoing", color.use = colors)
  ht_in <- make_cellchat_role_heatmap_original(cellchat_obj, pattern = "incoming", color.use = colors)
  ComplexHeatmap::draw(ht_out + ht_in, heatmap_legend_side = "right", annotation_legend_side = "bottom", merge_legends = TRUE, padding = unit(c(4, 3, 3, 3), "mm"))
  close_panel_pdf()
}

cellchat_pathway_strength <- function(cellchat_obj) {
  probs <- cellchat_obj@netP$prob
  pathways <- cellchat_obj@netP$pathways
  if (!is.null(probs) && length(dim(probs)) == 3) {
    strength <- apply(probs, 3, sum, na.rm = TRUE)
    strength <- strength[intersect(names(strength), pathways)]
  } else {
    strength <- setNames(rep(NA_real_, length(pathways)), pathways)
  }
  strength
}

top_cellchat_pathways <- function(cellchat_objects, n = 3) {
  scores <- lapply(cellchat_objects, cellchat_pathway_strength)
  common <- Reduce(intersect, lapply(scores, names))
  if (length(common) == 0) return(character())
  combined <- Reduce(`+`, lapply(scores, function(x) {
    x <- x[common]
    x[is.na(x)] <- 0
    x
  }))
  names(sort(combined, decreasing = TRUE))[seq_len(min(n, length(combined)))]
}

complete_named_colors <- function(values, preferred = NULL) {
  values <- unique(as.character(stats::na.omit(values)))
  if (length(values) == 0) return(character())
  preferred <- preferred %||% character()
  preferred <- preferred[names(preferred) %in% values]
  missing <- setdiff(values, names(preferred))
  c(preferred, nature_named_palette(missing))[values]
}

style_cellchat_heatmap <- function(ht, row_font = 6.0, column_font = 6.2, title_font = 10.8) {
  ht@row_names_param$gp <- gpar(fontsize = row_font, lineheight = 0.88, fontfamily = FONT_FAMILY)
  ht@column_names_param$gp <- gpar(fontsize = column_font, lineheight = 0.88, fontfamily = FONT_FAMILY)
  ht@column_title_param$gp <- gpar(fontsize = title_font, lineheight = 0.9, fontfamily = FONT_FAMILY, fontface = "plain")
  ht@matrix_legend_param$title_gp <- gpar(fontsize = 8.8, fontfamily = FONT_FAMILY, fontface = "plain")
  ht@matrix_legend_param$labels_gp <- gpar(fontsize = 8.0, fontfamily = FONT_FAMILY)
  ht
}

make_cellchat_role_heatmap_original <- function(cellchat_obj, pattern, title = NULL, display_title = NULL, color.use = NULL, signaling = NULL, width = 8.8, height = 14.2, title_font = 10.8) {
  ht <- CellChat::netAnalysis_signalingRole_heatmap(
    cellchat_obj,
    pattern = pattern,
    signaling = signaling,
    title = title,
    font.size = 5.6,
    font.size.title = 10,
    width = width,
    height = height,
    color.use = color.use,
    cluster.rows = FALSE,
    cluster.cols = FALSE
  )
  ht <- style_cellchat_heatmap(ht, title_font = title_font)
  if (!is.null(display_title)) {
    ht@column_title <- display_title
    ht@column_title_param$gp <- gpar(fontsize = title_font, lineheight = 0.9, fontfamily = FONT_FAMILY, fontface = "plain")
  }
  ht
}

make_cellchat_role_heatmap <- function(cellchat_obj, pattern, title, color.use = NULL, signaling = NULL, transpose = TRUE, heatmap_height = NULL, heatmap_width = NULL, sort_pathways = TRUE, show_celltype_legend = TRUE, title_position = c("column", "row")) {
  title_position <- match.arg(title_position)
  raw_ht <- CellChat::netAnalysis_signalingRole_heatmap(
    cellchat_obj,
    pattern = pattern,
    signaling = signaling,
    title = title,
    font.size = 7,
    color.use = color.use,
    cluster.rows = FALSE,
    cluster.cols = FALSE
  )
  mat <- raw_ht@matrix
  if (transpose) mat <- t(mat)
  mat[is.na(mat)] <- 0

  if (transpose && sort_pathways && ncol(mat) > 1) {
    pathway_order <- names(sort(colSums(mat, na.rm = TRUE), decreasing = TRUE))
    mat <- mat[, pathway_order, drop = FALSE]
  } else if (!transpose && sort_pathways && nrow(mat) > 1) {
    pathway_order <- names(sort(rowSums(mat, na.rm = TRUE), decreasing = TRUE))
    mat <- mat[pathway_order, , drop = FALSE]
  }

  max_val <- max(mat, na.rm = TRUE)
  if (!is.finite(max_val) || max_val <= 0) max_val <- 1
  column_strength <- colSums(mat, na.rm = TRUE)
  strength_anno <- anno_barplot(
    column_strength,
    gp = gpar(fill = "grey65", col = "grey35", lwd = 0.1),
    border = FALSE,
    height = unit(if (transpose) 10 else 8, "mm")
  )
  cell_type_colors <- if (is.null(color.use)) {
    NULL
  } else {
    palette_for(union(rownames(mat), colnames(mat)), color.use)
  }
  cell_type_legend <- list(
    title = "Cell type",
    direction = "horizontal",
    title_position = "leftcenter",
    nrow = 2,
    title_gp = gpar(fontsize = 8.2, fontfamily = FONT_FAMILY, fontface = "plain"),
    labels_gp = gpar(fontsize = 7.3, fontfamily = FONT_FAMILY),
    grid_height = unit(3, "mm"),
    grid_width = unit(3.5, "mm")
  )

  top_anno <- if (transpose) {
    HeatmapAnnotation(
      Spacer = anno_empty(height = unit(4, "mm"), border = FALSE),
      Strength = strength_anno,
      annotation_name_gp = gpar(fontsize = 8, fontfamily = FONT_FAMILY, fontface = "plain"),
      annotation_name_side = "left",
      show_annotation_name = c(Spacer = FALSE, Strength = TRUE)
    )
  } else {
    HeatmapAnnotation(
      `Cell type` = factor(colnames(mat), levels = colnames(mat)),
      Strength = strength_anno,
      col = list(`Cell type` = cell_type_colors[colnames(mat)]),
      annotation_name_gp = gpar(fontsize = 8, fontfamily = FONT_FAMILY, fontface = "plain"),
      annotation_name_side = "left",
      show_annotation_name = c(`Cell type` = TRUE, Strength = TRUE),
      show_legend = c(`Cell type` = show_celltype_legend, Strength = TRUE),
      annotation_legend_param = list(`Cell type` = cell_type_legend)
    )
  }
  left_anno <- NULL
  if (transpose) {
    left_anno <- rowAnnotation(
      `Cell type` = factor(rownames(mat), levels = rownames(mat)),
      col = list(`Cell type` = cell_type_colors[rownames(mat)]),
      simple_anno_size = unit(2.4, "mm"),
      show_annotation_name = FALSE,
      show_legend = c(`Cell type` = show_celltype_legend),
      annotation_legend_param = list(`Cell type` = cell_type_legend)
    )
  }

  Heatmap(
    mat,
    name = "Relative strength",
    col = colorRamp2(c(0, max_val * 0.5, max_val), c("white", "#C6DBEF", "#08519C")),
    top_annotation = top_anno,
    left_annotation = left_anno,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_names_side = "left",
    row_names_gp = gpar(fontsize = if (transpose) 7.4 else 6.4, fontfamily = FONT_FAMILY),
    column_names_gp = gpar(fontsize = if (transpose) 6.1 else 6.4, fontfamily = FONT_FAMILY),
    column_names_rot = if (transpose) 90 else 45,
    column_title = if (identical(title_position, "column")) title else NULL,
    column_title_gp = gpar(fontsize = 9, fontfamily = FONT_FAMILY, fontface = "plain"),
    row_title = if (identical(title_position, "row")) title else NULL,
    row_title_rot = 0,
    row_title_gp = gpar(fontsize = 9, fontfamily = FONT_FAMILY, fontface = "plain"),
    rect_gp = gpar(col = "white", lwd = 0.1),
    width = heatmap_width %||% unit(if (transpose) 6.6 else 2.9, "in"),
    height = heatmap_height %||% unit(if (transpose) 2.25 else 0.82, "in"),
    heatmap_legend_param = list(
      direction = "horizontal",
      title_gp = gpar(fontsize = 8, fontfamily = FONT_FAMILY, fontface = "plain"),
      labels_gp = gpar(fontsize = 7.3, fontfamily = FONT_FAMILY),
      title_position = "leftcenter",
      legend_width = unit(if (transpose) 36 else 28, "mm"),
      legend_height = unit(3, "mm")
    )
  )
}

`%||%` <- function(x, y) if (is.null(x)) y else x

message("Figure 2: CellChat communication panels")
if (requireNamespace("CellChat", quietly = TRUE)) {
  suppressPackageStartupMessages(library(CellChat))
  cellchat_all <- readRDS(file.path(script_dir, "cellchat_object.rds"))
  cellchat_pc <- readRDS(file.path(script_dir, "cellchat_PC.rds"))
  cellchat_cc <- readRDS(file.path(script_dir, "cellchat_CC.rds"))

  plot_cellchat_networks(cellchat_all, "fig2a_all_cellchat")
  plot_cellchat_networks(cellchat_pc, "fig2b_pc_cellchat")
  plot_cellchat_networks(cellchat_cc, "fig2c_cc_cellchat")

  all_pathways <- names(sort(cellchat_pathway_strength(cellchat_all), decreasing = TRUE))
  pathways <- intersect(all_pathways, intersect(cellchat_pc@netP$pathways, cellchat_cc@netP$pathways))
  pathways <- pathways[seq_len(min(3, length(pathways)))]
  if (length(pathways) > 0) {
    unlink(list.files(OUT_DIR, pattern = "^fig2h_.*_(pc|cc)_hierarchy\\.pdf$", full.names = TRUE))
    unlink(file.path(OUT_DIR, "fig2i_mif_fn1_app_signaling_comparison.pdf"))
    unlink(list.files(OUT_DIR, pattern = "^fig2i_top3_.*_signaling_comparison\\.pdf$", full.names = TRUE))
    for (pathway in pathways) {
      pc_cols <- palette_for(levels(cellchat_pc@idents), c(cell_type_palette, macrophage_palette, tcell_palette))
      cc_cols <- palette_for(levels(cellchat_cc@idents), c(cell_type_palette, macrophage_palette, tcell_palette))
      pc_group_size <- as.numeric(table(cellchat_pc@idents))
      names(pc_group_size) <- levels(cellchat_pc@idents)
      cc_group_size <- as.numeric(table(cellchat_cc@idents))
      names(cc_group_size) <- levels(cellchat_cc@idents)
      vertex_weight_max <- max(c(pc_group_size, cc_group_size), na.rm = TRUE)

      open_panel_pdf(file.path(OUT_DIR, paste0("fig2h_", pathway, "_pc_hierarchy.pdf")), 8.2, 5.2)
      par(family = FONT_FAMILY, mar = c(0.7, 3.5, 1.5, 3.5), xpd = NA)
      CellChat::netVisual_aggregate(cellchat_pc, signaling = pathway, vertex.receiver = seq_len(min(3, length(levels(cellchat_pc@idents)))), layout = "hierarchy", color.use = pc_cols, vertex.weight = pc_group_size, vertex.weight.max = vertex_weight_max, vertex.size.max = 20, vertex.label.cex = 0.58, edge.width.max = 2.35, alpha.edge = 0.45, arrow.size = 0.12, arrow.width = 0.8, pt.title = 10.5, title.space = 3.5)
      close_panel_pdf()

      open_panel_pdf(file.path(OUT_DIR, paste0("fig2h_", pathway, "_cc_hierarchy.pdf")), 8.2, 5.2)
      par(family = FONT_FAMILY, mar = c(0.7, 3.5, 1.5, 3.5), xpd = NA)
      CellChat::netVisual_aggregate(cellchat_cc, signaling = pathway, vertex.receiver = seq_len(min(3, length(levels(cellchat_cc@idents)))), layout = "hierarchy", color.use = cc_cols, vertex.weight = cc_group_size, vertex.weight.max = vertex_weight_max, vertex.size.max = 20, vertex.label.cex = 0.58, edge.width.max = 2.35, alpha.edge = 0.45, arrow.size = 0.12, arrow.width = 0.8, pt.title = 10.5, title.space = 3.5)
      close_panel_pdf()
    }

    comparison_file <- paste0("fig2i_top3_", paste(tolower(make.names(pathways)), collapse = "_"), "_signaling_comparison.pdf")
    open_panel_pdf(file.path(OUT_DIR, comparison_file), 6.2, 2.65)
    ht_pc <- make_cellchat_role_heatmap_original(
      cellchat_pc,
      pattern = "all",
      signaling = pathways,
      display_title = "PC",
      color.use = pc_cols,
      width = 3.1,
      height = 1.05,
      title_font = 9.2
    )
    ht_cc <- make_cellchat_role_heatmap_original(
      cellchat_cc,
      pattern = "all",
      signaling = pathways,
      display_title = "CC",
      color.use = cc_cols,
      width = 3.1,
      height = 1.05,
      title_font = 9.2
    )
    ComplexHeatmap::draw(ht_pc + ht_cc, heatmap_legend_side = "right", annotation_legend_side = "bottom", merge_legends = TRUE, padding = unit(c(4, 2, 2, 2), "mm"))
    close_panel_pdf()
  }
} else {
  message("  CellChat is not installed; skipped Figure 2 communication panels.")
}

message("Figure 2: sample composition clustering and survival panels")
sc <- readRDS(file.path(script_dir, "sc_combined_clean_clinical_full.rds"))
meta <- sc@meta.data
key_cells <- intersect(c("Chordoma", "Macrophage", "T cells"), unique(as.character(meta$cell_type)))
if (length(key_cells) >= 2) {
  prop_table <- table(meta$cell_type, meta$orig.ident) %>% prop.table(margin = 2)
  prop_key <- t(as.matrix(prop_table[key_cells, , drop = FALSE]))
  dist_mat <- dist(prop_key, method = "euclidean")
  hclust_res <- hclust(dist_mat, method = "ward.D2")
  sample_clusters <- data.frame(Sample = rownames(prop_key), Cluster = factor(cutree(hclust_res, k = 3)), stringsAsFactors = FALSE)

  anno_cols <- intersect(c("orig.ident", "clinical_group", "gender", "location", "status", "OS_event", "PFS_event"), colnames(meta))
  anno_df <- meta %>%
    as.data.frame() %>%
    remove_rownames() %>%
    select(all_of(anno_cols)) %>%
    distinct() %>%
    column_to_rownames("orig.ident")
  anno_df$Cluster <- sample_clusters$Cluster[match(rownames(anno_df), sample_clusters$Sample)]
  preferred_ann_colors <- list(
    clinical_group = group_palette[c("CC", "PC")],
    Cluster = setNames(c("#D55E00", "#0072B2", "#009E73"), c("1", "2", "3")),
    gender = c("female" = "#CC79A7", "male" = "#56B4E9"),
    location = c("clivus" = "#0072B2", "mobile spine" = "#E69F00", "sacrum" = "#009E73"),
    status = c("primary" = "#009E73", "recurrence" = "#D55E00"),
    OS_event = c("0" = "grey90", "1" = "black", "/" = "white"),
    PFS_event = c("0" = "grey90", "1" = "#D55E00", "/" = "white")
  )
  annotation_labels <- c(
    clinical_group = "Group",
    gender = "Gender",
    location = "Location",
    status = "Status",
    OS_event = "OS event",
    PFS_event = "PFS event",
    Cluster = "Cluster"
  )
  anno_df <- anno_df[, intersect(names(annotation_labels), colnames(anno_df)), drop = FALSE]
  anno_plot <- anno_df
  for (col_name in colnames(anno_plot)) {
    anno_plot[[col_name]] <- factor(as.character(anno_plot[[col_name]]), levels = unique(as.character(anno_plot[[col_name]])))
  }
  ann_colors <- lapply(colnames(anno_plot), function(col_name) {
    complete_named_colors(anno_plot[[col_name]], preferred_ann_colors[[col_name]])
  })
  names(ann_colors) <- colnames(anno_plot)
  annotation_legend_param <- lapply(colnames(anno_plot), function(col_name) {
    list(
      title = annotation_labels[[col_name]],
      direction = "horizontal",
      title_position = "leftcenter",
      nrow = 1,
      title_gp = gpar(fontsize = 8, fontfamily = FONT_FAMILY, fontface = "plain"),
      labels_gp = gpar(fontsize = 7.4, fontfamily = FONT_FAMILY),
      grid_height = unit(3.2, "mm"),
      grid_width = unit(4, "mm")
    )
  })
  names(annotation_legend_param) <- colnames(anno_plot)
  heatmap_pdf <- file.path(OUT_DIR, "fig2d_patient_composition_cluster_heatmap.pdf")
  open_panel_pdf(heatmap_pdf, 5.8, 2.8)
  heat_mat <- t(prop_key)
  heat_mat <- t(scale(t(heat_mat)))
  heat_mat[is.na(heat_mat)] <- 0
  top_anno <- HeatmapAnnotation(
    df = anno_plot,
    col = ann_colors,
    simple_anno_size = unit(1.45, "mm"),
    annotation_name_gp = gpar(fontsize = 6.4, fontfamily = FONT_FAMILY, fontface = "plain"),
    annotation_name_side = "left",
    annotation_label = unname(annotation_labels[colnames(anno_plot)]),
    show_annotation_name = TRUE,
    show_legend = setNames(rep(TRUE, ncol(anno_plot)), colnames(anno_plot)),
    annotation_legend_param = annotation_legend_param
  )
  ht <- Heatmap(
    heat_mat,
    name = "Scaled fraction",
    col = colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B")),
    top_annotation = top_anno,
    cluster_rows = FALSE,
    cluster_columns = hclust_res,
    show_column_names = TRUE,
    row_names_side = "left",
    column_names_gp = gpar(fontsize = 7.2, fontfamily = FONT_FAMILY),
    row_names_gp = gpar(fontsize = 7.2, fontfamily = FONT_FAMILY),
    row_names_max_width = unit(16, "mm"),
    column_names_rot = 90,
    rect_gp = gpar(col = NA),
    heatmap_legend_param = list(
      direction = "horizontal",
      title_position = "leftcenter",
      title_gp = gpar(fontsize = 8, fontfamily = FONT_FAMILY),
      labels_gp = gpar(fontsize = 7.4, fontfamily = FONT_FAMILY),
      legend_width = unit(20, "mm"),
      legend_height = unit(2.5, "mm")
    )
  )
  draw(ht, heatmap_legend_side = "bottom", annotation_legend_side = "bottom", merge_legends = TRUE)
  close_panel_pdf()

  find_col <- function(candidates) {
    hit <- intersect(candidates, colnames(meta))
    if (length(hit) == 0) NA_character_ else hit[[1]]
  }
  os_time_source <- find_col(c("OS.Month.", "OS.Month", "OS_month", "OS_months", "OS.time", "OS_time"))
  pfs_time_source <- find_col(c("PFS.Month.", "PFS.Month", "PFS_month", "PFS_months", "PFS.time", "PFS_time"))
  age_source <- find_col(c("age", "Age"))

  if (!is.na(os_time_source) && !is.na(pfs_time_source) && !is.na(age_source)) {
    sample_clinical <- meta %>%
      as.data.frame() %>%
      remove_rownames() %>%
      select(any_of(c("orig.ident", os_time_source, "OS_event", pfs_time_source, "PFS_event", age_source, "gender"))) %>%
      distinct() %>%
      rename(Sample = orig.ident) %>%
      left_join(sample_clusters, by = "Sample") %>%
      mutate(
        OS_time = suppressWarnings(as.numeric(.data[[os_time_source]])),
        OS_status = suppressWarnings(as.numeric(as.character(OS_event))),
        PFS_time = suppressWarnings(as.numeric(.data[[pfs_time_source]])),
        PFS_status = suppressWarnings(as.numeric(as.character(PFS_event))),
        age = suppressWarnings(as.numeric(.data[[age_source]]))
      )

    plot_adjusted <- function(df, time_col, event_col, out_file, ylab) {
      df <- df %>% filter(!is.na(.data[[time_col]]), !is.na(.data[[event_col]]), !is.na(age), !is.na(Cluster))
      if (nrow(df) == 0 || length(unique(df$Cluster)) < 2) return(invisible(NULL))
      df$Cluster <- droplevels(factor(df$Cluster))
      cox_fit <- coxph(as.formula(paste0("Surv(", time_col, ", ", event_col, ") ~ Cluster + age")), data = df)
      cox_null <- coxph(as.formula(paste0("Surv(", time_col, ", ", event_col, ") ~ age")), data = df)
      compare <- anova(cox_null, cox_fit, test = "Chisq")
      p_col <- grep("P", colnames(compare), value = TRUE)[1]
      p_val <- suppressWarnings(as.numeric(compare[2, p_col]))
      p_label <- if (is.na(p_val)) "Adjusted Cox P = NA" else if (p_val < 0.001) "Adjusted Cox P < 0.001" else paste0("Adjusted Cox P = ", formatC(p_val, format = "f", digits = 3))
      newdata <- data.frame(age = rep(median(df$age, na.rm = TRUE), length(levels(df$Cluster))), Cluster = factor(levels(df$Cluster), levels = levels(df$Cluster)))
      fit <- survfit(cox_fit, newdata = newdata)
      pal <- setNames(c("#D55E00", "#0072B2", "#009E73")[seq_along(levels(df$Cluster))], levels(df$Cluster))
      p <- ggsurvplot(
        fit,
        data = newdata,
        conf.int = FALSE,
        risk.table = FALSE,
        censor = FALSE,
        break.time.by = 12,
        palette = pal,
        xlab = "Time (months)",
        ylab = ylab,
        legend.title = "Cluster",
        legend.labs = levels(df$Cluster),
        size = 0.85,
        ggtheme = theme_nature()
      )
      p$plot <- p$plot +
        annotate("text", x = max(df[[time_col]], na.rm = TRUE) * 0.58, y = 0.14, label = p_label, size = 2.6, family = FONT_FAMILY) +
        theme(legend.position = c(0.82, 0.80), legend.background = element_rect(fill = "white", color = "grey85", linewidth = 0.25))
      save_panel_pdf(p$plot, file.path(OUT_DIR, out_file), 4.5, 3.5)
    }

    plot_adjusted(sample_clinical, "OS_time", "OS_status", "fig2e_overall_survival_by_patient_cluster.pdf", "Overall survival probability")
    plot_adjusted(sample_clinical, "PFS_time", "PFS_status", "fig2f_progression_free_survival_by_patient_cluster.pdf", "Progression-free survival probability")
  } else {
    message("  Skipped Figure 2 survival panels: OS/PFS time and age columns are not present in the local clinical Seurat metadata.")
  }
}

export_large_panel_pngs(OUT_DIR)
message("Figure 2 nature panels written to: ", OUT_DIR)
