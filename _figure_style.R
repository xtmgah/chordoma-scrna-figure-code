options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
})

FONT_FAMILY <- "Roboto Condensed"

pdfhr2 <- function(...) {
  if (!requireNamespace("showtext", quietly = TRUE)) {
    warning("Package showtext is not installed; using the active PDF device font fallback.", call. = FALSE)
    return(invisible(FALSE))
  }
  if (!requireNamespace("sysfonts", quietly = TRUE)) {
    warning("Package sysfonts is not installed; using the active PDF device font fallback.", call. = FALSE)
    return(invisible(FALSE))
  }
  font_files <- c(
    regular = Sys.getenv("ROBOTO_CONDENSED_REGULAR", "/Users/zhangt8/Library/Fonts/RobotoCondensed-Regular.ttf"),
    bold = Sys.getenv("ROBOTO_CONDENSED_BOLD", "/Users/zhangt8/Library/Fonts/RobotoCondensed-Bold.ttf"),
    italic = Sys.getenv("ROBOTO_CONDENSED_ITALIC", "/Users/zhangt8/Library/Fonts/RobotoCondensed-Italic.ttf"),
    bolditalic = Sys.getenv("ROBOTO_CONDENSED_BOLDITALIC", "/Users/zhangt8/Library/Fonts/RobotoCondensed-BlackItalic.ttf")
  )
  if (!all(file.exists(font_files))) {
    missing <- paste(names(font_files)[!file.exists(font_files)], collapse = ", ")
    warning(
      "Roboto Condensed font files were not found for: ", missing,
      ". Install Roboto Condensed or set ROBOTO_CONDENSED_* environment variables.",
      call. = FALSE
    )
    return(invisible(FALSE))
  }
  sysfonts::font_add(
    family = FONT_FAMILY,
    regular = font_files[["regular"]],
    bold = font_files[["bold"]],
    italic = font_files[["italic"]],
    bolditalic = font_files[["bolditalic"]]
  )
  showtext::showtext_auto()
  showtext::showtext_opts(dpi = 600)
  invisible(TRUE)
}

pdfhr2()

nature_palette <- c(
  "#3B5BA5", "#D55E00", "#009E73", "#CC79A7", "#0072B2", "#E69F00",
  "#6A3D9A", "#56B4E9", "#B15928", "#2F7F6F", "#8C564B", "#E15759",
  "#4E79A7", "#59A14F", "#F28E2B", "#9C755F", "#76B7B2", "#EDC948",
  "#B07AA1", "#FF9DA7", "#BAB0AC", "#1F77B4", "#FF7F0E", "#2CA02C",
  "#9467BD", "#8C6D31", "#17BECF", "#A55194", "#393B79", "#637939",
  "#8C6D31", "#843C39", "#7B4173", "#3182BD", "#E6550D", "#31A354",
  "#756BB1", "#636363", "#DE2D26", "#08519C", "#006D2C", "#54278F"
)

ncicolor <- {
  pal <- if (requireNamespace("ggsci", quietly = TRUE)) {
    c(
      ggsci::pal_npg(alpha = 1)(10),
      ggsci::pal_jco(alpha = 1)(10),
      ggsci::pal_lancet(alpha = 1)(9),
      ggsci::pal_nejm(alpha = 1)(8)
    )
  } else {
    nature_palette
  }
  unique(stats::na.omit(pal))
}

nci_at <- function(i, fallback) {
  if (length(ncicolor) >= i && !is.na(ncicolor[[i]])) ncicolor[[i]] else fallback
}

subgroup_palette_values <- function(n) {
  pal <- if (requireNamespace("ggsci", quietly = TRUE)) {
    ggsci::pal_primer(palette = c("mark17"), alpha = 1)(min(max(n, 1), 17))
  } else {
    c(
      "#006EDB", "#EB670F", "#DF0C24", "#179B9B", "#30A147", "#894CEB",
      "#B88700", "#CE2C85", "#856D4C", "#527A29", "#D43511", "#167E53",
      "#9D615C", "#64762D", "#A830E8", "#866E04", "#808FA3"
    )
  }
  pal <- unique(stats::na.omit(pal))
  if (n <= length(pal)) {
    pal[seq_len(n)]
  } else {
    grDevices::colorRampPalette(pal)(n)
  }
}

group_palette <- c(
  "CC" = "#4DBBD5",
  "PC" = "#7E5A9B",
  "DC" = "#F28E2B",
  "CC-PC_enriched" = "#D55E00",
  "CC-enriched" = "#4DBBD5",
  "PC-enriched" = "#7E5A9B"
)

cell_type_palette <- c(
  "Chordoma" = nci_at(1, "#E64B35"),
  "Macrophage" = nci_at(2, "#4DBBD5"),
  "T cells" = nci_at(3, "#00A087"),
  "Fibroblast" = nci_at(4, "#3C5488"),
  "Neutrophil" = nci_at(5, "#F39B7F"),
  "Cycling" = nci_at(6, "#8491B4"),
  "B cells" = nci_at(11, "#0073C2"),
  "plasmablasts" = nci_at(12, "#EFC000"),
  "Endothelial" = nci_at(7, "#91D1C2"),
  "Mast cells" = nci_at(9, "#7E6148"),
  "Osteoclast" = nci_at(10, "#B09C85"),
  "DCs" = nci_at(14, "#CD534C"),
  "Mural cells" = nci_at(15, "#7AA6DC"),
  "Plasma cell" = nci_at(18, "#3B3B3B")
)

cc4_subcell_levels <- c(
  "CC_tumor1", "CC_tumor2", "CC_tumor3", "CC_tumor4", "Tumor_PC_enriched",
  "FTL+ Mac", "CCL3+ Mac", "FN1+ Mac", "CD1C+ Mac", "CD3E+ Mac", "NAMPT+ Mac",
  "CD4 memory/helper T", "Effector-memory T", "Activated cytotoxic T",
  "NK-like cytotoxic T", "Treg", "MAIT-like T"
)

cc4_subcell_palette <- setNames(subgroup_palette_values(length(cc4_subcell_levels)), cc4_subcell_levels)

macrophage_palette <- cc4_subcell_palette[c("FTL+ Mac", "CCL3+ Mac", "FN1+ Mac", "CD1C+ Mac", "CD3E+ Mac", "NAMPT+ Mac")]

tcell_palette <- cc4_subcell_palette[c("CD4 memory/helper T", "Effector-memory T", "Activated cytotoxic T", "NK-like cytotoxic T", "Treg", "MAIT-like T")]

tumor_subcell_palette <- cc4_subcell_palette[c("CC_tumor1", "CC_tumor2", "CC_tumor3", "CC_tumor4", "Tumor_PC_enriched")]

lineage_palette_shared <- c(
  "Tumor" = cell_type_palette[["Chordoma"]],
  "Macrophage" = cell_type_palette[["Macrophage"]],
  "T cell" = cell_type_palette[["T cells"]],
  "T cells" = cell_type_palette[["T cells"]]
)

nature_discrete <- function(n) {
  if (n <= length(nature_palette)) {
    return(nature_palette[seq_len(n)])
  }
  grDevices::colorRampPalette(nature_palette)(n)
}

unused_palette_values <- function(n, avoid = character()) {
  pool <- unique(c(ncicolor, nature_palette, subgroup_palette_values(max(n, 17))))
  pool <- pool[!is.na(pool)]
  pool <- setdiff(pool, avoid)
  if (n <= length(pool)) {
    pool[seq_len(n)]
  } else {
    grDevices::colorRampPalette(pool)(n)
  }
}

nature_named_palette <- function(values) {
  values <- unique(as.character(values))
  setNames(unused_palette_values(length(values)), values)
}

canonical_cell_type <- function(values) {
  values <- as.character(values)
  key <- tolower(trimws(values))
  aliases <- c(
    "t cell" = "T cells",
    "t cells" = "T cells",
    "t-cell" = "T cells",
    "t-cells" = "T cells",
    "b cell" = "B cells",
    "b cells" = "B cells",
    "plasmablast" = "plasmablasts",
    "plasmablasts" = "plasmablasts",
    "plasma cells" = "Plasma cell",
    "plasma cell" = "Plasma cell",
    "dc" = "DCs",
    "dcs" = "DCs",
    "dendritic cells" = "DCs",
    "mural cell" = "Mural cells",
    "mural cells" = "Mural cells",
    "mast cell" = "Mast cells",
    "mast cells" = "Mast cells",
    "endothelial cells" = "Endothelial",
    "endothelial" = "Endothelial",
    "tumor" = "Chordoma",
    "tumour" = "Chordoma",
    "chordoma" = "Chordoma"
  )
  out <- values
  hit <- key %in% names(aliases)
  out[hit] <- aliases[key[hit]]
  out
}

palette_for <- function(values, base_palette = nature_palette) {
  values <- unique(as.character(values))
  if (is.null(names(base_palette))) {
    palette_values <- if (length(values) <= length(base_palette)) {
      base_palette[seq_len(length(values))]
    } else {
      grDevices::colorRampPalette(base_palette)(length(values))
    }
    return(setNames(palette_values, values))
  }

  colors <- setNames(rep(NA_character_, length(values)), values)
  canonical <- canonical_cell_type(values)
  for (i in seq_along(values)) {
    value <- values[[i]]
    canonical_value <- canonical[[i]]
    if (value %in% names(base_palette)) {
      colors[[i]] <- base_palette[[value]]
    } else if (canonical_value %in% names(base_palette)) {
      colors[[i]] <- base_palette[[canonical_value]]
    }
  }

  missing <- names(colors)[is.na(colors)]
  if (length(missing) > 0) {
    colors[missing] <- unused_palette_values(length(missing), avoid = colors[!is.na(colors)])
  }
  colors[values]
}

sort_numeric_labels <- function(values) {
  values <- unique(as.character(values))
  numeric_values <- suppressWarnings(as.numeric(values))
  if (all(!is.na(numeric_values))) {
    values[order(numeric_values)]
  } else {
    sort(values)
  }
}

theme_nature <- function(base_size = 8.5, base_family = FONT_FAMILY) {
  theme_classic(base_size = base_size, base_family = base_family) %+replace%
    theme(
      text = element_text(family = base_family, color = "black"),
      plot.title = element_text(size = base_size + 1.5, face = "bold", hjust = 0, margin = margin(b = 4)),
      plot.subtitle = element_text(size = base_size, color = "grey25", margin = margin(b = 5)),
      axis.title = element_text(size = base_size + 0.5, face = "plain"),
      axis.text = element_text(size = base_size, color = "black"),
      axis.line = element_line(linewidth = 0.35, color = "black"),
      axis.ticks = element_line(linewidth = 0.3, color = "black"),
      axis.ticks.length = unit(1.5, "mm"),
      legend.title = element_text(size = base_size, face = "plain"),
      legend.text = element_text(size = base_size - 0.3),
      legend.key.size = unit(2.6, "mm"),
      legend.box.spacing = unit(0.8, "mm"),
      legend.spacing.x = unit(0.8, "mm"),
      legend.spacing.y = unit(0.6, "mm"),
      legend.margin = margin(0, 0, 0, 0),
      legend.background = element_blank(),
      strip.background = element_rect(fill = "grey94", color = NA),
      strip.text = element_text(size = base_size, face = "bold", color = "black"),
      panel.grid.major = element_line(color = "grey92", linewidth = 0.25),
      panel.grid.minor = element_blank(),
      plot.margin = margin(4, 5, 4, 5)
    )
}

theme_umap_nature <- function(base_size = 8.5, base_family = FONT_FAMILY) {
  theme_void(base_size = base_size, base_family = base_family) %+replace%
    theme(
      text = element_text(family = base_family, color = "black"),
      plot.title = element_text(size = base_size + 1.5, face = "bold", hjust = 0, margin = margin(b = 3)),
      legend.title = element_text(size = base_size, face = "plain"),
      legend.text = element_text(size = base_size - 0.4),
      legend.key.size = unit(2.6, "mm"),
      legend.box.spacing = unit(0.8, "mm"),
      legend.spacing.x = unit(0.8, "mm"),
      legend.spacing.y = unit(0.6, "mm"),
      legend.margin = margin(0, 0, 0, 0),
      plot.margin = margin(4, 5, 4, 5)
    )
}

theme_set(theme_nature())

panel_pdf_device <- function(filename, width, height, ...) {
  grDevices::cairo_pdf(
    filename = filename,
    width = width,
    height = height,
    family = FONT_FAMILY,
    onefile = FALSE,
    ...
  )
}

save_panel_pdf <- function(plot, filename, width, height, dpi = 450, limitsize = FALSE) {
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  ggsave(
    filename = filename,
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = dpi,
    device = panel_pdf_device,
    bg = "white",
    limitsize = limitsize
  )
  invisible(filename)
}

open_panel_pdf <- function(filename, width, height) {
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  grDevices::cairo_pdf(filename = filename, width = width, height = height, family = FONT_FAMILY, onefile = FALSE)
}

close_panel_pdf <- function() {
  if (grDevices::dev.cur() > 1) grDevices::dev.off()
}

percent_axis <- function(accuracy = 1) {
  scale_y_continuous(labels = label_percent(accuracy = accuracy), expand = expansion(mult = c(0, 0.04)))
}

clean_label <- function(x) {
  gsub("_", " ", x, fixed = TRUE)
}

current_script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) == 0) return(normalizePath(getwd()))
  dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE))
}

message_panel <- function(path) {
  message("  wrote ", normalizePath(path, mustWork = FALSE))
}

pdf_panel_to_png <- function(pdf, dest = sub("\\.pdf$", ".png", pdf), size = 1800) {
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  if (nzchar(Sys.which("qlmanage"))) {
    tmp_dir <- tempfile("panel_preview_")
    dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
    on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)
    status <- suppressWarnings(system2("qlmanage", c("-t", "-s", as.character(size), "-o", tmp_dir, pdf), stdout = TRUE, stderr = TRUE))
    invisible(status)
    produced <- file.path(tmp_dir, paste0(basename(pdf), ".png"))
    if (file.exists(produced)) {
      file.copy(produced, dest, overwrite = TRUE)
      return(dest)
    }
  }
  if (nzchar(Sys.which("sips"))) {
    status <- suppressWarnings(system2("sips", c("-s", "format", "png", pdf, "--out", dest), stdout = TRUE, stderr = TRUE))
    invisible(status)
    if (file.exists(dest)) return(dest)
  }
  warning("No local PDF-to-PNG preview tool produced output for: ", pdf, call. = FALSE)
  invisible(NA_character_)
}

export_large_panel_pngs <- function(panel_dir, threshold_mb = 1, size = 1800) {
  if (!dir.exists(panel_dir)) return(invisible(character()))
  pdfs <- list.files(panel_dir, pattern = "\\.pdf$", full.names = TRUE)
  pdfs <- pdfs[file.info(pdfs)$size > threshold_mb * 1024^2]
  if (!length(pdfs)) return(invisible(character()))
  out <- vapply(pdfs, function(pdf) {
    dest <- sub("\\.pdf$", ".png", pdf)
    if (!file.exists(dest) || file.info(dest)$mtime < file.info(pdf)$mtime) {
      pdf_panel_to_png(pdf, dest = dest, size = size)
    } else {
      dest
    }
  }, character(1))
  invisible(out)
}
