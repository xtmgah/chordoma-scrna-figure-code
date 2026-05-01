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
  showtext::font_add(
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

group_palette <- c(
  "CC" = "#4DBBD5",
  "PC" = "#7E5A9B",
  "DC" = "#F28E2B",
  "CC-PC_enriched" = "#D55E00",
  "CC-enriched" = "#4DBBD5",
  "PC-enriched" = "#7E5A9B"
)

cell_type_palette <- c(
  "Chordoma" = "#D55E00",
  "Macrophage" = "#4E79A7",
  "T cells" = "#009E73",
  "Fibroblast" = "#E69F00",
  "Neutrophil" = "#56B4E9",
  "Cycling" = "#CC79A7",
  "B cells" = "#0072B2",
  "plasmablasts" = "#B07AA1",
  "Endothelial" = "#59A14F",
  "Mast cells" = "#8C564B",
  "Osteoclast" = "#B15928",
  "DCs" = "#6A3D9A",
  "Mural cells" = "#9C755F",
  "Plasma cell" = "#E15759"
)

macrophage_palette <- c(
  "FTL+ Mac" = "#3B5BA5",
  "CCL3+ Mac" = "#D55E00",
  "FN1+ Mac" = "#009E73",
  "CD1C+ Mac" = "#CC79A7",
  "CD3E+ Mac" = "#E69F00",
  "NAMPT+ Mac" = "#7E5A9B"
)

tcell_palette <- c(
  "CD4 memory/helper T" = "#3B5BA5",
  "Effector-memory T" = "#009E73",
  "Activated cytotoxic T" = "#D55E00",
  "NK-like cytotoxic T" = "#0072B2",
  "Treg" = "#CC79A7",
  "MAIT-like T" = "#E69F00"
)

lineage_palette_shared <- c(
  "Tumor" = "#D55E00",
  "Macrophage" = "#4E79A7",
  "T cell" = "#009E73",
  "T cells" = "#009E73"
)

nature_discrete <- function(n) {
  if (n <= length(nature_palette)) {
    return(nature_palette[seq_len(n)])
  }
  grDevices::colorRampPalette(nature_palette)(n)
}

nature_named_palette <- function(values) {
  values <- unique(as.character(values))
  setNames(nature_discrete(length(values)), values)
}

palette_for <- function(values, base_palette = nature_palette) {
  values <- unique(as.character(values))
  named <- base_palette[intersect(names(base_palette), values)]
  missing <- setdiff(values, names(named))
  if (length(missing) > 0) {
    named <- c(named, setNames(nature_discrete(length(missing)), missing))
  }
  named[values]
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
