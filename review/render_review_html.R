options(stringsAsFactors = FALSE)

find_script_path <- function() {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg)) {
    return(normalizePath(sub("^--file=", "", file_arg[1]), mustWork = TRUE))
  }

  frame_paths <- vapply(
    sys.frames(),
    function(frame) {
      if (!is.null(frame$ofile)) {
        return(frame$ofile)
      }
      NA_character_
    },
    character(1)
  )
  frame_paths <- frame_paths[!is.na(frame_paths)]
  if (length(frame_paths)) {
    return(normalizePath(frame_paths[length(frame_paths)], mustWork = TRUE))
  }

  stop("Could not determine the path to render_review_html.R. Run with Rscript or source() from R.")
}

script_path <- find_script_path()
review_dir <- dirname(script_path)
repo_root <- normalizePath(file.path(review_dir, ".."), mustWork = TRUE)

panel_source_root <- Sys.getenv(
  "CHORDOMA_PANEL_SOURCE_ROOT",
  unset = file.path(dirname(repo_root), "Chrodoma")
)
panel_source_root <- normalizePath(panel_source_root, mustWork = FALSE)

preview_root <- file.path(review_dir, "panel_previews")
dir.create(preview_root, recursive = TRUE, showWarnings = FALSE)

figures <- paste0("Figure", 1:5)
thumbnail_size <- Sys.getenv("CHORDOMA_REVIEW_THUMBNAIL_SIZE", unset = "1400")

convert_pdf_with_qlmanage <- function(pdf, out_dir, size = thumbnail_size) {
  tmp_dir <- tempfile("ql_preview_")
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  status <- system2(
    "qlmanage",
    args = c("-t", "-s", size, "-o", tmp_dir, pdf),
    stdout = TRUE,
    stderr = TRUE
  )
  invisible(status)

  produced <- file.path(tmp_dir, paste0(basename(pdf), ".png"))
  if (!file.exists(produced)) {
    stop("Quick Look did not produce a preview for: ", pdf)
  }

  dest <- file.path(out_dir, paste0(tools::file_path_sans_ext(basename(pdf)), ".png"))
  file.copy(produced, dest, overwrite = TRUE)
  dest
}

convert_pdf_with_sips <- function(pdf, out_dir) {
  dest <- file.path(out_dir, paste0(tools::file_path_sans_ext(basename(pdf)), ".png"))
  status <- system2(
    "sips",
    args = c("-s", "format", "png", pdf, "--out", dest),
    stdout = TRUE,
    stderr = TRUE
  )
  invisible(status)
  if (!file.exists(dest)) {
    stop("sips did not produce a preview for: ", pdf)
  }
  dest
}

make_preview <- function(pdf, out_dir) {
  dest <- file.path(out_dir, paste0(tools::file_path_sans_ext(basename(pdf)), ".png"))
  if (file.exists(dest) && file.info(dest)$mtime >= file.info(pdf)$mtime) {
    return(dest)
  }

  if (nzchar(Sys.which("qlmanage"))) {
    return(convert_pdf_with_qlmanage(pdf, out_dir))
  }
  if (nzchar(Sys.which("sips"))) {
    return(convert_pdf_with_sips(pdf, out_dir))
  }
  stop("No PDF preview tool found. Install pdftools, or render this report on macOS with qlmanage or sips available.")
}

message("Panel source root: ", panel_source_root)
for (figure in figures) {
  source_dir <- file.path(panel_source_root, figure, "nature_panels")
  out_dir <- file.path(preview_root, figure)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  pdfs <- sort(list.files(source_dir, pattern = "\\.pdf$", full.names = TRUE))
  if (!length(pdfs)) {
    warning("No PDF panels found for ", figure, " in ", source_dir)
    next
  }

  message("Creating previews for ", figure, ": ", length(pdfs), " panels")
  for (pdf in pdfs) {
    make_preview(pdf, out_dir)
  }
}

if (!requireNamespace("rmarkdown", quietly = TRUE)) {
  stop("The rmarkdown package is required to render the review HTML.")
}

configure_pandoc <- function() {
  if (rmarkdown::pandoc_available()) {
    return(invisible(TRUE))
  }

  arch <- Sys.info()[["machine"]]
  candidate_dirs <- c(
    if (identical(arch, "arm64")) {
      "/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/aarch64"
    },
    if (!identical(arch, "arm64")) {
      "/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/x86_64"
    },
    "/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/aarch64",
    "/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/x86_64"
  )
  candidate_dirs <- unique(candidate_dirs[!is.na(candidate_dirs)])
  candidate_dirs <- candidate_dirs[file.exists(file.path(candidate_dirs, "pandoc"))]

  if (length(candidate_dirs)) {
    Sys.setenv(RSTUDIO_PANDOC = candidate_dirs[1])
    return(invisible(rmarkdown::pandoc_available()))
  }

  invisible(FALSE)
}

if (!configure_pandoc()) {
  stop("Pandoc was not found. Install Pandoc or open this project in RStudio/Quarto before rendering.")
}

message("Rendering HTML review report")
rmarkdown::render(
  input = file.path(review_dir, "chordoma_figure_reproducibility.Rmd"),
  output_file = "chordoma_figure_reproducibility.html",
  output_dir = review_dir,
  clean = TRUE,
  quiet = FALSE
)

message("Wrote: ", file.path(review_dir, "chordoma_figure_reproducibility.html"))
