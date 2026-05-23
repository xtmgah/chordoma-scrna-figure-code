options(stringsAsFactors = FALSE)

root_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), mustWork = TRUE))
source(file.path(root_dir, "_figure_style.R"))

scripts <- c(
  file.path(root_dir, "Figure1", "Figure1_Nature_Panels.R"),
  file.path(root_dir, "Figure2", "Figure2_Nature_Panels.R"),
  file.path(root_dir, "Figure3", "Figure3_Nature_Panels.R"),
  file.path(root_dir, "Figure4", "Figure4_Nature_Panels.R"),
  file.path(root_dir, "Figure5", "Figure5_Nature_Panels.R")
)

for (script in scripts) {
  message("\n=== Running ", basename(script), " ===")
  tryCatch(
    {
      env <- new.env(parent = globalenv())
      env$SCRIPT_DIR <- dirname(script)
      source(script, local = env)
      export_large_panel_pngs(file.path(dirname(script), "nature_panels"))
    },
    error = function(e) {
      message("FAILED: ", script)
      message(conditionMessage(e))
    }
  )
}
