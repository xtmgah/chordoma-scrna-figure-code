# Chordoma scRNA-seq Figure Code

This repository contains the R code used to generate the revised publication-style figure panels for a chordoma single-cell RNA-seq study.

Manuscript information is provisional and can be updated before release. The intended use is to provide a stable code link for the manuscript data/code availability section.

## Repository Contents

- `_figure_style.R`: shared Nature-style plotting theme, palettes, PDF device, and font setup.
- `run_all_nature_figures.R`: runs the cleaned Figure 1-5 plotting scripts in order.
- `Figure1/Figure1_Nature_Panels.R`: final Figure 1 panel-generation code.
- `Figure2/Figure2_Nature_Panels.R`: final Figure 2 panel-generation code.
- `Figure3/Figure3_Nature_Panels.R`: final Figure 3 panel-generation code.
- `Figure4/Figure4_Nature_Panels.R`: final Figure 4 panel-generation code.
- `Figure5/Figure5_Nature_Panels.R`: final Figure 5 panel-generation code.
- `Figure*/Figure*.R` and helper scripts: earlier/original analysis plotting scripts retained for transparency.
- `data_manifest.tsv`: required and optional input files expected by the final plotting scripts.
- `check_required_inputs.R`: quick check for required local input files.

Large input data objects and generated figures are intentionally not included in this repository.

## Data and Input Files

The final plotting scripts expect local `.rds` and `.csv` objects in the same relative folder structure shown in `data_manifest.tsv`. These files include processed Seurat objects, CellChat objects, AUC/modeling outputs, and pseudotime/network inputs.

Before running the scripts, place the required input files into the listed locations, then run:

```r
source("check_required_inputs.R")
```

If any required file is missing, the checker prints the missing path and exits with a non-zero status.

## R Environment

The code was developed in R with these main packages:

- Seurat
- scCustomize
- CellChat
- ggplot2
- dplyr
- tidyr
- patchwork
- ggpubr
- rstatix
- ComplexHeatmap
- circlize
- ggraph
- igraph
- ggrepel
- survminer
- survival
- fmsb
- ggsci
- viridis
- showtext
- sysfonts
- scales
- RColorBrewer
- pheatmap

Install missing CRAN packages with `install.packages()`. `ComplexHeatmap` may be installed through Bioconductor, and CellChat installation may depend on the version used in the analysis environment.

## Font Setup

The figures use Roboto Condensed as the primary font. If the font files are not installed at the default local path, set these environment variables before running:

```r
Sys.setenv(
  ROBOTO_CONDENSED_REGULAR = "/path/to/RobotoCondensed-Regular.ttf",
  ROBOTO_CONDENSED_BOLD = "/path/to/RobotoCondensed-Bold.ttf",
  ROBOTO_CONDENSED_ITALIC = "/path/to/RobotoCondensed-Italic.ttf",
  ROBOTO_CONDENSED_BOLDITALIC = "/path/to/RobotoCondensed-BlackItalic.ttf"
)
```

If the font files are not available, the scripts warn and continue with the active PDF device fallback.

## Reproducing the Figure Panels

After placing the required input files, run all final panels with:

```r
source("run_all_nature_figures.R")
```

Or run a single figure:

```r
source("Figure1/Figure1_Nature_Panels.R")
source("Figure2/Figure2_Nature_Panels.R")
source("Figure3/Figure3_Nature_Panels.R")
source("Figure4/Figure4_Nature_Panels.R")
source("Figure5/Figure5_Nature_Panels.R")
```

Each final script writes PDF panels to:

```text
Figure*/nature_panels/
```

Generated figure files are ignored by Git and should not be committed.

## Suggested Manuscript Code Availability Text

The R scripts used to generate the single-cell RNA-seq figure panels are available at:

```text
https://github.com/xtmgah/chordoma-scrna-figure-code
```

Processed input objects required to reproduce the panels are described in `data_manifest.tsv` and are available through the study data-access mechanism described in the manuscript.

## Notes

Some legacy scripts retain absolute paths from the original analysis workstation. The recommended reproducible entry points are the `*_Nature_Panels.R` scripts and `run_all_nature_figures.R`, which use paths relative to this repository.
