dir.create(
  path = Sys.getenv("R_LIBS_USER"),
  showWarnings = FALSE,
  recursive = TRUE
)

packages <- c(
  "tidyverse",
  "phangorn",
  "ape",
  "scales",
  "ggforce",
  "cowplot",
  "magrittr",
  "gridExtra",
  "patchwork",
  "plotly",
  "finalfit",
  "shiny"
)

# Install packages that are not already installed
install.packages(
  setdiff(packages, rownames(installed.packages())),
  lib = Sys.getenv("R_LIBS_USER"),
  repos = "https://cran.rstudio.com/"
)
