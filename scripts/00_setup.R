args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = TRUE)) else "scripts"
root_dir <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = TRUE)
setwd(root_dir)

cran_packages <- c(
  "data.table",
  "yaml",
  "readr",
  "jsonlite",
  "ggplot2",
  "rmarkdown",
  "knitr",
  "coloc",
  "remotes"
)

install_if_missing <- function(pkgs) {
  installed <- rownames(installed.packages())
  missing <- pkgs[!pkgs %in% installed]
  if (length(missing) == 0) {
    return(invisible(TRUE))
  }
  install.packages(missing, repos = "https://cloud.r-project.org")
}

install_if_missing(cran_packages)

if (!("TwoSampleMR" %in% rownames(installed.packages()))) {
  message("Installing TwoSampleMR from MRCIEU r-universe.")
  install.packages(
    "TwoSampleMR",
    repos = c("https://mrcieu.r-universe.dev", "https://cloud.r-project.org")
  )
}

if (!("ieugwasr" %in% rownames(installed.packages()))) {
  message("Installing optional ieugwasr package for OpenGWAS validation.")
  try(
    install.packages(
      "ieugwasr",
      repos = c("https://mrcieu.r-universe.dev", "https://cloud.r-project.org")
    ),
    silent = TRUE
  )
}

attach_packages <- c("data.table", "yaml", "TwoSampleMR")
invisible(lapply(attach_packages, function(pkg) {
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}))

if (!("coloc" %in% rownames(installed.packages()))) {
  stop("Package coloc was not installed successfully.", call. = FALSE)
}

source(file.path("scripts", "utils.R"))
cfg <- read_config()
invisible(lapply(unlist(cfg$project, use.names = FALSE), function(path) {
  if (is.character(path) && grepl("/", path, fixed = TRUE)) {
    ensure_dir(path)
  }
}))

status_message("Setup complete. Installed and loaded required packages where available.")
sessionInfo()
