# CLI: Rscript gudenaa_pipeline/run_pipeline.R
# RStudio: source("gudenaa_pipeline/run_pipeline.R")
#          run_gudenaa_pipeline("gudenaa_pipeline/config.R")
gp_entry_file <- function() {
  frames <- sys.frames()
  for (i in rev(seq_along(frames))) {
    f <- frames[[i]]$ofile
    if (!is.null(f)) return(normalizePath(f, winslash = "/", mustWork = TRUE))
  }
  arg <- grep("^--file=", commandArgs(), value = TRUE)
  if (length(arg)) return(normalizePath(sub("^--file=", "", arg[1]), winslash = "/"))
  stop("Source run_pipeline.R or use Rscript.")
}
gp_package_dir <- dirname(gp_entry_file())
source(file.path(gp_package_dir, "R", "pipeline.R"))
source(file.path(gp_package_dir, "R", "network_exports.R"))
source(file.path(gp_package_dir, "R", "compare_footprints.R"))
if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  run_gudenaa_pipeline(if (length(args)) args[1] else file.path(gp_package_dir, "config.R"))
}
