# Source from the project root; nothing is solved merely by sourcing this file.
.dp_file <- function() {
  for (e in rev(sys.frames())) if (!is.null(e$ofile)) return(normalizePath(e$ofile))
  stop("Use source('denmark_pipeline/run.R').")
}
dp_home <- dirname(.dp_file())
# terra supplies raster-aware operators such as %in%; attach its S4 generics.
suppressPackageStartupMessages(library(terra))
for (.dp_mod in c("config", "spatial", "graph", "solver", "projection"))
  sys.source(file.path(dp_home, "R", paste0(.dp_mod, ".R")), envir = globalenv())
rm(.dp_file, .dp_mod)

launch_config_app <- function(project_root = getwd(), launch.browser = TRUE, ...) {
  for (pkg in c("shiny", "jsonlite")) if (!requireNamespace(pkg, quietly = TRUE))
    stop("Install the app dependency first: install.packages('", pkg, "')")
  root <- normalizePath(project_root, winslash = "/", mustWork = TRUE)
  old <- options(denmark.project_root = root)
  on.exit(options(old), add = TRUE)
  shiny::runApp(file.path(dp_home, "shiny"), launch.browser = launch.browser, ...)
}
