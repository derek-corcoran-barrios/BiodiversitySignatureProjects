# Separate process: input is produced by the local app, never an uploaded R script.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Expected one preview request file.")
request <- readRDS(args[1])
result <- tryCatch({
  source(file.path(request$package, "run.R"))
  source(file.path(request$package, "shiny", "helpers.R"))
  list(ok = TRUE, result = dp_app_make_preview(request$config, request$folder))
}, error = function(e) list(ok = FALSE, message = conditionMessage(e)))
path <- file.path(request$folder, "result.rds")
saveRDS(result, paste0(path, ".building"))
if (!file.rename(paste0(path, ".building"), path)) stop("Cannot publish preview result.")
