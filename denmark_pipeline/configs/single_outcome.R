source("denmark_pipeline/configs/area.R", local = TRUE)
config$run_id <- "single_v3"
config$model$mode <- "single"
config$model$outcomes <- list(restore = list(class = NULL,
  scores = c(Score = "Results/independent_model/restoration_score_10m.tif"), allowed_mask = NULL))
config$model$metric_weights <- c(Score = 1)
config$model$metric_scales <- NULL
config$eligibility$forest_only <- FALSE
# The sole score is supplied by your separate scientific model; the pipeline
# does not derive it by maximizing or averaging the eight pathway scores.
