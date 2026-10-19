source("denmark_pipeline/configs/area.R", local = TRUE)
config$projection$id <- "denmark_pathways_v3"
config$projection$boundary <- "Datasets/Denmark.gpkg"
config$projection$workers <- 8L
# run_projection() only projects saved models; it does not optimize Denmark.
# Keep the complete national occurrences, thresholds, buffers and range weights.
