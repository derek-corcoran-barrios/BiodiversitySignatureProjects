# Projection-only targets pipeline: no occurrence downloads or model refits.
library(targets)
source("denmark_pipeline/run.R")
config_path <- Sys.getenv("DENMARK_CONFIG", "denmark_pipeline/configs/area.R")
initial <- dp_config(config_path)
tar_option_set(packages = c("terra", "SpeciesPoolR", "digest", "dplyr", "readr", "ape", "maxnet", "Matrix"),
  controller = crew::crew_controller_local(workers = initial$projection$workers),
  error = "stop", memory = "transient", storage = "main", retrieval = "main")
list(
  tar_target(config_file, config_path, format = "file"),
  tar_target(code_files, list.files("denmark_pipeline/R", full.names = TRUE), format = "file"),
  tar_target(adapter_file, "denmark_pipeline/adapters/speciespoolr_saved_models.R", format = "file"),
  tar_target(cfg, { code_files; adapter_file; dp_config(config_file) }),
  tar_target(input_files, {
    r <- dp_projection_registry(cfg)
    training <- unname(unlist(readRDS(cfg$projection$training_manifest)))
    climate_files <- if (is.null(cfg$projection$climate)) c(cfg$projection$climate_manifest,
      unname(unlist(readRDS(cfg$projection$climate_manifest)))[2]) else cfg$projection$climate
    dp_files(c(dp_region_inputs(cfg), cfg$projection$boundary, cfg$projection$training_manifest,
               training, climate_files, r$files))
  }, format = "file"),
  tar_target(plan, { input_files; dp_projection_plan(cfg) }),
  tar_target(tile_job, plan$jobs, iteration = "list"),
  tar_target(tile_metrics, dp_project_tile(plan, tile_job), pattern = map(tile_job), format = "file"),
  tar_target(metric_mosaics, dp_mosaic_metrics(plan, tile_metrics), format = "file")
)
