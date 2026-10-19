# Source denmark_pipeline/run.R first. Copy this config for each area/model.
config <- denmark_config()
config$area_id <- "gudenaa"
config$run_id <- "pathways_v3"
config$region$boundary <- "Datasets/gudenaaplanen_afgreansning_15sep_delområder.shp"
config$region$subareas <- config$region$boundary
config$connectivity$subarea_min_ha <- 1000
config$projection$id <- "gudenaa_pathways_v3"
config$model$outcomes <- dp_eight_outcomes("Results/denmark_projection/gudenaa_pathways_v3/metrics", "vrt")

# Warm-up and final spatial preferences can be controlled independently.
config$objective$warmup_contiguity <- 1.2
config$objective$quality_contiguity <- 1.2 # Set to 0 for a biodiversity-only final objective.
config$repair <- list(enabled = TRUE, seconds = 1800, gap = 0.005)
config$quality <- list(footprint = "free", max_change_ha = NULL)
# Alternatives: list(footprint = "capped", max_change_ha = 100), or
# list(footprint = "fixed", max_change_ha = NULL).

# Choose any combination of sources. The default uses all current nature.
config$anchors$sources <- list(nature = list(type = "nature"))
# config$anchors$sources$political <- list(type = "vector", path = "Datasets/political_anchors.gpkg")
# Lakes anchor the adjoining land by a one-pixel shoreline buffer. This default
# does not make open water a terrestrial corridor or credit it as restored nature.
# config$anchors$sources$lakes <- list(type = "vector", path = "Datasets/lakes.gpkg", buffer_m = 10)
# To use a hand-built 0/1 10 m raster as the ONLY anchor definition:
# config$anchors$sources <- list(custom = list(type = "raster", path = "Datasets/anchors_10m.tif", values = 1))

# NULL subareas is fully supported:
# config$region$boundary <- "Datasets/another_area.gpkg"
# config$region$subareas <- NULL
# config$connectivity$subarea_min_ha <- NULL

# Budget alternatives (values here illustrate syntax, not a recommendation):
# config$budget$target_nature_fraction <- 0.25
# Or an absolute restoration area interval:
# config$budget$target_nature_fraction <- NULL
# config$budget$min_restoration_ha <- 2000
# config$budget$max_restoration_ha <- 2016
# Or a monetary cap, with optional area/connectivity requirements:
# config$budget$target_nature_fraction <- NULL
# config$budget$cost_per_ha <- "Datasets/restoration_DKK_per_ha.tif"
# config$budget$cost_cap <- 50000000

# Context is an explicit modelling choice. Default connectivity is within the
# planning polygon. Supply a wider polygon and choose within_context to credit
# existing nature there. A clipped context never claims full national patches.
# config$connectivity$context_boundary <- "Datasets/wider_connectivity_context.gpkg"
# config$connectivity$global_area_scope <- "within_context"

# Reuse the original climate source recorded in Predictors/current_predictor_paths.rds.
# Set config$projection$climate to its actual TIFF if that manifest has moved.
config$projection$workers <- 8L
config$projection$tile_cells <- 1000L  # 10 x 10 km at 10 m resolution
config$projection$memmax_gb_per_worker <- 4
