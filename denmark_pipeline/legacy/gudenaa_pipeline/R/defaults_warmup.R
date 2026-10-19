# ---- Spatial inputs ------------------------------------------------------
fine_landuse_path <- "Datasets/FinalLanduse.tif"

# Reuse the established 200 m grid so cell alignment remains comparable with
# the earlier Gudenå run. The notebook crops and masks it to the new polygon.
coarse_landuse_path <- file.path(
  "Predictors",
  "gudenaa_optimization",
  "Landuse_Gudenaa_200m_hierarchical.tif"
)

subareas_path <- file.path(
  "Datasets",
  "gudenaaplanen_afgreansning_15sep_delområder.shp"
)

# NULL assigns stable IDs 1:4 in shapefile row order. To preserve meaningful
# names in output tables, replace NULL with an existing attribute name.
subarea_label_field <- NULL
expected_number_of_subareas <- 4L

# These values audit that the intended shapefile was loaded. The target itself
# is always calculated from the file, not hard-coded from these rounded values.
expected_subarea_areas_ha <- c(7991.14551899414, 8653.25670292969, 13009.8300014893, 14526.539438208)
expected_total_area_ha <- 44180.77
polygon_area_audit_tolerance_ha <- 5

# ---- Biodiversity inputs -------------------------------------------------
metrics_dir <- file.path(
  "Results",
  "gudenaa_plant_biodiversity",
  "full",
  "metrics"
)

metric_weights <- c(
  Richness = 1 / 3,
  PD = 1 / 3,
  RangeRarity = 1 / 3
)
biodiversity_scale <- 1000

# Candidate land must have complete support in the chosen biodiversity raster.
restrict_eligibility_to_biodiversity_support <- TRUE
biodiversity_support_metric <- "Richness"
biodiversity_support_scenario <- "current"

# ---- Policy and 200 m rules ---------------------------------------------
target_nature_proportion <- 0.30

# At most one nominal 200 m cell of overshoot is allowed in each subarea.
restoration_tolerance_ha <- 4

# Ignore coarse cells containing only a very small eligible sliver. The
# capacity audit stops if this makes a subarea target impossible.
minimum_eligible_area_ha <- 1

# Biodiversity and the combined spatial score are normalized. A value of 1
# gives the combined spatial score the same headline weight as biodiversity.
# A value of 10 is an intentionally extreme spatial stress test.
spatial_contiguity_bonus <- 1.2

# These partition the combined spatial score and must sum to one. The default
# gives equal influence to new-to-new and existing-to-new adjacency without
# silently doubling the spatial weight used by the completed baseline.
new_to_new_adjacency_weight <- 0.5
existing_to_new_adjacency_weight <- 0.5

# Any positive spatial weight activates all restoration actions. A positive
# new-to-new weight also activates adjacency binaries.
prune_dominated_landuse_actions <- TRUE

# Proposed phase-two network goals. They are diagnosed after this solve but
# are deliberately not enforced in this warm-up model.
proposed_global_network_ha <- 5000
proposed_subarea_network_ha <- 1000

# ---- Outputs and solver --------------------------------------------------
problem_dir <- file.path(
  "Optimization",
  "gudenaa_new_area_200m_30pct_existing_adjacency_thr1"
)
problem_name <- "gudenaa_new_area_200m_30pct_existing_adjacency_thr1"
model_path <- "Gudenaa_30pct_by_subarea_200m_existing_adjacency.mod"

# Matching expensive checkpoints from the completed baseline are read-only
# fallbacks. New or changed checkpoints are written only to problem_dir.
baseline_checkpoint_dir <- file.path(
  "Optimization",
  "gudenaa_new_area_200m_30pct"
)

solver_name <- "gurobi"
solver_options <- paste(
  "lim:time=21600",
  "mip:gap=0.005",
  "tech:outlev=1",
  paste0(
    "tech:logfile=",
    file.path(problem_dir, paste0(problem_name, "_gurobi.log"))
  ),
  collapse = " "
)

# Expensive tables are reused only when their input/configuration signature
# matches. Re-knit safely after interruption or after the AMPL run.
reuse_matching_checkpoints <- TRUE
