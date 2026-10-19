# ---- Spatial inputs ------------------------------------------------------
fine_landuse_path <- "Datasets/FinalLanduse.tif"
coarse_landuse_path <- file.path(
  "Predictors",
  "gudenaa_optimization",
  "Landuse_Gudenaa_200m_hierarchical.tif"
)
subareas_path <- file.path(
  "Datasets",
  "gudenaaplanen_afgreansning_15sep_delområder.shp"
)

# National current-nature patches preserve connectivity and full patch area
# beyond the new planning boundary. These are the same prepared inputs used by
# the earlier large-polygon anchored model.
national_patch_path <- file.path(
  "Predictors",
  "gudenaa_optimization",
  "Current_nature_patches_Denmark_10m.gpkg"
)
national_patch_layer <- "current_nature_patches"
national_patch_dictionary_path <- file.path(
  "Predictors",
  "gudenaa_optimization",
  "Current_nature_patch_dictionary_Denmark.csv"
)

# Geometry contact settings inherited from the earlier patch-supernode model.
contact_buffer_m <- 0.10
minimum_contact_proxy_m2 <- 0.10

subarea_label_field <- NULL
expected_number_of_subareas <- 4L
expected_subarea_areas_ha <- c(
  7991.14551899414,
  8653.25670292969,
  13009.8300014893,
  14526.539438208
)
expected_total_area_ha <- 44180.77
polygon_area_audit_tolerance_ha <- 5

# Needed to recreate exactly the biodiversity-supported eligibility mask used
# by notebook 07.
metrics_dir <- file.path(
  "Results",
  "gudenaa_plant_biodiversity",
  "full",
  "metrics"
)
biodiversity_support_metric <- "Richness"
biodiversity_support_scenario <- "current"

# ---- Upstream successful warm-up ----------------------------------------
warmup_dir <- file.path(
  "Optimization",
  "gudenaa_new_area_200m_30pct_existing_adjacency_thr1_forest_rule"
)

warmup_name <-
  "gudenaa_new_area_200m_30pct_existing_adjacency_thr1_forest_rule"

# Notebook 07 can legitimately reuse this matching baseline checkpoint
# without copying it into warmup_dir. Notebook 08 therefore resolves the
# biodiversity RDS from the warm-up folder first and this fallback second.
baseline_checkpoint_dir <- file.path(
  "Optimization",
  "gudenaa_new_area_200m_30pct"
)

# ---- Policy --------------------------------------------------------------
global_network_minimum_ha <- 5000
subarea_network_minimum_ha <- 1000
global_nature_target_fraction <- 0.30

# A narrow global upper tolerance prevents an unconstrained overshoot (for
# example, 50% nature) while retaining the successful notebook-07 footprint.
# Sixteen hectares equals the former 4 ha tolerance in each of four subareas.
global_restoration_tolerance_ha <- 16

# Keep the same area tolerance and spatial weights used by the warm-up.
spatial_contiguity_bonus <- 1.2
new_to_new_adjacency_weight <- 0.5
existing_to_new_adjacency_weight <- 0.5

# Strict lexicographic tolerances. Stage 3 may use only this much additional
# footprint-change area above the stage-2 incumbent; 0.01 ha is numerical
# slack, not permission for a materially larger repair.
network_shortfall_tolerance_ha <- 1e-6
repair_change_tolerance_ha <- 0.01

# Additional preference for joining several distinct existing-nature patches
# inside the certified global network. This is a normalized soft reward, not a
# hard requirement to connect a particular number of patches.
connected_existing_patch_bonus <- 0.25

# Established coarse rule: when at least half of a cell's eligible area is
# ProductionForest, the entire 200 m decision receives a forest destination.
production_forest_share_threshold <- 0.50

# ---- Outputs and solver --------------------------------------------------
problem_dir <- file.path(
  "Optimization",
  "gudenaa_new_area_200m_global30_feasible_network_repair"
)

problem_name <-
  "gudenaa_new_area_200m_global30_feasible_network_repair"

model_path <- "Gudenaa_global30_feasible_network_repair_200m.mod"

solver_name <- "gurobi"
repair_feasibility_time_limit_seconds <- 10800
repair_minimum_change_time_limit_seconds <- 7200
repair_quality_time_limit_seconds <- 3600
repair_mip_gap <- 0.005

solver_options_base <- paste(
  "tech:outlev=1",
  collapse = " "
)

reuse_matching_checkpoints <- TRUE
