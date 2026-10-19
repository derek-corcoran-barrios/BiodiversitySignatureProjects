# Run from the BiodiversitySignatureProjects root. Paths are relative to it.
# Change run_id for changes to data or scientific settings.
config <- list(
  project_root = ".",
  run_id = "global30_eight_habitats_v1",
  output_root = "Optimization/gudenaa_pipeline",
  ampl_executable = "ampl", # Or the full path to ampl.exe / ampl.
  solver = "gurobi",
  solver_options = "tech:outlev=1",
  mip_gap = 0.005,
  fine_landuse_path = "Datasets/FinalLanduse.tif",
  coarse_landuse_path = "Predictors/gudenaa_optimization/Landuse_Gudenaa_200m_hierarchical.tif",
  subareas_path = "Datasets/gudenaaplanen_afgreansning_15sep_delområder.shp",
  national_patch_path = "Predictors/gudenaa_optimization/Current_nature_patches_Denmark_10m.gpkg",
  national_patch_layer = "current_nature_patches",
  national_patch_dictionary_path = "Predictors/gudenaa_optimization/Current_nature_patch_dictionary_Denmark.csv",
  metrics_dir = "Results/gudenaa_plant_biodiversity/full/metrics",
  baseline_checkpoint_dir = "Optimization/gudenaa_new_area_200m_30pct",
  global_nature_target_fraction = 0.30,
  global_restoration_tolerance_ha = 16,
  global_network_minimum_ha = 5000,
  subarea_network_minimum_ha = 1000,
  minimum_eligible_area_ha = 1,
  production_forest_share_threshold = 0.50,
  spatial_contiguity_bonus = 1.2,
  new_to_new_adjacency_weight = 0.5,
  existing_to_new_adjacency_weight = 0.5,
  connected_existing_patch_bonus = 0.25,
  metric_weights = c(Richness = 1/3, PD = 1/3, RangeRarity = 1/3),
  biodiversity_scale = 1000,
  contact_buffer_m = 0.10,
  minimum_contact_proxy_m2 = 0.10,
  subarea_label_field = NULL,
  # Reference checks copied from the supplied notebooks; targets use geometry.
  expected_subarea_areas_ha = c(7991.14551899414, 8653.25670292969,
                              13009.8300014893, 14526.539438208),
  expected_total_area_ha = 44180.77,
  polygon_area_audit_tolerance_ha = 5,
  # "minimum_change" reproduces the existing repair priority.
  # "free" uses the repaired solution as a start without a footprint-change cap.
  quality_mode = "minimum_change",
  repair_change_tolerance_ha = 0.01,
  network_shortfall_tolerance_ha = 1e-6,
  warmup_seconds = 21600,
  feasibility_seconds = 10800,
  feasibility_continuation_seconds = c(10800),
  repair_verification_seconds = 600,
  minimum_change_seconds = 7200,
  quality_seconds = c(28800, 14400), # Separate checkpoints after 8h and 4h.
  # TRUE starts AMPL; FALSE prepares inputs/run files and pauses before a solve.
  execute_solver = TRUE,
  make_maps = TRUE,
  sld_path = "gudenaa_incumbent_maps/NaturetypeStyle.sld",
  raster_maxcell = 1500000,
  png_dpi = 600
)
