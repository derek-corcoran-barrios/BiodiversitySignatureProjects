# Pure configuration functions, shared by the app and its tests.
dp_app_defaults <- function(root = getwd()) list(
  project_root = root, area_id = "new_area", run_id = "v3_01", output_root = "Optimization/denmark",
  boundary = "Datasets/area.gpkg", landuse = "Predictors/gudenaa_landuse/Landuse_Denmark_current.tif",
  planning_m = 200, min_unit_ha = 1, directions = "4", origin_x = 0, origin_y = 0,
  subareas_on = FALSE, subareas = "", subarea_id = "", subarea_goal_on = TRUE,
  subarea_min = 1000, subarea_goals = "", context = "", area_scope = "within_region",
  agriculture = TRUE, production_forest = TRUE, eligibility_mask = "", forest_only = TRUE, forest_share = 0.5,
  nature_on = TRUE, political_on = FALSE, lakes_on = FALSE, custom_on = FALSE,
  political_path = "", political_buffer = 0, political_existing = FALSE, political_budget = FALSE, political_network = FALSE,
  lakes_path = "", lakes_buffer = 10, lakes_existing = FALSE, lakes_budget = FALSE, lakes_network = FALSE,
  custom_path = "", custom_buffer = 0, custom_values = "1", custom_existing = FALSE, custom_budget = FALSE, custom_network = FALSE,
  anchor_combine = "union", extra_anchors = "{}",
  mode = "pathways", metrics_dir = "Results/new_area/metrics", extension = "vrt", single_score = "",
  richness_weight = 1, pd_weight = 1, rarity_weight = 1,
  score_type = "endpoint", baseline_dir = "Results/new_area/metrics", single_baseline = "",
  outcome_overrides = "{}", metric_scales = "",
  target_on = TRUE, target_percent = 30, tolerance_on = TRUE, tolerance_ha = 16,
  min_ha = 0, max_on = FALSE, max_ha = 5000, cost_on = FALSE, cost_cap = 50000000, cost_raster = "",
  global_on = TRUE, global_min = 5000, require_new = TRUE,
  warmup_weight = 1.2, quality_weight = 1.2, new_weight = 0.5, anchor_weight = 0.5,
  repair_on = TRUE, repair_minutes = 30, footprint = "free", change_cap = 100,
  solver = "gurobi", ampl = "ampl", solver_options = "tech:outlev=1", gap_percent = 0.5,
  warmup_minutes = 60, feasibility_minutes = "180, 180", quality_minutes = "480, 240",
  projection_on = FALSE, projection_id = "new_area", projection_boundary = "",
  projection_output = "Results/denmark_projection", saved_fit_root = "Results/gudenaa_plant_biodiversity/full",
  training_manifest = "Predictors/gudenaa_scenarios/training_current_predictor_paths_Denmark.rds",
  climate = "", climate_manifest = "Predictors/current_predictor_paths.rds",
  tile_cells = 1000, workers = 4, mem_gb = 4, blocks = 512, chunk = 20, pd_cells = 10000, keep_species = FALSE
)
dp_app_optional <- function(x) if (is.null(x) || !nzchar(trimws(x))) NULL else trimws(x)
dp_app_numbers <- function(x, label) {
  z <- suppressWarnings(as.numeric(trimws(strsplit(x, ",", fixed = TRUE)[[1]])))
  if (!length(z) || any(!is.finite(z))) stop(label, ": supply comma-separated numbers.")
  z
}
dp_app_object <- function(x, label) {
  if (!nzchar(trimws(x))) return(list())
  z <- tryCatch(jsonlite::fromJSON(x, simplifyVector = TRUE), error = function(e) stop(label, ": invalid JSON."))
  if (!is.list(z) || (length(z) && (is.null(names(z)) || any(!nzchar(names(z)))))) stop(label, ": use a named JSON object.")
  z
}
dp_app_named_numbers <- function(x, label) {
  z <- dp_app_object(x, label)
  if (!length(z) || any(vapply(z, length, integer(1)) != 1L) || !all(vapply(z, is.numeric, logical(1))))
    stop(label, ": use a JSON object mapping names to single numbers.")
  setNames(as.numeric(unlist(z, use.names = FALSE)), names(z))
}
dp_app_config <- function(state) {
  s <- utils::modifyList(dp_app_defaults(), state)
  c <- denmark_config()
  c$project_root <- s$project_root; c$area_id <- s$area_id; c$run_id <- s$run_id; c$output_root <- s$output_root
  c$region$boundary <- s$boundary; c$grid$landuse <- s$landuse
  c$grid$planning_unit_m <- s$planning_m; c$grid$directions <- as.integer(s$directions)
  c$grid$origin <- c(s$origin_x, s$origin_y); c$eligibility$min_unit_ha <- s$min_unit_ha
  c$eligibility$classes <- c(if (isTRUE(s$agriculture)) "Agriculture", if (isTRUE(s$production_forest)) "ProductionForest")
  if (!length(c$eligibility$classes)) stop("Choose at least one selectable land-use class.")
  c$eligibility$mask <- dp_app_optional(s$eligibility_mask)
  c$eligibility$forest_only <- s$mode == "pathways" && isTRUE(s$forest_only)
  c$eligibility$forest_share <- s$forest_share
  if (isTRUE(s$subareas_on)) {
    c$region$subareas <- dp_app_optional(s$subareas)
    if (is.null(c$region$subareas)) stop("Supply the subarea file or switch subareas off.")
    c$region$subarea_id <- dp_app_optional(s$subarea_id)
    if (isTRUE(s$subarea_goal_on)) c$connectivity$subarea_min_ha <-
      if (nzchar(trimws(s$subarea_goals))) dp_app_named_numbers(s$subarea_goals, "Subarea goals") else s$subarea_min
  }
  c$connectivity$context_boundary <- dp_app_optional(s$context)
  c$connectivity$global_area_scope <- s$area_scope
  c$connectivity$global_min_ha <- if (isTRUE(s$global_on)) s$global_min else NULL
  c$connectivity$require_new_restoration <- s$require_new
  sources <- if (isTRUE(s$nature_on)) list(nature = list(type = "nature")) else list()
  for (nm in c("political", "lakes", "custom")) if (isTRUE(s[[paste0(nm, "_on")]])) {
    sources[[nm]] <- list(type = if (nm == "custom") "raster" else "vector",
      path = s[[paste0(nm, "_path")]], buffer_m = s[[paste0(nm, "_buffer")]],
      include_as_existing = s[[paste0(nm, "_existing")]], credit_budget = s[[paste0(nm, "_budget")]],
      credit_network = s[[paste0(nm, "_network")]])
    if (nm == "custom" && nzchar(trimws(s$custom_values))) sources[[nm]]$values <- dp_app_numbers(s$custom_values, "Custom anchor values")
  }
  extras <- dp_app_object(s$extra_anchors, "Additional anchors")
  if (length(intersect(names(sources), names(extras)))) stop("Additional anchor names duplicate an enabled source.")
  for (nm in names(extras)) {
    if (!is.null(extras[[nm]]$values)) extras[[nm]]$values <- unlist(extras[[nm]]$values)
    if (!is.null(extras[[nm]]$filter_values)) extras[[nm]]$filter_values <- unlist(extras[[nm]]$filter_values)
  }
  c$anchors$sources <- c(sources, extras); c$anchors$combine <- s$anchor_combine
  c$model$mode <- s$mode; c$model$score_type <- s$score_type
  if (s$mode == "single") {
    if (!nzchar(trimws(s$single_score))) stop("Supply the independent model's score raster.")
    c$model$outcomes <- list(restore = list(class = NULL, scores = c(Score = s$single_score), allowed_mask = NULL))
    c$model$metric_weights <- c(Score = 1)
    if (s$score_type == "gain") c$model$current_scores <- c(Score = s$single_baseline)
  } else {
    if (isTRUE(s$projection_on)) c$model$outcomes <- dp_eight_outcomes(file.path(s$projection_output, s$projection_id, "metrics"), "vrt")
    else c$model$outcomes <- dp_eight_outcomes(s$metrics_dir, s$extension)
    w <- c(Richness = s$richness_weight, PD = s$pd_weight, RangeRarity = s$rarity_weight)
    if (any(!is.finite(w) | w < 0) || sum(w) <= 0) stop("Metric weights must be nonnegative, with at least one positive value.")
    c$model$metric_weights <- w/sum(w)
    if (s$score_type == "gain") {
      bdir <- if (isTRUE(s$projection_on)) file.path(s$projection_output, s$projection_id, "metrics") else s$baseline_dir
      ext <- if (isTRUE(s$projection_on)) "vrt" else s$extension
      c$model$current_scores <- setNames(file.path(bdir, paste0(names(w), "_current.", ext)), names(w))
    }
  }
  overrides <- dp_app_object(s$outcome_overrides, "Outcome overrides")
  if (length(setdiff(names(overrides), names(c$model$outcomes)))) stop("Outcome overrides contain unknown outcome names.")
  for (nm in names(overrides)) {
    if (length(setdiff(names(overrides[[nm]]), c("scores", "allowed_mask")))) stop("Override scores or allowed_mask only.")
    spec <- c$model$outcomes[[nm]]; override <- overrides[[nm]]
    if (!is.null(override$scores)) {
      vals <- unlist(override$scores, use.names = TRUE)
      if (any(!names(vals) %in% names(spec$scores)) || anyDuplicated(names(vals))) stop("Unknown overridden score metric.")
      spec$scores[names(vals)] <- vals
    }
    if (!is.null(override$allowed_mask)) spec$allowed_mask <- override$allowed_mask
    c$model$outcomes[[nm]] <- spec
  }
  if (nzchar(trimws(s$metric_scales))) c$model$metric_scales <- dp_app_named_numbers(s$metric_scales, "Metric scales")
  c$budget$target_nature_fraction <- if (isTRUE(s$target_on)) s$target_percent/100 else NULL
  c$budget$target_tolerance_ha <- if (isTRUE(s$tolerance_on)) s$tolerance_ha else NULL
  c$budget$min_restoration_ha <- s$min_ha; c$budget$max_restoration_ha <- if (isTRUE(s$max_on)) s$max_ha else NULL
  c$budget$cost_cap <- if (isTRUE(s$cost_on)) s$cost_cap else NULL
  c$budget$cost_per_ha <- if (isTRUE(s$cost_on)) dp_app_optional(s$cost_raster) else NULL
  c$objective$warmup_contiguity <- s$warmup_weight; c$objective$quality_contiguity <- s$quality_weight
  c$objective$new_to_new <- s$new_weight; c$objective$to_anchors <- s$anchor_weight
  c$repair$enabled <- s$repair_on; c$repair$seconds <- s$repair_minutes * 60
  c$repair$gap <- s$gap_percent/100
  c$quality$footprint <- s$footprint; c$quality$max_change_ha <- if (s$footprint == "capped") s$change_cap else NULL
  c$solver$ampl <- s$ampl; c$solver$name <- s$solver; c$solver$options <- s$solver_options; c$solver$gap <- s$gap_percent/100
  c$solver$warmup_seconds <- s$warmup_minutes * 60
  c$solver$feasibility_seconds <- dp_app_numbers(s$feasibility_minutes, "Feasibility minutes") * 60
  c$solver$quality_seconds <- dp_app_numbers(s$quality_minutes, "Quality minutes") * 60
  c$projection$id <- s$projection_id; c$projection$boundary <- dp_app_optional(s$projection_boundary)
  c$projection$output_root <- s$projection_output; c$projection$saved_fit_root <- s$saved_fit_root
  c$projection$training_manifest <- s$training_manifest; c$projection$climate <- dp_app_optional(s$climate)
  c$projection$climate_manifest <- dp_app_optional(s$climate_manifest)
  for (pair in list(c("tile_cells", "tile_cells"), c("workers", "workers"), c("memmax_gb_per_worker", "mem_gb"),
                    c("blocks_in_memory", "blocks"), c("metric_chunk_size", "chunk"), c("pd_max_cells", "pd_cells")))
    c$projection[[pair[1]]] <- s[[pair[2]]]
  c$projection$keep_species_rasters <- s$keep_species
  if (isTRUE(s$projection_on) && s$mode == "pathways") {
    if (!grepl("^[A-Za-z][A-Za-z0-9_]*$", c$projection$id)) stop("Projection ID must start with a letter and use letters, numbers or underscores.")
    for (nm in c("tile_cells", "workers", "blocks_in_memory", "metric_chunk_size", "pd_max_cells"))
      if (length(c$projection[[nm]]) != 1 || !is.finite(c$projection[[nm]]) || c$projection[[nm]] < 1 || c$projection[[nm]] %% 1 != 0)
        stop("Projection ", nm, " must be a positive integer.")
    if (!is.finite(s$mem_gb) || s$mem_gb <= 0) stop("Memory per worker must be positive.")
  }
  dp_config(c)
}
dp_app_export <- function(cfg) c(
  "# Denmark restoration pipeline v3. Generated by the configuration app.",
  "# Source denmark_pipeline/run.R before using this configuration.",
  "# File paths are interpreted relative to project_root. Review before running.",
  paste0("config <- ", paste(capture.output(dput(cfg,
    control = c("keepNA", "keepInteger", "niceNames", "showAttributes", "digits17"))), collapse = "\n"))
)
dp_app_summary <- function(cfg) {
  b <- cfg$budget; n <- cfg$connectivity
  fmt <- function(x) format(round(x, 3), scientific = FALSE, trim = TRUE, big.mark = ",")
  c(if (is.null(b$target_nature_fraction)) "No whole-polygon percentage target." else
      sprintf("Nature target: %.1f%% of the whole planning polygon, including credited current nature.", b$target_nature_fraction*100),
    if (!is.null(b$target_nature_fraction) && !is.null(b$target_tolerance_ha))
      sprintf("Restoration is capped at the area needed for that target plus %s ha.", fmt(b$target_tolerance_ha)),
    sprintf("New restoration: at least %s ha%s.", fmt(b$min_restoration_ha),
      if (is.null(b$max_restoration_ha)) "" else sprintf("; explicit maximum %s ha", fmt(b$max_restoration_ha))),
    if (!is.null(b$cost_cap)) sprintf("Restoration cost must be at most %s in the cost layer's currency.", format(b$cost_cap, big.mark = ",", scientific = FALSE)),
    if (is.null(n$global_min_ha) || n$global_min_ha == 0) "No Global network target." else sprintf("At least one anchored Global network of %s ha.", fmt(n$global_min_ha)),
    if (is.null(cfg$region$subareas)) "No subareas." else if (is.null(n$subarea_min_ha)) "Subareas are displayed; no subarea network targets." else
      if (length(n$subarea_min_ha) == 1 && is.null(names(n$subarea_min_ha))) sprintf("Each subarea needs an anchored network of %s ha; no per-subarea percentage target.", fmt(n$subarea_min_ha)) else
        "Subareas have individually specified network-area targets; no per-subarea percentage target.",
    paste0("Eligible anchor sources: ", if (length(cfg$anchors$sources)) paste(names(cfg$anchors$sources), collapse = ", ") else "none", ". Designation does not force selection."),
    if (cfg$model$mode == "single") "One independent model supplies the restoration score." else "Eight possible habitat pathways; at most one per selected planning unit.",
    if (cfg$quality$footprint == "free") "Final optimization can move restoration and change pathways while maintaining all requirements." else
      if (cfg$quality$footprint == "fixed") "Final optimization keeps the feasible footprint; pathway choices can still change." else
        sprintf("Final optimization may change at most %s ha from the feasible footprint (additions + removals).", fmt(cfg$quality$max_change_ha)))
}
dp_app_file_check <- function(cfg, projection = FALSE) {
  old <- setwd(cfg$project_root); on.exit(setwd(old), add = TRUE)
  paths <- unique(c(dp_region_inputs(cfg), cfg$budget$cost_per_ha,
    cfg$model$current_scores, unlist(lapply(cfg$model$outcomes, function(o) c(o$scores, o$allowed_mask)))))
  generated <- if (projection) {
    dir <- file.path(cfg$projection$output_root, cfg$projection$id, "metrics")
    c(unlist(lapply(dp_eight_outcomes(dir,"vrt"), `[[`, "scores")),
      file.path(dir, paste0(c("Richness", "PD", "RangeRarity"), "_current.vrt")))
  } else character()
  if (projection) paths <- unique(c(paths, cfg$projection$saved_fit_root, cfg$projection$training_manifest,
    cfg$projection$climate %or% cfg$projection$climate_manifest, cfg$projection$boundary))
  paths <- paths[!is.na(paths) & nzchar(paths)]
  status <- vapply(paths, function(p) {
    if (file.exists(p)) {
      if (tolower(tools::file_ext(p)) %in% c("shp", "vrt")) {
        err <- try(dp_files(p), silent = TRUE)
        if (inherits(err, "try-error")) return("Missing dependencies")
      }
      "Present"
    } else if (p %in% generated) "To be projected" else "Missing"
  }, character(1))
  data.frame(path = paths, status = status, row.names = NULL)
}
dp_app_read_settings <- function(path) {
  if (file.info(path)$size > 2*1024^2) stop("Settings file is too large.")
  x <- jsonlite::fromJSON(path, simplifyVector = TRUE)
  if (!is.list(x) || !identical(as.integer(x$schema_version), 1L) || !is.list(x$settings)) stop("Use a settings JSON exported by this app.")
  defaults <- dp_app_defaults(); unknown <- setdiff(names(x$settings), names(defaults))
  if (length(unknown)) stop("Unknown setting: ", paste(unknown, collapse = ", "))
  for (nm in names(x$settings)) {
    z <- x$settings[[nm]]; expected <- defaults[[nm]]
    if (length(z) != 1 || is.na(z) || (is.character(expected) && !is.character(z)) ||
        (is.logical(expected) && !is.logical(z)) || (is.numeric(expected) && (!is.numeric(z) || !is.finite(z))))
      stop("Invalid setting type for ", nm)
    if (is.numeric(expected)) x$settings[[nm]] <- as.numeric(z)
  }
  utils::modifyList(defaults, x$settings)
}

# Called in a separate R process. Generates exact anchor layers once, then
# samples only the display rasters. The optimizer always retains 10 m inputs.
dp_app_make_preview <- function(cfg, folder, max_cells = 150000L) {
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  old <- setwd(cfg$project_root); on.exit(setwd(old), add = TRUE)
  build_anchors(cfg, file.path(folder, "layers"))
  t <- dp_template(cfg, cfg$connectivity$context_boundary %or% cfg$region$boundary)
  b <- dp_vector(cfg$region$boundary, terra::crs(t$lu)); region <- dp_mask_polygon(b, t$lu)
  mask <- dp_mask_polygon(t$boundary, t$lu)
  fields <- c("nature", "anchors", "eligible", "fixed", "budget", "credit")
  paths <- setNames(file.path(folder, "layers", paste0(fields, "_10m.tif")), fields)
  counts <- vapply(paths, function(f) dp_sum(terra::rast(f) & region)*0.01, numeric(1))
  source <- file.path(folder, "layers", "anchor_sources_10m.tif")
  layers <- lapply(paths, terra::rast)
  if (file.exists(source)) {
    z <- terra::rast(source)
    for (i in seq_len(terra::nlyr(z))) layers[[paste0("source_", names(z)[i])]] <- z[[i]]
  }
  preview_paths <- vapply(names(layers), function(nm) {
    r <- terra::ifel(mask, layers[[nm]], NA)
    if (terra::ncell(r) > max_cells) r <- terra::spatSample(r, size = max_cells, method = "regular", as.raster = TRUE)
    f <- file.path(folder, paste0("preview_", nm, ".tif")); dp_raster(r, f, "INT1U"); f
  }, character(1))
  area <- sum(terra::expanse(terra::aggregate(b), unit = "ha", transform = FALSE))
  needed <- if (is.null(cfg$budget$target_nature_fraction)) 0 else max(0, cfg$budget$target_nature_fraction*area-counts[["budget"]])
  list(paths = preview_paths, counts = counts, polygon_ha = area, needed_ha = needed,
    minimum_ha = max(needed, cfg$budget$min_restoration_ha), boundary = cfg$region$boundary,
    subareas = cfg$region$subareas, project_root = cfg$project_root,
    folder = folder, exact_anchors = file.path(folder, "layers", "anchors_10m.tif"))
}
