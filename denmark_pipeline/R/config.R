`%or%` <- function(x, y) if (is.null(x)) y else x
dp_pathways <- function() data.frame(
  outcome = c("forest_dry_poor", "forest_wet_poor", "forest_dry_rich", "forest_wet_rich",
              "open_dry_poor", "open_wet_poor", "open_dry_rich", "open_wet_rich"),
  class = c("ForestDryPoor", "ForestWetPoor", "ForestDryRich", "ForestWetRich",
            "OpenDryPoor", "OpenWetPoor", "OpenDryRich", "OpenWetRich"),
  stringsAsFactors = FALSE)
dp_eight_outcomes <- function(metrics_dir, extension = "tif") {
  p <- dp_pathways()
  setNames(lapply(seq_len(nrow(p)), function(i) list(class = p$class[i],
    scores = setNames(file.path(metrics_dir, paste0(c("Richness", "PD", "RangeRarity"),
      "_all_to_", p$outcome[i], ".", extension)), c("Richness", "PD", "RangeRarity")),
    allowed_mask = NULL)), p$outcome)
}
denmark_config <- function() list(
  project_root = ".", area_id = "new_area", run_id = "v3_01", output_root = "Optimization/denmark",
  region = list(boundary = "Datasets/area.gpkg", subareas = NULL, subarea_id = NULL),
  grid = list(landuse = "Predictors/gudenaa_landuse/Landuse_Denmark_current.tif",
    resolution_m = 10, planning_unit_m = 200, origin = c(0, 0), directions = 4L),
  eligibility = list(classes = c("Agriculture", "ProductionForest"), mask = NULL,
    min_unit_ha = 1, forest_only = TRUE, forest_share = 0.5),
  anchors = list(sources = list(nature = list(type = "nature")), combine = "union"),
  connectivity = list(context_boundary = NULL, global_area_scope = "within_region",
    global_min_ha = 5000, subarea_min_ha = NULL, require_new_restoration = TRUE),
  budget = list(target_nature_fraction = 0.30, target_tolerance_ha = 16,
    min_restoration_ha = 0, max_restoration_ha = NULL, cost_cap = NULL, cost_per_ha = NULL),
  model = list(mode = "pathways", outcomes = dp_eight_outcomes("Results/new_area/metrics"),
    metric_weights = c(Richness = 1/3, PD = 1/3, RangeRarity = 1/3),
    metric_scales = NULL, score_type = "endpoint", current_scores = NULL),
  # NULL phase weights inherit contiguity, preserving existing configurations.
  objective = list(contiguity = 1.2, warmup_contiguity = NULL, quality_contiguity = NULL,
    new_to_new = 0.5, to_anchors = 0.5),
  repair = list(enabled = TRUE, seconds = 1800, gap = 0.005),
  quality = list(footprint = "free", max_change_ha = NULL),
  solver = list(ampl = "ampl", name = "gurobi", options = "tech:outlev=1", gap = 0.005,
    warmup_seconds = 3600, feasibility_seconds = c(10800, 10800),
    quality_seconds = c(28800, 14400), execute = TRUE),
  projection = list(id = "new_area", boundary = NULL, output_root = "Results/denmark_projection",
    saved_fit_root = "Results/gudenaa_plant_biodiversity/full",
    training_manifest = "Predictors/gudenaa_scenarios/training_current_predictor_paths_Denmark.rds",
    climate = NULL, climate_manifest = "Predictors/current_predictor_paths.rds",
    tile_cells = 1000L, workers = 4L, memmax_gb_per_worker = 4,
    blocks_in_memory = 512L, metric_chunk_size = 20L, pd_max_cells = 10000L,
    keep_species_rasters = FALSE)
)
dp_config <- function(path) {
  if (is.list(path)) x <- path else {
    e <- new.env(parent = globalenv()); sys.source(path, e); x <- e$config
  }
  x$repair <- utils::modifyList(list(enabled = TRUE, seconds = 1800, gap = 0.005), x$repair %or% list(), keep.null = TRUE)
  x$quality <- utils::modifyList(list(footprint = "free", max_change_ha = NULL), x$quality %or% list(), keep.null = TRUE)
  stopifnot(is.list(x), x$model$mode %in% c("pathways", "single"),
    x$anchors$combine %in% c("union", "intersection"),
    x$connectivity$global_area_scope %in% c("within_region", "within_context"),
    x$grid$directions %in% c(4L, 8L), x$grid$resolution_m == 10,
    x$grid$planning_unit_m >= 10, x$grid$planning_unit_m %% 10 == 0,
    x$model$score_type %in% c("endpoint", "gain"))
  for (id in c(x$area_id, x$run_id, names(x$model$outcomes), names(x$anchors$sources)))
    if (!grepl("^[A-Za-z][A-Za-z0-9_]*$", id)) stop("Use safe identifiers starting with a letter: ", id)
  if (!length(x$model$outcomes) || anyDuplicated(names(x$model$outcomes))) stop("Invalid outcomes.")
  if (x$model$mode == "single" && length(x$model$outcomes) != 1L) stop("Single mode requires one supplied outcome.")
  if (x$model$mode == "single" && isTRUE(x$eligibility$forest_only))
    stop("For the independent single model, explicitly set eligibility$forest_only = FALSE, or supply a forest-specific model.")
  if (x$model$score_type == "gain" && is.null(x$model$current_scores)) stop("Gain scoring needs current_scores.")
  w <- x$model$metric_weights
  if (!length(w) || any(!is.finite(w) | w < 0) || abs(sum(w) - 1) > 1e-9) stop("Metric weights must sum to one.")
  for (o in x$model$outcomes) {
    if (!setequal(names(o$scores), names(w))) stop("Each outcome needs the configured metrics.")
    if (!is.character(o$scores) || anyNA(o$scores) || any(!nzchar(trimws(o$scores)))) stop("Supply a raster path for each outcome metric.")
  }
  if (x$model$score_type == "gain" && !setequal(names(x$model$current_scores), names(w)))
    stop("Gain scoring needs a current score for every configured metric.")
  if (x$model$score_type == "gain" && (!is.character(x$model$current_scores) || anyNA(x$model$current_scores) || any(!nzchar(trimws(x$model$current_scores)))))
    stop("Supply current score raster paths for gain scoring.")
  if (!is.null(x$model$metric_scales) && (!setequal(names(x$model$metric_scales), names(w)) ||
      any(!is.finite(x$model$metric_scales) | x$model$metric_scales <= 0))) stop("Metric scales must be positive and match the configured metrics.")
  if (!is.null(x$budget$target_nature_fraction) &&
      (!is.finite(x$budget$target_nature_fraction) || x$budget$target_nature_fraction < 0 || x$budget$target_nature_fraction > 1))
    stop("Nature fraction must be in [0,1], or NULL.")
  if (!is.null(x$connectivity$subarea_min_ha) && is.null(x$region$subareas)) stop("Subarea goals require a subarea file.")
  if (!is.null(x$budget$cost_cap) && is.null(x$budget$cost_per_ha)) stop("A cost cap needs a cost-per-ha raster.")
  positive <- function(v, name, zero = TRUE, nullable = FALSE) {
    if (is.null(v) && nullable) return(invisible(NULL))
    if (!is.numeric(v) || !length(v) || any(!is.finite(v) | if (zero) v < 0 else v <= 0))
      stop(name, if (zero) " must be finite and nonnegative." else " must be finite and positive.")
  }
  for (nm in c("min_restoration_ha", "max_restoration_ha", "target_tolerance_ha", "cost_cap"))
    positive(x$budget[[nm]], paste0("budget$", nm), nullable = nm != "min_restoration_ha")
  if (!is.null(x$budget$max_restoration_ha) && x$budget$max_restoration_ha < x$budget$min_restoration_ha)
    stop("Maximum restoration area is below the minimum.")
  positive(x$eligibility$min_unit_ha, "Minimum planning-unit area")
  positive(x$eligibility$forest_share, "Production-forest share")
  if (x$eligibility$forest_share > 1) stop("Production-forest share must be in [0,1].")
  for (nm in c("global_min_ha", "subarea_min_ha")) positive(x$connectivity[[nm]], nm, nullable = TRUE)
  if (any(c(x$connectivity$global_min_ha, x$connectivity$subarea_min_ha) > 0) &&
      !length(x$anchors$sources))
    stop("Network requirements need at least one enabled anchor source.")
  for (nm in c("contiguity", "warmup_contiguity", "quality_contiguity", "new_to_new", "to_anchors"))
    positive(x$objective[[nm]], paste0("objective$", nm), nullable = nm %in% c("warmup_contiguity", "quality_contiguity"))
  if (!x$quality$footprint %in% c("free", "capped", "fixed")) stop("Quality footprint must be free, capped, or fixed.")
  if (x$quality$footprint == "capped") positive(x$quality$max_change_ha, "Quality footprint change cap")
  if (!is.logical(x$repair$enabled) || length(x$repair$enabled) != 1L || is.na(x$repair$enabled)) stop("Invalid repair switch.")
  positive(x$repair$seconds, "Repair time", zero = FALSE)
  if (length(x$repair$seconds) != 1L) stop("Repair time must be a single number.")
  for (v in list(x$solver$gap, x$repair$gap))
    if (length(v) != 1L || !is.finite(v) || v < 0 || v >= 1) stop("Solver gaps must be in [0,1).")
  if (!x$solver$name %in% c("gurobi", "highs")) stop("Supported solver names: gurobi, highs.")
  if (length(x$solver$options) != 1L || !is.character(x$solver$options) || is.na(x$solver$options) ||
      grepl("[\"'\n\r]", x$solver$options)) stop("Solver options must be one string without quotes or newlines.")
  if (length(x$grid$origin) != 2L || any(!is.finite(x$grid$origin))) stop("Supply a finite x/y grid origin.")
  for (nm in names(x$anchors$sources)) {
    s <- x$anchors$sources[[nm]]
    if (!s$type %in% c("nature", "vector", "raster")) stop("Invalid anchor type: ", nm)
    if (s$type != "nature" && (is.null(s$path) || !nzchar(s$path))) stop("Anchor source needs a path: ", nm)
    positive(s$buffer_m %or% 0, paste0("Anchor buffer for ", nm))
    if ((s$buffer_m %or% 0) %% 10 != 0) stop("Anchor buffers must be multiples of 10 m.")
    if ((isTRUE(s$credit_budget) || isTRUE(s$credit_network)) && !isTRUE(s$include_as_existing))
      stop("Anchor area credit requires include_as_existing: ", nm)
  }
  times <- unlist(x$solver[c("warmup_seconds", "feasibility_seconds", "quality_seconds")])
  if (any(!is.finite(times) | times < 1) || length(x$solver$warmup_seconds) != 1L ||
      !length(x$solver$feasibility_seconds) || !length(x$solver$quality_seconds) || x$repair$seconds < 1)
    stop("Supply time budgets of at least one second, including feasibility and quality stages.")
  x
}
dp_phase_weight <- function(cfg, phase) cfg$objective[[paste0(phase, "_contiguity")]] %or% cfg$objective$contiguity
dp_sum <- function(x) as.numeric(terra::global(x, "sum", na.rm = TRUE)[1, 1])
dp_max <- function(x) as.numeric(terra::global(x, "max", na.rm = TRUE)[1, 1])
dp_any <- function(x) dp_sum(terra::ifel(!is.na(x) & x, 1, 0)) > 0
dp_write <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(x, path, useBytes = TRUE); invisible(path)
}
dp_csv <- function(x, path) { utils::write.csv(x, path, row.names = FALSE); invisible(path) }
dp_save <- function(x, path) {
  tmp <- paste0(path, ".building"); saveRDS(x, tmp)
  if (!file.copy(tmp, path, overwrite = TRUE)) stop("Cannot save ", path)
  unlink(tmp); invisible(path)
}
dp_files <- function(paths, seen = character()) {
  paths <- unique(unname(unlist(paths))); paths <- paths[!is.na(paths) & nzchar(paths)]
  out <- unique(unlist(lapply(paths, function(p) {
    if (tolower(tools::file_ext(p)) == "shp") return(paste0(tools::file_path_sans_ext(p), c(".shp", ".shx", ".dbf", ".prj")))
    c(p, paste0(p, ".aux.xml")[file.exists(paste0(p, ".aux.xml"))])
  })))
  if (any(!file.exists(out))) stop("Missing files:\n", paste(out[!file.exists(out)], collapse = "\n"))
  # terra::sources() can return the VRT itself. Read its explicit file references
  # so changes in a backing tile invalidate cached optimization inputs too.
  for (p in paths[tolower(tools::file_ext(paths)) == "vrt"]) {
    p <- normalizePath(p, winslash = "/", mustWork = TRUE)
    if (p %in% seen) next
    xml <- paste(readLines(p, warn = FALSE), collapse = "\n")
    tags <- regmatches(xml, gregexpr("(?s)<Source(?:Filename|Dataset)\\b[^>]*>.*?</Source(?:Filename|Dataset)>", xml, perl = TRUE))[[1]]
    if (!length(tags)) stop("Cannot fingerprint the sources of VRT: ", p)
    refs <- sub("(?s)^<[^>]+>(.*?)</[^>]+>$", "\\1", tags, perl = TRUE)
    for (entry in list(c("&quot;", '"'), c("&apos;", "'"), c("&lt;", "<"), c("&gt;", ">"), c("&amp;", "&")))
      refs <- gsub(entry[1], entry[2], refs, fixed = TRUE)
    if (any(grepl("&#[0-9xXa-fA-F]+;", refs))) stop("Use plain filename references in VRT: ", p)
    relative <- grepl("relativeToVRT=[\"']1[\"']", tags)
    refs[relative] <- file.path(dirname(p), refs[relative])
    out <- unique(c(out, dp_files(refs, c(seen, p))))
  }
  out
}
dp_signature <- function(cfg, paths) {
  src <- list.files(file.path(dp_home, "R"), full.names = TRUE)
  src <- c(src, list.files(file.path(dp_home, "models"), full.names = TRUE),
           list.files(file.path(dp_home, "adapters"), full.names = TRUE))
  digest::digest(list(cfg, tools::md5sum(dp_files(paths)), tools::md5sum(src)))
}
dp_manifest <- function(folder, signature, config) {
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  f <- file.path(folder, "manifest.rds")
  if (file.exists(f) && !identical(readRDS(f)$signature, signature))
    stop("Inputs/settings/code changed. Use a new run_id or projection id: ", folder)
  if (!file.exists(f)) dp_save(list(signature = signature, config = config, created = Sys.time()), f)
}
dp_lock <- function(folder) {
  lock <- file.path(folder, "RUNNING.lock")
  if (!dir.create(lock, showWarnings = FALSE)) stop("Run is locked: ", lock,
    ". Remove a stale lock only after confirming no process is active.")
  lock
}
