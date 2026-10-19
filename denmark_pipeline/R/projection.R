dp_legacy <- function() {
  e <- new.env(parent = globalenv())
  sys.source(file.path(dp_home, "adapters", "speciespoolr_saved_models.R"), e)
  e
}
dp_projection_registry <- function(cfg) {
  h <- dp_legacy(); root <- cfg$projection$saved_fit_root
  status_files <- list.files(file.path(root, "national_baseline"), pattern = "_status[.]csv$", full.names = TRUE)
  if (!length(status_files)) stop("No saved national status files found in ", root)
  statuses <- h$read_csv_files(status_files)
  statuses <- statuses[statuses$usable_for_metrics %in% TRUE, , drop = FALSE]
  if (!nrow(statuses) || anyDuplicated(statuses$species)) stop("Invalid national species set.")
  statuses <- statuses[order(statuses$species), , drop = FALSE]
  for (col in c("threshold_path", "buffer_path")) {
    moved <- file.path(root, "national_baseline", basename(statuses[[col]]))
    statuses[[col]][!file.exists(statuses[[col]])] <- moved[!file.exists(statuses[[col]])]
  }
  statuses$model_path <- file.path(root, "models", paste0(vapply(statuses$species, h$species_id, character(1)), "_model.rds"))
  tree <- file.path(root, "phylogeny", "plant_tree.rds")
  coverage <- file.path(root, "audits", "tree_coverage.csv")
  files <- dp_files(c(status_files, statuses$model_path, statuses$threshold_path, statuses$buffer_path, tree, coverage))
  if (any(!is.finite(statuses$weight) | statuses$weight <= 0) ||
      any(!is.finite(statuses$range_size) | statuses$range_size <= 0) ||
      any(abs(statuses$weight*statuses$range_size - 1) > 1e-10) || any(statuses$range_unit != "km2"))
    stop("National rarity weights must be inverse range area in km2.")
  cov <- utils::read.csv(coverage, stringsAsFactors = FALSE)
  if (!all(statuses$species %in% cov$species) || any(!cov$matched[match(statuses$species, cov$species)]))
    stop("The saved phylogeny does not cover the same national species set.")
  list(status = statuses, tree = tree, coverage = cov[match(statuses$species, cov$species), , drop = FALSE], files = files)
}
dp_projection_plan <- function(config, validate_models = TRUE) {
  cfg <- dp_config(config); p <- cfg$projection
  if (!grepl("^[A-Za-z][A-Za-z0-9_]*$", p$id) || p$tile_cells < 1 || p$tile_cells != as.integer(p$tile_cells) ||
      p$workers < 1 || p$workers != as.integer(p$workers)) stop("Invalid projection ID, tile size, or worker count.")
  train <- unname(unlist(readRDS(p$training_manifest)))
  if (length(train) != 2) stop("Training manifest must contain Landuse and climate paths.")
  dp_files(c(p$training_manifest, train))
  train_lu <- terra::rast(train[1]); climate <- p$climate
  if (is.null(climate)) {
    if (is.null(p$climate_manifest) || !file.exists(p$climate_manifest))
      stop("Set projection$climate to the original climate source used in document 02, or supply its climate_manifest.")
    original <- unname(unlist(readRDS(p$climate_manifest)))
    if (length(original) != 2L) stop("Original predictor manifest must contain Landuse and climate paths.")
    climate <- original[2]
  }
  lu <- terra::rast(cfg$grid$landuse)
  if (!identical(terra::cats(train_lu)[[1]], terra::cats(lu)[[1]])) stop("Projection and training Landuse category dictionaries differ.")
  if (!identical(names(terra::rast(climate)), names(terra::rast(train[2])))) stop("Climate predictor names differ from training.")
  dp_codes(lu, c(dp_pathways()$class, cfg$eligibility$classes))
  boundary <- p$boundary %or% cfg$region$boundary
  v <- dp_vector(boundary, terra::crs(lu))
  grid <- terra::crop(terra::rast(lu), terra::ext(v), snap = "out")
  whole <- terra::vect(terra::ext(lu), crs = terra::crs(lu))
  if (sum(terra::expanse(terra::erase(v, whole), unit = "ha", transform = FALSE)) > 1e-6) stop("Projection boundary exceeds land-use raster coverage.")
  jobs <- list(); j <- 0L
  for (row in seq.int(1L, nrow(grid), by = p$tile_cells)) for (col in seq.int(1L, ncol(grid), by = p$tile_cells)) {
    nr <- min(p$tile_cells, nrow(grid)-row+1L); nc <- min(p$tile_cells, ncol(grid)-col+1L)
    e <- as.vector(terra::ext(grid)); resolution <- terra::res(grid)
    ex <- terra::ext(e[1]+(col-1L)*resolution[1], e[1]+(col+nc-1L)*resolution[1],
                     e[4]-(row+nr-1L)*resolution[2], e[4]-(row-1L)*resolution[2])
    tile <- terra::crop(grid, ex)
    if (!dp_any(dp_mask_polygon(v, tile))) next
    j <- j+1L
    jobs[[j]] <- list(id = sprintf("tile_%06d", j), extent = as.vector(ex), nrow = nr, ncol = nc)
  }
  if (!length(jobs)) stop("No projection tiles intersect the boundary.")
  registry <- if (validate_models) dp_projection_registry(cfg) else NULL
  build <- NULL
  if (validate_models) {
    api <- c("PredictSuitability", "threshold_suitability", "sum_binary_files_tree", "calc_range_rarity_files", "calc_pd_files")
    if (!requireNamespace("SpeciesPoolR", quietly = TRUE)) stop("Install the tested SpeciesPoolR version before projection.")
    if (any(!api %in% getNamespaceExports("SpeciesPoolR"))) stop("SpeciesPoolR lacks the saved-model/file-backed metric API.")
    db <- system.file("R", "SpeciesPoolR.rdb", package = "SpeciesPoolR")
    build <- list(version = as.character(utils::packageVersion("SpeciesPoolR")), db_md5 = tools::md5sum(db))
  }
  inputs <- dp_files(c(cfg$grid$landuse, cfg$eligibility$mask, boundary, p$training_manifest,
    train, climate, if (is.null(p$climate)) p$climate_manifest,
    unlist(lapply(cfg$anchors$sources, `[[`, "path")), registry$files))
  # Runtime concurrency and tile memory limits do not change model semantics.
  scientific <- list(grid = cfg$grid, eligibility = cfg$eligibility, anchors = cfg$anchors,
    projection = p[setdiff(names(p), c("workers", "memmax_gb_per_worker", "keep_species_rasters"))],
    boundary = boundary, scenarios = dp_pathways(), speciespoolr = build)
  signature <- dp_signature(scientific, inputs)
  root <- file.path(p$output_root, p$id)
  dp_manifest(root, signature, cfg); root <- normalizePath(root, winslash = "/")
  plan <- list(cfg = cfg, jobs = jobs, climate = normalizePath(climate, winslash = "/"),
    boundary = normalizePath(boundary, winslash = "/"), root = root, registry = registry, signature = signature)
  # Paths only cross process boundaries; terra external pointers are never sent.
  dp_save(plan, file.path(root, "projection_plan.rds"))
  dp_csv(data.frame(tile = vapply(jobs, `[[`, character(1), "id"),
    cells = vapply(jobs, function(j) j$nrow*j$ncol, numeric(1))), file.path(root, "tiles.csv"))
  plan
}
dp_tile_predictors <- function(plan, job, folder) {
  cfg <- plan$cfg; full <- terra::rast(cfg$grid$landuse)
  lu <- terra::crop(full, terra::ext(job$extent)); names(lu) <- "Landuse"
  boundary <- dp_vector(plan$boundary, terra::crs(lu)); region <- dp_mask_polygon(boundary, lu)
  layers <- dp_layers(cfg, lu, region, file.path(folder, "layers"))
  # Project directly from the unchanged full climate source. No tile-clipped
  # climate interpolation is used, avoiding artificial edge effects.
  climate <- terra::project(terra::rast(plan$climate), lu, method = "bilinear")
  complete <- terra::app(is.finite(as.numeric(c(lu, climate))), "sum") == (1+terra::nlyr(climate))
  if (dp_any(region & !complete)) stop("Incomplete predictors in tile ", job$id, "; fix input coverage.")
  domain <- terra::ifel(region & complete, 1, NA)
  lu <- terra::mask(lu, domain); lu <- terra::categories(lu, value = terra::cats(full)[[1]]); names(lu) <- "Landuse"
  climate <- terra::mask(climate, domain)
  paths <- list(current = file.path(folder, "Landuse_current.tif"))
  dp_raster(lu, paths$current, "INT2U")
  cp <- file.path(folder, "climate.tif"); dp_raster(climate, cp)
  dom <- file.path(folder, "domain.tif"); dp_raster(domain, dom, "INT1U")
  area <- file.path(folder, "cell_area_km2.tif"); dp_raster(terra::cellSize(domain, unit = "km"), area)
  selectable <- layers$eligible & !is.na(domain)
  p <- dp_pathways(); changes <- list()
  for (i in seq_len(nrow(p))) {
    scenario <- paste0("all_to_", p$outcome[i]); code <- dp_codes(full, p$class[i])
    z <- terra::ifel(selectable, code, as.numeric(lu)); z <- terra::categories(z, value = terra::cats(full)[[1]]); names(z) <- "Landuse"
    if (dp_any(!selectable & !is.na(lu) & as.numeric(z) != as.numeric(lu))) stop("A scenario changed protected cells.")
    if (dp_any(is.na(z) != is.na(lu))) stop("A scenario changed the missing-cell domain.")
    paths[[scenario]] <- file.path(folder, paste0("Landuse_", scenario, ".tif")); dp_raster(z, paths[[scenario]], "INT2U")
    changes[[i]] <- data.frame(scenario = scenario, selectable_cells = dp_sum(selectable),
      changed_cells = dp_sum(selectable & as.numeric(z) != as.numeric(lu)))
  }
  dp_csv(do.call(rbind, changes), file.path(folder, "scenario_audit.csv"))
  list(landuse = paths, climate = cp, domain = dom, area = area)
}
dp_project_tile <- function(plan, job) {
  cfg <- plan$cfg; p <- cfg$projection
  folder <- file.path(plan$root, "tiles", job$id); dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  done <- file.path(folder, "complete.rds")
  if (file.exists(done)) {
    d <- readRDS(done)
    if (identical(d$signature, plan$signature) && all(file.exists(d$metrics)) &&
        identical(unname(tools::md5sum(d$metrics)), d$hashes)) return(d$metrics)
    stop("Completed tile changed: ", folder)
  }
  lock <- dp_lock(folder); on.exit(unlink(lock, recursive = TRUE), add = TRUE)
  Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1", GDAL_NUM_THREADS = "1")
  needed <- c("SpeciesPoolR", "terra", "dplyr", "readr", "digest", "maxnet", "ape", "Matrix")
  missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Projection requires: ", paste(missing, collapse = ", "))
  h <- dp_legacy(); h$terra_memmax_gb <- p$memmax_gb_per_worker
  h$configure_terra_worker(p$memmax_gb_per_worker)
  inputs <- dp_tile_predictors(plan, job, folder)
  reg <- plan$registry; weights <- reg$status[c("species", "range_size", "weight", "range_unit")]
  dir.create(file.path(folder, "metrics"), showWarnings = FALSE)
  statusdir <- file.path(folder, "national_status"); dir.create(statusdir, showWarnings = FALSE)
  statuses <- vapply(seq_len(nrow(reg$status)), function(i) {
    f <- file.path(statusdir, paste0(h$species_id(reg$status$species[i]), "_status.csv"))
    dp_csv(reg$status[i, , drop = FALSE], f)
  }, character(1))
  metrics <- character()
  for (scenario in names(inputs$landuse)) {
    scenario_done <- file.path(folder, "metrics", paste0(scenario, "_complete.rds"))
    cached_scenario <- if (file.exists(scenario_done)) readRDS(scenario_done) else NULL
    if (!is.null(cached_scenario) && identical(cached_scenario$signature, plan$signature) &&
        all(file.exists(cached_scenario$metrics)) &&
        identical(unname(tools::md5sum(cached_scenario$metrics)), cached_scenario$hashes)) {
      metrics <- c(metrics, cached_scenario$metrics); next
    }
    products <- vector("list", nrow(reg$status))
    for (i in seq_len(nrow(reg$status))) {
      species <- reg$status$species[i]
      fit <- readRDS(reg$status$model_path[i])
      if (!identical(fit$species, species) || !identical(fit$status, "fitted") || is.null(fit$models)) stop("Invalid saved fit: ", species)
      rm(fit)
      # Reuse individual, fully completed species jobs after an interrupted tile.
      mark <- file.path(folder, "local_projections", scenario, paste0(h$species_id(species), "_done.rds"))
      cached <- if (file.exists(mark)) readRDS(mark) else NULL
      if (!is.null(cached) && identical(cached$signature, plan$signature) && all(file.exists(cached$products))) {
        products[[i]] <- cached$products
      } else {
        products[[i]] <- h$project_local_plant(reg$status$model_path[i], c(status = statuses[i]),
          c(inputs$landuse[[scenario]], inputs$climate), scenario, inputs$domain, inputs$area,
          folder, p$blocks_in_memory, p$memmax_gb_per_worker, FALSE)
        dp_save(list(signature = plan$signature, products = products[[i]]), mark)
      }
    }
    summary <- h$read_csv_files(h$collect_files(products, "_summary.csv"))
    binary <- h$select_local_paths(summary, weights, scenario)
    out <- setNames(file.path(folder, "metrics", paste0(c("Richness", "PD", "RangeRarity"), "_", scenario, ".tif")), c("Richness", "PD", "RangeRarity"))
    h$richness_file(binary, inputs$domain, out[["Richness"]], p$metric_chunk_size, p$memmax_gb_per_worker)
    h$rarity_file(binary, weights, inputs$domain, out[["RangeRarity"]], p$metric_chunk_size, p$memmax_gb_per_worker)
    h$pd_file(binary, reg$tree, reg$coverage, inputs$domain, out[["PD"]], p$pd_max_cells, p$memmax_gb_per_worker)
    audit <- h$metric_audit(scenario, out[["Richness"]], out[["RangeRarity"]], out[["PD"]],
      summary, weights, inputs$domain)
    dp_csv(audit, file.path(folder, "metrics", paste0("audit_", scenario, ".csv")))
    metrics <- c(metrics, out)
    dp_save(list(signature = plan$signature, metrics = unname(out), hashes = unname(tools::md5sum(out))), scenario_done)
    # Only remove this tile's generated species products after metric audits.
    if (!isTRUE(p$keep_species_rasters)) unlink(file.path(folder, "local_projections", scenario), recursive = TRUE)
  }
  dp_save(list(signature = plan$signature, metrics = unname(metrics), hashes = unname(tools::md5sum(metrics))), done)
  unname(metrics)
}
dp_mosaic_metrics <- function(plan, tile_products) {
  paths <- unname(unlist(tile_products)); folder <- file.path(plan$root, "metrics")
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  outputs <- vapply(unique(basename(paths)), function(nm) {
    tiles <- paths[basename(paths) == nm]
    if (length(tiles) != length(plan$jobs)) stop("Incomplete tile set for ", nm)
    dest <- file.path(folder, sub("[.]tif$", ".vrt", nm))
    terra::vrt(tiles, filename = dest, overwrite = TRUE)
    dest
  }, character(1))
  if (length(outputs) != 27L) stop("Expected current + eight pathways for three metrics.")
  dp_csv(data.frame(metric_file = names(outputs), path = outputs), file.path(plan$root, "metric_manifest.csv"))
  dp_write("COMPLETE", file.path(plan$root, "status.txt")); unname(outputs)
}
run_projection <- function(config) {
  cfg <- dp_config(config); old <- setwd(cfg$project_root); on.exit(setwd(old), add = TRUE)
  plan <- dp_projection_plan(cfg)
  lock <- dp_lock(plan$root); on.exit(unlink(lock, recursive = TRUE), add = TRUE)
  dp_write("RUNNING", file.path(plan$root, "status.txt"))
  if (cfg$projection$workers == 1L) products <- lapply(plan$jobs, function(job) dp_project_tile(plan, job)) else {
    cl <- parallel::makePSOCKcluster(min(cfg$projection$workers, length(plan$jobs)))
    on.exit(parallel::stopCluster(cl), add = TRUE)
    entry <- normalizePath(file.path(dp_home, "run.R"), winslash = "/")
    wd <- getwd()
    parallel::clusterCall(cl, function(entry, wd) { setwd(wd); source(entry); NULL }, entry, wd)
    products <- parallel::parLapplyLB(cl, plan$jobs, function(job, plan) dp_project_tile(plan, job), plan = plan)
  }
  outputs <- dp_mosaic_metrics(plan, products)
  message("Projection metrics: ", file.path(plan$root, "metrics")); invisible(outputs)
}
