# Plant biodiversity: fit in Denmark; project current + 8 scenarios in Gudenå.
#
# Put this file in the project root. Run 02_build_gudenaa_scenarios.Rmd first.
# Do not source this file to run the models: use targets::tar_make().
#
# 1. Tiny offline tests (no GBIF or modelling):
# targets::tar_make(
#   names = Workflow_smoke_test,
#   script = "targets_plant_biodiversity.R", store = "_targets_plant_biodiversity"
# )
# 2. Check predictors and the initial 30-species selection:
# targets::tar_make(
#   names = c(Input_audit_file, Taxon_scope_audit_file, Plant_selection_file),
#   script = "targets_plant_biodiversity.R", store = "_targets_plant_biodiversity"
# )
# 3. Fit nationally, calibrate nationally, calculate CURRENT local diversity:
# targets::tar_make(
#   names = c(Presence_audit_file, Model_audit_file, National_audit_file,
#             Range_weights_Denmark_file, Tree_coverage_file,
#             Richness_current_file, RangeRarity_current_file,
#             PD_current_file, Current_metric_audit_file),
#   script = "targets_plant_biodiversity.R", store = "_targets_plant_biodiversity"
# )
# 4. Only then run all eight local scenarios:
# targets::tar_make(
#   script = "targets_plant_biodiversity.R", store = "_targets_plant_biodiversity"
# )
#
# Always pass BOTH script and store to tar_outdated()/tar_visnetwork().
# tar_read() needs the store, not the script.
#
# Scope/defaults:
# - Vascular plants (Tracheophyta), not all Plantae: the CSV also includes algae.
# - All Denmark is the geographic sampling domain, NOT just Gudenå.
# - GBIF retrieval is CAPPED at gbif_download_limit records per species.
#   This is not a complete GBIF bulk download. Capped/filtered counts are audited.
# - Same SampleEnv() + ModelSpecies() sequence used by FitSpeciesModels().
#   Training tables are retained in each model RDS for inspection/reassembly.
# - Thresholds: national current suitability at Danish occurrence points.
# - Range sizes: national current binary area in km2 on the training grid.
# - Local richness, fixed-baseline range rarity and plant PD use one species set.
# - A 500 m occurrence-buffer restriction matches the earlier workflow. It is
#   NOT a fitted dispersal model. Change it only as a deliberate sensitivity run.
# - Binary masking happens AFTER thresholding, so a threshold of zero cannot
#   turn cells outside the occurrence buffer into presences.
# - Strict prediction checks: unknown categories/missing predictions INSIDE
#   accessible valid habitat stop that branch. They are never filled with zero.
# - PD requires complete tree coverage; unmatched taxa are never silently dropped.
# - No old _targets_climate or _targets_diversity_pilot stores are modified.
#
# Full run: set maximum_species <- Inf and run_label <- "full"; use a NEW store
# such as "_targets_plant_biodiversity_full" so the pilot remains recoverable.
# Do not update packages while a run is active. Package-code changes are tracked.
# Keep Results/ AND the selected targets store; neither is a substitute for the other.
# Cached GBIF failures are documented outcomes, not automatic infinite retries.
# To retry only those downloads after inspecting Presence_audit:
# failed <- subset(targets::tar_read(Presence_audit,
#   store = "_targets_plant_biodiversity"), status != "downloaded")$target_branch
# if (length(failed)) targets::tar_invalidate(
#   tidyselect::all_of(failed), store = "_targets_plant_biodiversity")
# Then re-run the same stage. Do not delete or invalidate the entire store.
#
# Before a full run, inspect model warnings, per-class support and spatial
# validation separately: the metric audits test arithmetic, not predictive skill.
# Full PD opens one binary raster per retained species; check the operating
# system's open-file limit before scaling up (the 30-species pilot is small).

# ---- Editable configuration ------------------------------------------------
maximum_species <- 30L
run_label <- "pilot30"
plant_phyla <- "Tracheophyta"
excluded_species <- character() # Explicit, documented taxonomic/ecological exclusions.
minimum_gbif_count <- 7L
minimum_complete_records <- 7L
gbif_years <- c(1999L, 2026L) # Fixed, reproducible query window; not Sys.Date().
gbif_download_limit <- 5000L
gbif_retries <- 5L
background_points <- 10000L
buffer_distance_m <- 500
retain_local_suitability <- FALSE

# 24 compute + 2 GBIF + 2 reduction + 1 tree = 29 workers (+ main R process).
# These are concurrency limits, not a hard RAM reservation. Enforce a real
# 30-CPU / ~300-GB cap with your server scheduler if required.
compute_workers <- 3L#24L
gbif_workers <- 2L
reduce_workers <- 2L
tree_workers <- 1L
terra_memmax_gb <- 6
projection_blocks_in_memory <- 512L # A block-sizing hint, NOT a hard cell cap.
metric_chunk_size <- 20L
pd_max_cells_per_block <- 10000L # calc_pd_files processes at least one raster row.

taxa_path <- "Datasets/Clean_Taxa.csv"
predictor_root <- "Predictors/gudenaa_scenarios"
output_root <- file.path("Results", "gudenaa_plant_biodiversity", run_label)

# ---- Environment and API checks --------------------------------------------
stopifnot(
  is.numeric(maximum_species), length(maximum_species) == 1L,
  !is.na(maximum_species), maximum_species >= 2,
  is.infinite(maximum_species) || maximum_species == floor(maximum_species),
  grepl("^[A-Za-z0-9_-]+$", run_label), length(gbif_years) == 2L,
  gbif_years[1] <= gbif_years[2], minimum_complete_records >= 7L,
  gbif_download_limit >= minimum_complete_records, gbif_download_limit <= 100000L,
  background_points >= 1L, buffer_distance_m > 0,
  compute_workers + gbif_workers + reduce_workers + tree_workers <= 29L
)
Sys.setenv(
  OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1", GDAL_NUM_THREADS = "1"
)
required_packages <- c(
  "targets", "crew", "SpeciesPoolR", "terra", "dplyr", "readr", "digest",
  "rgbif", "maxnet", "ape", "Matrix", "rtrees", "megatrees", "piggyback"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Install the missing packages before running: ",
       paste(missing_packages, collapse = ", "), call. = FALSE)
}
required_api <- c(
  "SampleEnv", "ModelSpecies", "PredictSuitability", "create_thresholds",
  "threshold_suitability", "make_buffer", "count_presences_simple",
  "sum_binary_files_tree", "calc_range_weights_files", "calc_range_rarity_files",
  "build_rtrees_tree", "audit_tree_coverage", "calc_pd_files"
)
missing_api <- setdiff(required_api, getNamespaceExports("SpeciesPoolR"))
if (length(missing_api)) {
  stop("The installed SpeciesPoolR is missing: ", paste(missing_api, collapse = ", "),
       ". Install the tested Geotargets version and restart R.", call. = FALSE)
}
# A local devtools installation may retain the same package version. Hash its
# lazy-load database as well, so changes to internal functions are not invisible.
package_db <- system.file("R", "SpeciesPoolR.rdb", package = "SpeciesPoolR")
SpeciesPoolR_build <- list(
  version = as.character(utils::packageVersion("SpeciesPoolR")),
  database_md5 = if (nzchar(package_db) && file.exists(package_db)) {
    unname(tools::md5sum(package_db))
  } else {
    digest::digest(lapply(required_api, function(nm) {
      paste(deparse(body(getExportedValue("SpeciesPoolR", nm))), collapse = "\n")
    }))
  },
  terra = as.character(utils::packageVersion("terra")),
  maxnet = as.character(utils::packageVersion("maxnet"))
)
library(targets)
library(terra)

controller_compute <- crew::crew_controller_local(
  name = "compute", workers = compute_workers, tasks_max = 50L
)
controller_gbif <- crew::crew_controller_local(
  name = "gbif", workers = gbif_workers, tasks_max = 50L
)
controller_reduce <- crew::crew_controller_local(
  name = "reduce", workers = reduce_workers, tasks_max = 20L
)
controller_tree <- crew::crew_controller_local(
  name = "tree", workers = tree_workers, tasks_max = 5L
)
targets::tar_option_set(
  packages = c("SpeciesPoolR", "terra", "dplyr", "readr", "digest", "rgbif", "maxnet"),
  seed = 20260902L, garbage_collection = TRUE,
  error = "continue", # Finish independent work; never hide an actual projection error.
  controller = crew::crew_controller_group(
    controller_compute, controller_gbif, controller_reduce, controller_tree
  ),
  resources = targets::tar_resources(
    crew = targets::tar_resources_crew(controller = "compute")
  )
)
gbif_resources <- targets::tar_resources(
  crew = targets::tar_resources_crew(controller = "gbif"))
reduce_resources <- targets::tar_resources(
  crew = targets::tar_resources_crew(controller = "reduce"))
tree_resources <- targets::tar_resources(
  crew = targets::tar_resources_crew(controller = "tree"))

# ---- Small I/O and validation helpers --------------------------------------
configure_terra_worker <- function(memmax_gb) {
  settings <- list(progress = 0L, memmax = memmax_gb, memmin = 0, memfrac = 0.02)
  available <- names(terra::terraOptions(print = FALSE))
  if ("threads" %in% available) settings$threads <- 1L
  do.call(terra::terraOptions, settings)
  invisible(NULL)
}
assert_paths <- function(paths) {
  if (!is.character(paths) || !length(paths) || anyNA(paths) ||
      any(!nzchar(paths)) || any(!file.exists(paths))) {
    stop("Missing or invalid file paths: ", paste(paths, collapse = ", "), call. = FALSE)
  }
  invisible(paths)
}
track_raster_files <- function(paths) {
  assert_paths(paths)
  sidecars <- unlist(lapply(paths, function(x) paste0(x, c(".aux.xml", ".vat.dbf"))),
                     use.names = FALSE)
  unique(c(unname(paths), sidecars[file.exists(sidecars)]))
}
read_pair <- function(manifest) {
  paths <- unname(unlist(readRDS(manifest), use.names = FALSE))
  assert_paths(paths)
  if (length(paths) != 2L) stop("Each predictor manifest must contain exactly two TIFFs.")
  paths
}
check_pair <- function(paths) {
  assert_paths(paths)
  lu <- terra::rast(paths[[1]])
  clim <- terra::rast(paths[[2]])
  if (terra::nlyr(lu) != 1L || !identical(names(lu), "Landuse") ||
      !isTRUE(terra::is.factor(lu)) || !nzchar(terra::crs(lu)) ||
      !terra::compareGeom(lu, clim, stopOnError = FALSE)) {
    stop("Predictor geometry / Landuse category metadata is invalid.")
  }
  rat <- levels(lu)[[1]]
  if (!is.data.frame(rat) || ncol(rat) < 2L ||
      anyDuplicated(rat[[1]]) || anyDuplicated(rat[[2]])) {
    stop("Invalid Landuse category table.")
  }
  list(
    names = names(c(lu, clim)),
    ids = as.integer(rat[[1]]), classes = as.character(rat[[2]]),
    rows = terra::nrow(lu), cols = terra::ncol(lu),
    resolution = terra::res(lu), crs = terra::crs(lu)
  )
}
check_input_set <- function(training_paths, local_paths, specs,
                            national_domain_path, local_domain_path) {
  training <- check_pair(training_paths)
  local <- check_pair(local_paths)
  if (!identical(training$names, local$names) ||
      !identical(training$ids, local$ids) ||
      !identical(training$classes, local$classes)) {
    stop("Training and local predictors must have identical names and class dictionaries.")
  }
  expected <- c(
    "all_to_forest_dry_poor", "all_to_forest_wet_poor",
    "all_to_forest_dry_rich", "all_to_forest_wet_rich",
    "all_to_open_dry_poor", "all_to_open_wet_poor",
    "all_to_open_dry_rich", "all_to_open_wet_rich"
  )
  if (!all(c("scenario", "path") %in% names(specs)) ||
      nrow(specs) != 8L || anyDuplicated(specs$scenario) ||
      !setequal(as.character(specs$scenario), expected)) {
    stop("The scenario manifest must contain the eight expected restoration scenarios.")
  }
  current <- terra::rast(local_paths[[1]])
  climate_normal <- normalizePath(local_paths[[2]], winslash = "/", mustWork = TRUE)
  for (paths in specs$path) {
    paths <- unname(unlist(paths, use.names = FALSE))
    if (length(paths) != 2L) stop("Invalid scenario path vector.")
    got <- check_pair(paths)
    if (!identical(got$names, local$names) || !identical(got$ids, local$ids) ||
        !identical(got$classes, local$classes) ||
        !terra::compareGeom(terra::rast(paths[[1]]), current, stopOnError = FALSE) ||
        !identical(normalizePath(paths[[2]], winslash = "/", mustWork = TRUE),
                   climate_normal)) {
      stop("A scenario differs in geometry, predictor names, classes or climate source.")
    }
    if (count_true(is.na(terra::rast(paths[[1]])) != is.na(current)) > 0) {
      stop("A scenario changed the baseline land-use missing-cell mask.")
    }
  }
  assert_paths(c(national_domain_path, local_domain_path))
  stopifnot(
    terra::compareGeom(terra::rast(national_domain_path), terra::rast(training_paths[[1]])),
    terra::compareGeom(terra::rast(local_domain_path), current)
  )
  for (item in list(list(paths = training_paths, domain = national_domain_path),
                    list(paths = local_paths, domain = local_domain_path))) {
    env <- terra::rast(item$paths)
    complete <- terra::app(is.finite(as.numeric(env)), fun = "sum") == terra::nlyr(env)
    domain <- terra::rast(item$domain)
    if (count_true(!is.na(domain)) == 0 ||
        count_true(!is.na(domain) != complete) > 0) {
      stop("A supplied domain is empty or differs from the complete predictor cells. ",
           "Re-run document 02; do not reuse a stale domain raster.")
    }
  }
  data.frame(
    training_region = "Denmark", projection_region = "Gudenå",
    training_rows = training$rows, training_cols = training$cols,
    local_rows = local$rows, local_cols = local$cols,
    predictors = length(training$names), local_conditions = 9L,
    category_dictionary_matches = TRUE, local_geometry_matches = TRUE,
    climate_held_constant = TRUE
  )
}
species_id <- function(species) {
  stopifnot(length(species) == 1L, !is.na(species), nzchar(species))
  stem <- substr(gsub("[^A-Za-z0-9]+", "_", species), 1L, 100L)
  paste0(stem, "_", digest::digest(species, algo = "xxhash64", serialize = FALSE))
}
save_rds_file <- function(object, filename) {
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  saveRDS(object, filename)
  filename
}
save_csv_file <- function(object, filename) {
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(object, filename, na = "")
  filename
}
collect_files <- function(products, suffix) {
  paths <- unname(unlist(products, use.names = FALSE))
  unique(paths[!is.na(paths) & nzchar(paths) & endsWith(paths, suffix)])
}
read_csv_files <- function(paths) {
  if (!length(paths)) return(data.frame())
  # Per-species status files can have entirely missing columns. Avoid guessing
  # a different type for the same column in successful and excluded branches.
  out <- dplyr::bind_rows(lapply(paths, function(path) readr::read_csv(
    path, col_types = readr::cols(.default = readr::col_character()),
    show_col_types = FALSE, progress = FALSE)))
  numeric_fields <- intersect(names(out), c(
    "range_size", "weight", "complete_presence_records", "occupied_cells",
    "occupied_area_km2", "accessible_cells", "national_range_km2",
    "fixed_national_weight"))
  logical_fields <- intersect(names(out), c(
    "usable_for_metrics", "projected", "no_buffer_overlap"))
  for (nm in numeric_fields) {
    value <- suppressWarnings(readr::parse_double(out[[nm]]))
    if (any(!is.na(out[[nm]]) & is.na(value))) {
      stop("Invalid numeric value in saved projection column: ", nm)
    }
    out[[nm]] <- value
  }
  for (nm in logical_fields) {
    value <- suppressWarnings(readr::parse_logical(out[[nm]]))
    if (any(!is.na(out[[nm]]) & is.na(value))) {
      stop("Invalid logical value in saved projection column: ", nm)
    }
    out[[nm]] <- value
  }
  out
}
count_true <- function(x) {
  as.numeric(terra::global(terra::ifel(!is.na(x) & x, 1, 0),
                           "sum", na.rm = TRUE)[1, 1])
}
sum_raster <- function(x) {
  as.numeric(terra::global(x, "sum", na.rm = TRUE)[1, 1])
}
write_raster_file <- function(x, filename, datatype = "FLT8S", layer_name = NULL) {
  if (!is.null(layer_name)) names(x) <- layer_name
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  wopt <- list(datatype = datatype,
               gdal = c("COMPRESS=DEFLATE", "TILED=YES", "BIGTIFF=IF_SAFER"))
  if (datatype == "INT1U") wopt$NAflag <- 255
  if (datatype == "INT2U") wopt$NAflag <- 65535
  terra::writeRaster(x, filename, overwrite = TRUE, wopt = wopt)
  filename
}
assert_real_prediction <- function(x, species) {
  if (any(grepl("SpeciesPoolR_zero_suitability_",
                basename(terra::sources(x)), fixed = TRUE))) {
    stop("Prediction failed for ", species,
         ": SpeciesPoolR returned its zero-raster error fallback.", call. = FALSE)
  }
  invisible(x)
}
# Delete only transient raster files beneath this worker's actual temp folders.
# Explicit model, binary, threshold, buffer and status files are never in this list.
temporary_sources <- function(...) {
  paths <- unique(unlist(lapply(list(...), function(x) {
    if (inherits(x, "SpatRaster")) as.character(terra::sources(x)) else character()
  }), use.names = FALSE))
  paths <- paths[!is.na(paths) & nzchar(paths)]
  roots <- unique(c(tempdir(), terra::terraOptions(print = FALSE)$tempdir))
  roots <- roots[!is.na(roots) & nzchar(roots)]
  norm <- function(x) {
    out <- normalizePath(x, winslash = "/", mustWork = FALSE)
    if (.Platform$OS.type == "windows") tolower(out) else out
  }
  keep <- vapply(norm(paths), function(path) {
    any(startsWith(path, paste0(norm(roots), "/")))
  }, logical(1))
  paths[keep]
}
remove_transient_sources <- function(paths) {
  paths <- unique(paths[!is.na(paths) & nzchar(paths)])
  if (length(paths)) {
    invisible(gc())
    unlink(paths[file.exists(paths)], force = TRUE)
  }
  invisible(NULL)
}
# A reachability mask is shared across current and future habitat scenarios.
# Unknown predictions inside that mask are errors; outside it genuine zeros.
binary_with_access <- function(suitability, thresholds, buffer, domain, species) {
  stopifnot(terra::compareGeom(suitability, domain))
  access <- terra::mask(domain, buffer)
  missing <- count_true(!is.na(access) & !is.finite(suitability))
  if (missing > 0) {
    stop(species, ": ", missing, " accessible cells have no finite prediction. ",
         "Check unsupported Landuse levels and predictor coverage; these are not absences.",
         call. = FALSE)
  }
  th <- thresholds$Thres_95[match(species, thresholds$species)]
  if (length(th) != 1L || !is.finite(th)) stop("No finite national threshold for ", species)
  names(suitability) <- species
  binary <- SpeciesPoolR::threshold_suitability(
    Model = suitability, Thresholds = thresholds, threshold = "Thres_95")
  binary <- terra::mask(binary, buffer, updatevalue = 0)
  binary <- terra::mask(binary, domain)
  if (count_true(!is.na(domain) & is.na(binary)) > 0) {
    stop("Binary projection has holes inside the valid domain for ", species)
  }
  names(binary) <- species
  binary
}

# ---- Taxonomy and reproducible occurrence selection ------------------------
read_plant_taxa <- function(path, phyla, excluded) {
  # Read explicitly as character to avoid acceptedUsageKey being inferred as
  # logical, the parsing issue seen in the earlier Clean_Taxa.csv run.
  raw <- readr::read_csv(
    path, col_types = readr::cols(.default = readr::col_character()),
    na = c("", "NA"), show_col_types = FALSE, progress = FALSE
  )
  required <- c("kingdom", "phylum", "class", "family", "genus", "species",
                "gbif_speciesKey")
  if (!all(required %in% names(raw))) stop("Clean_Taxa.csv is missing required columns.")
  x <- as.data.frame(raw[required], stringsAsFactors = FALSE)
  for (nm in required) x[[nm]] <- trimws(x[[nm]])
  in_scope <- x$kingdom %in% "Plantae" & x$phylum %in% phyla
  good <- in_scope & !is.na(x$species) & nzchar(x$species) &
    !is.na(x$genus) & nzchar(x$genus) & !is.na(x$family) & nzchar(x$family) &
    !x$species %in% excluded
  numeric_keys <- suppressWarnings(as.numeric(x$gbif_speciesKey))
  good <- good & is.finite(numeric_keys) & numeric_keys > 0
  x$gbif_speciesKey <- numeric_keys
  kept <- x[good, , drop = FALSE]
  duplicates <- split(kept[c("family", "genus", "gbif_speciesKey")], kept$species)
  conflict <- names(duplicates)[vapply(duplicates,
    function(z) nrow(unique(z)) > 1L, logical(1))]
  if (length(conflict)) {
    stop("Conflicting taxonomy/GBIF keys for: ", paste(conflict, collapse = ", "))
  }
  kept <- kept[!duplicated(kept$species), , drop = FALSE]
  kept <- kept[order(kept$species), , drop = FALSE]
  norm <- gsub("[ .]+", "_", kept$species)
  if (anyDuplicated(norm)) stop("Plant names collide after normalisation.")
  if (nrow(kept) < 2L) stop("Fewer than two eligible plants in the selected scope.")
  list(
    taxa = kept,
    audit = data.frame(
      input_rows = nrow(raw), plantae_rows = sum(x$kingdom %in% "Plantae"),
      scope_rows = sum(in_scope), valid_scope_rows = sum(good),
      unique_scope_species = nrow(kept),
      selected_phyla = paste(phyla, collapse = ";"),
      explicit_exclusions = paste(excluded, collapse = ";")
    )
  )
}
select_plants <- function(taxa, counts, minimum_count, maximum) {
  counts <- as.data.frame(counts)
  if (anyDuplicated(counts$species)) stop("Duplicate species in the GBIF count result.")
  x <- taxa
  x$N <- counts$N[match(x$species, counts$species)]
  x <- x[!is.na(x$N) & x$N >= minimum_count, , drop = FALSE]
  # Rank high-count species within families, then take a round-robin across
  # families. This is a reproducible technical pilot, not a probability sample.
  x <- x[order(x$family, -x$N, x$species), , drop = FALSE]
  if (!nrow(x)) stop("No plants passed the GBIF count threshold.")
  rank_in_family <- ave(seq_len(nrow(x)), x$family, FUN = seq_along)
  x <- x[order(rank_in_family, -x$N, x$species), , drop = FALSE]
  if (is.finite(maximum)) x <- utils::head(x, as.integer(maximum))
  x$selection_rank <- seq_len(nrow(x))
  x <- x[order(x$species), , drop = FALSE]
  row.names(x) <- NULL
  x
}
download_plant_presences <- function(species, taxon_key, genus, family, root,
                                    years, limit, retries) {
  filename <- file.path(root, "occurrences", paste0(species_id(species), ".rds"))
  query <- list(taxonKey = as.numeric(taxon_key), country = "DK",
                year = paste(years, collapse = ","), hasCoordinate = TRUE,
                hasGeospatialIssue = FALSE, limit = limit)
  last_error <- NULL
  response <- NULL
  for (attempt in seq_len(retries)) {
    response <- tryCatch(do.call(rgbif::occ_data, query), error = function(e) {
      last_error <<- conditionMessage(e)
      NULL
    })
    if (is.list(response) && is.data.frame(response$data) && nrow(response$data) > 0) break
    if (attempt < retries) Sys.sleep(min(30, 2^(attempt - 1L)))
  }
  ok <- is.list(response) && is.data.frame(response$data) && nrow(response$data) > 0
  records <- if (ok) response$data else data.frame()
  # Preserve source identifiers and quality fields as well as modelling columns.
  keep_columns <- intersect(c(
    "key", "gbifID", "species", "speciesKey", "scientificName", "taxonKey",
    "acceptedTaxonKey", "family", "genus", "countryCode", "decimalLongitude",
    "decimalLatitude", "coordinateUncertaintyInMeters", "year", "basisOfRecord",
    "occurrenceStatus", "datasetKey", "license"
  ), names(records))
  records <- records[keep_columns]
  matched_count <- if (is.list(response) && length(response$meta$count)) {
    as.numeric(response$meta$count[[1L]])
  } else NA_real_
  out <- list(
    species = species, taxon_key = as.numeric(taxon_key), genus = genus, family = family,
    target_branch = targets::tar_name(),
    query = query, retrieved_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
    status = if (ok) "downloaded" else "download_failed_or_empty",
    reason = if (ok) NA_character_ else if (is.null(last_error)) {
      "No usable GBIF records returned after retries."
    } else last_error,
    matched_count = matched_count, download_limit = limit, records = records,
    download_capped = isTRUE(matched_count > limit) || nrow(records) >= limit
  )
  if (!ok) warning("Skipping ", species, ": ", out$reason, call. = FALSE)
  save_rds_file(out, filename)
}
clean_occurrences <- function(download) {
  raw <- download$records
  required <- c("decimalLongitude", "decimalLatitude")
  if (!is.data.frame(raw) || !all(required %in% names(raw)) || !nrow(raw)) {
    return(data.frame(species = character(), genus = character(), family = character(),
                      decimalLongitude = double(), decimalLatitude = double()))
  }
  lon <- suppressWarnings(as.numeric(raw$decimalLongitude))
  lat <- suppressWarnings(as.numeric(raw$decimalLatitude))
  keep <- is.finite(lon) & is.finite(lat) & lon >= -180 & lon <= 180 &
    lat >= -90 & lat <= 90
  if ("speciesKey" %in% names(raw)) {
    keep <- keep & suppressWarnings(as.numeric(raw$speciesKey)) %in% download$taxon_key
  } else if ("species" %in% names(raw)) {
    keep <- keep & raw$species %in% download$species
  } else {
    stop("The GBIF response cannot be matched to the requested species.")
  }
  if ("countryCode" %in% names(raw)) keep <- keep & raw$countryCode %in% "DK"
  if ("occurrenceStatus" %in% names(raw)) {
    keep <- keep & !toupper(raw$occurrenceStatus) %in% "ABSENT"
  }
  x <- data.frame(
    species = rep(download$species, sum(keep)),
    genus = rep(download$genus, sum(keep)), family = rep(download$family, sum(keep)),
    decimalLongitude = lon[keep], decimalLatitude = lat[keep]
  )
  unique(x)
}
complete_occurrences <- function(x, predictors) {
  if (!nrow(x)) return(x)
  env <- terra::rast(predictors)
  pts <- terra::vect(x, geom = c("decimalLongitude", "decimalLatitude"), crs = "EPSG:4326")
  pts <- terra::project(pts, terra::crs(env))
  v <- terra::extract(env, pts, ID = FALSE)
  keep <- stats::complete.cases(v)
  numeric_columns <- vapply(v, is.numeric, logical(1))
  if (any(numeric_columns)) {
    keep <- keep & rowSums(!is.finite(as.matrix(v[numeric_columns]))) == 0
  }
  x[keep, , drop = FALSE]
}

# ---- Model and national baseline products ----------------------------------
fit_plant <- function(occurrence_file, predictors, root, minimum_records,
                      n_background, memmax_gb) {
  configure_terra_worker(memmax_gb)
  download <- readRDS(occurrence_file)
  species <- download$species
  presences <- clean_occurrences(download)
  before <- nrow(presences)
  presences <- complete_occurrences(presences, predictors)
  status <- if (download$status != "downloaded") download$status else if (
    nrow(presences) < minimum_records) "insufficient_complete_presences" else "ready"
  reason <- if (status == "insufficient_complete_presences") {
    paste("Only", nrow(presences), "unique records have complete national predictors;",
          "at least", minimum_records, "are required.")
  } else download$reason
  training <- NULL
  model <- NULL
  warnings_seen <- character()
  if (status == "ready") {
    # This is the single-species internals of FitSpeciesModels(), made explicit
    # solely so the exact presence/background table can be retained.
    sampling_error <- NULL
    training <- tryCatch(withCallingHandlers({
      pres <- SpeciesPoolR::SampleEnv(
        presences, file = predictors, categorical = "Landuse", type = "pres")
      bg <- SpeciesPoolR::SampleEnv(
        presences, file = predictors, categorical = "Landuse",
        type = "bg", n_bg = n_background)
      dplyr::bind_rows(pres, bg)
    }, warning = function(w) {
      warnings_seen <<- c(warnings_seen, conditionMessage(w))
    }), error = function(e) {
      sampling_error <<- conditionMessage(e)
      NULL
    })
    if (!is.null(sampling_error)) {
      status <- "model_failed"
      reason <- paste("Training-data sampling failed:", sampling_error)
    } else if (!any(training$Pres == 0L)) {
      status <- "no_background_records"
      reason <- "No complete background records were available in the sampling extent."
    } else {
      model <- withCallingHandlers(SpeciesPoolR::ModelSpecies(training),
        warning = function(w) {
          warnings_seen <<- c(warnings_seen, conditionMessage(w))
        })
      status <- if (is.null(model)) "model_failed" else "fitted"
      if (is.null(model)) reason <- "ModelSpecies() returned NULL; inspect model warnings."
    }
  }
  category_support <- if (is.data.frame(training) && nrow(training)) {
    as.data.frame(table(Landuse = as.character(training$Landuse), Pres = training$Pres),
                  stringsAsFactors = FALSE)
  } else data.frame()
  out <- list(
    species = species, status = status, reason = reason,
    model = model, models = stats::setNames(list(model), species),
    presences = presences, training_data = training, category_support = category_support,
    audit = data.frame(
      species = species, status = status, reason = reason,
      downloaded_records = nrow(download$records), unique_valid_records = before,
      complete_presence_records = nrow(presences),
      model_ok = !is.null(model),
      model_class = if (is.null(model)) NA_character_ else paste(class(model), collapse = "/"),
      download_capped = download$download_capped,
      warnings = paste(unique(warnings_seen), collapse = " | ")
    )
  )
  save_rds_file(out, file.path(root, "models", paste0(species_id(species), "_model.rds")))
}
national_baseline <- function(model_file, predictors, domain_file, cell_area_file,
                              root, buffer_m, blocks, memmax_gb) {
  configure_terra_worker(memmax_gb)
  fit <- readRDS(model_file)
  species <- fit$species
  id <- species_id(species)
  folder <- file.path(root, "national_baseline")
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  status_path <- file.path(folder, paste0(id, "_status.csv"))
  if (fit$status != "fitted") {
    record <- data.frame(
      species = species, status = fit$status, reason = fit$reason,
      usable_for_metrics = FALSE, range_size = NA_real_, weight = NA_real_,
      range_unit = "km2", binary_path = NA_character_, threshold_path = NA_character_,
      buffer_path = NA_character_, complete_presence_records = nrow(fit$presences))
    return(c(status = save_csv_file(record, status_path)))
  }
  scratch <- character()
  on.exit(remove_transient_sources(scratch), add = TRUE)
  domain <- terra::rast(domain_file)
  suitability <- SpeciesPoolR::PredictSuitability(
    Models = fit$models, file = predictors, categorical = "Landuse",
    blocks_in_memory = blocks)
  scratch <- unique(c(scratch, temporary_sources(suitability)))
  assert_real_prediction(suitability, species)
  thresholds <- SpeciesPoolR::create_thresholds(Model = suitability, reference = fit$presences)
  th <- if (all(c("species", "Thres_95") %in% names(thresholds))) {
    thresholds$Thres_95[match(species, thresholds$species)]
  } else NA_real_
  if (length(th) != 1L || !is.finite(th)) {
    stop("No finite national threshold was calibrated for ", species, call. = FALSE)
  }
  buffer <- SpeciesPoolR::make_buffer(fit$presences, dist = buffer_m)
  buffer_for_grid <- terra::project(buffer, terra::crs(domain))
  binary <- binary_with_access(suitability, thresholds, buffer_for_grid, domain, species)
  scratch <- unique(c(scratch, temporary_sources(binary)))
  binary_path <- file.path(folder, paste0(id, "_binary.tif"))
  threshold_path <- file.path(folder, paste0(id, "_threshold.rds"))
  buffer_path <- file.path(folder, paste0(id, "_buffer.gpkg"))
  write_raster_file(binary, binary_path, datatype = "INT1U", layer_name = species)
  save_rds_file(thresholds, threshold_path)
  terra::writeVector(buffer, buffer_path, filetype = "GPKG", overwrite = TRUE)
  # Same area-based definition as calc_range_weights_files(), but cell areas
  # are precomputed once for the national grid instead of once per species.
  quantity <- binary * terra::rast(cell_area_file)
  scratch <- unique(c(scratch, temporary_sources(quantity)))
  range_size <- sum_raster(quantity)
  if (!is.finite(range_size) || range_size < 0) stop("Invalid national range for ", species)
  positive <- range_size > 0
  record <- data.frame(
    species = species, status = if (positive) "baseline_ready" else "zero_national_range",
    reason = if (positive) NA_character_ else
      "No occupied national cells after thresholding and buffer masking; inverse range undefined.",
    usable_for_metrics = positive, range_size = range_size,
    weight = if (positive) 1 / range_size else NA_real_, range_unit = "km2",
    binary_path = binary_path, threshold_path = threshold_path, buffer_path = buffer_path,
    complete_presence_records = nrow(fit$presences)
  )
  save_csv_file(record, status_path)
  rm(suitability, binary, quantity, buffer, buffer_for_grid, domain)
  invisible(gc())
  c(binary = binary_path, threshold = threshold_path, buffer = buffer_path,
    status = status_path)
}
read_national_record <- function(products) {
  status <- collect_files(products, "_status.csv")
  if (length(status) != 1L) stop("Expected exactly one national species status file.")
  read_csv_files(status)
}
make_range_weights <- function(national_status) {
  x <- national_status[national_status$usable_for_metrics %in% TRUE, , drop = FALSE]
  if (nrow(x) < 2L || anyDuplicated(x$species) ||
      any(!is.finite(x$range_size)) || any(x$range_size <= 0) ||
      any(!is.finite(x$weight)) || any(!file.exists(x$binary_path)) ||
      max(abs(x$range_size * x$weight - 1)) > 1e-12) {
    stop("The national range-weight table failed validation (need >=2 usable species).")
  }
  x[c("species", "range_size", "weight", "range_unit", "binary_path")]
}

# ---- Local current and scenario projections --------------------------------
project_local_plant <- function(model_file, national_products, predictors, scenario,
                                domain_file, cell_area_file, root, blocks,
                                memmax_gb, save_suitability = FALSE) {
  configure_terra_worker(memmax_gb)
  baseline <- read_national_record(national_products)
  species <- baseline$species[[1L]]
  folder <- file.path(root, "local_projections", scenario)
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  id <- species_id(species)
  status_path <- file.path(folder, paste0(id, "_status.csv"))
  if (!isTRUE(baseline$usable_for_metrics[[1L]])) {
    record <- data.frame(
      scenario = scenario, species = species, status = "excluded_from_national_baseline",
      reason = paste(baseline$status[[1]], baseline$reason[[1]], sep = ": "),
      projected = FALSE, binary_path = NA_character_)
    return(c(status = save_csv_file(record, status_path)))
  }
  fit <- readRDS(model_file)
  thresholds <- readRDS(baseline$threshold_path[[1]])
  buffer <- terra::vect(baseline$buffer_path[[1]])
  domain <- terra::rast(domain_file)
  buffer <- terra::project(buffer, terra::crs(domain))
  access <- terra::mask(domain, buffer)
  reachable_cells <- count_true(!is.na(access))
  scratch <- temporary_sources(access)
  on.exit(remove_transient_sources(scratch), add = TRUE)
  suitability <- NULL
  if (reachable_cells == 0) {
    # This is a documented buffer restriction, not a substituted error raster.
    binary <- terra::ifel(!is.na(domain), 0, NA)
    names(binary) <- species
  } else {
    # Predict only inside the fixed occurrence buffer. This is equivalent to
    # applying the same buffer afterwards, but avoids expensive model-matrix
    # evaluation over inaccessible cells in the 10 m local grid.
    prediction_env <- terra::mask(terra::rast(predictors), buffer)
    scratch <- unique(c(scratch, temporary_sources(prediction_env)))
    suitability <- SpeciesPoolR::PredictSuitability(
      Models = fit$models, file = prediction_env, categorical = "Landuse",
      blocks_in_memory = blocks)
    scratch <- unique(c(scratch, temporary_sources(suitability)))
    assert_real_prediction(suitability, species)
    binary <- binary_with_access(suitability, thresholds, buffer, domain, species)
    rm(prediction_env)
  }
  scratch <- unique(c(scratch, temporary_sources(binary)))
  binary_path <- file.path(folder, paste0(id, "_binary.tif"))
  summary_path <- file.path(folder, paste0(id, "_summary.csv"))
  write_raster_file(binary, binary_path, "INT1U", species)
  area <- binary * terra::rast(cell_area_file)
  scratch <- unique(c(scratch, temporary_sources(area)))
  summary <- data.frame(
    scenario = scenario, species = species, occupied_cells = sum_raster(binary),
    occupied_area_km2 = sum_raster(area), accessible_cells = reachable_cells,
    national_range_km2 = baseline$range_size[[1]],
    fixed_national_weight = baseline$weight[[1]], binary_path = binary_path,
    no_buffer_overlap = reachable_cells == 0)
  if (any(!is.finite(summary$occupied_area_km2)) || summary$occupied_area_km2 < 0) {
    stop("Invalid local occupied area for ", species, " / ", scenario)
  }
  save_csv_file(summary, summary_path)
  record <- data.frame(
    scenario = scenario, species = species, status = "projected",
    reason = if (reachable_cells == 0) "Zero inside domain: no overlap with fixed occurrence buffer."
             else NA_character_,
    projected = TRUE, binary_path = binary_path)
  save_csv_file(record, status_path)
  products <- c(binary = binary_path, summary = summary_path, status = status_path)
  if (isTRUE(save_suitability) && !is.null(suitability)) {
    path <- file.path(folder, paste0(id, "_suitability.tif"))
    write_raster_file(suitability, path, "FLT4S", species)
    products <- c(products, suitability = path)
  }
  rm(suitability, binary, area, domain, buffer, access)
  invisible(gc())
  products
}
select_local_paths <- function(summary, weights, scenario) {
  if (!all(c("scenario", "species", "binary_path") %in% names(summary))) {
    stop("No usable local summaries for ", scenario)
  }
  selected <- summary[summary$scenario == scenario, , drop = FALSE]
  expected <- as.character(weights$species)
  if (anyDuplicated(selected$species) ||
      !setequal(selected$species, expected)) {
    stop("Species set is incomplete/different for ", scenario,
         ". Do not calculate a metric from a partial projection set.")
  }
  paths <- as.character(selected$binary_path[match(expected, selected$species)])
  assert_paths(paths)
  stats::setNames(paths, expected)
}

# ---- Metric wrappers and audits --------------------------------------------
richness_file <- function(paths, domain_file, filename, chunk_size, memmax_gb) {
  configure_terra_worker(memmax_gb)
  x <- SpeciesPoolR::sum_binary_files_tree(
    paths = paths, predictors = domain_file, filename = filename,
    name = "Richness", chunk_size = chunk_size, overwrite = TRUE)
  rm(x)
  assert_paths(filename)
  filename
}
rarity_file <- function(paths, weights, domain_file, filename, chunk_size, memmax_gb) {
  configure_terra_worker(memmax_gb)
  x <- SpeciesPoolR::calc_range_rarity_files(
    paths = paths, Weights = weights, template = domain_file, filename = filename,
    name = "RangeRarity", output = "cell_contribution", chunk_size = chunk_size,
    align = "error", overwrite = TRUE)
  rm(x)
  assert_paths(filename)
  filename
}
pd_file <- function(paths, tree_file, coverage, domain_file, filename,
                    max_cells, memmax_gb) {
  configure_terra_worker(memmax_gb)
  if (!nrow(coverage) || any(!coverage$matched) ||
      !setequal(coverage$species, names(paths))) {
    stop("Complete tree coverage of the SAME metric species set is required for PD.")
  }
  tree <- readRDS(tree_file)
  x <- SpeciesPoolR::calc_pd_files(
    paths = paths, Tree = tree, template = domain_file, filename = filename,
    name = "PD", include_root = FALSE, unmatched = "error", align = "error",
    max_cells_per_block = max_cells, overwrite = TRUE)
  rm(x)
  assert_paths(filename)
  filename
}
metric_audit <- function(scenario, richness_path, rarity_path, pd_path,
                         summary, weights, domain_file) {
  configure_terra_worker(terra_memmax_gb)
  domain <- terra::rast(domain_file)
  r <- terra::rast(richness_path)
  rr <- terra::rast(rarity_path)
  pd <- terra::rast(pd_path)
  for (x in list(r, rr, pd)) {
    if (!terra::compareGeom(x, domain, stopOnError = FALSE) ||
        count_true(is.na(x) != is.na(domain)) > 0) {
      stop("A metric has wrong geometry or a different domain in ", scenario)
    }
  }
  s <- summary[summary$scenario == scenario, , drop = FALSE]
  if (anyDuplicated(s$species) || !setequal(s$species, weights$species)) {
    stop("Summary and weight species differ in the metric audit.")
  }
  w <- weights$weight[match(s$species, weights$species)]
  expected_rarity <- sum(s$occupied_area_km2 * w)
  observed_rarity <- sum_raster(rr)
  richness_stats <- terra::global(r, c("min", "mean", "max", "sum"), na.rm = TRUE)
  pd_stats <- terra::global(pd, c("min", "mean", "max"), na.rm = TRUE)
  rr_min <- as.numeric(terra::global(rr, "min", na.rm = TRUE)[1, 1])
  difference <- observed_rarity - expected_rarity
  ok <- is.finite(expected_rarity) && is.finite(observed_rarity) &&
    abs(difference) <= max(1e-8, abs(expected_rarity) * 1e-7) &&
    richness_stats$sum == sum(s$occupied_cells) &&
    richness_stats$min >= 0 && richness_stats$max <= nrow(weights) &&
    is.finite(rr_min) && rr_min >= 0 &&
    all(is.finite(as.matrix(pd_stats))) && pd_stats$min >= 0
  if (!isTRUE(ok)) stop("Metric identity/bounds audit failed for ", scenario)
  data.frame(
    scenario = scenario, modelled_species = nrow(weights),
    locally_present_species = sum(s$occupied_area_km2 > 0),
    richness_min = richness_stats$min, richness_mean = richness_stats$mean,
    richness_max = richness_stats$max,
    national_positive_range_species = nrow(weights),
    expected_local_range_rarity = expected_rarity,
    total_local_range_rarity = observed_rarity,
    range_rarity_difference = difference, metric_checks_pass = ok,
    pd_min = pd_stats$min, pd_mean = pd_stats$mean, pd_max = pd_stats$max,
    richness_file = richness_path, range_rarity_file = rarity_path, pd_file = pd_path)
}
workflow_smoke_test <- function() {
  # Deliberately independent of the real inputs: catches API/math errors in seconds.
  folder <- tempfile("gudenaa_plant_test_")
  dir.create(folder)
  template <- terra::rast(nrows = 2, ncols = 2, xmin = 560000, xmax = 560020,
                          ymin = 6220000, ymax = 6220020, crs = "EPSG:25832")
  values(template) <- c(1, 1, 1, NA_real_)
  a <- b <- template
  values(a) <- c(1, 1, 0, NA_real_)
  values(b) <- c(1, 0, 1, NA_real_)
  paths <- c(Species_a = file.path(folder, "a.tif"),
             Species_b = file.path(folder, "b.tif"))
  write_raster_file(a, paths[[1]], "INT1U")
  write_raster_file(b, paths[[2]], "INT1U")
  domain_path <- write_raster_file(template, file.path(folder, "domain.tif"), "INT1U")
  weights <- SpeciesPoolR::calc_range_weights_files(
    paths, template = domain_path, unit = "km2", align = "error", verbose = FALSE)
  rich_path <- richness_file(paths, domain_path, file.path(folder, "rich.tif"), 2L, 1)
  rr_path <- rarity_file(paths, weights, domain_path, file.path(folder, "rr.tif"), 2L, 1)
  tree <- ape::read.tree(text = "(Species_a:1,Species_b:1);")
  tree_path <- save_rds_file(tree, file.path(folder, "tree.rds"))
  coverage <- SpeciesPoolR::audit_tree_coverage(names(paths), tree, group = "plant")
  pd_path <- pd_file(paths, tree_path, coverage, domain_path,
                     file.path(folder, "pd.tif"), 2L, 1)
  rich_values <- as.numeric(terra::values(terra::rast(rich_path)))
  pd_values <- as.numeric(terra::values(terra::rast(pd_path)))
  stopifnot(all(rich_values[1:3] == c(2, 1, 1)), is.na(rich_values[4]),
            all(pd_values[1:3] == c(2, 0, 0)), is.na(pd_values[4]),
            abs(sum_raster(terra::rast(rr_path)) - 2) < 1e-8)
  # Threshold 0 must not create a presence outside the permitted geometry.
  buffer <- terra::as.polygons(terra::ext(560000, 560010, 6220000, 6220020),
                               crs = terra::crs(template))
  suitability <- template
  values(suitability) <- c(0.5, 0.5, 0.5, NA_real_)
  names(suitability) <- "Species_a"
  th <- data.frame(species = "Species_a", Thres_95 = 0)
  result <- binary_with_access(suitability, th, buffer, template, "Species_a")
  got <- as.numeric(terra::values(result))
  stopifnot(all(got[1:3] == c(1, 0, 1)), is.na(got[4]))
  values(suitability) <- c(NA_real_, 0.5, 0.5, NA_real_)
  rejected <- tryCatch({
    binary_with_access(suitability, th, buffer, template, "Species_a")
    FALSE
  }, error = function(e) TRUE)
  stopifnot(rejected)
  # One-row successful/excluded status files must retain a common schema.
  save_csv_file(data.frame(species = "Species_a", range_size = 1,
    usable_for_metrics = TRUE, binary_path = paths[[1]], reason = NA_character_),
    file.path(folder, "success_status.csv"))
  save_csv_file(data.frame(species = "Species_b", range_size = NA_real_,
    usable_for_metrics = FALSE, binary_path = NA_character_, reason = "Excluded"),
    file.path(folder, "excluded_status.csv"))
  status <- read_csv_files(file.path(folder, c("success_status.csv", "excluded_status.csv")))
  stopifnot(is.double(status$range_size), is.logical(status$usable_for_metrics),
            is.character(status$binary_path), identical(status$usable_for_metrics, c(TRUE, FALSE)),
            is.na(status$binary_path[[2]]))
  data.frame(test = c("richness", "area_based_rarity", "PD", "buffer_threshold_zero",
                      "missing_prediction_rejected", "status_file_schema"), passed = TRUE)
}

# ---- Targets ---------------------------------------------------------------
list(
  tar_target(Workflow_smoke_test, { SpeciesPoolR_build; workflow_smoke_test() },
             resources = reduce_resources),

  tar_target(Taxa_file, taxa_path, format = "file"),
  tar_target(Training_manifest_file, file.path(
    predictor_root, "training_current_predictor_paths_Denmark.rds"), format = "file"),
  tar_target(Projection_manifest_file, file.path(
    predictor_root, "projection_current_predictor_paths_gudenaa.rds"), format = "file"),
  tar_target(Scenario_manifest_file, file.path(
    predictor_root, "mixed_scenario_specs_gudenaa.rds"), format = "file"),
  tar_target(Training_paths, read_pair(Training_manifest_file)),
  tar_target(Projection_paths, read_pair(Projection_manifest_file)),
  tar_target(Scenario_specs, readRDS(Scenario_manifest_file)),
  tar_target(Training_files, track_raster_files(Training_paths), format = "file"),
  tar_target(Projection_files, track_raster_files(Projection_paths), format = "file"),
  tar_target(Scenario_files, track_raster_files(
    unique(unlist(Scenario_specs$path, use.names = FALSE))), format = "file"),
  tar_target(Training_domain_file, file.path(
    predictor_root, "Domain_Denmark_training.tif"), format = "file"),
  tar_target(Projection_domain_file, file.path(
    predictor_root, "Domain_gudenaa_projection.tif"), format = "file"),

  tar_target(Input_audit, {
    Workflow_smoke_test; Training_files; Projection_files; Scenario_files
    configure_terra_worker(terra_memmax_gb)
    check_input_set(Training_paths, Projection_paths, Scenario_specs,
                    Training_domain_file, Projection_domain_file)
  }, resources = reduce_resources),
  tar_target(Input_audit_file, save_csv_file(
    Input_audit, file.path(output_root, "audits", "input_audit.csv")), format = "file"),

  tar_target(Taxa_import, read_plant_taxa(Taxa_file, plant_phyla, excluded_species)),
  tar_target(Plant_taxa, Taxa_import$taxa),
  tar_target(Taxon_scope_audit_file, save_csv_file(
    Taxa_import$audit, file.path(output_root, "audits", "taxon_scope.csv")), format = "file"),
  tar_target(Count_Presences, {
    Input_audit; SpeciesPoolR_build
    SpeciesPoolR::count_presences_simple(
      Plant_taxa, country = "DK", year = gbif_years,
      retries = gbif_retries, verbose = TRUE)
  }, resources = gbif_resources),
  tar_target(Plant_selection, select_plants(
    Plant_taxa, Count_Presences, minimum_gbif_count, maximum_species)),
  tar_target(Plant_selection_file, save_csv_file(
    Plant_selection, file.path(output_root, "audits", "plant_selection.csv")), format = "file"),
  tar_target(Species, Plant_selection$species, iteration = "vector"),
  tar_target(Taxon_key, Plant_selection$gbif_speciesKey, iteration = "vector"),
  tar_target(Taxon_genus, Plant_selection$genus, iteration = "vector"),
  tar_target(Taxon_family, Plant_selection$family, iteration = "vector"),

  tar_target(Presence_file, download_plant_presences(
    Species, Taxon_key, Taxon_genus, Taxon_family, output_root,
    gbif_years, gbif_download_limit, gbif_retries),
    pattern = map(Species, Taxon_key, Taxon_genus, Taxon_family),
    iteration = "vector", format = "file", resources = gbif_resources),
  tar_target(Presence_audit, dplyr::bind_rows(lapply(Presence_file, function(path) {
    x <- readRDS(path)
    data.frame(species = x$species, status = x$status, reason = x$reason,
               downloaded_records = nrow(x$records), matched_count = x$matched_count,
               download_limit = x$download_limit, download_capped = x$download_capped,
               target_branch = x$target_branch,
               file = path)
  }))),
  tar_target(Presence_audit_file, save_csv_file(
    Presence_audit, file.path(output_root, "audits", "presence_audit.csv")), format = "file"),

  tar_target(Model_file, {
    Training_files; SpeciesPoolR_build
    fit_plant(Presence_file, Training_paths, output_root, minimum_complete_records,
               background_points, terra_memmax_gb)
  }, pattern = map(Presence_file), iteration = "vector", format = "file"),
  tar_target(Model_audit, dplyr::bind_rows(lapply(
    Model_file, function(path) readRDS(path)$audit))),
  tar_target(Model_audit_file, save_csv_file(
    Model_audit, file.path(output_root, "audits", "model_audit.csv")), format = "file"),
  tar_target(Model_category_support_file, {
    rows <- dplyr::bind_rows(lapply(Model_file, function(path) {
      x <- readRDS(path)
      out <- x$category_support
      if (!nrow(out)) return(NULL)
      out$species <- x$species
      out
    }))
    save_csv_file(rows, file.path(output_root, "audits", "model_landuse_support.csv"))
  }, format = "file"),

  tar_target(Training_cell_area_file, {
    configure_terra_worker(terra_memmax_gb)
    domain <- terra::rast(Training_domain_file)
    write_raster_file(terra::mask(terra::cellSize(domain, unit = "km"), domain),
      file.path(output_root, "intermediates", "cell_area_Denmark_km2.tif"),
      "FLT8S", "cell_area_km2")
  }, format = "file", resources = reduce_resources),
  tar_target(Projection_cell_area_file, {
    configure_terra_worker(terra_memmax_gb)
    domain <- terra::rast(Projection_domain_file)
    write_raster_file(terra::mask(terra::cellSize(domain, unit = "km"), domain),
      file.path(output_root, "intermediates", "cell_area_gudenaa_km2.tif"),
      "FLT8S", "cell_area_km2")
  }, format = "file", resources = reduce_resources),

  tar_target(National_species_products, {
    Training_files; SpeciesPoolR_build
    national_baseline(Model_file, Training_paths, Training_domain_file,
      Training_cell_area_file, output_root, buffer_distance_m,
      projection_blocks_in_memory, terra_memmax_gb)
  }, pattern = map(Model_file), iteration = "list", format = "file"),
  tar_target(National_status, read_csv_files(collect_files(
    National_species_products, "_status.csv"))),
  tar_target(National_audit_file, save_csv_file(
    National_status, file.path(output_root, "audits", "national_baseline_audit.csv")),
    format = "file"),
  tar_target(Range_weights_Denmark, make_range_weights(National_status)),
  tar_target(Range_weights_Denmark_file, save_csv_file(
    Range_weights_Denmark, file.path(output_root, "audits", "range_weights_Denmark.csv")),
    format = "file"),

  tar_target(Current_species_products, {
    Projection_files; SpeciesPoolR_build
    project_local_plant(Model_file, National_species_products, Projection_paths,
      "current", Projection_domain_file, Projection_cell_area_file, output_root,
      projection_blocks_in_memory, terra_memmax_gb, retain_local_suitability)
  }, pattern = map(Model_file, National_species_products), iteration = "list", format = "file"),
  tar_target(Current_summary, read_csv_files(collect_files(
    Current_species_products, "_summary.csv"))),
  tar_target(Current_projection_audit_file, save_csv_file(
    read_csv_files(collect_files(Current_species_products, "_status.csv")),
    file.path(output_root, "audits", "current_projection_audit.csv")), format = "file"),
  tar_target(Current_paths, select_local_paths(
    Current_summary, Range_weights_Denmark, "current")),

  tar_target(Tree_taxa, Plant_selection[
    match(Range_weights_Denmark$species, Plant_selection$species),
    c("species", "genus", "family"), drop = FALSE]),
  tar_target(Plant_tree_file, {
    SpeciesPoolR_build
    tree <- SpeciesPoolR::build_rtrees_tree(
      species = Tree_taxa, group = "plant", tree_index = 1L,
      scenario = "at_basal_node", show_grafted = FALSE)
    save_rds_file(tree, file.path(output_root, "phylogeny", "plant_tree.rds"))
  }, format = "file", resources = tree_resources),
  # Keep the coverage table readable even if a later PD calculation is blocked.
  tar_target(Tree_coverage, {
    SpeciesPoolR_build
    SpeciesPoolR::audit_tree_coverage(
      species = Tree_taxa, Tree = readRDS(Plant_tree_file), group = "plant")
  }),
  tar_target(Tree_coverage_file, save_csv_file(
    Tree_coverage, file.path(output_root, "audits", "tree_coverage.csv")), format = "file"),

  tar_target(Richness_current_file, {
    SpeciesPoolR_build
    richness_file(Current_paths, Projection_domain_file,
      file.path(output_root, "metrics", "Richness_current.tif"),
      metric_chunk_size, terra_memmax_gb)
  }, format = "file", resources = reduce_resources),
  tar_target(RangeRarity_current_file, {
    SpeciesPoolR_build
    rarity_file(Current_paths, Range_weights_Denmark, Projection_domain_file,
      file.path(output_root, "metrics", "RangeRarity_current.tif"),
      metric_chunk_size, terra_memmax_gb)
  }, format = "file", resources = reduce_resources),
  tar_target(PD_current_file, {
    SpeciesPoolR_build
    pd_file(Current_paths, Plant_tree_file, Tree_coverage, Projection_domain_file,
      file.path(output_root, "metrics", "PD_current.tif"),
      pd_max_cells_per_block, terra_memmax_gb)
  }, format = "file", resources = reduce_resources),
  tar_target(Current_metric_audit, metric_audit(
    "current", Richness_current_file, RangeRarity_current_file, PD_current_file,
    Current_summary, Range_weights_Denmark, Projection_domain_file),
    resources = reduce_resources),
  tar_target(Current_metric_audit_file, save_csv_file(
    Current_metric_audit, file.path(output_root, "audits", "current_metric_audit.csv")),
    format = "file"),

  # Branch each scenario's file dependencies BEFORE crossing with species.
  tar_target(Scenario_pair_files, track_raster_files(
    unname(unlist(Scenario_specs$path, use.names = FALSE))),
    pattern = map(Scenario_specs), iteration = "list", format = "file"),
  tar_target(Scenario_species_products, {
    Scenario_pair_files; SpeciesPoolR_build
    project_local_plant(Model_file, National_species_products,
      unname(unlist(Scenario_specs$path, use.names = FALSE)),
      as.character(Scenario_specs$scenario[[1L]]), Projection_domain_file,
      Projection_cell_area_file, output_root, projection_blocks_in_memory,
      terra_memmax_gb, retain_local_suitability)
  }, pattern = cross(map(Model_file, National_species_products),
                      map(Scenario_specs, Scenario_pair_files)),
    iteration = "list", format = "file"),
  tar_target(Scenario_summary, read_csv_files(collect_files(
    Scenario_species_products, "_summary.csv"))),
  tar_target(Scenario_projection_audit_file, save_csv_file(
    read_csv_files(collect_files(Scenario_species_products, "_status.csv")),
    file.path(output_root, "audits", "scenario_projection_audit.csv")), format = "file"),
  tar_target(Scenario_paths, select_local_paths(
    Scenario_summary, Range_weights_Denmark, as.character(Scenario_specs$scenario[[1L]])),
    pattern = map(Scenario_specs), iteration = "list"),

  tar_target(Richness_scenario_file, {
    SpeciesPoolR_build
    richness_file(Scenario_paths, Projection_domain_file,
      file.path(output_root, "metrics", paste0("Richness_", Scenario_specs$scenario[[1]], ".tif")),
      metric_chunk_size, terra_memmax_gb)
  }, pattern = map(Scenario_paths, Scenario_specs), format = "file",
    iteration = "vector", resources = reduce_resources),
  tar_target(RangeRarity_scenario_file, {
    SpeciesPoolR_build
    rarity_file(Scenario_paths, Range_weights_Denmark, Projection_domain_file,
      file.path(output_root, "metrics", paste0("RangeRarity_", Scenario_specs$scenario[[1]], ".tif")),
      metric_chunk_size, terra_memmax_gb)
  }, pattern = map(Scenario_paths, Scenario_specs), format = "file",
    iteration = "vector", resources = reduce_resources),
  tar_target(PD_scenario_file, {
    SpeciesPoolR_build
    pd_file(Scenario_paths, Plant_tree_file, Tree_coverage, Projection_domain_file,
      file.path(output_root, "metrics", paste0("PD_", Scenario_specs$scenario[[1]], ".tif")),
      pd_max_cells_per_block, terra_memmax_gb)
  }, pattern = map(Scenario_paths, Scenario_specs), format = "file",
    iteration = "vector", resources = reduce_resources),
  tar_target(Scenario_metric_audit, metric_audit(
    as.character(Scenario_specs$scenario[[1]]),
    Richness_scenario_file, RangeRarity_scenario_file, PD_scenario_file,
    Scenario_summary, Range_weights_Denmark, Projection_domain_file),
    pattern = map(Scenario_specs, Richness_scenario_file, RangeRarity_scenario_file,
                  PD_scenario_file), iteration = "list", resources = reduce_resources),

  tar_target(Metric_audit_file, save_csv_file(
    dplyr::bind_rows(c(list(Current_metric_audit), Scenario_metric_audit)),
    file.path(output_root, "audits", "all_metric_audits.csv")), format = "file"),
  tar_target(Projection_summary_file, save_csv_file(
    dplyr::bind_rows(Current_summary, Scenario_summary),
    file.path(output_root, "audits", "all_projection_summaries.csv")), format = "file"),
  tar_target(Run_configuration_file, save_rds_file(list(
    scope = plant_phyla, maximum_species = maximum_species, run_label = run_label,
    excluded_species = excluded_species, minimum_complete_records = minimum_complete_records,
    gbif_years = gbif_years, gbif_download_limit = gbif_download_limit,
    background_points = background_points, buffer_distance_m = buffer_distance_m,
    threshold = "Thres_95", range_reference = "national current binary area, km2",
    rarity_output = "cell_contribution", PD_include_root = FALSE,
    scientific_note = paste("Coarse CHELSA-grid training, fine local projection;",
                            "interpolation does not create fine climate information.",
                            "Tree branch-length units are not assumed to be years.",
                            "Grafted tips do not represent independently observed phylogenies."),
    package_build = SpeciesPoolR_build,
    workers = c(compute = compute_workers, gbif = gbif_workers,
                reduce = reduce_workers, tree = tree_workers),
    training_paths = Training_paths, projection_paths = Projection_paths,
    scenario_specs = Scenario_specs
  ), file.path(output_root, "run_configuration.rds")), format = "file")
)
