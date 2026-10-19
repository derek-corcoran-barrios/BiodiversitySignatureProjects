# Projection functions extracted unchanged from targets_plant_biodiversity.R, 2026-09-02.
configure_terra_worker <- function(memmax_gb) {
    settings <- list(progress = 0L, memmax = memmax_gb, memmin = 0, memfrac = 0.02)
    available <- names(terra::terraOptions(print = FALSE))
    if ("threads" %in% available) 
        settings$threads <- 1L
    do.call(terra::terraOptions, settings)
    invisible(NULL)
}

assert_paths <- function(paths) {
    if (!is.character(paths) || !length(paths) || anyNA(paths) || any(!nzchar(paths)) || any(!file.exists(paths))) {
        stop("Missing or invalid file paths: ", paste(paths, collapse = ", "), call. = FALSE)
    }
    invisible(paths)
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
    if (!length(paths)) 
        return(data.frame())
    out <- dplyr::bind_rows(lapply(paths, function(path) readr::read_csv(path, col_types = readr::cols(.default = readr::col_character()), 
        show_col_types = FALSE, progress = FALSE)))
    numeric_fields <- intersect(names(out), c("range_size", "weight", "complete_presence_records", "occupied_cells", 
        "occupied_area_km2", "accessible_cells", "national_range_km2", "fixed_national_weight"))
    logical_fields <- intersect(names(out), c("usable_for_metrics", "projected", "no_buffer_overlap"))
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
    as.numeric(terra::global(terra::ifel(!is.na(x) & x, 1, 0), "sum", na.rm = TRUE)[1, 1])
}

sum_raster <- function(x) {
    as.numeric(terra::global(x, "sum", na.rm = TRUE)[1, 1])
}

write_raster_file <- function(x, filename, datatype = "FLT8S", layer_name = NULL) {
    if (!is.null(layer_name)) 
        names(x) <- layer_name
    dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
    wopt <- list(datatype = datatype, gdal = c("COMPRESS=DEFLATE", "TILED=YES", "BIGTIFF=IF_SAFER"))
    if (datatype == "INT1U") 
        wopt$NAflag <- 255
    if (datatype == "INT2U") 
        wopt$NAflag <- 65535
    terra::writeRaster(x, filename, overwrite = TRUE, wopt = wopt)
    filename
}

assert_real_prediction <- function(x, species) {
    if (any(grepl("SpeciesPoolR_zero_suitability_", basename(terra::sources(x)), fixed = TRUE))) {
        stop("Prediction failed for ", species, ": SpeciesPoolR returned its zero-raster error fallback.", 
            call. = FALSE)
    }
    invisible(x)
}

temporary_sources <- function(...) {
    paths <- unique(unlist(lapply(list(...), function(x) {
        if (inherits(x, "SpatRaster")) 
            as.character(terra::sources(x))
        else character()
    }), use.names = FALSE))
    paths <- paths[!is.na(paths) & nzchar(paths)]
    roots <- unique(c(tempdir(), terra::terraOptions(print = FALSE)$tempdir))
    roots <- roots[!is.na(roots) & nzchar(roots)]
    norm <- function(x) {
        out <- normalizePath(x, winslash = "/", mustWork = FALSE)
        if (.Platform$OS.type == "windows") 
            tolower(out)
        else out
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

binary_with_access <- function(suitability, thresholds, buffer, domain, species) {
    stopifnot(terra::compareGeom(suitability, domain))
    access <- terra::mask(domain, buffer)
    missing <- count_true(!is.na(access) & !is.finite(suitability))
    if (missing > 0) {
        stop(species, ": ", missing, " accessible cells have no finite prediction. ", "Check unsupported Landuse levels and predictor coverage; these are not absences.", 
            call. = FALSE)
    }
    th <- thresholds$Thres_95[match(species, thresholds$species)]
    if (length(th) != 1L || !is.finite(th)) 
        stop("No finite national threshold for ", species)
    names(suitability) <- species
    binary <- SpeciesPoolR::threshold_suitability(Model = suitability, Thresholds = thresholds, threshold = "Thres_95")
    binary <- terra::mask(binary, buffer, updatevalue = 0)
    binary <- terra::mask(binary, domain)
    if (count_true(!is.na(domain) & is.na(binary)) > 0) {
        stop("Binary projection has holes inside the valid domain for ", species)
    }
    names(binary) <- species
    binary
}

read_national_record <- function(products) {
    status <- collect_files(products, "_status.csv")
    if (length(status) != 1L) 
        stop("Expected exactly one national species status file.")
    read_csv_files(status)
}

project_local_plant <- function(model_file, national_products, predictors, scenario, domain_file, cell_area_file, 
    root, blocks, memmax_gb, save_suitability = FALSE) {
    configure_terra_worker(memmax_gb)
    baseline <- read_national_record(national_products)
    species <- baseline$species[[1L]]
    folder <- file.path(root, "local_projections", scenario)
    dir.create(folder, recursive = TRUE, showWarnings = FALSE)
    id <- species_id(species)
    status_path <- file.path(folder, paste0(id, "_status.csv"))
    if (!isTRUE(baseline$usable_for_metrics[[1L]])) {
        record <- data.frame(scenario = scenario, species = species, status = "excluded_from_national_baseline", 
            reason = paste(baseline$status[[1]], baseline$reason[[1]], sep = ": "), projected = FALSE, 
            binary_path = NA_character_)
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
        binary <- terra::ifel(!is.na(domain), 0, NA)
        names(binary) <- species
    }
    else {
        prediction_env <- terra::mask(terra::rast(predictors), buffer)
        scratch <- unique(c(scratch, temporary_sources(prediction_env)))
        suitability <- SpeciesPoolR::PredictSuitability(Models = fit$models, file = prediction_env, categorical = "Landuse", 
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
    summary <- data.frame(scenario = scenario, species = species, occupied_cells = sum_raster(binary), 
        occupied_area_km2 = sum_raster(area), accessible_cells = reachable_cells, national_range_km2 = baseline$range_size[[1]], 
        fixed_national_weight = baseline$weight[[1]], binary_path = binary_path, no_buffer_overlap = reachable_cells == 
            0)
    if (any(!is.finite(summary$occupied_area_km2)) || summary$occupied_area_km2 < 0) {
        stop("Invalid local occupied area for ", species, " / ", scenario)
    }
    save_csv_file(summary, summary_path)
    record <- data.frame(scenario = scenario, species = species, status = "projected", reason = if (reachable_cells == 
        0) 
        "Zero inside domain: no overlap with fixed occurrence buffer."
    else NA_character_, projected = TRUE, binary_path = binary_path)
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
    if (anyDuplicated(selected$species) || !setequal(selected$species, expected)) {
        stop("Species set is incomplete/different for ", scenario, ". Do not calculate a metric from a partial projection set.")
    }
    paths <- as.character(selected$binary_path[match(expected, selected$species)])
    assert_paths(paths)
    stats::setNames(paths, expected)
}

richness_file <- function(paths, domain_file, filename, chunk_size, memmax_gb) {
    configure_terra_worker(memmax_gb)
    x <- SpeciesPoolR::sum_binary_files_tree(paths = paths, predictors = domain_file, filename = filename, 
        name = "Richness", chunk_size = chunk_size, overwrite = TRUE)
    rm(x)
    assert_paths(filename)
    filename
}

rarity_file <- function(paths, weights, domain_file, filename, chunk_size, memmax_gb) {
    configure_terra_worker(memmax_gb)
    x <- SpeciesPoolR::calc_range_rarity_files(paths = paths, Weights = weights, template = domain_file, 
        filename = filename, name = "RangeRarity", output = "cell_contribution", chunk_size = chunk_size, 
        align = "error", overwrite = TRUE)
    rm(x)
    assert_paths(filename)
    filename
}

pd_file <- function(paths, tree_file, coverage, domain_file, filename, max_cells, memmax_gb) {
    configure_terra_worker(memmax_gb)
    if (!nrow(coverage) || any(!coverage$matched) || !setequal(coverage$species, names(paths))) {
        stop("Complete tree coverage of the SAME metric species set is required for PD.")
    }
    tree <- readRDS(tree_file)
    x <- SpeciesPoolR::calc_pd_files(paths = paths, Tree = tree, template = domain_file, filename = filename, 
        name = "PD", include_root = FALSE, unmatched = "error", align = "error", max_cells_per_block = max_cells, 
        overwrite = TRUE)
    rm(x)
    assert_paths(filename)
    filename
}

metric_audit <- function(scenario, richness_path, rarity_path, pd_path, summary, weights, domain_file) {
    configure_terra_worker(terra_memmax_gb)
    domain <- terra::rast(domain_file)
    r <- terra::rast(richness_path)
    rr <- terra::rast(rarity_path)
    pd <- terra::rast(pd_path)
    for (x in list(r, rr, pd)) {
        if (!terra::compareGeom(x, domain, stopOnError = FALSE) || count_true(is.na(x) != is.na(domain)) > 
            0) {
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
    ok <- is.finite(expected_rarity) && is.finite(observed_rarity) && abs(difference) <= max(1e-08, abs(expected_rarity) * 
        1e-07) && richness_stats$sum == sum(s$occupied_cells) && richness_stats$min >= 0 && richness_stats$max <= 
        nrow(weights) && is.finite(rr_min) && rr_min >= 0 && all(is.finite(as.matrix(pd_stats))) && pd_stats$min >= 
        0
    if (!isTRUE(ok)) 
        stop("Metric identity/bounds audit failed for ", scenario)
    data.frame(scenario = scenario, modelled_species = nrow(weights), locally_present_species = sum(s$occupied_area_km2 > 
        0), richness_min = richness_stats$min, richness_mean = richness_stats$mean, richness_max = richness_stats$max, 
        national_positive_range_species = nrow(weights), expected_local_range_rarity = expected_rarity, 
        total_local_range_rarity = observed_rarity, range_rarity_difference = difference, metric_checks_pass = ok, 
        pd_min = pd_stats$min, pd_mean = pd_stats$mean, pd_max = pd_stats$max, richness_file = richness_path, 
        range_rarity_file = rarity_path, pd_file = pd_path)
}

