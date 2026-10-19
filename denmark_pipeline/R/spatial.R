dp_raster <- function(x, filename, datatype = "FLT8S") {
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  terra::writeRaster(x, filename, overwrite = TRUE, datatype = datatype,
    gdal = c("COMPRESS=DEFLATE", "TILED=YES", "BIGTIFF=IF_SAFER"))
}
dp_vector <- function(path, crs) {
  v <- terra::makeValid(terra::vect(path))
  if (!nzchar(terra::crs(v))) stop("Missing vector CRS: ", path)
  if (!terra::same.crs(v, crs)) v <- terra::project(v, crs)
  v
}
dp_mask_polygon <- function(v, template) {
  terra::rasterize(v, template, field = 1, background = 0, touches = FALSE) > 0
}
dp_aligned <- function(path, template) {
  x <- terra::rast(path)
  if (terra::nlyr(x) != 1L || !terra::same.crs(x, template) ||
      any(abs(terra::res(x) - terra::res(template)) > 1e-8) ||
      any(abs((terra::origin(x) - terra::origin(template)) / terra::res(template) -
              round((terra::origin(x) - terra::origin(template)) / terra::res(template))) > 1e-8))
    stop("Raster must use the same 10 m CRS and alignment: ", path)
  a <- as.vector(terra::ext(x)); b <- as.vector(terra::ext(template))
  if (a[2] <= b[1] || a[1] >= b[2] || a[4] <= b[3] || a[3] >= b[4])
    return(terra::init(template, NA_real_))
  x <- terra::crop(x, terra::ext(template), snap = "near")
  x <- terra::extend(x, template)
  terra::compareGeom(x, template)
  x
}
dp_codes <- function(lu, classes) {
  rat <- terra::cats(lu)[[1]]
  if (is.null(rat) || ncol(rat) < 2L || anyDuplicated(rat[[1]]) || anyDuplicated(rat[[2]]))
    stop("Landuse needs the fitted model's id/class category table.")
  if (any(!classes %in% rat[[2]])) stop("Unknown land-use classes: ", paste(setdiff(classes, rat[[2]]), collapse = ", "))
  as.numeric(rat[[1]][match(classes, rat[[2]])])
}
dp_indicator <- function(lu, codes) {
  raw <- as.numeric(lu)
  terra::ifel(!is.na(raw) & raw %in% codes, 1, 0)
}
dp_source <- function(s, lu, nature) {
  if (s$type == "nature") return(nature)
  if (s$type == "raster") {
    a <- as.numeric(dp_aligned(s$path, lu))
    a <- if (!is.null(s$values)) a %in% s$values else a > 0
    a <- terra::ifel(!is.na(a) & a, 1, 0)
  } else if (s$type == "vector") {
    v <- dp_vector(s$path, terra::crs(lu))
    if (!is.null(s$filter_field)) v <- v[as.data.frame(v)[[s$filter_field]] %in% s$filter_values, ]
    a <- if (nrow(v)) dp_mask_polygon(v, lu) else terra::init(lu, 0)
  } else stop("Anchor type must be nature, raster, or vector.")
  a
}
dp_layers <- function(cfg, lu, region, folder) {
  nature <- dp_indicator(lu, dp_codes(lu, dp_pathways()$class))
  fixed <- nature; budget <- nature; credit <- nature
  roots <- list()
  for (nm in names(cfg$anchors$sources)) {
    s <- cfg$anchors$sources[[nm]]
    a <- dp_source(s, lu, nature)
    # A shoreline buffer marks land eligible to root a network; it does not add
    # the buffered area to existing habitat or to an area ledger.
    b <- s$buffer_m %or% 0
    if (b < 0 || b %% cfg$grid$resolution_m != 0) stop("Anchor buffer must be a nonnegative multiple of 10 m.")
    if (b > 0) {
      offsets <- seq.int(-as.integer(b/cfg$grid$resolution_m), as.integer(b/cfg$grid$resolution_m)) * cfg$grid$resolution_m
      window <- ifelse(outer(offsets^2, offsets^2, "+") <= b^2, 1, NA_real_)
      roots[[nm]] <- terra::focal(a, w = window, fun = "max", na.rm = TRUE, fillvalue = 0)
    } else roots[[nm]] <- a
    if (isTRUE(s$include_as_existing)) {
      fixed <- fixed | a
      if (isTRUE(s$credit_budget)) budget <- budget | a
      if (isTRUE(s$credit_network)) credit <- credit | a
    } else if (isTRUE(s$credit_budget) || isTRUE(s$credit_network))
      stop("Area credit requires include_as_existing = TRUE for source ", nm)
  }
  anchor <- if (!length(roots)) terra::init(lu, 0) else {
    z <- roots[[1]]
    if (length(roots) > 1L) for (i in 2:length(roots))
      z <- if (cfg$anchors$combine == "union") z | roots[[i]] else z & roots[[i]]
    z
  }
  eligible <- dp_indicator(lu, dp_codes(lu, cfg$eligibility$classes)) & !fixed & region
  if (!is.null(cfg$eligibility$mask)) {
    a <- dp_aligned(cfg$eligibility$mask, lu)
    eligible <- eligible & terra::ifel(!is.na(a) & a > 0, 1, 0)
  }
  out <- list(nature = nature, fixed = fixed, budget = budget, credit = credit,
              anchors = anchor, eligible = eligible)
  for (nm in names(out)) out[[nm]] <- dp_raster(out[[nm]], file.path(folder, paste0(nm, "_10m.tif")), "INT1U")
  if (length(roots)) {
    bands <- terra::rast(roots); names(bands) <- names(roots)
    dp_raster(bands, file.path(folder, "anchor_sources_10m.tif"), "INT1U")
  }
  out
}
build_anchors <- function(config, output_dir = NULL) {
  cfg <- dp_config(config); old <- setwd(cfg$project_root); on.exit(setwd(old), add = TRUE)
  folder <- output_dir %or% file.path(cfg$output_root, cfg$area_id, "anchor_layers")
  t <- dp_template(cfg, cfg$connectivity$context_boundary %or% cfg$region$boundary)
  b <- dp_vector(cfg$region$boundary, terra::crs(t$lu))
  dp_layers(cfg, t$lu, dp_mask_polygon(b, t$lu), folder)
  invisible(normalizePath(file.path(folder, "anchors_10m.tif"), winslash = "/"))
}
dp_zonal <- function(x, zones, fun = "sum") {
  a <- as.data.frame(terra::zonal(x, zones, fun, na.rm = TRUE))
  names(a)[1] <- "id"; a
}
dp_region_inputs <- function(cfg) c(cfg$grid$landuse, cfg$region$boundary,
  cfg$region$subareas, cfg$connectivity$context_boundary, cfg$eligibility$mask,
  unlist(lapply(cfg$anchors$sources, `[[`, "path")))
dp_input_files <- function(cfg) dp_files(c(dp_region_inputs(cfg), cfg$budget$cost_per_ha,
  cfg$model$current_scores, unlist(lapply(cfg$model$outcomes, function(o) c(o$scores, o$allowed_mask)))))
dp_template <- function(cfg, boundary) {
  full <- terra::rast(cfg$grid$landuse)
  if (terra::is.lonlat(full) || any(abs(terra::res(full) - cfg$grid$resolution_m) > 1e-8))
    stop("Use a projected 10 m land-use raster with metre units (e.g. EPSG:25832).")
  if (!grepl('LENGTHUNIT["metre",1', terra::crs(full), fixed = TRUE) &&
      !grepl('+units=m', terra::crs(full, proj = TRUE), fixed = TRUE)) stop("Grid units must be metres.")
  b <- dp_vector(boundary, terra::crs(full))
  e <- as.vector(terra::ext(b)); r <- cfg$grid$planning_unit_m; o <- cfg$grid$origin
  e <- c(floor((e[1]-o[1])/r)*r+o[1], ceiling((e[2]-o[1])/r)*r+o[1],
         floor((e[3]-o[2])/r)*r+o[2], ceiling((e[4]-o[2])/r)*r+o[2])
  coarse <- terra::rast(terra::ext(e), resolution = r, crs = terra::crs(full))
  fine <- terra::rast(terra::ext(e), resolution = cfg$grid$resolution_m, crs = terra::crs(full))
  lu <- dp_aligned(cfg$grid$landuse, fine)
  list(lu = lu, coarse = coarse, boundary = b)
}
prepare_area <- function(config, folder) {
  cfg <- dp_config(config); dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  t <- dp_template(cfg, cfg$connectivity$context_boundary %or% cfg$region$boundary)
  lu <- t$lu; coarse <- t$coarse
  boundary <- dp_vector(cfg$region$boundary, terra::crs(lu))
  region <- dp_mask_polygon(boundary, lu)
  if (!is.null(cfg$connectivity$context_boundary)) {
    if (sum(terra::expanse(terra::erase(boundary, t$boundary), unit = "ha", transform = FALSE)) > 1e-6)
      stop("Connectivity context must contain the entire planning region.")
  }
  context <- dp_mask_polygon(t$boundary, lu)
  if (dp_any(region & is.na(lu))) stop("The planning region has missing Landuse pixels; fix coverage first.")
  dp_raster(region, file.path(folder, "region_10m.tif"), "INT1U")
  dp_raster(context, file.path(folder, "context_10m.tif"), "INT1U")
  a <- dp_layers(cfg, lu, region, folder)
  a$fixed <- a$fixed & context; a$credit <- a$credit & context
  unit <- terra::resample(terra::init(coarse, "cell"), lu, method = "near")
  unit <- terra::ifel(a$eligible, unit, NA)
  ha <- prod(terra::res(lu)) / 10000
  areas <- dp_zonal(a$eligible * ha, unit); names(areas) <- c("unit", "area_ha")
  units <- areas[areas$area_ha >= cfg$eligibility$min_unit_ha & areas$area_ha > 0, , drop = FALSE]
  if (!nrow(units)) stop("No selectable planning units remain.")
  unit <- terra::ifel(unit %in% units$unit, unit, NA)
  eligible <- !is.na(unit)
  unit_path <- file.path(folder, "unit_ids_10m.tif"); dp_raster(unit, unit_path)
  dp_raster(eligible, file.path(folder, "selectable_10m.tif"), "INT1U")
  forest <- dp_zonal(dp_indicator(lu, dp_codes(lu, "ProductionForest")) * ha, unit)
  units$forest_share <- forest[[2]][match(units$unit, forest$id)] / units$area_ha
  units$cost <- units$area_ha
  if (!is.null(cfg$budget$cost_per_ha)) {
    costs <- dp_aligned(cfg$budget$cost_per_ha, lu)
    if (dp_any(eligible & (!is.finite(costs) | costs < 0))) stop("Missing or negative candidate restoration costs.")
    z <- dp_zonal(costs * ha, unit); units$cost <- z[[2]][match(units$unit, z$id)]
  }
  scales <- cfg$model$metric_scales
  if (is.null(scales)) {
    scales <- setNames(vapply(names(cfg$model$metric_weights), function(m) {
      paths <- c(vapply(cfg$model$outcomes, function(o) o$scores[[m]], character(1)), cfg$model$current_scores[m])
      paths <- paths[!is.na(paths)]
      max(vapply(paths, function(p) dp_max(abs(dp_aligned(p, lu))), numeric(1)), 1e-12)
    }, numeric(1)), names(cfg$model$metric_weights))
  }
  if (!setequal(names(scales), names(cfg$model$metric_weights)) || any(!is.finite(scales) | scales <= 0)) stop("Invalid metric scales.")
  actions <- list()
  for (nm in names(cfg$model$outcomes)) {
    spec <- cfg$model$outcomes[[nm]]; score <- terra::init(lu, 0)
    for (m in names(cfg$model$metric_weights)) {
      z <- dp_aligned(spec$scores[[m]], lu)
      if (dp_any(eligible & !is.finite(z))) stop("Missing candidate score: ", nm, " / ", m)
      if (cfg$model$score_type == "gain") {
        base <- dp_aligned(cfg$model$current_scores[[m]], lu)
        if (dp_any(eligible & !is.finite(base))) stop("Missing baseline scores.")
        z <- z - base
      }
      score <- score + z * cfg$model$metric_weights[[m]] / scales[[m]]
    }
    z <- dp_zonal(score * ha, unit)
    keep <- rep(TRUE, nrow(units))
    if (isTRUE(cfg$eligibility$forest_only) && !grepl("^Forest", spec$class %or% ""))
      keep <- units$forest_share < cfg$eligibility$forest_share
    if (!is.null(spec$allowed_mask)) {
      allowed <- dp_aligned(spec$allowed_mask, lu)
      bad <- dp_zonal(terra::ifel(is.na(allowed) | allowed <= 0, 1, 0), unit)
      keep <- keep & bad[[2]][match(units$unit, bad$id)] == 0
    }
    actions[[nm]] <- data.frame(outcome = nm, unit = units$unit[keep],
      score = z[[2]][match(units$unit[keep], z$id)])
  }
  actions <- do.call(rbind, actions); rownames(actions) <- NULL
  # Units with no scientifically allowed outcome are removed from the domain.
  units <- units[units$unit %in% actions$unit, , drop = FALSE]
  if (!nrow(units)) stop("All units were excluded by the outcome feasibility rules.")
  unit <- terra::ifel(unit %in% units$unit, unit, NA); eligible <- !is.na(unit)
  dp_raster(unit, unit_path); dp_raster(eligible, file.path(folder, "selectable_10m.tif"), "INT1U")
  polygon_ha <- sum(terra::expanse(terra::aggregate(boundary), unit = "ha", transform = FALSE))
  current_ha <- dp_sum(a$budget & region) * ha
  target <- if (is.null(cfg$budget$target_nature_fraction)) 0 else cfg$budget$target_nature_fraction * polygon_ha
  required <- max(0, target - current_ha)
  minimum <- max(cfg$budget$min_restoration_ha, required)
  maximum <- cfg$budget$max_restoration_ha %or% sum(units$area_ha)
  if (!is.null(cfg$budget$target_nature_fraction) && !is.null(cfg$budget$target_tolerance_ha))
    maximum <- min(maximum, required + cfg$budget$target_tolerance_ha)
  if (minimum > maximum + 1e-8 || minimum > sum(units$area_ha) + 1e-8) stop("Insufficient eligible area or conflicting budget bounds.")
  ledger <- data.frame(polygon_ha, represented_ha = dp_sum(region) * ha, current_ha, target,
    minimum_restoration_ha = minimum, maximum_restoration_ha = maximum,
    eligible_ha = sum(units$area_ha), cost_cap = cfg$budget$cost_cap %or% sum(units$cost))
  goals <- list()
  if (!is.null(cfg$connectivity$global_min_ha) && cfg$connectivity$global_min_ha > 0)
    goals$Global <- list(mask = context,
      credit_mask = if (cfg$connectivity$global_area_scope == "within_region") region else context,
      minimum = cfg$connectivity$global_min_ha)
  if (!is.null(cfg$region$subareas)) {
    s <- dp_vector(cfg$region$subareas, terra::crs(lu))
    labels <- if (is.null(cfg$region$subarea_id)) as.character(seq_len(nrow(s))) else as.character(as.data.frame(s)[[cfg$region$subarea_id]])
    if (length(labels) != nrow(s) || anyNA(labels) || anyDuplicated(labels)) stop("Subarea IDs must be unique.")
    if (sum(terra::expanse(s, unit = "ha", transform = FALSE)) -
        sum(terra::expanse(terra::aggregate(s), unit = "ha", transform = FALSE)) > 1e-5) stop("Subareas overlap.")
    if (sum(terra::expanse(terra::erase(s, boundary), unit = "ha", transform = FALSE)) > 1e-5) stop("Subareas extend outside region.")
    sid <- terra::rasterize(s, lu, field = seq_len(nrow(s)), background = NA, touches = FALSE)
    dp_raster(sid, file.path(folder, "subareas_10m.tif"))
    dp_csv(data.frame(subarea = seq_along(labels), label = labels), file.path(folder, "subareas.csv"))
    mins <- cfg$connectivity$subarea_min_ha
    if (!is.null(mins)) for (i in seq_len(nrow(s))) {
      m <- if (length(mins) == 1L && is.null(names(mins))) mins else mins[labels[i]]
      if (length(m) != 1L || !is.finite(m) || m < 0) stop("Supply a minimum for each subarea ID.")
      if (m > 0) goals[[paste0("Subarea_", i)]] <- list(mask = dp_mask_polygon(s[i, ], lu),
        credit_mask = dp_mask_polygon(s[i, ], lu), minimum = unname(m))
    }
  }
  graphs <- lapply(names(goals), function(k) dp_graph(k, goals[[k]], unit, a, cfg, folder))
  names(graphs) <- names(goals)
  # The same full-context graph is also used for physically touching unit pairs.
  spatial <- dp_graph("Spatial", list(mask = context, credit_mask = region, minimum = 0), unit, a, cfg, folder)
  dp_csv(units, file.path(folder, "candidate_units.csv")); dp_csv(actions, file.path(folder, "actions.csv"))
  dp_csv(ledger, file.path(folder, "area_budget_audit.csv")); dp_save(scales, file.path(folder, "metric_scales.rds"))
  spec <- list(crs = terra::crs(lu), resolution = terra::res(lu), origin = cfg$grid$origin,
    planning_unit_m = cfg$grid$planning_unit_m, extent = as.vector(terra::ext(lu)),
    selectable_md5 = unname(tools::md5sum(file.path(folder, "selectable_10m.tif"))), units = units[c("unit", "area_ha")])
  dp_write(digest::digest(spec), file.path(folder, "grid_id.txt")); dp_save(spec, file.path(folder, "grid_specification.rds"))
  prepared <- list(units = units, actions = actions, ledger = ledger, graphs = graphs, spatial = spatial,
    folder = normalizePath(folder, winslash = "/"), unit_path = normalizePath(unit_path, winslash = "/"),
    metric_scales = scales)
  dp_save(prepared, file.path(folder, "prepared.rds")); prepared
}
