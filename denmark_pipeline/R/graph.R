# Graph nodes are connected 10 m fragments. Several fragments may share one
# 200 m decision, but that decision never creates a connection between them.
dp_edges <- function(r, directions = 4L, block_rows = 128L) {
  terra::readStart(r); on.exit(terra::readStop(r), add = TRUE)
  parts <- list(); nc <- ncol(r); nr <- nrow(r); j <- 0L
  add <- function(a, b) {
    ok <- !is.na(a) & !is.na(b) & a != b
    if (any(ok)) {
      j <<- j + 1L
      parts[[j]] <<- unique(data.frame(from = pmin(a[ok], b[ok]), to = pmax(a[ok], b[ok])))
    }
  }
  for (row in seq.int(1L, nr, by = block_rows)) {
    n <- min(block_rows + 1L, nr - row + 1L)
    m <- matrix(terra::readValues(r, row = row, nrows = n, mat = FALSE), nrow = n, ncol = nc, byrow = TRUE)
    take <- seq_len(min(block_rows, n))
    if (nc > 1) add(as.vector(m[take, -nc, drop = FALSE]), as.vector(m[take, -1, drop = FALSE]))
    if (n > 1) {
      add(as.vector(m[-n, , drop = FALSE]), as.vector(m[-1, , drop = FALSE]))
      if (directions == 8L && nc > 1) {
        add(as.vector(m[-n, -nc, drop = FALSE]), as.vector(m[-1, -1, drop = FALSE]))
        add(as.vector(m[-n, -1, drop = FALSE]), as.vector(m[-1, -nc, drop = FALSE]))
      }
    }
  }
  if (!length(parts)) return(data.frame(from = numeric(), to = numeric()))
  unique(do.call(rbind, parts))
}
dp_value_patches <- function(x, directions) {
  if (!dp_any(!is.na(x))) return(x)
  result <- tryCatch(terra::patches(x, directions = directions, values = TRUE, allowGaps = TRUE),
    error = function(e) {
      if (!grepl("unknown option.*values", conditionMessage(e))) stop(e)
      if (directions != 4L) stop("Eight-neighbour value patches require a current terra version.")
      # Exact rook-connectivity fallback for terra 1.7. New terra versions use
      # the faster file-backed values=TRUE implementation above.
      v <- terra::disagg(terra::as.polygons(x, dissolve = TRUE, values = TRUE, na.rm = TRUE))
      v$node_id <- seq_len(nrow(v))
      terra::rasterize(v, x, field = "node_id", background = NA, touches = FALSE)
    })
  result
}
dp_graph <- function(key, goal, unit, layers, cfg, folder) {
  zones <- terra::ifel(goal$mask, terra::ifel(layers$fixed, 1, unit + 2), NA)
  p <- dp_value_patches(zones, cfg$grid$directions)
  path <- file.path(folder, paste0("nodes_", key, "_10m.tif")); p <- dp_raster(p, path)
  owner <- terra::ifel(layers$fixed | is.na(unit), 0, unit)
  credit <- terra::ifel(layers$fixed, layers$credit & goal$credit_mask, !is.na(unit)) *
    prod(terra::res(unit)) / 10000
  stats <- c(owner, layers$anchors, credit)
  names(stats) <- c("unit", "anchor", "area_ha")
  z <- dp_zonal(stats[[1:2]], p, "max")
  q <- dp_zonal(stats[[3]], p)
  z$area_ha <- q[[2]][match(z$id, q$id)]
  z$unit <- as.numeric(z$unit); z$anchor <- as.integer(z$anchor > 0)
  edges <- dp_edges(p, cfg$grid$directions)
  if (nrow(z)) {
    g <- igraph::graph_from_data_frame(edges, directed = FALSE, vertices = as.character(z$id))
    groups <- igraph::components(g)$membership
    z$component <- unname(groups[as.character(z$id)])
    z$flow_cap <- as.numeric(table(z$component)[as.character(z$component)])
    potential <- aggregate(z$area_ha, list(component = z$component), sum); names(potential)[2] <- "area_ha"
    potential$has_anchor <- vapply(potential$component, function(k) any(z$anchor[z$component == k] > 0), logical(1))
    potential$has_restoration <- vapply(potential$component, function(k) any(z$unit[z$component == k] > 0), logical(1))
    dp_csv(potential, file.path(folder, paste0("potential_", key, ".csv")))
    if (key != "Spatial") {
      keep <- potential$component[potential$has_anchor &
        (!cfg$connectivity$require_new_restoration | potential$has_restoration)]
      z <- z[z$component %in% keep, , drop = FALSE]
      edges <- edges[edges$from %in% z$id & edges$to %in% z$id, , drop = FALSE]
    }
  } else { z$component <- numeric(); z$flow_cap <- numeric() }
  list(key = key, minimum = goal$minimum, nodes = z, edges = edges, raster = normalizePath(path, winslash = "/"))
}
dp_spatial_terms <- function(prepared) {
  g <- prepared$spatial; n <- g$nodes; e <- g$edges
  u <- n$unit[match(e$from, n$id)]; v <- n$unit[match(e$to, n$id)]
  valid <- u > 0 & v > 0 & u != v
  pairs <- unique(data.frame(from = pmin(u[valid], v[valid]), to = pmax(u[valid], v[valid])))
  contact <- setNames(rep(0, nrow(prepared$units)), prepared$units$unit)
  # Root-capable selected fragments also receive a soft anchoring preference.
  hit <- unique(c(n$unit[n$unit > 0 & n$anchor > 0],
    u[u > 0 & v == 0 & n$anchor[match(e$to, n$id)] > 0],
    v[v > 0 & u == 0 & n$anchor[match(e$from, n$id)] > 0]))
  contact[as.character(hit)] <- 1
  list(pairs = pairs, contact = contact)
}
audit_solution <- function(cfg, prepared, selected, output_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  if (anyDuplicated(selected$unit) || any(!selected$unit %in% prepared$units$unit)) stop("Invalid selected-unit export.")
  keys <- paste(prepared$actions$outcome, prepared$actions$unit)
  if (any(!paste(selected$outcome, selected$unit) %in% keys)) stop("Export includes an illegal outcome.")
  idx <- match(selected$unit, prepared$units$unit)
  area <- sum(prepared$units$area_ha[idx]); cost <- sum(prepared$units$cost[idx]); ledger <- prepared$ledger
  if (area < ledger$minimum_restoration_ha - 1e-5 || area > ledger$maximum_restoration_ha + 1e-5 || cost > ledger$cost_cap + 1e-5)
    stop("Export fails the area or cost budget.")
  unit <- terra::rast(prepared$unit_path)
  cfg_template <- dp_template(cfg, cfg$connectivity$context_boundary %or% cfg$region$boundary)
  lu <- cfg_template$lu
  region <- terra::rast(file.path(prepared$folder, "region_10m.tif"))
  context <- terra::rast(file.path(prepared$folder, "context_10m.tif"))
  fixed <- terra::rast(file.path(prepared$folder, "fixed_10m.tif")) & context
  roots <- terra::rast(file.path(prepared$folder, "anchors_10m.tif"))
  credit <- terra::rast(file.path(prepared$folder, "credit_10m.tif"))
  chosen <- terra::ifel(!is.na(unit) & unit %in% selected$unit, 1, 0)
  dp_raster(chosen, file.path(output_dir, "selected_10m.tif"), "INT1U")
  # Rebuild connected components directly from the final fine-resolution map.
  reports <- list()
  for (k in names(prepared$graphs)) {
    g <- prepared$graphs[[k]]
    if (k == "Global") {
      mask <- context; credit_mask <- if (cfg$connectivity$global_area_scope == "within_region") region else context
    } else {
      sid <- as.integer(sub("Subarea_", "", k))
      mask <- terra::rast(file.path(prepared$folder, "subareas_10m.tif")) == sid
      mask <- terra::ifel(!is.na(mask) & mask, 1, 0); credit_mask <- mask
    }
    p <- terra::patches(terra::ifel(mask & (fixed | chosen), 1, NA), directions = cfg$grid$directions)
    v <- c(terra::ifel(fixed, credit & credit_mask, chosen) * prod(terra::res(unit))/10000, roots, chosen)
    names(v) <- c("area_ha", "roots", "restored")
    z <- dp_zonal(v, p)
    valid <- z$roots > 0 & (!cfg$connectivity$require_new_restoration | z$restored > 0)
    largest <- if (any(valid)) max(z$area_ha[valid]) else 0
    reports[[k]] <- data.frame(network = k, required_ha = g$minimum, largest_anchored_ha = largest,
      passes = largest + 1e-5 >= g$minimum)
    dp_csv(z, file.path(output_dir, paste0("components_", k, ".csv")))
  }
  report <- if (length(reports)) do.call(rbind, reports) else data.frame(network = character(), required_ha = numeric(), largest_anchored_ha = numeric(), passes = logical())
  dp_csv(report, file.path(output_dir, "exact_10m_network_audit.csv"))
  dp_csv(data.frame(restored_ha = area, cost = cost, current_nature_ha = ledger$current_ha,
    final_nature_ha = ledger$current_ha + area, polygon_ha = ledger$polygon_ha), file.path(output_dir, "solution_summary.csv"))
  if (!all(report$passes)) stop("Fine-resolution network audit failed; no validated pointer published.")
  if (cfg$model$mode == "pathways") {
    final <- as.numeric(lu)
    for (nm in unique(selected$outcome)) final <- terra::ifel(chosen & unit %in% selected$unit[selected$outcome == nm],
      dp_codes(lu, cfg$model$outcomes[[nm]]$class), final)
    final <- terra::categories(final, value = terra::cats(lu)[[1]]); names(final) <- "Landuse"
    dp_raster(final, file.path(output_dir, "restored_landuse_10m.tif"), "INT2U")
  }
  invisible(report)
}
