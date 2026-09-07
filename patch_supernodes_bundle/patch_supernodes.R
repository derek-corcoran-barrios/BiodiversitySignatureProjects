# Utilities for exporting an optimized 10 m restoration graph to AMPL.
#
# These functions are intentionally standalone so they can be sourced from the
# BiodiversityandCarbonOpt workflow now, then moved into TroublemakeR later.
# They do not calculate biodiversity or redefine the user's biodiversity index.

.ps_require_single_layer <- function(x, arg) {
  if (!inherits(x, "SpatRaster")) {
    stop(sprintf("`%s` must be a terra SpatRaster.", arg), call. = FALSE)
  }

  if (terra::nlyr(x) != 1L) {
    stop(sprintf("`%s` must have exactly one layer.", arg), call. = FALSE)
  }

  invisible(TRUE)
}

.ps_same_geometry <- function(x, y, x_name, y_name) {
  ok <- terra::compareGeom(
    x,
    y,
    lyrs = FALSE,
    stopOnError = FALSE
  )

  if (!isTRUE(ok)) {
    stop(
      sprintf("`%s` and `%s` must have identical geometry and CRS.", x_name, y_name),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

.ps_extract_cells <- function(x, cells) {
  if (!length(cells)) {
    return(numeric())
  }

  # The SpatRaster,numeric method interprets `cells` as cell numbers and does
  # not use the SpatVector-only `ID` argument.
  out <- terra::extract(x, cells)

  if (is.atomic(out) && is.null(dim(out))) {
    return(as.vector(out))
  }

  out <- as.data.frame(out)
  as.vector(out[[ncol(out)]])
}

.ps_integer_ids <- function(x, label) {
  x <- suppressWarnings(as.numeric(x))

  if (any(!is.finite(x)) || any(abs(x - round(x)) > 1e-8) || any(x <= 0)) {
    stop(
      sprintf("`%s` must contain finite positive integer IDs.", label),
      call. = FALSE
    )
  }

  as.integer(round(x))
}

.ps_number <- function(x) {
  format(x, scientific = FALSE, trim = TRUE, digits = 15)
}

#' Make a candidate-only raster domain
#'
#' @param current_landuse One-layer categorical or integer SpatRaster.
#' @param eligible_values Values representing restorable land, normally
#'   Agriculture (9) and ProductionForest (11).
#' @param filename Optional checkpoint filename.
#' @param overwrite Passed to terra::writeRaster().
#'
#' @return A one-layer raster named `candidate`, with value 1 on eligible
#'   cells and NA elsewhere.
#' @export
make_candidate_domain <- function(
    current_landuse,
    eligible_values = c(9L, 11L),
    filename = "",
    overwrite = FALSE) {
  .ps_require_single_layer(current_landuse, "current_landuse")

  eligible_values <- unique(as.numeric(eligible_values))
  if (!length(eligible_values) || any(!is.finite(eligible_values))) {
    stop("`eligible_values` must contain finite raster values.", call. = FALSE)
  }

  out <- terra::ifel(current_landuse %in% eligible_values, 1, NA)
  names(out) <- "candidate"

  if (nzchar(filename)) {
    dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
    out <- terra::writeRaster(
      out,
      filename,
      overwrite = overwrite,
      datatype = "INT1U",
      gdal = c("COMPRESS=DEFLATE")
    )
  }

  out
}

#' Rasterize stable, precomputed natural-patch IDs on the planning grid
#'
#' @param patches A SpatVector (or path readable by terra::vect) whose ID field
#'   identifies true connected natural patches.
#' @param template The 10 m raster template used by the optimization.
#' @param patch_id_col Name of the positive-integer patch ID column.
#' @param filename Optional checkpoint filename.
#' @param overwrite Passed to terra::writeRaster().
#'
#' @details One ID must represent one rook-connected natural patch. Do not use
#'   polygon area itself as an ID: two patches may have equal areas. Rasterize
#'   the stable ID, and keep hectares in a separate lookup table.
#'
#' @return An integer patch-ID SpatRaster aligned with `template`.
#' @export
rasterize_precomputed_patch_ids <- function(
    patches,
    template,
    patch_id_col = "patch_id",
    filename = "",
    overwrite = FALSE) {
  .ps_require_single_layer(template, "template")

  if (is.character(patches) && length(patches) == 1L) {
    patches <- terra::vect(patches)
  }

  if (!inherits(patches, "SpatVector")) {
    stop("`patches` must be a terra SpatVector or a vector-file path.", call. = FALSE)
  }

  if (!patch_id_col %in% names(patches)) {
    stop(sprintf("Patch polygons have no `%s` field.", patch_id_col), call. = FALSE)
  }

  patch_ids <- .ps_integer_ids(patches[[patch_id_col]][, 1], patch_id_col)
  patches[[patch_id_col]] <- patch_ids

  if (!nzchar(terra::crs(patches))) {
    stop("Patch polygons have no CRS.", call. = FALSE)
  }
  if (!nzchar(terra::crs(template))) {
    stop("The raster template has no CRS.", call. = FALSE)
  }

  if (!identical(terra::crs(patches), terra::crs(template))) {
    patches <- terra::project(patches, terra::crs(template))
  }

  out <- terra::rasterize(
    patches,
    template,
    field = patch_id_col,
    background = NA,
    touches = FALSE
  )
  names(out) <- "patch_id"

  if (nzchar(filename)) {
    dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
    out <- terra::writeRaster(
      out,
      filename,
      overwrite = overwrite,
      datatype = "INT4S",
      gdal = c("COMPRESS=DEFLATE")
    )
  }

  out
}

#' Build the compressed candidate-cell/current-patch graph
#'
#' @param candidate_domain One-layer raster with non-NA values only on cells
#'   eligible for restoration.
#' @param patch_id One-layer positive-integer ID raster for connected current
#'   natural patches, aligned with `candidate_domain`.
#' @param patch_areas Data frame with one row per connected patch and its
#'   precomputed area in hectares.
#' @param patch_id_col ID column in `patch_areas`.
#' @param patch_area_col Hectare column in `patch_areas`.
#' @param current_landuse Optional aligned current-land-use raster. If supplied,
#'   same-type current-nature contact counts are also prepared.
#' @param class_lookup Optional data frame mapping current raster values to AMPL
#'   nature-type names. Required with `current_landuse`.
#' @param class_value_col Value column in `class_lookup`.
#' @param class_name_col AMPL nature-type name column in `class_lookup`.
#' @param directions Must be 4. The hard rule uses shared-edge connectivity.
#'
#' @return A `patch_supernodes` list ready for write_patch_supernodes_dat().
#' @export
prepare_patch_supernodes <- function(
    candidate_domain,
    patch_id,
    patch_areas,
    patch_id_col = "patch_id",
    patch_area_col = "area_ha",
    current_landuse = NULL,
    class_lookup = NULL,
    class_value_col = "value",
    class_name_col = "landuse",
    directions = 4L) {
  .ps_require_single_layer(candidate_domain, "candidate_domain")
  .ps_require_single_layer(patch_id, "patch_id")
  .ps_same_geometry(candidate_domain, patch_id, "candidate_domain", "patch_id")

  if (!identical(as.integer(directions), 4L)) {
    stop("Only rook connectivity (`directions = 4`) is supported.", call. = FALSE)
  }

  if (!is.data.frame(patch_areas)) {
    patch_areas <- as.data.frame(patch_areas)
  }

  needed <- c(patch_id_col, patch_area_col)
  missing_cols <- setdiff(needed, names(patch_areas))
  if (length(missing_cols)) {
    stop(
      sprintf("`patch_areas` is missing: %s.", paste(missing_cols, collapse = ", ")),
      call. = FALSE
    )
  }

  patch_areas <- patch_areas[, needed, drop = FALSE]
  names(patch_areas) <- c("patch_id", "area_ha")
  patch_areas$patch_id <- .ps_integer_ids(patch_areas$patch_id, patch_id_col)
  patch_areas$area_ha <- suppressWarnings(as.numeric(patch_areas$area_ha))

  if (anyDuplicated(patch_areas$patch_id)) {
    stop("`patch_areas` must have exactly one row per patch ID.", call. = FALSE)
  }
  if (any(!is.finite(patch_areas$area_ha)) || any(patch_areas$area_ha <= 0)) {
    stop("Every precomputed patch area must be finite and positive.", call. = FALSE)
  }

  candidate_values <- terra::values(candidate_domain, mat = FALSE)
  candidate_cells <- which(!is.na(candidate_values))

  if (!length(candidate_cells)) {
    stop("The candidate domain contains no non-missing cells.", call. = FALSE)
  }

  overlap_ids <- .ps_extract_cells(patch_id, candidate_cells)
  if (any(!is.na(overlap_ids))) {
    stop(
      "Candidate cells overlap current natural patches; the two domains must be disjoint.",
      call. = FALSE
    )
  }

  adjacent_all <- terra::adjacent(
    candidate_domain,
    cells = candidate_cells,
    directions = 4,
    pairs = TRUE,
    include = FALSE
  )

  adjacent_all <- as.matrix(adjacent_all)
  if (!nrow(adjacent_all)) {
    adjacent_all <- matrix(integer(), ncol = 2L)
  }
  if (ncol(adjacent_all) != 2L) {
    stop("terra::adjacent() did not return a two-column pair matrix.", call. = FALSE)
  }
  storage.mode(adjacent_all) <- "integer"
  colnames(adjacent_all) <- c("cell", "neighbor")
  adjacent_all <- unique(adjacent_all)

  neighbor_is_candidate <-
    match(adjacent_all[, "neighbor"], candidate_cells, nomatch = 0L) > 0L

  candidate_pairs <- adjacent_all[neighbor_is_candidate, , drop = FALSE]
  if (nrow(candidate_pairs)) {
    edges <- data.frame(
      i = pmin(candidate_pairs[, 1], candidate_pairs[, 2]),
      j = pmax(candidate_pairs[, 1], candidate_pairs[, 2])
    )
    edges <- unique(edges[edges$i != edges$j, , drop = FALSE])
    edges <- edges[order(edges$i, edges$j), , drop = FALSE]
    rownames(edges) <- NULL
  } else {
    edges <- data.frame(i = integer(), j = integer())
  }

  neighbor_patch <- .ps_extract_cells(patch_id, adjacent_all[, "neighbor"])
  has_patch <- !is.na(neighbor_patch)

  if (any(has_patch)) {
    patch_edges <- data.frame(
      cell = adjacent_all[has_patch, "cell"],
      patch_id = .ps_integer_ids(neighbor_patch[has_patch], "patch_id raster")
    )
    patch_edges <- unique(patch_edges)
    patch_edges <- patch_edges[order(patch_edges$cell, patch_edges$patch_id), , drop = FALSE]
    rownames(patch_edges) <- NULL
  } else {
    patch_edges <- data.frame(cell = integer(), patch_id = integer())
  }

  touching_ids <- sort(unique(patch_edges$patch_id))
  missing_area_ids <- setdiff(touching_ids, patch_areas$patch_id)
  if (length(missing_area_ids)) {
    stop(
      sprintf(
        "No precomputed area for touching patch ID(s): %s.",
        paste(utils::head(missing_area_ids, 20L), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  touching_patches <- patch_areas[
    match(touching_ids, patch_areas$patch_id),
    ,
    drop = FALSE
  ]
  rownames(touching_patches) <- NULL

  contacts <- data.frame(
    landuse = character(),
    cell = integer(),
    contacts = integer()
  )

  if (!is.null(current_landuse)) {
    .ps_require_single_layer(current_landuse, "current_landuse")
    .ps_same_geometry(
      candidate_domain,
      current_landuse,
      "candidate_domain",
      "current_landuse"
    )

    if (is.null(class_lookup)) {
      stop("Supply `class_lookup` with `current_landuse`.", call. = FALSE)
    }
    class_lookup <- as.data.frame(class_lookup)
    lookup_needed <- c(class_value_col, class_name_col)
    if (!all(lookup_needed %in% names(class_lookup))) {
      stop(
        sprintf("`class_lookup` must contain %s.", paste(lookup_needed, collapse = " and ")),
        call. = FALSE
      )
    }

    lookup <- class_lookup[, lookup_needed, drop = FALSE]
    names(lookup) <- c("value", "landuse")
    if (anyDuplicated(lookup$value)) {
      stop("`class_lookup` values must be unique.", call. = FALSE)
    }

    neighbor_class <- .ps_extract_cells(
      current_landuse,
      adjacent_all[, "neighbor"]
    )
    lookup_pos <- match(neighbor_class, lookup$value)
    keep <- !is.na(lookup_pos)

    if (any(keep)) {
      raw_contacts <- data.frame(
        landuse = as.character(lookup$landuse[lookup_pos[keep]]),
        cell = adjacent_all[keep, "cell"],
        contacts = 1L
      )
      contacts <- stats::aggregate(
        contacts ~ landuse + cell,
        data = raw_contacts,
        FUN = sum
      )
      contacts$contacts <- as.integer(contacts$contacts)
      contacts <- contacts[order(contacts$landuse, contacts$cell), , drop = FALSE]
      rownames(contacts) <- NULL
    }
  }

  out <- list(
    cells = as.integer(candidate_cells),
    edges = edges,
    patches = touching_patches,
    cell_patch_edges = patch_edges,
    existing_nature_contacts = contacts,
    resolution = terra::res(candidate_domain),
    crs = terra::crs(candidate_domain)
  )
  class(out) <- c("patch_supernodes", class(out))
  out
}

#' Append compressed graph sections to an AMPL data file
#'
#' @param x Output of prepare_patch_supernodes().
#' @param dat_path Exact `.dat` path to create or append.
#' @param append Append to an existing file made by TroublemakeR.
#' @param write_cells Also write `set Cells`. Set FALSE when
#'   TroublemakeR::define_cells() has already written it.
#'
#' @return The normalized `.dat` path, invisibly.
#' @export
write_patch_supernodes_dat <- function(
    x,
    dat_path,
    append = TRUE,
    write_cells = FALSE) {
  if (!inherits(x, "patch_supernodes")) {
    stop("`x` must come from prepare_patch_supernodes().", call. = FALSE)
  }

  dir.create(dirname(dat_path), recursive = TRUE, showWarnings = FALSE)
  con <- file(dat_path, open = if (append) "a" else "w", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)

  put <- function(...) writeLines(c(...), con = con, sep = "\n", useBytes = TRUE)

  put("", "# ---- Compressed current-nature patch network --------------------")

  if (isTRUE(write_cells)) {
    put("set Cells :=")
    if (length(x$cells)) put(paste0("  ", x$cells))
    put(";")
  }

  put("set E :=")
  if (nrow(x$edges)) {
    put(sprintf("  (%d,%d)", x$edges$i, x$edges$j))
  }
  put(";")

  put("", "set ExistingPatches :=")
  if (nrow(x$patches)) {
    put(paste0("  ", x$patches$patch_id))
  }
  put(";")

  put("", "set CellPatchEdges :=")
  if (nrow(x$cell_patch_edges)) {
    put(sprintf(
      "  (%d,%d)",
      x$cell_patch_edges$cell,
      x$cell_patch_edges$patch_id
    ))
  }
  put(";")

  put("", "param ExistingPatchAreaHa :=")
  if (nrow(x$patches)) {
    put(sprintf(
      "  %d %s",
      x$patches$patch_id,
      .ps_number(x$patches$area_ha)
    ))
  }
  put(";")

  if (nrow(x$existing_nature_contacts)) {
    put("", "param ExistingNatureContacts default 0 :=")
    put(sprintf(
      "  %s %d %d",
      x$existing_nature_contacts$landuse,
      x$existing_nature_contacts$cell,
      x$existing_nature_contacts$contacts
    ))
    put(";")
  }

  invisible(normalizePath(dat_path, mustWork = TRUE))
}

#' Summarize a compressed network before writing or solving it
#'
#' @param x Output of prepare_patch_supernodes().
#' @param min_patch_area_ha Intended hard threshold.
#'
#' @return A one-row data frame suitable for an audit CSV.
#' @export
audit_patch_supernodes <- function(x, min_patch_area_ha = 5000) {
  if (!inherits(x, "patch_supernodes")) {
    stop("`x` must come from prepare_patch_supernodes().", call. = FALSE)
  }

  linked_cells <- unique(c(
    x$edges$i,
    x$edges$j,
    x$cell_patch_edges$cell
  ))

  data.frame(
    candidate_cells = length(x$cells),
    candidate_edges = nrow(x$edges),
    touching_existing_patches = nrow(x$patches),
    cell_patch_edges = nrow(x$cell_patch_edges),
    isolated_candidate_cells = sum(!x$cells %in% linked_cells),
    touching_patch_area_ha = sum(x$patches$area_ha),
    largest_touching_patch_ha = if (nrow(x$patches)) {
      max(x$patches$area_ha)
    } else {
      0
    },
    min_patch_area_ha = min_patch_area_ha,
    cell_area_ha_from_resolution = prod(x$resolution) / 10000,
    stringsAsFactors = FALSE
  )
}
