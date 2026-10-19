# Source chunk: packages
library(terra)
library(dplyr)
library(tidyr)
library(purrr)
library(readr)
library(stringr)
library(digest)
library(igraph)
library(knitr)


# Source chunk: names-and-validation
nature_landuses <- c(
  "ForestDryPoor",
  "ForestDryRich",
  "ForestWetPoor",
  "ForestWetRich",
  "OpenDryPoor",
  "OpenDryRich",
  "OpenWetPoor",
  "OpenWetRich"
)
forest_landuses <- c(
  "ForestDryPoor",
  "ForestDryRich",
  "ForestWetPoor",
  "ForestWetRich"
)

required_networks <- c(
  "Global",
  paste0("Subarea", seq_len(expected_number_of_subareas))
)
network_minimums <- c(
  Global = global_network_minimum_ha,
  stats::setNames(
    rep(subarea_network_minimum_ha, expected_number_of_subareas),
    paste0("Subarea", seq_len(expected_number_of_subareas))
  )
)

stopifnot(
  length(nature_landuses) == 8L,
  length(forest_landuses) == 4L,
  all(forest_landuses %in% nature_landuses),
  global_network_minimum_ha > 0,
  subarea_network_minimum_ha > 0,
  global_nature_target_fraction > 0,
  global_nature_target_fraction <= 1,
  global_restoration_tolerance_ha >= 0,
  spatial_contiguity_bonus >= 0,
  connected_existing_patch_bonus >= 0,
  new_to_new_adjacency_weight >= 0,
  existing_to_new_adjacency_weight >= 0,
  abs(
    new_to_new_adjacency_weight +
      existing_to_new_adjacency_weight - 1
  ) < 1e-10,
  production_forest_share_threshold > 0,
  production_forest_share_threshold <= 1,
  network_shortfall_tolerance_ha >= 0,
  repair_change_tolerance_ha >= 0,
  repair_feasibility_time_limit_seconds > 0,
  repair_minimum_change_time_limit_seconds > 0,
  repair_quality_time_limit_seconds > 0,
  repair_mip_gap >= 0,
  repair_mip_gap < 1
)


# Source chunk: helpers
shapefile_components <- function(path) {
  if (tolower(tools::file_ext(path)) != "shp") return(path)
  stem <- tools::file_path_sans_ext(path)
  paste0(stem, c(".shp", ".shx", ".dbf", ".prj"))
}

stop_for_missing_files <- function(paths, description) {
  paths <- unique(paths)
  missing <- paths[!file.exists(paths)]
  if (length(missing)) {
    stop(
      paste0(
        "Missing ", description, ":\n- ",
        paste(missing, collapse = "\n- ")
      ),
      call. = FALSE
    )
  }
  invisible(paths)
}

same_crs <- function(x, y) isTRUE(terra::same.crs(x, y))

align_to_template <- function(x, template, method) {
  if (!same_crs(x, template)) {
    x <- terra::project(x, template, method = method)
  } else if (!isTRUE(terra::compareGeom(
    x,
    template,
    lyrs = FALSE,
    stopOnError = FALSE
  ))) {
    x <- terra::resample(x, template, method = method)
  }
  terra::mask(x, template)
}

category_codes <- function(x, requested_labels) {
  lev <- terra::levels(x)[[1]]
  if (is.null(lev) || !nrow(lev)) {
    stop("The fine land-use raster has no category table.", call. = FALSE)
  }

  numeric_columns <- names(lev)[vapply(lev, is.numeric, logical(1))]
  label_columns <- names(lev)[vapply(
    lev,
    function(z) is.character(z) || is.factor(z),
    logical(1)
  )]
  if (!length(numeric_columns) || !length(label_columns)) {
    stop("Could not identify the land-use category columns.", call. = FALSE)
  }

  labels <- as.character(lev[[label_columns[[1]]]])
  missing_labels <- setdiff(requested_labels, labels)
  if (length(missing_labels)) {
    stop(
      "Missing land-use classes: ",
      paste(missing_labels, collapse = ", "),
      call. = FALSE
    )
  }

  stats::setNames(
    lev[[numeric_columns[[1]]]][
      match(requested_labels, labels)
    ],
    requested_labels
  )
}

indicator_for_codes <- function(x, codes) {
  terra::ifel(
    is.na(x),
    NA,
    terra::ifel(x %in% codes, 1, 0)
  )
}

zonal_sum_table <- function(x, zones, id_name) {
  out <- terra::zonal(
    x,
    zones,
    fun = "sum",
    na.rm = TRUE,
    wide = TRUE
  ) |>
    as.data.frame()
  names(out)[1] <- id_name
  out[[id_name]] <- as.numeric(out[[id_name]])
  out
}

global_number <- function(x, fun) {
  out <- if (fun %in% c("notNA", "isNA")) {
    terra::global(x, fun)
  } else {
    terra::global(x, fun, na.rm = TRUE)
  }
  as.numeric(out[1, 1])
}

unwrap_checkpoint <- function(path) {
  x <- readRDS(path)
  if (is.list(x) && "value" %in% names(x)) x$value else x
}

file_fingerprint <- function(paths) {
  paths <- normalizePath(paths, winslash = "/", mustWork = FALSE)
  info <- file.info(paths)
  tibble(path = paths, size = info$size, modified = as.numeric(info$mtime))
}

load_or_build_checkpoint <- function(
    path,
    signature,
    builder,
    required_sidecars = character()
) {
  if (
    reuse_matching_checkpoints &&
    file.exists(path) &&
    all(file.exists(required_sidecars))
  ) {
    checkpoint <- readRDS(path)
    if (
      is.list(checkpoint) &&
      identical(checkpoint$signature, signature) &&
      !is.null(checkpoint$value)
    ) {
      return(checkpoint$value)
    }
  }
  value <- builder()
  saveRDS(list(signature = signature, value = value), path)
  value
}

ampl_number <- function(x) {
  format(x, scientific = FALSE, trim = TRUE, digits = 15)
}

append_dat <- function(path, lines) {
  cat(
    paste(lines, collapse = "\n"),
    "\n",
    file = path,
    append = TRUE,
    sep = ""
  )
}

chunk_text <- function(x, chunk_size = 20L) {
  if (!length(x)) return(character())
  groups <- ceiling(seq_along(x) / chunk_size)
  vapply(
    split(x, groups),
    paste,
    collapse = " ",
    FUN.VALUE = character(1)
  )
}

write_ampl_set <- function(path, name, values) {
  append_dat(
    path,
    c(paste("set", name, ":="), chunk_text(values), ";", "")
  )
}

write_ampl_pair_set <- function(path, name, first, second) {
  tuples <- paste0("(", first, ",", second, ")")
  append_dat(
    path,
    c(
      paste("set", name, ":="),
      chunk_text(tuples, chunk_size = 8L),
      ";",
      ""
    )
  )
}

write_ampl_triple_set <- function(path, name, first, second, third) {
  tuples <- paste0("(", first, ",", second, ",", third, ")")
  append_dat(
    path,
    c(
      paste("set", name, ":="),
      chunk_text(tuples, chunk_size = 5L),
      ";",
      ""
    )
  )
}

write_ampl_scalar <- function(path, name, value) {
  append_dat(
    path,
    c(paste0("param ", name, " := ", ampl_number(value), ";"), "")
  )
}

write_ampl_vector <- function(path, name, key, value) {
  rows <- paste(key, ampl_number(value))
  append_dat(path, c(paste("param", name, ":="), rows, ";", ""))
}

write_ampl_matrix <- function(path, name, row_key, column_key, value) {
  rows <- paste(row_key, column_key, ampl_number(value))
  append_dat(path, c(paste("param", name, ":="), rows, ";", ""))
}

commit_generated_file <- function(temporary_path, final_path) {
  unchanged <- file.exists(final_path) &&
    unname(tools::md5sum(temporary_path)) ==
      unname(tools::md5sum(final_path))
  if (unchanged) {
    unlink(temporary_path)
    return(invisible(FALSE))
  }
  copied <- file.copy(temporary_path, final_path, overwrite = TRUE)
  unlink(temporary_path)
  if (!copied) stop("Could not write: ", final_path, call. = FALSE)
  invisible(TRUE)
}

write_lines_if_changed <- function(lines, path) {
  temporary_path <- paste0(path, ".building")
  writeLines(lines, temporary_path)
  commit_generated_file(temporary_path, path)
}

ampl_path <- function(path) chartr(intToUtf8(92L), "/", path)

append_start <- function(path, lines) {
  if (!length(lines)) return(invisible(NULL))
  cat(
    paste(lines, collapse = "\n"),
    "\n",
    file = path,
    append = TRUE,
    sep = ""
  )
}


# Source chunk: read-upstream
dir.create(problem_dir, recursive = TRUE, showWarnings = FALSE)

warmup_selected_path <- file.path(
  warmup_dir,
  paste0(warmup_name, "_selected_cells.tsv")
)
warmup_summary_path <- file.path(
  warmup_dir,
  paste0(warmup_name, "_solution_summary.tsv")
)
candidate_table_path <- file.path(
  warmup_dir,
  "candidate_cell_dictionary.csv"
)
candidate_subarea_path <- file.path(
  warmup_dir,
  "candidate_area_by_subarea.csv"
)
candidate_edges_path <- file.path(
  warmup_dir,
  "candidate_rook_edges.csv"
)
subarea_budget_path <- file.path(
  warmup_dir,
  "subarea_30pct_budget.csv"
)
existing_contact_path <- file.path(
  warmup_dir,
  "existing_nature_contact_by_action.csv"
)
biodiversity_checkpoint_candidates <- c(
  file.path(warmup_dir, "coarse_biodiversity_benefits.rds"),
  file.path(
    baseline_checkpoint_dir,
    "coarse_biodiversity_benefits.rds"
  )
)
biodiversity_checkpoint_paths_found <-
  biodiversity_checkpoint_candidates[
    file.exists(biodiversity_checkpoint_candidates)
  ]
if (!length(biodiversity_checkpoint_paths_found)) {
  stop(
    paste0(
      "No biodiversity checkpoint was found. Expected either:\n- ",
      paste(biodiversity_checkpoint_candidates, collapse = "\n- ")
    ),
    call. = FALSE
  )
}
biodiversity_checkpoint_path <-
  biodiversity_checkpoint_paths_found[[1]]

biodiversity_support_path <- file.path(
  metrics_dir,
  paste0(
    biodiversity_support_metric,
    "_",
    biodiversity_support_scenario,
    ".tif"
  )
)

required_paths <- c(
  fine_landuse_path,
  coarse_landuse_path,
  shapefile_components(subareas_path),
  national_patch_path,
  national_patch_dictionary_path,
  biodiversity_support_path,
  model_path,
  warmup_selected_path,
  warmup_summary_path,
  candidate_table_path,
  candidate_subarea_path,
  candidate_edges_path,
  subarea_budget_path,
  existing_contact_path,
  biodiversity_checkpoint_path
)
stop_for_missing_files(required_paths, "08a repair input file(s)")

warmup_summary <- readr::read_tsv(
  warmup_summary_path,
  show_col_types = FALSE,
  col_types = cols(.default = col_character())
)
if (!gp_accepted_result(as.numeric(warmup_summary$value[warmup_summary$item == "solve_result_num"]))) {
  stop(
    "The upstream summary does not report a solved/optimal warm-up.",
    call. = FALSE
  )
}

# Make the repaired output folder a self-contained upstream folder for the
# subsequent 08b refinement. These are immutable support inputs; solution
# tables are written later by AMPL under problem_name.
downstream_support_sources <- c(
  candidate_table_path,
  candidate_subarea_path,
  candidate_edges_path,
  subarea_budget_path,
  existing_contact_path,
  biodiversity_checkpoint_path
)
downstream_support_names <- c(
  "candidate_cell_dictionary.csv",
  "candidate_area_by_subarea.csv",
  "candidate_rook_edges.csv",
  "subarea_30pct_budget.csv",
  "existing_nature_contact_by_action.csv",
  "coarse_biodiversity_benefits.rds"
)
downstream_support_destinations <- file.path(
  problem_dir,
  downstream_support_names
)
support_copy_ok <- map2_lgl(
  downstream_support_sources,
  downstream_support_destinations,
  function(source, destination) {
    if (
      file.exists(destination) &&
        identical(
          unname(tools::md5sum(source)),
          unname(tools::md5sum(destination))
        )
    ) {
      return(TRUE)
    }
    isTRUE(file.copy(source, destination, overwrite = TRUE))
  }
)
if (!all(support_copy_ok)) {
  stop(
    "Could not copy all downstream warm-start support files into problem_dir.",
    call. = FALSE
  )
}

warmup_selected <- readr::read_tsv(
  warmup_selected_path,
  show_col_types = FALSE
) |>
  transmute(
    coarse_cell_id = as.integer(coarse_cell_id),
    landuse = as.character(landuse)
  )

candidate_table <- readr::read_csv(
  candidate_table_path,
  show_col_types = FALSE
) |>
  mutate(
    coarse_cell_id = as.integer(coarse_cell_id),
    ProductionForestShare = if_else(
      EligibleAreaHa > 0,
      ProductionForestHa / EligibleAreaHa,
      0
    ),
    IsProductionForestCell = as.integer(
      ProductionForestShare >= production_forest_share_threshold
    )
  ) |>
  filter(IsCandidate) |>
  arrange(coarse_cell_id)

candidate_cells <- candidate_table$coarse_cell_id

restorable_area_by_subarea_cell <- readr::read_csv(
  candidate_subarea_path,
  show_col_types = FALSE
) |>
  transmute(
    subarea_id = as.integer(subarea_id),
    coarse_cell_id = as.integer(coarse_cell_id),
    RestorableAreaHaBySubarea = as.numeric(
      RestorableAreaHaBySubarea
    )
  ) |>
  filter(coarse_cell_id %in% candidate_cells)

candidate_edges <- readr::read_csv(
  candidate_edges_path,
  show_col_types = FALSE
) |>
  transmute(i = as.integer(i), j = as.integer(j)) |>
  filter(i %in% candidate_cells, j %in% candidate_cells) |>
  distinct() |>
  arrange(i, j)

subarea_budget <- readr::read_csv(
  subarea_budget_path,
  show_col_types = FALSE
) |>
  mutate(subarea_id = as.integer(subarea_id)) |>
  arrange(subarea_id)

global_budget <- subarea_budget |>
  summarise(
    GlobalPolygonAreaHa = sum(PolygonAreaHa),
    GlobalCurrentNatureAreaHa = sum(CurrentNatureAreaHa)
  ) |>
  mutate(
    GlobalNatureTargetFraction = global_nature_target_fraction,
    GlobalFinalNatureTargetHa =
      GlobalNatureTargetFraction * GlobalPolygonAreaHa,
    GlobalRestorationRequirementHa = pmax(
      0,
      GlobalFinalNatureTargetHa - GlobalCurrentNatureAreaHa
    ),
    GlobalRestorationToleranceHa = global_restoration_tolerance_ha,
    GlobalMaximumRestorationHa =
      GlobalRestorationRequirementHa + GlobalRestorationToleranceHa,
    GlobalMaximumFinalNatureHa =
      GlobalFinalNatureTargetHa + GlobalRestorationToleranceHa,
    GlobalMaximumFinalNaturePercent =
      100 * GlobalMaximumFinalNatureHa / GlobalPolygonAreaHa
  )

readr::write_csv(
  global_budget,
  file.path(problem_dir, "whole_polygon_30pct_budget.csv")
)

existing_nature_contact <- readr::read_csv(
  existing_contact_path,
  show_col_types = FALSE
) |>
  transmute(
    landuse = as.character(landuse),
    coarse_cell_id = as.integer(coarse_cell_id),
    ExistingNatureContact = as.numeric(ExistingNatureContact)
  ) |>
  filter(coarse_cell_id %in% candidate_cells)

biodiversity_checkpoint <- unwrap_checkpoint(
  biodiversity_checkpoint_path
)
if (
  !is.list(biodiversity_checkpoint) ||
  !"benefits" %in% names(biodiversity_checkpoint)
) {
  stop("Unexpected biodiversity checkpoint structure.", call. = FALSE)
}
biodiversity_table <- biodiversity_checkpoint$benefits |>
  transmute(
    coarse_cell_id = as.integer(coarse_cell_id),
    landuse = as.character(landuse),
    BioDiversity = as.numeric(BioDiversity)
  ) |>
  filter(coarse_cell_id %in% candidate_cells)

stopifnot(
  nrow(subarea_budget) == expected_number_of_subareas,
  nrow(warmup_selected) > 0L,
  !anyDuplicated(warmup_selected$coarse_cell_id),
  all(warmup_selected$coarse_cell_id %in% candidate_cells),
  all(warmup_selected$landuse %in% nature_landuses),
  length(candidate_cells) > 0L,
  nrow(candidate_edges) > 0L,
  nrow(biodiversity_table) ==
    length(candidate_cells) * length(nature_landuses),
  nrow(existing_nature_contact) ==
    length(candidate_cells) * length(nature_landuses),
  !anyNA(biodiversity_table),
  !anyNA(existing_nature_contact)
)

input_audit <- file_fingerprint(required_paths)
readr::write_csv(
  input_audit,
  file.path(problem_dir, "input_file_audit.csv")
)
readr::write_csv(
  tibble(
    candidate_path = biodiversity_checkpoint_candidates,
    exists = file.exists(biodiversity_checkpoint_candidates),
    selected = biodiversity_checkpoint_candidates ==
      biodiversity_checkpoint_path
  ),
  file.path(problem_dir, "biodiversity_checkpoint_resolution.csv")
)
knitr::kable(input_audit)


# Source chunk: aligned-domain
fine_landuse_full <- terra::rast(fine_landuse_path)
coarse_landuse_full <- terra::rast(coarse_landuse_path)

stopifnot(
  terra::nlyr(fine_landuse_full) == 1L,
  terra::nlyr(coarse_landuse_full) == 1L,
  all(abs(terra::res(fine_landuse_full) - 10) < 1e-8),
  all(abs(terra::res(coarse_landuse_full) - 200) < 1e-8),
  same_crs(fine_landuse_full, coarse_landuse_full)
)

fine_cell_area_ha <- prod(terra::res(fine_landuse_full)) / 10000

subareas <- terra::vect(subareas_path) |>
  terra::makeValid()
if (nrow(subareas) != expected_number_of_subareas) {
  stop(
    "Expected four subareas; found ", nrow(subareas), ".",
    call. = FALSE
  )
}
if (!same_crs(subareas, fine_landuse_full)) {
  subareas <- terra::project(subareas, terra::crs(fine_landuse_full))
}

if (is.null(subarea_label_field)) {
  subarea_labels <- paste0("subarea_", seq_len(nrow(subareas)))
} else {
  if (!subarea_label_field %in% names(subareas)) {
    stop("Invalid subarea_label_field.", call. = FALSE)
  }
  subarea_labels <- make.unique(as.character(
    as.data.frame(subareas)[[subarea_label_field]]
  ))
}

subareas$subarea_id <- seq_len(nrow(subareas))
subareas$subarea_label <- subarea_labels
polygon_areas_ha <- as.numeric(
  terra::expanse(subareas, unit = "ha", transform = FALSE)
)
planning_domain <- terra::aggregate(subareas)

if (
  any(
    abs(sort(polygon_areas_ha) - sort(expected_subarea_areas_ha)) >
      polygon_area_audit_tolerance_ha
  ) ||
  abs(sum(polygon_areas_ha) - expected_total_area_ha) >
    polygon_area_audit_tolerance_ha
) {
  stop(
    "The polygon areas do not match the intended smaller Gudenå polygon.",
    call. = FALSE
  )
}

fine_landuse <- terra::crop(
  fine_landuse_full,
  planning_domain,
  snap = "out"
) |>
  terra::mask(planning_domain, touches = FALSE)

coarse_landuse <- terra::crop(
  coarse_landuse_full,
  planning_domain,
  snap = "out"
) |>
  terra::mask(planning_domain, touches = TRUE)

coarse_cell_id <- terra::init(coarse_landuse, "cell")
coarse_cell_id <- terra::ifel(
  is.na(coarse_landuse),
  NA,
  coarse_cell_id
)
names(coarse_cell_id) <- "coarse_cell_id"

coarse_id_on_fine <- align_to_template(
  coarse_cell_id,
  fine_landuse,
  method = "near"
)
names(coarse_id_on_fine) <- "coarse_cell_id"

subarea_id_on_fine <- terra::rasterize(
  subareas[, "subarea_id"],
  fine_landuse,
  field = "subarea_id",
  background = NA,
  touches = FALSE
)
subarea_id_on_fine <- terra::mask(
  subarea_id_on_fine,
  fine_landuse
)
names(subarea_id_on_fine) <- "subarea_id"

all_required_classes <- c(
  nature_landuses,
  "Agriculture",
  "ProductionForest",
  "Urban",
  "Other"
)
landuse_codes <- category_codes(fine_landuse, all_required_classes)
nature_codes <- unname(landuse_codes[nature_landuses])
eligible_codes <- unname(
  landuse_codes[c("Agriculture", "ProductionForest")]
)

fine_nature <- indicator_for_codes(fine_landuse, nature_codes)
names(fine_nature) <- "CurrentNatureCells"

biodiversity_support <- terra::rast(biodiversity_support_path)
biodiversity_support <- align_to_template(
  biodiversity_support,
  fine_landuse,
  method = "bilinear"
)
fine_biodiversity_support <- terra::ifel(
  is.na(fine_landuse),
  NA,
  terra::ifel(!is.na(biodiversity_support), 1, 0)
)

fine_eligible <-
  indicator_for_codes(fine_landuse, eligible_codes) *
  fine_biodiversity_support
fine_eligible_mask <- terra::ifel(fine_eligible == 1, 1, NA)
names(fine_eligible_mask) <- "EligibleCells"

candidate_mask_on_fine <- terra::ifel(
  coarse_id_on_fine %in% candidate_cells,
  1,
  NA
)
fine_candidate_eligible_mask <- terra::mask(
  fine_eligible_mask,
  candidate_mask_on_fine
)

if (
  global_number(
    terra::ifel(
      !is.na(fine_landuse) & is.na(coarse_id_on_fine),
      1,
      0
    ),
    "sum"
  ) > 0
) {
  stop("Some valid 10 m cells have no 200 m ID.", call. = FALSE)
}

coarse_ids_reconstructed <- sort(unique(na.omit(as.integer(
  terra::values(coarse_cell_id, mat = FALSE)
))))
if (!all(candidate_cells %in% coarse_ids_reconstructed)) {
  stop(
    "Reconstructed coarse IDs do not match notebook 07 outputs.",
    call. = FALSE
  )
}


# Source chunk: patch-network
national_patches <- terra::vect(
  national_patch_path,
  layer = national_patch_layer
) |>
  terra::makeValid()

if (!same_crs(national_patches, coarse_landuse)) {
  national_patches <- terra::project(
    national_patches,
    terra::crs(coarse_landuse)
  )
}
if (!"patch_id" %in% names(national_patches)) {
  stop("The national patch layer has no patch_id field.", call. = FALSE)
}
national_patches$patch_id <- as.integer(national_patches$patch_id)

patch_dictionary_full <- readr::read_csv(
  national_patch_dictionary_path,
  show_col_types = FALSE
) |>
  transmute(
    patch_id = as.integer(patch_id),
    ExistingPatchAreaHa = as.numeric(area_network_ha)
  )

stopifnot(
  !anyDuplicated(patch_dictionary_full$patch_id),
  all(is.finite(patch_dictionary_full$ExistingPatchAreaHa)),
  all(patch_dictionary_full$ExistingPatchAreaHa > 0)
)

patches_near_domain <- terra::crop(
  national_patches,
  terra::buffer(planning_domain, width = contact_buffer_m)
)
if (!nrow(patches_near_domain)) {
  stop("No national current-nature patches occur near the domain.")
}

# Candidate planning units become individual 200 m polygons. This preserves
# the established coarse contact convention from notebook 04.
candidate_id_raster <- coarse_cell_id
candidate_id_values <- rep(NA_real_, terra::ncell(candidate_id_raster))
candidate_id_values[candidate_cells] <- candidate_cells
terra::values(candidate_id_raster) <- candidate_id_values
names(candidate_id_raster) <- "coarse_cell_id"

candidate_polygons <- terra::as.polygons(
  candidate_id_raster,
  aggregate = FALSE,
  values = TRUE,
  na.rm = TRUE
)
candidate_polygons$coarse_cell_id <- as.integer(
  candidate_polygons$coarse_cell_id
)

# Exact patch area inside each delområde. These coefficients are distinct from
# the full national area used by the Global certificate.
patch_fragments_by_subarea <- terra::intersect(
  patches_near_domain[, "patch_id"],
  subareas[, "subarea_id"]
)
if (!nrow(patch_fragments_by_subarea)) {
  stop("No national patch fragments intersect the four subareas.")
}
patch_fragments_by_subarea$fragment_area_ha <- terra::expanse(
  patch_fragments_by_subarea,
  unit = "ha",
  transform = FALSE
)

patch_area_by_subarea_all <- patch_fragments_by_subarea |>
  as.data.frame() |>
  transmute(
    patch_id = as.integer(patch_id),
    subarea_id = as.integer(subarea_id),
    fragment_area_ha = as.numeric(fragment_area_ha)
  ) |>
  group_by(patch_id, subarea_id) |>
  summarise(
    ExistingPatchAreaHaBySubarea = sum(fragment_area_ha),
    .groups = "drop"
  ) |>
  filter(ExistingPatchAreaHaBySubarea > 0)

national_current_by_subarea <- patch_area_by_subarea_all |>
  group_by(subarea_id) |>
  summarise(
    NationalPatchNatureHa = sum(ExistingPatchAreaHaBySubarea),
    .groups = "drop"
  ) |>
  right_join(
    subarea_budget |>
      select(subarea_id, LanduseRasterNatureHa = CurrentNatureAreaHa),
    by = "subarea_id"
  ) |>
  mutate(
    NationalPatchNatureHa = coalesce(NationalPatchNatureHa, 0),
    DifferenceHa = NationalPatchNatureHa - LanduseRasterNatureHa,
    DifferencePercent =
      100 * DifferenceHa / pmax(LanduseRasterNatureHa, 1e-9)
  ) |>
  arrange(subarea_id)

patch_signature <- digest::digest(list(
  files = file_fingerprint(c(
    coarse_landuse_path,
    shapefile_components(subareas_path),
    national_patch_path,
    national_patch_dictionary_path,
    candidate_table_path,
    candidate_edges_path
  )),
  national_patch_layer = national_patch_layer,
  candidate_cells = candidate_cells,
  contact_buffer_m = contact_buffer_m,
  minimum_contact_proxy_m2 = minimum_contact_proxy_m2
))

patch_network <- load_or_build_checkpoint(
  file.path(problem_dir, "existing_patch_network.rds"),
  patch_signature,
  function() {
    candidate_contact_buffers <- terra::buffer(
      candidate_polygons,
      width = contact_buffer_m
    )
    contacts <- terra::intersect(
      candidate_contact_buffers,
      patches_near_domain[, "patch_id"]
    )
    if (!nrow(contacts)) {
      stop("No candidate-to-current-patch contacts were found.")
    }

    contacts$contact_proxy_m2 <- terra::expanse(
      contacts,
      unit = "m",
      transform = FALSE
    )
    links <- contacts |>
      as.data.frame() |>
      transmute(
        coarse_cell_id = as.integer(coarse_cell_id),
        patch_id = as.integer(patch_id),
        contact_proxy_m2 = as.numeric(contact_proxy_m2)
      ) |>
      group_by(coarse_cell_id, patch_id) |>
      summarise(
        contact_proxy_m2 = sum(contact_proxy_m2),
        .groups = "drop"
      ) |>
      filter(
        contact_proxy_m2 >= minimum_contact_proxy_m2,
        coarse_cell_id %in% candidate_cells
      ) |>
      arrange(coarse_cell_id, patch_id)

    linked_patch_areas <- links |>
      distinct(patch_id) |>
      left_join(patch_dictionary_full, by = "patch_id")
    if (anyNA(linked_patch_areas$ExistingPatchAreaHa)) {
      missing_ids <- linked_patch_areas$patch_id[
        is.na(linked_patch_areas$ExistingPatchAreaHa)
      ]
      stop(
        "Missing full national area for patch IDs: ",
        paste(head(missing_ids, 20), collapse = ", "),
        call. = FALSE
      )
    }

    list(
      patches = linked_patch_areas |> arrange(patch_id),
      links = links
    )
  }
)

patch_table <- patch_network$patches
cell_patch_links <- patch_network$links |>
  select(coarse_cell_id, patch_id) |>
  distinct()
existing_patches <- patch_table$patch_id
patch_area_by_subarea <- patch_area_by_subarea_all |>
  filter(patch_id %in% existing_patches) |>
  arrange(subarea_id, patch_id)

stopifnot(
  nrow(patch_table) > 0L,
  nrow(cell_patch_links) > 0L,
  all(cell_patch_links$coarse_cell_id %in% candidate_cells),
  all(cell_patch_links$patch_id %in% existing_patches),
  all(patch_table$ExistingPatchAreaHa > 0),
  !anyDuplicated(patch_table$patch_id),
  !anyDuplicated(cell_patch_links)
)

patch_area_crosscheck <- patch_area_by_subarea |>
  group_by(patch_id) |>
  summarise(
    ExistingPatchAreaInsidePlanningHa =
      sum(ExistingPatchAreaHaBySubarea),
    .groups = "drop"
  ) |>
  right_join(patch_table, by = "patch_id") |>
  mutate(
    ExistingPatchAreaInsidePlanningHa = coalesce(
      ExistingPatchAreaInsidePlanningHa,
      0
    ),
    ExistingPatchAreaOutsidePlanningHa = pmax(
      0,
      ExistingPatchAreaHa - ExistingPatchAreaInsidePlanningHa
    ),
    OutsideSharePercent =
      100 * ExistingPatchAreaOutsidePlanningHa / ExistingPatchAreaHa
  )

if (any(
  patch_area_crosscheck$ExistingPatchAreaInsidePlanningHa >
    patch_area_crosscheck$ExistingPatchAreaHa + 0.1
)) {
  stop(
    "A clipped patch area exceeds its full national dictionary area.",
    call. = FALSE
  )
}

# Keep the scopes explicit. national_current_by_subarea contains every
# national existing-nature patch fragment inside each delområde, whereas the
# certificate ledger below contains only patches that can actually be reached
# from at least one candidate restoration cell under the contact convention.
# The difference is expected and must never be silently counted as reachable
# certificate area.
linked_patch_area_by_subarea <- patch_area_by_subarea |>
  group_by(subarea_id) |>
  summarise(
    LinkedCertificatePatchNatureHa =
      sum(ExistingPatchAreaHaBySubarea),
    .groups = "drop"
  )

patch_ledger_scope_audit <- national_current_by_subarea |>
  left_join(linked_patch_area_by_subarea, by = "subarea_id") |>
  mutate(
    LinkedCertificatePatchNatureHa = coalesce(
      LinkedCertificatePatchNatureHa,
      0
    ),
    UnlinkedNationalPatchNatureHa =
      NationalPatchNatureHa - LinkedCertificatePatchNatureHa,
    LinkedDoesNotExceedAll =
      LinkedCertificatePatchNatureHa <= NationalPatchNatureHa + 1e-6,
    MaterialPatchVsLanduseDifference =
      abs(DifferencePercent) > 5
  ) |>
  arrange(subarea_id)

if (!all(patch_ledger_scope_audit$LinkedDoesNotExceedAll)) {
  stop(
    "Linked certificate-patch area exceeds the all-patch ledger. ",
    "Delete the 08a output folder and rebuild its checkpoints.",
    call. = FALSE
  )
}

linked_inside_from_fragments <-
  sum(patch_area_by_subarea$ExistingPatchAreaHaBySubarea)
linked_inside_from_patch_audit <-
  sum(patch_area_crosscheck$ExistingPatchAreaInsidePlanningHa)
if (abs(
  linked_inside_from_fragments - linked_inside_from_patch_audit
) > 1e-6) {
  stop(
    "Linked patch fragments do not reconcile to the inside/outside audit.",
    call. = FALSE
  )
}

readr::write_csv(
  patch_table,
  file.path(problem_dir, "existing_patch_dictionary_full_national.csv")
)
readr::write_csv(
  patch_area_by_subarea,
  file.path(problem_dir, "existing_patch_area_by_subarea.csv")
)
readr::write_csv(
  patch_area_crosscheck,
  file.path(problem_dir, "existing_patch_inside_outside_audit.csv")
)
readr::write_csv(
  national_current_by_subarea,
  file.path(problem_dir, "national_patches_vs_landuse_budget_audit.csv")
)
readr::write_csv(
  patch_ledger_scope_audit,
  file.path(problem_dir, "patch_ledger_scope_audit.csv")
)
readr::write_csv(
  cell_patch_links,
  file.path(problem_dir, "candidate_to_existing_patch_links.csv")
)

print(knitr::kable(national_current_by_subarea, digits = 3))
print(knitr::kable(patch_ledger_scope_audit, digits = 3))
knitr::kable(tibble(
  retained_existing_patches = nrow(patch_table),
  full_national_area_of_retained_patches_ha =
    sum(patch_table$ExistingPatchAreaHa),
  area_of_retained_patches_inside_planning_ha =
    sum(patch_area_crosscheck$ExistingPatchAreaInsidePlanningHa),
  area_of_retained_patches_outside_planning_ha =
    sum(patch_area_crosscheck$ExistingPatchAreaOutsidePlanningHa),
  linked_full_patches_already_at_least_5000_ha = sum(
    patch_table$ExistingPatchAreaHa >= global_network_minimum_ha
  ),
  candidate_patch_contacts = nrow(cell_patch_links),
  candidate_cells_touching_current_nature =
    n_distinct(cell_patch_links$coarse_cell_id)
), digits = 2)


# Source chunk: potential-graph
potential_vertices <- bind_rows(
  candidate_table |>
    transmute(
      vertex = paste0("c_", coarse_cell_id),
      node_type = "cell",
      node_id = coarse_cell_id,
      FlowAreaHa = EligibleAreaHa
    ),
  patch_table |>
    transmute(
      vertex = paste0("p_", patch_id),
      node_type = "patch",
      node_id = patch_id,
      FlowAreaHa = ExistingPatchAreaHa
    )
)

potential_edges <- bind_rows(
  candidate_edges |>
    transmute(from = paste0("c_", i), to = paste0("c_", j)),
  cell_patch_links |>
    transmute(
      from = paste0("c_", coarse_cell_id),
      to = paste0("p_", patch_id)
    )
) |>
  distinct()

potential_graph <- igraph::graph_from_data_frame(
  potential_edges,
  directed = FALSE,
  vertices = potential_vertices
)
potential_membership_raw <- igraph::components(
  potential_graph
)$membership
potential_membership <- tibble(
  vertex = names(potential_membership_raw),
  potential_component = as.integer(potential_membership_raw)
)

potential_vertices <- potential_vertices |>
  left_join(potential_membership, by = "vertex")

potential_component_summary <- potential_vertices |>
  group_by(potential_component) |>
  summarise(
    ComponentMaxAreaHa = sum(FlowAreaHa),
    candidate_cells = sum(node_type == "cell"),
    existing_patches = sum(node_type == "patch"),
    .groups = "drop"
  )

cell_component_bounds <- potential_vertices |>
  filter(node_type == "cell") |>
  select(coarse_cell_id = node_id, potential_component) |>
  left_join(potential_component_summary, by = "potential_component") |>
  select(coarse_cell_id, ComponentMaxAreaHa)

patch_component_bounds <- potential_vertices |>
  filter(node_type == "patch") |>
  select(patch_id = node_id, potential_component) |>
  left_join(potential_component_summary, by = "potential_component") |>
  select(patch_id, ComponentMaxAreaHa)

# Area coefficients for the Global and four Subarea certificates.
network_cell_area <- bind_rows(
  candidate_table |>
    transmute(
      network_name = "Global",
      coarse_cell_id,
      NetworkCellAreaHa = EligibleAreaHa
    ),
  restorable_area_by_subarea_cell |>
    transmute(
      network_name = paste0("Subarea", subarea_id),
      coarse_cell_id,
      NetworkCellAreaHa = RestorableAreaHaBySubarea
    )
) |>
  filter(NetworkCellAreaHa > 0)

network_patch_area <- bind_rows(
  patch_table |>
    transmute(
      network_name = "Global",
      patch_id,
      NetworkPatchAreaHa = ExistingPatchAreaHa
    ),
  patch_area_by_subarea |>
    transmute(
      network_name = paste0("Subarea", subarea_id),
      patch_id,
      NetworkPatchAreaHa = ExistingPatchAreaHaBySubarea
    )
) |>
  filter(NetworkPatchAreaHa > 0)

potential_component_network_area <- bind_rows(
  network_cell_area |>
    left_join(
      potential_vertices |>
        filter(node_type == "cell") |>
        transmute(
          coarse_cell_id = node_id,
          potential_component
        ),
      by = "coarse_cell_id"
    ) |>
    transmute(
      network_name,
      potential_component,
      CertificateAreaHa = NetworkCellAreaHa
    ),
  network_patch_area |>
    left_join(
      potential_vertices |>
        filter(node_type == "patch") |>
        transmute(patch_id = node_id, potential_component),
      by = "patch_id"
    ) |>
    transmute(
      network_name,
      potential_component,
      CertificateAreaHa = NetworkPatchAreaHa
    )
) |>
  group_by(network_name, potential_component) |>
  summarise(
    PotentialCertificateAreaHa = sum(CertificateAreaHa),
    .groups = "drop"
  ) |>
  left_join(
    potential_component_summary,
    by = "potential_component"
  ) |>
  mutate(
    RequiredAreaHa = unname(network_minimums[network_name]),
    CanAnchor = candidate_cells > 0 & existing_patches > 0,
    PassesPotentialArea =
      PotentialCertificateAreaHa >= RequiredAreaHa - 1e-8,
    FeasiblePotentialComponent = CanAnchor & PassesPotentialArea
  )

potential_network_audit <- tibble(
  network_name = required_networks,
  RequiredAreaHa = unname(network_minimums[required_networks])
) |>
  left_join(
    potential_component_network_area |>
      group_by(network_name) |>
      summarise(
        RequiredAreaHa = first(RequiredAreaHa),
        feasible_potential_components =
          sum(FeasiblePotentialComponent),
        largest_potential_certificate_ha = max(
          PotentialCertificateAreaHa[CanAnchor],
          na.rm = TRUE
        ),
        .groups = "drop"
      ),
    by = c("network_name", "RequiredAreaHa")
  )

if (
  anyNA(potential_network_audit$feasible_potential_components) ||
  any(potential_network_audit$feasible_potential_components < 1L)
) {
  stop(
    "At least one requested network is impossible in the potential graph. ",
    "Inspect potential_network_feasibility_audit.csv.",
    call. = FALSE
  )
}

# Exact network-specific reduction: retain only potential components whose
# counted area can satisfy the corresponding certificate and which contain at
# least one candidate cell and one existing patch.
feasible_potential_components <- potential_component_network_area |>
  filter(FeasiblePotentialComponent) |>
  distinct(network_name, potential_component)

cell_potential_membership <- potential_vertices |>
  filter(node_type == "cell") |>
  transmute(
    coarse_cell_id = as.integer(node_id),
    potential_component
  )
patch_potential_membership <- potential_vertices |>
  filter(node_type == "patch") |>
  transmute(
    patch_id = as.integer(node_id),
    potential_component
  )

eligible_network_cells <- feasible_potential_components |>
  inner_join(cell_potential_membership, by = "potential_component") |>
  select(network_name, coarse_cell_id, potential_component) |>
  distinct() |>
  arrange(network_name, coarse_cell_id)

eligible_network_patches <- feasible_potential_components |>
  inner_join(patch_potential_membership, by = "potential_component") |>
  select(network_name, patch_id, potential_component) |>
  distinct() |>
  arrange(network_name, patch_id)

eligible_i <- eligible_network_cells |>
  transmute(network_name, i = coarse_cell_id)
eligible_j <- eligible_network_cells |>
  transmute(network_name, j = coarse_cell_id)

eligible_network_edges <- candidate_edges |>
  inner_join(eligible_i, by = "i") |>
  inner_join(eligible_j, by = c("network_name", "j")) |>
  select(network_name, i, j) |>
  distinct()

eligible_network_arcs <- eligible_network_edges |>
  bind_rows(
    eligible_network_edges |>
      transmute(network_name, i = j, j = i)
  ) |>
  distinct() |>
  arrange(network_name, i, j)

eligible_network_cell_patch_edges <- cell_patch_links |>
  inner_join(
    eligible_network_cells |>
      select(network_name, coarse_cell_id),
    by = "coarse_cell_id"
  ) |>
  inner_join(
    eligible_network_patches |>
      select(network_name, patch_id),
    by = c("network_name", "patch_id")
  ) |>
  select(network_name, coarse_cell_id, patch_id) |>
  distinct() |>
  arrange(network_name, coarse_cell_id, patch_id)

network_reduction_audit <- tibble(
  network_name = required_networks
) |>
  left_join(
    eligible_network_cells |>
      count(network_name, name = "eligible_cells"),
    by = "network_name"
  ) |>
  left_join(
    eligible_network_patches |>
      count(network_name, name = "eligible_patches"),
    by = "network_name"
  ) |>
  left_join(
    eligible_network_arcs |>
      count(network_name, name = "directed_cell_arcs"),
    by = "network_name"
  ) |>
  left_join(
    eligible_network_cell_patch_edges |>
      count(network_name, name = "cell_patch_edges"),
    by = "network_name"
  ) |>
  mutate(
    across(
      c(
        eligible_cells,
        eligible_patches,
        directed_cell_arcs,
        cell_patch_edges
      ),
      ~ coalesce(.x, 0L)
    )
  )

if (any(
  network_reduction_audit$eligible_cells < 1L |
    network_reduction_audit$eligible_patches < 1L |
    network_reduction_audit$cell_patch_edges < 1L
)) {
  stop(
    "Internal error while reducing network variables to feasible components.",
    call. = FALSE
  )
}

patch_connection_normalizer <- max(
  1,
  max(potential_component_summary$existing_patches) - 1
)

readr::write_csv(
  potential_component_summary,
  file.path(problem_dir, "potential_graph_components.csv")
)
readr::write_csv(
  potential_component_network_area,
  file.path(problem_dir, "potential_component_network_area.csv")
)
readr::write_csv(
  potential_network_audit,
  file.path(problem_dir, "potential_network_feasibility_audit.csv")
)
readr::write_csv(
  feasible_potential_components,
  file.path(problem_dir, "feasible_potential_components.csv")
)
readr::write_csv(
  network_reduction_audit,
  file.path(problem_dir, "network_variable_reduction_audit.csv")
)

print(knitr::kable(potential_network_audit, digits = 2))
knitr::kable(network_reduction_audit, digits = 0)


# Source chunk: model-actions
allowed_actions <- tidyr::crossing(
  landuse = nature_landuses,
  candidate_table |>
    select(coarse_cell_id, IsProductionForestCell)
) |>
  filter(
    IsProductionForestCell == 0 |
      landuse %in% forest_landuses
  ) |>
  select(landuse, coarse_cell_id) |>
  arrange(landuse, coarse_cell_id)

objective_table <- allowed_actions |>
  left_join(
    biodiversity_table,
    by = c("landuse", "coarse_cell_id")
  ) |>
  left_join(
    existing_nature_contact,
    by = c("landuse", "coarse_cell_id")
  )

stopifnot(
  n_distinct(allowed_actions$coarse_cell_id) ==
    length(candidate_cells),
  !anyDuplicated(allowed_actions),
  !anyNA(objective_table),
  all(is.finite(objective_table$BioDiversity)),
  all(is.finite(objective_table$ExistingNatureContact))
)

biodiversity_normalizer <- objective_table |>
  group_by(coarse_cell_id) |>
  summarise(best = max(BioDiversity), .groups = "drop") |>
  summarise(total = sum(best)) |>
  pull(total) |>
  max(1)

contiguity_normalizer <- max(1, nrow(candidate_edges))

existing_nature_contact_normalizer <- objective_table |>
  group_by(coarse_cell_id) |>
  summarise(
    best = max(ExistingNatureContact),
    .groups = "drop"
  ) |>
  summarise(total = sum(best)) |>
  pull(total) |>
  max(1)

allowed_i <- allowed_actions |>
  transmute(landuse, i = coarse_cell_id)
allowed_j <- allowed_actions |>
  transmute(landuse, j = coarse_cell_id)

same_habitat_edges <- tidyr::crossing(
  landuse = nature_landuses,
  candidate_edges
) |>
  semi_join(allowed_i, by = c("landuse", "i")) |>
  semi_join(allowed_j, by = c("landuse", "j")) |>
  arrange(landuse, i, j)

best_forest_action <- objective_table |>
  filter(landuse %in% forest_landuses) |>
  mutate(
    local_start_score =
      BioDiversity / biodiversity_normalizer +
      spatial_contiguity_bonus *
        existing_to_new_adjacency_weight *
        ExistingNatureContact /
          existing_nature_contact_normalizer
  ) |>
  group_by(coarse_cell_id) |>
  arrange(
    desc(local_start_score),
    desc(BioDiversity),
    landuse,
    .by_group = TRUE
  ) |>
  slice_head(n = 1L) |>
  ungroup() |>
  select(
    coarse_cell_id,
    repaired_forest_landuse = landuse
  )

warm_start_assignment <- warmup_selected |>
  rename(original_landuse = landuse) |>
  left_join(
    candidate_table |>
      select(coarse_cell_id, IsProductionForestCell),
    by = "coarse_cell_id"
  ) |>
  left_join(best_forest_action, by = "coarse_cell_id") |>
  mutate(
    forest_rule_repair =
      IsProductionForestCell == 1L &
      !original_landuse %in% forest_landuses,
    landuse = if_else(
      forest_rule_repair,
      repaired_forest_landuse,
      original_landuse
    )
  ) |>
  select(
    coarse_cell_id,
    original_landuse,
    landuse,
    IsProductionForestCell,
    forest_rule_repair
  ) |>
  arrange(coarse_cell_id)

invalid_start_actions <- warm_start_assignment |>
  anti_join(
    allowed_actions,
    by = c("landuse", "coarse_cell_id")
  )
if (nrow(invalid_start_actions)) {
  stop("Warm-start action repair failed.", call. = FALSE)
}

warm_start_subarea_reference_audit <- restorable_area_by_subarea_cell |>
  filter(
    coarse_cell_id %in% warm_start_assignment$coarse_cell_id
  ) |>
  group_by(subarea_id) |>
  summarise(
    WarmStartRestorationHa = sum(RestorableAreaHaBySubarea),
    .groups = "drop"
  ) |>
  right_join(subarea_budget, by = "subarea_id") |>
  mutate(
    WarmStartRestorationHa = coalesce(WarmStartRestorationHa, 0),
    WarmStartFinalNatureHa =
      CurrentNatureAreaHa + WarmStartRestorationHa,
    WarmStartFinalNaturePercent =
      100 * WarmStartFinalNatureHa / PolygonAreaHa,
    ReferenceThirtyPercentDifferenceHa =
      WarmStartFinalNatureHa - FinalNatureTargetHa
  ) |>
  arrange(subarea_id)

warm_start_global_area_audit <-
  warm_start_subarea_reference_audit |>
  summarise(
    WarmStartRestorationHa = sum(WarmStartRestorationHa)
  ) |>
  bind_cols(global_budget) |>
  mutate(
    WarmStartFinalNatureHa =
      GlobalCurrentNatureAreaHa + WarmStartRestorationHa,
    WarmStartFinalNaturePercent =
      100 * WarmStartFinalNatureHa / GlobalPolygonAreaHa,
    ExcessAboveThirtyPercentHa =
      WarmStartFinalNatureHa - GlobalFinalNatureTargetHa,
    PassesGlobalThirtyPercentMinimum =
      WarmStartFinalNatureHa >= GlobalFinalNatureTargetHa - 1e-6,
    PassesGlobalThirtyPercentMaximum =
      WarmStartFinalNatureHa <= GlobalMaximumFinalNatureHa + 1e-6,
    PassesGlobalThirtyPercentBudget =
      PassesGlobalThirtyPercentMinimum &
        PassesGlobalThirtyPercentMaximum
  )

if (!warm_start_global_area_audit$PassesGlobalThirtyPercentBudget[[1]]) {
  stop(
    paste0(
      "The upstream footprint does not satisfy the whole-polygon 30% ",
      "budget and its configured upper tolerance."
    ),
    call. = FALSE
  )
}

readr::write_csv(
  warm_start_assignment,
  file.path(problem_dir, "warm_start_assignment_and_forest_repairs.csv")
)
readr::write_csv(
  warm_start_global_area_audit,
  file.path(problem_dir, "warm_start_30pct_audit.csv")
)
readr::write_csv(
  warm_start_global_area_audit,
  file.path(problem_dir, "warm_start_global_30pct_audit.csv")
)
readr::write_csv(
  warm_start_subarea_reference_audit,
  file.path(problem_dir, "warm_start_subarea_reference_audit.csv")
)

knitr::kable(tibble(
  warm_start_selected_cells = nrow(warm_start_assignment),
  production_forest_dominated_selected_cells = sum(
    warm_start_assignment$IsProductionForestCell == 1L
  ),
  destination_labels_repaired = sum(
    warm_start_assignment$forest_rule_repair
  ),
  warm_start_global_nature_percent =
    warm_start_global_area_audit$WarmStartFinalNaturePercent,
  passes_global_30pct_budget =
    warm_start_global_area_audit$PassesGlobalThirtyPercentBudget
))


# Source chunk: warm-start-components
warm_start_ids <- warm_start_assignment$coarse_cell_id
warm_start_patch_ids <- cell_patch_links |>
  filter(coarse_cell_id %in% warm_start_ids) |>
  distinct(patch_id) |>
  pull(patch_id) |>
  sort()

warm_start_vertices <- bind_rows(
  candidate_table |>
    filter(coarse_cell_id %in% warm_start_ids) |>
    transmute(
      vertex = paste0("c_", coarse_cell_id),
      node_type = "cell",
      node_id = coarse_cell_id,
      FlowAreaHa = EligibleAreaHa
    ),
  patch_table |>
    filter(patch_id %in% warm_start_patch_ids) |>
    transmute(
      vertex = paste0("p_", patch_id),
      node_type = "patch",
      node_id = patch_id,
      FlowAreaHa = ExistingPatchAreaHa
    )
)

warm_start_edges <- bind_rows(
  candidate_edges |>
    filter(i %in% warm_start_ids, j %in% warm_start_ids) |>
    transmute(from = paste0("c_", i), to = paste0("c_", j)),
  cell_patch_links |>
    filter(
      coarse_cell_id %in% warm_start_ids,
      patch_id %in% warm_start_patch_ids
    ) |>
    transmute(
      from = paste0("c_", coarse_cell_id),
      to = paste0("p_", patch_id)
    )
) |>
  distinct()

warm_start_graph <- igraph::graph_from_data_frame(
  warm_start_edges,
  directed = FALSE,
  vertices = warm_start_vertices
)
warm_membership_raw <- igraph::components(warm_start_graph)$membership
warm_membership <- tibble(
  vertex = names(warm_membership_raw),
  warm_component = as.integer(warm_membership_raw)
)
warm_start_vertices <- warm_start_vertices |>
  left_join(warm_membership, by = "vertex")

warm_component_summary <- warm_start_vertices |>
  group_by(warm_component) |>
  summarise(
    selected_cells = sum(node_type == "cell"),
    existing_patches = sum(node_type == "patch"),
    FlowAreaHa = sum(FlowAreaHa),
    .groups = "drop"
  )

warm_cell_membership <- warm_start_vertices |>
  filter(node_type == "cell") |>
  transmute(
    coarse_cell_id = node_id,
    warm_component
  )
warm_patch_membership <- warm_start_vertices |>
  filter(node_type == "patch") |>
  transmute(patch_id = node_id, warm_component)

warm_component_network_area <- bind_rows(
  network_cell_area |>
    filter(coarse_cell_id %in% warm_start_ids) |>
    left_join(warm_cell_membership, by = "coarse_cell_id") |>
    transmute(
      network_name,
      warm_component,
      CertificateAreaHa = NetworkCellAreaHa
    ),
  network_patch_area |>
    filter(patch_id %in% warm_start_patch_ids) |>
    left_join(warm_patch_membership, by = "patch_id") |>
    transmute(
      network_name,
      warm_component,
      CertificateAreaHa = NetworkPatchAreaHa
    )
) |>
  group_by(network_name, warm_component) |>
  summarise(
    CertifiedAreaHa = sum(CertificateAreaHa),
    .groups = "drop"
  ) |>
  left_join(warm_component_summary, by = "warm_component") |>
  mutate(
    RequiredAreaHa = unname(network_minimums[network_name]),
    HasExistingRoot = existing_patches > 0,
    HasNewRestoration = selected_cells > 0,
    PassesArea = CertifiedAreaHa >= RequiredAreaHa - 1e-8,
    FeasibleWarmCertificate =
      HasExistingRoot & HasNewRestoration & PassesArea
  )

# Preserve the diagnostic, but do not stop when the warm footprint is not yet
# feasible. That is precisely the condition notebook 08a is designed to repair.
readr::write_csv(
  warm_component_network_area,
  file.path(problem_dir, "warm_start_component_network_audit.csv")
)

warm_start_certificate_availability <- tibble(
  network_name = required_networks
) |>
  left_join(
    warm_component_network_area |>
      group_by(network_name) |>
      summarise(
        feasible_warm_components = sum(FeasibleWarmCertificate),
        largest_anchored_warm_certificate_ha = max(c(
          0,
          CertifiedAreaHa[HasExistingRoot & HasNewRestoration]
        )),
        can_seed_complete_start = any(FeasibleWarmCertificate),
        .groups = "drop"
      ),
    by = "network_name"
  ) |>
  mutate(
    feasible_warm_components = coalesce(
      feasible_warm_components,
      0L
    ),
    largest_anchored_warm_certificate_ha = coalesce(
      largest_anchored_warm_certificate_ha,
      0
    ),
    can_seed_complete_start = coalesce(
      can_seed_complete_start,
      FALSE
    ),
    RequiredAreaHa = unname(network_minimums[network_name]),
    shortfall_ha = pmax(
      0,
      RequiredAreaHa - largest_anchored_warm_certificate_ha
    )
  )

complete_feasible_mip_start_available <- all(
  warm_start_certificate_availability$can_seed_complete_start
)

effective_warm_start_strategy <-
  "complete_slack_start_then_lexicographic_repair"

warm_start_strategy_audit <- warm_start_certificate_availability |>
  mutate(
    configured_mode = "08a_feasible_network_repair",
    effective_strategy = effective_warm_start_strategy,
    potential_graph_passed = TRUE
  )

readr::write_csv(
  warm_start_strategy_audit,
  file.path(problem_dir, "warm_start_strategy_audit.csv")
)

chosen_certificate_components <- if (
  complete_feasible_mip_start_available
) {
  map_dfr(
    required_networks,
    function(k) {
      warm_component_network_area |>
        filter(
          network_name == k,
          FeasibleWarmCertificate
        ) |>
        arrange(desc(CertifiedAreaHa), warm_component) |>
        slice_head(n = 1L)
    }
  )
} else {
  warm_component_network_area |>
    slice(0)
}

readr::write_csv(
  chosen_certificate_components,
  file.path(problem_dir, "warm_start_chosen_certificates.csv")
)

print(knitr::kable(warm_start_strategy_audit, digits = 2))

if (complete_feasible_mip_start_available) {
  knitr::kable(
    chosen_certificate_components |>
      select(
        network_name,
        warm_component,
        RequiredAreaHa,
        CertifiedAreaHa,
        selected_cells,
        existing_patches
      ),
    digits = 2
  )
} else {
  cat(
    "Notebook 07 is not yet a complete feasible network start. ",
    "Notebook 08a will start from its valid global-30% footprint with explicit ",
    "certificate shortfalls, repair all five networks, and only then ",
    "minimize footprint changes and optimize ecological quality.\n"
  )
}
