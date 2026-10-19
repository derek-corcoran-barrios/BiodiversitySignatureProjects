# Source chunk: packages
library(terra)
library(dplyr)
library(tidyr)
library(purrr)
library(readr)
library(digest)
library(igraph)
library(knitr)


# Source chunk: names-and-scenarios
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

scenario_index <- tibble::tribble(
  ~landuse,          ~scenario,
  "ForestDryPoor",  "all_to_forest_dry_poor",
  "ForestDryRich",  "all_to_forest_dry_rich",
  "ForestWetPoor",  "all_to_forest_wet_poor",
  "ForestWetRich",  "all_to_forest_wet_rich",
  "OpenDryPoor",    "all_to_open_dry_poor",
  "OpenDryRich",    "all_to_open_dry_rich",
  "OpenWetPoor",    "all_to_open_wet_poor",
  "OpenWetRich",    "all_to_open_wet_rich"
)

stopifnot(
  setequal(nature_landuses, scenario_index$landuse),
  setequal(names(metric_weights), c("Richness", "PD", "RangeRarity")),
  abs(sum(metric_weights) - 1) < 1e-10,
  target_nature_proportion > 0,
  target_nature_proportion < 1,
  restoration_tolerance_ha >= 0,
  minimum_eligible_area_ha > 0,
  spatial_contiguity_bonus >= 0,
  new_to_new_adjacency_weight >= 0,
  existing_to_new_adjacency_weight >= 0,
  abs(
    new_to_new_adjacency_weight +
      existing_to_new_adjacency_weight - 1
  ) < 1e-10
)


# Source chunk: helpers
shapefile_components <- function(path) {
  if (tolower(tools::file_ext(path)) != "shp") {
    return(path)
  }

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

file_fingerprint <- function(paths) {
  paths <- normalizePath(paths, winslash = "/", mustWork = FALSE)
  info <- file.info(paths)

  tibble(
    path = paths,
    size = info$size,
    modified = as.numeric(info$mtime)
  )
}

load_or_build_checkpoint <- function(
    path,
    signature,
    builder,
    fallback_paths = character()
) {
  if (reuse_matching_checkpoints) {
    candidate_paths <- unique(c(path, fallback_paths))

    for (candidate_path in candidate_paths[file.exists(candidate_paths)]) {
      checkpoint <- readRDS(candidate_path)

      if (
        is.list(checkpoint) &&
        identical(checkpoint$signature, signature) &&
        !is.null(checkpoint$value)
      ) {
        return(checkpoint$value)
      }
    }
  }

  value <- builder()
  saveRDS(list(signature = signature, value = value), path)
  value
}

same_crs <- function(x, y) {
  isTRUE(terra::same.crs(x, y))
}

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

  if (is.null(lev) || nrow(lev) == 0L) {
    stop("The fine land-use raster has no category table.", call. = FALSE)
  }

  numeric_columns <- names(lev)[vapply(lev, is.numeric, logical(1))]
  label_columns <- names(lev)[vapply(
    lev,
    function(z) is.character(z) || is.factor(z),
    logical(1)
  )]

  if (!length(numeric_columns) || !length(label_columns)) {
    stop(
      "Could not identify value and label columns in the category table.",
      call. = FALSE
    )
  }

  value_column <- numeric_columns[[1]]
  label_column <- label_columns[[1]]
  labels <- as.character(lev[[label_column]])
  missing_labels <- setdiff(requested_labels, labels)

  if (length(missing_labels)) {
    stop(
      paste0(
        "Missing land-use classes: ",
        paste(missing_labels, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  stats::setNames(
    lev[[value_column]][match(requested_labels, labels)],
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
  if (!length(x)) {
    return(character())
  }

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
    c(
      paste("set", name, ":="),
      chunk_text(as.character(values)),
      ";",
      ""
    )
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
  append_dat(
    path,
    c(paste("param", name, ":="), rows, ";", "")
  )
}

write_ampl_matrix <- function(path, name, row_key, column_key, value) {
  rows <- paste(row_key, column_key, ampl_number(value))
  append_dat(
    path,
    c(paste("param", name, ":="), rows, ";", "")
  )
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

  if (!copied) {
    stop("Could not write generated file: ", final_path, call. = FALSE)
  }

  invisible(TRUE)
}

write_lines_if_changed <- function(lines, path) {
  temporary_path <- paste0(path, ".building")
  writeLines(lines, temporary_path)
  commit_generated_file(temporary_path, path)
}

ampl_path <- function(path) {
  chartr(intToUtf8(92L), "/", path)
}


# Source chunk: input-audit
dir.create(problem_dir, recursive = TRUE, showWarnings = FALSE)

metric_file_index <- tidyr::crossing(
  metric = names(metric_weights),
  scenario_index
) |>
  mutate(
    path = file.path(
      metrics_dir,
      paste0(metric, "_", scenario, ".tif")
    )
  )

metric_current_index <- tibble(
  metric = names(metric_weights),
  scenario = "current",
  landuse = NA_character_,
  path = file.path(
    metrics_dir,
    paste0(names(metric_weights), "_current.tif")
  )
)

required_core_paths <- c(
  fine_landuse_path,
  coarse_landuse_path,
  shapefile_components(subareas_path),
  model_path,
  metric_file_index$path,
  metric_current_index$path
)

stop_for_missing_files(required_core_paths, "optimization input file(s)")

input_audit <- file_fingerprint(required_core_paths)
knitr::kable(input_audit)


# Source chunk: read-subareas
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
coarse_nominal_area_ha <- prod(terra::res(coarse_landuse_full)) / 10000

stopifnot(
  abs(fine_cell_area_ha - 0.01) < 1e-8,
  abs(coarse_nominal_area_ha - 4) < 1e-8
)

subareas <- terra::vect(subareas_path) |>
  terra::makeValid()

if (nrow(subareas) != expected_number_of_subareas) {
  stop(
    "Expected ", expected_number_of_subareas,
    " subarea features but found ", nrow(subareas), ".",
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
    stop(
      "subarea_label_field is not present in the shapefile: ",
      subarea_label_field,
      call. = FALSE
    )
  }

  subarea_labels <- as.character(
    as.data.frame(subareas)[[subarea_label_field]]
  )
  if (anyNA(subarea_labels) || any(!nzchar(subarea_labels))) {
    stop("Subarea labels contain missing or empty values.", call. = FALSE)
  }
  subarea_labels <- make.unique(subarea_labels)
}

subareas$subarea_id <- seq_len(nrow(subareas))
subareas$subarea_label <- subarea_labels

polygon_areas_ha <- as.numeric(
  terra::expanse(subareas, unit = "ha", transform = FALSE)
)
planning_domain <- terra::aggregate(subareas)
union_area_ha <- sum(
  terra::expanse(planning_domain, unit = "ha", transform = FALSE)
)
overlap_area_ha <- sum(polygon_areas_ha) - union_area_ha

if (overlap_area_ha > 0.1) {
  stop(
    "The four subareas overlap by approximately ",
    round(overlap_area_ha, 3),
    " ha. Targets would be double-counted; repair the polygons first.",
    call. = FALSE
  )
}

if (
  length(expected_subarea_areas_ha) == expected_number_of_subareas &&
  any(
    abs(sort(polygon_areas_ha) - sort(expected_subarea_areas_ha)) >
      polygon_area_audit_tolerance_ha
  )
) {
  stop(
    "The four polygon areas do not match the supplied Gudenå values within ",
    polygon_area_audit_tolerance_ha,
    " ha. Confirm subareas_path and its sidecar files.",
    call. = FALSE
  )
}

if (
  abs(sum(polygon_areas_ha) - expected_total_area_ha) >
    polygon_area_audit_tolerance_ha
) {
  stop(
    "Total polygon area is ", round(sum(polygon_areas_ha), 2),
    " ha rather than approximately ", expected_total_area_ha, " ha.",
    call. = FALSE
  )
}

subarea_dictionary <- tibble(
  subarea_id = seq_len(nrow(subareas)),
  subarea_label = subarea_labels,
  PolygonAreaHa = polygon_areas_ha,
  ThirtyPercentHa = target_nature_proportion * polygon_areas_ha
)

readr::write_csv(
  subarea_dictionary,
  file.path(problem_dir, "subarea_dictionary.csv")
)

knitr::kable(subarea_dictionary, digits = 2)


# Source chunk: align-grids
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

# IDs are local to this cropped 200 m grid. This makes them valid both as
# raster values and as cell indices for terra::adjacent().
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

fine_domain <- terra::ifel(is.na(fine_landuse), NA, 1)
names(fine_domain) <- "DomainCells"

fine_without_coarse_id <- global_number(
  terra::ifel(
    !is.na(fine_landuse) & is.na(coarse_id_on_fine),
    1,
    0
  ),
  "sum"
)

fine_without_subarea_id <- global_number(
  terra::ifel(
    !is.na(fine_landuse) & is.na(subarea_id_on_fine),
    1,
    0
  ),
  "sum"
)

if (fine_without_coarse_id > 0) {
  stop(
    fine_without_coarse_id,
    " valid 10 m cells have no 200 m ID. Rebuild/extend ",
    coarse_landuse_path,
    " before continuing.",
    call. = FALSE
  )
}

if (fine_without_subarea_id > 0) {
  stop(
    fine_without_subarea_id,
    " valid 10 m cells have no subarea ID.",
    call. = FALSE
  )
}

modelled_area_by_subarea <- zonal_sum_table(
  fine_domain,
  subarea_id_on_fine,
  "subarea_id"
) |>
  transmute(
    subarea_id = as.integer(subarea_id),
    ModelledRasterAreaHa = DomainCells * fine_cell_area_ha
  )

area_representation_audit <- subarea_dictionary |>
  left_join(modelled_area_by_subarea, by = "subarea_id") |>
  mutate(
    RasterMinusPolygonHa = ModelledRasterAreaHa - PolygonAreaHa,
    RasterMinusPolygonPercent =
      100 * RasterMinusPolygonHa / PolygonAreaHa
  )

readr::write_csv(
  area_representation_audit,
  file.path(problem_dir, "polygon_vs_10m_area_audit.csv")
)

knitr::kable(area_representation_audit, digits = 3)


# Source chunk: cell-composition
all_required_classes <- c(
  nature_landuses,
  "Agriculture",
  "ProductionForest",
  "Urban",
  "Other"
)

landuse_codes <- category_codes(fine_landuse, all_required_classes)
nature_codes <- unname(landuse_codes[nature_landuses])
agriculture_code <- unname(landuse_codes["Agriculture"])
production_forest_code <- unname(landuse_codes["ProductionForest"])
eligible_codes <- c(agriculture_code, production_forest_code)

fine_nature <- indicator_for_codes(fine_landuse, nature_codes)
fine_agriculture_landuse <- indicator_for_codes(
  fine_landuse,
  agriculture_code
)
fine_production_forest_landuse <- indicator_for_codes(
  fine_landuse,
  production_forest_code
)
fine_eligible_landuse <- indicator_for_codes(
  fine_landuse,
  eligible_codes
)

names(fine_nature) <- "NatureCells"
names(fine_agriculture_landuse) <- "AgricultureLanduseCells"
names(fine_production_forest_landuse) <-
  "ProductionForestLanduseCells"

biodiversity_support_path <- bind_rows(
  metric_current_index,
  metric_file_index
) |>
  filter(
    .data$metric == biodiversity_support_metric,
    .data$scenario == biodiversity_support_scenario
  ) |>
  pull(path)

if (length(biodiversity_support_path) != 1L) {
  stop(
    "Expected exactly one raster for the biodiversity support mask.",
    call. = FALSE
  )
}

if (restrict_eligibility_to_biodiversity_support) {
  biodiversity_support_raster <- terra::rast(
    biodiversity_support_path
  )
  biodiversity_support_raster <- align_to_template(
    biodiversity_support_raster,
    fine_landuse,
    method = "bilinear"
  )

  fine_biodiversity_support <- terra::ifel(
    is.na(fine_landuse),
    NA,
    terra::ifel(!is.na(biodiversity_support_raster), 1, 0)
  )
} else {
  fine_biodiversity_support <- terra::ifel(
    is.na(fine_landuse),
    NA,
    1
  )
}

fine_uncovered_eligible <-
  fine_eligible_landuse * (1 - fine_biodiversity_support)
fine_agriculture <-
  fine_agriculture_landuse * fine_biodiversity_support
fine_production_forest <-
  fine_production_forest_landuse * fine_biodiversity_support
fine_eligible <- fine_agriculture + fine_production_forest
fine_eligible_mask <- terra::ifel(fine_eligible == 1, 1, NA)

names(fine_agriculture) <- "AgricultureCells"
names(fine_production_forest) <- "ProductionForestCells"
names(fine_eligible) <- "EligibleCells"
names(fine_uncovered_eligible) <- "UncoveredEligibleCells"

composition_stack <- c(
  fine_nature,
  fine_agriculture,
  fine_production_forest,
  fine_eligible,
  fine_uncovered_eligible
)

composition_signature <- digest::digest(
  list(
    files = file_fingerprint(c(
      fine_landuse_path,
      coarse_landuse_path,
      shapefile_components(subareas_path),
      biodiversity_support_path
    )),
    subarea_areas = polygon_areas_ha,
    restrict_eligibility_to_biodiversity_support =
      restrict_eligibility_to_biodiversity_support,
    biodiversity_support_metric = biodiversity_support_metric,
    biodiversity_support_scenario = biodiversity_support_scenario,
    minimum_eligible_area_ha = minimum_eligible_area_ha
  )
)

cell_table <- load_or_build_checkpoint(
  file.path(problem_dir, "coarse_cell_composition.rds"),
  composition_signature,
  function() {
    zonal_sum_table(
      composition_stack,
      coarse_id_on_fine,
      "coarse_cell_id"
    ) |>
      mutate(
        coarse_cell_id = as.integer(round(coarse_cell_id)),
        NatureHa = NatureCells * fine_cell_area_ha,
        AgricultureHa = AgricultureCells * fine_cell_area_ha,
        ProductionForestHa =
          ProductionForestCells * fine_cell_area_ha,
        EligibleAreaHa = EligibleCells * fine_cell_area_ha,
        UncoveredEligibleHa =
          UncoveredEligibleCells * fine_cell_area_ha,
        IsCandidate = EligibleAreaHa >= minimum_eligible_area_ha
      ) |>
      select(
        coarse_cell_id,
        NatureHa,
        AgricultureHa,
        ProductionForestHa,
        EligibleAreaHa,
        UncoveredEligibleHa,
        IsCandidate
      )
  },
  fallback_paths = file.path(
    baseline_checkpoint_dir,
    "coarse_cell_composition.rds"
  )
)

candidate_table <- cell_table |>
  filter(IsCandidate) |>
  arrange(coarse_cell_id)
candidate_table <- candidate_table |> mutate(
  ProductionForestShare = ProductionForestHa / EligibleAreaHa,
  IsProductionForestCell = as.integer(ProductionForestShare >= production_forest_share_threshold)
)
forest_landuses <- nature_landuses[grepl("^Forest", nature_landuses)]
candidate_cells <- candidate_table$coarse_cell_id

if (!length(candidate_cells)) {
  stop("No 200 m candidate cells passed the eligibility rule.")
}

composition_audit <- tibble(
  current_nature_ha = sum(cell_table$NatureHa),
  agriculture_ha = sum(cell_table$AgricultureHa),
  production_forest_ha = sum(cell_table$ProductionForestHa),
  supported_eligible_ha = sum(cell_table$EligibleAreaHa),
  uncovered_eligible_ha = sum(cell_table$UncoveredEligibleHa),
  retained_candidate_eligible_ha =
    sum(candidate_table$EligibleAreaHa),
  candidate_cells = nrow(candidate_table)
)

knitr::kable(t(composition_audit), digits = 2)


# Source chunk: subarea-accounting
current_nature_by_subarea <- zonal_sum_table(
  fine_nature,
  subarea_id_on_fine,
  "subarea_id"
) |>
  transmute(
    subarea_id = as.integer(round(subarea_id)),
    CurrentNatureAreaHa = NatureCells * fine_cell_area_ha
  )

# One numeric zone encodes (coarse cell, subarea), avoiding a multi-million-row
# extraction table. With four subareas, multiplier 5 is collision-free.
pair_multiplier <- expected_number_of_subareas + 1L
pair_zone <- terra::ifel(
  is.na(coarse_id_on_fine) | is.na(subarea_id_on_fine),
  NA,
  coarse_id_on_fine * pair_multiplier + subarea_id_on_fine
)
names(pair_zone) <- "cell_subarea_pair"

candidate_mask_on_fine <- terra::ifel(
  coarse_id_on_fine %in% candidate_cells,
  1,
  NA
)
fine_candidate_eligible <- terra::mask(
  fine_eligible,
  candidate_mask_on_fine
)
names(fine_candidate_eligible) <- "EligibleCells"

fine_candidate_eligible_mask <- terra::mask(
  fine_eligible_mask,
  candidate_mask_on_fine
)
names(fine_candidate_eligible_mask) <- "ModelEligibleCells"

restorable_area_by_subarea_cell <- zonal_sum_table(
  fine_candidate_eligible,
  pair_zone,
  "pair_id"
) |>
  mutate(
    coarse_cell_id = as.integer(floor(pair_id / pair_multiplier)),
    subarea_id = as.integer(
      round(pair_id - coarse_cell_id * pair_multiplier)
    ),
    RestorableAreaHaBySubarea =
      EligibleCells * fine_cell_area_ha
  ) |>
  filter(
    coarse_cell_id %in% candidate_cells,
    subarea_id %in% subarea_dictionary$subarea_id,
    RestorableAreaHaBySubarea > 0
  ) |>
  select(
    subarea_id,
    coarse_cell_id,
    RestorableAreaHaBySubarea
  ) |>
  arrange(subarea_id, coarse_cell_id)

cell_area_crosscheck <- restorable_area_by_subarea_cell |>
  group_by(coarse_cell_id) |>
  summarise(
    PairAreaHa = sum(RestorableAreaHaBySubarea),
    .groups = "drop"
  ) |>
  right_join(
    candidate_table |>
      select(coarse_cell_id, EligibleAreaHa),
    by = "coarse_cell_id"
  )

if (
  anyNA(cell_area_crosscheck$PairAreaHa) ||
  any(
    abs(
      cell_area_crosscheck$PairAreaHa -
        cell_area_crosscheck$EligibleAreaHa
    ) > fine_cell_area_ha + 1e-8
  )
) {
  stop(
    "Candidate area does not reconcile with its subarea split.",
    call. = FALSE
  )
}

candidate_capacity_by_subarea <-
  restorable_area_by_subarea_cell |>
  group_by(subarea_id) |>
  summarise(
    CandidateCapacityHa = sum(RestorableAreaHaBySubarea),
    .groups = "drop"
  )

subarea_budget <- subarea_dictionary |>
  left_join(current_nature_by_subarea, by = "subarea_id") |>
  left_join(candidate_capacity_by_subarea, by = "subarea_id") |>
  mutate(
    CurrentNatureAreaHa = coalesce(CurrentNatureAreaHa, 0),
    CandidateCapacityHa = coalesce(CandidateCapacityHa, 0),
    FinalNatureTargetHa = target_nature_proportion * PolygonAreaHa,
    RestorationRequirementHa = pmax(
      0,
      FinalNatureTargetHa - CurrentNatureAreaHa
    ),
    RestorationToleranceHa = restoration_tolerance_ha,
    MaximumFinalNatureHa =
      CurrentNatureAreaHa + CandidateCapacityHa,
    CapacityMarginHa =
      CandidateCapacityHa - RestorationRequirementHa
  )

# Subarea 30% columns are reporting references, never constraints.
global_required_ha <- max(0, sum(subarea_budget$FinalNatureTargetHa) -
                               sum(subarea_budget$CurrentNatureAreaHa))
if (sum(candidate_table$EligibleAreaHa) + 1e-8 < global_required_ha) {
  stop("Insufficient candidate capacity for the whole-polygon target.")
}
if (sum(subarea_budget$CurrentNatureAreaHa) > sum(subarea_budget$FinalNatureTargetHa)) {
  stop("Current nature exceeds 30%; review the configured area policy.")
}
readr::write_csv(
  subarea_budget,
  file.path(problem_dir, "subarea_30pct_budget.csv")
)
readr::write_csv(
  candidate_table,
  file.path(problem_dir, "candidate_cell_dictionary.csv")
)
readr::write_csv(
  restorable_area_by_subarea_cell,
  file.path(problem_dir, "candidate_area_by_subarea.csv")
)

knitr::kable(
  subarea_budget |>
    select(
      subarea_id,
      subarea_label,
      PolygonAreaHa,
      FinalNatureTargetHa,
      CurrentNatureAreaHa,
      RestorationRequirementHa,
      CandidateCapacityHa,
      CapacityMarginHa
    ),
  digits = 2
)


# Source chunk: candidate-adjacency
adjacency_raw <- terra::adjacent(
  coarse_cell_id,
  cells = candidate_cells,
  directions = 4,
  pairs = TRUE,
  include = FALSE
) |>
  as.matrix()

candidate_edges <- tibble(
  i = pmin(adjacency_raw[, 1], adjacency_raw[, 2]),
  j = pmax(adjacency_raw[, 1], adjacency_raw[, 2])
) |>
  filter(
    i != j,
    i %in% candidate_cells,
    j %in% candidate_cells
  ) |>
  distinct() |>
  arrange(i, j)

if (!nrow(candidate_edges)) {
  stop("The 200 m candidate graph contains no rook edges.")
}

readr::write_csv(
  candidate_edges,
  file.path(problem_dir, "candidate_rook_edges.csv")
)


# Source chunk: existing-nature-contact
rook_kernel <- matrix(
  c(
    NA, 1, NA,
     1, 0,  1,
    NA, 1, NA
  ),
  nrow = 3,
  byrow = TRUE
)

existing_contact_signature <- digest::digest(
  list(
    files = file_fingerprint(c(
      fine_landuse_path,
      coarse_landuse_path,
      shapefile_components(subareas_path),
      biodiversity_support_path
    )),
    candidate_cells = candidate_cells,
    candidate_eligible_area_ha = candidate_table$EligibleAreaHa,
    landuse_codes = landuse_codes[nature_landuses],
    directions = 4L,
    fine_resolution = terra::res(fine_landuse)
  )
)

empty_existing_contact <- tidyr::crossing(
  landuse = nature_landuses,
  coarse_cell_id = candidate_cells
) |>
  mutate(ExistingNatureContact = 0)

if (
  spatial_contiguity_bonus > 0 &&
  existing_to_new_adjacency_weight > 0
) {
  existing_nature_contact <- load_or_build_checkpoint(
    file.path(problem_dir, "existing_nature_contact_by_action.rds"),
    existing_contact_signature,
    function() {
      contact_values <- purrr::map_dfr(
        nature_landuses,
        function(landuse_name) {
          existing_habitat <- indicator_for_codes(
            fine_landuse,
            unname(landuse_codes[[landuse_name]])
          )

          neighboring_existing_cells <- terra::focal(
            existing_habitat,
            w = rook_kernel,
            fun = "sum",
            na.rm = TRUE,
            fillvalue = 0
          )
          neighboring_existing_cells <- terra::mask(
            neighboring_existing_cells,
            fine_candidate_eligible_mask
          )
          names(neighboring_existing_cells) <-
            "ExistingNatureContact"

          zonal_sum_table(
            neighboring_existing_cells,
            coarse_id_on_fine,
            "coarse_cell_id"
          ) |>
            transmute(
              landuse = landuse_name,
              coarse_cell_id =
                as.integer(round(coarse_cell_id)),
              ExistingNatureContact =
                as.numeric(ExistingNatureContact)
            ) |>
            filter(coarse_cell_id %in% candidate_cells)
        }
      )

      empty_existing_contact |>
        select(landuse, coarse_cell_id) |>
        left_join(
          contact_values,
          by = c("landuse", "coarse_cell_id")
        ) |>
        mutate(
          ExistingNatureContact = coalesce(
            ExistingNatureContact,
            0
          )
        ) |>
        arrange(landuse, coarse_cell_id)
    }
  )
} else {
  existing_nature_contact <- empty_existing_contact
}

stopifnot(
  nrow(existing_nature_contact) ==
    length(candidate_cells) * length(nature_landuses),
  !anyNA(existing_nature_contact),
  all(is.finite(existing_nature_contact$ExistingNatureContact)),
  all(existing_nature_contact$ExistingNatureContact >= 0)
)

if (
  spatial_contiguity_bonus > 0 &&
  existing_to_new_adjacency_weight > 0 &&
  sum(existing_nature_contact$ExistingNatureContact) == 0
) {
  warning(
    "Existing-to-new adjacency is active, but no matching 10 m contacts ",
    "were found. Inspect the land-use categories and candidate mask."
  )
}

fine_edge_length_m <- min(terra::res(fine_landuse))
existing_nature_contact_audit <- existing_nature_contact |>
  mutate(
    ExistingNatureContactM =
      ExistingNatureContact * fine_edge_length_m
  )

readr::write_csv(
  existing_nature_contact_audit,
  file.path(problem_dir, "existing_nature_contact_by_action.csv")
)

knitr::kable(
  existing_nature_contact_audit |>
    group_by(landuse) |>
    summarise(
      candidate_actions_with_contact =
        sum(ExistingNatureContact > 0),
      total_available_contact_km =
        sum(ExistingNatureContactM) / 1000,
      .groups = "drop"
    ),
  digits = 3
)


# Source chunk: biodiversity-benefits
read_metric_on_fine_grid <- function(path) {
  x <- terra::rast(path)

  if (terra::nlyr(x) != 1L) {
    stop("Expected one metric layer in ", path, call. = FALSE)
  }

  align_to_template(x, fine_landuse, method = "bilinear")
}

model_eligible_cells <- global_number(
  fine_candidate_eligible_mask,
  "notNA"
)

biodiversity_signature <- digest::digest(
  list(
    files = file_fingerprint(c(
      metric_file_index$path,
      metric_current_index$path,
      fine_landuse_path,
      coarse_landuse_path,
      shapefile_components(subareas_path)
    )),
    metric_weights = metric_weights,
    biodiversity_scale = biodiversity_scale,
    candidate_cells = candidate_cells,
    candidate_eligible_area_ha = candidate_table$EligibleAreaHa,
    support_metric = biodiversity_support_metric,
    support_scenario = biodiversity_support_scenario
  )
)

biodiversity_checkpoint <- load_or_build_checkpoint(
  file.path(problem_dir, "coarse_biodiversity_benefits.rds"),
  biodiversity_signature,
  function() {
    all_metric_files <- bind_rows(
      metric_current_index,
      metric_file_index
    )

    metric_normalization <- map_dfr(
      names(metric_weights),
      function(metric_name) {
        paths <- all_metric_files$path[
          all_metric_files$metric == metric_name
        ]

        maxima <- map_dbl(
          paths,
          function(path) {
            x <- terra::rast(path)
            global_number(x, "max")
          }
        )
        normalizer <- max(maxima)

        if (!is.finite(normalizer) || normalizer <= 0) {
          stop("Invalid normalizer for ", metric_name, call. = FALSE)
        }

        tibble(
          metric = metric_name,
          normalizer = normalizer,
          weight = unname(metric_weights[[metric_name]])
        )
      }
    )

    biodiversity_table <- map_dfr(
      seq_len(nrow(scenario_index)),
      function(index) {
        landuse_name <- scenario_index$landuse[[index]]
        scenario_name <- scenario_index$scenario[[index]]
        combined_score <- NULL

        for (metric_name in names(metric_weights)) {
          metric_path <- metric_file_index |>
            filter(
              .data$metric == metric_name,
              .data$scenario == scenario_name
            ) |>
            pull(path)

          stopifnot(length(metric_path) == 1L)
          metric_raster <- read_metric_on_fine_grid(metric_path)
          metric_raster <- terra::mask(
            metric_raster,
            fine_candidate_eligible_mask
          )

          covered_cells <- global_number(metric_raster, "notNA")
          if (covered_cells != model_eligible_cells) {
            stop(
              "Incomplete biodiversity coverage for ",
              scenario_name, " / ", metric_name, ": ",
              model_eligible_cells - covered_cells,
              " eligible 10 m cells are missing.",
              call. = FALSE
            )
          }

          normalizer <- metric_normalization$normalizer[
            metric_normalization$metric == metric_name
          ]
          contribution <-
            metric_raster / normalizer *
            unname(metric_weights[[metric_name]])

          if (is.null(combined_score)) {
            combined_score <- contribution
          } else {
            combined_score <- combined_score + contribution
          }
        }

        integrated_score <- terra::mask(
          combined_score,
          fine_candidate_eligible_mask
        ) * biodiversity_scale * fine_cell_area_ha
        names(integrated_score) <- "BioDiversity"

        zonal_sum_table(
          integrated_score,
          coarse_id_on_fine,
          "coarse_cell_id"
        ) |>
          transmute(
            coarse_cell_id = as.integer(round(coarse_cell_id)),
            landuse = landuse_name,
            BioDiversity = as.numeric(BioDiversity)
          ) |>
          filter(coarse_cell_id %in% candidate_cells)
      }
    )

    list(
      normalization = metric_normalization,
      benefits = biodiversity_table
    )
  },
  fallback_paths = file.path(
    baseline_checkpoint_dir,
    "coarse_biodiversity_benefits.rds"
  )
)

metric_normalization <- biodiversity_checkpoint$normalization
biodiversity_table <- biodiversity_checkpoint$benefits

expected_biodiversity_rows <-
  length(candidate_cells) * length(nature_landuses)

if (
  nrow(biodiversity_table) != expected_biodiversity_rows ||
  any(!is.finite(biodiversity_table$BioDiversity))
) {
  stop(
    "The biodiversity table is incomplete or non-finite.",
    call. = FALSE
  )
}

readr::write_csv(
  metric_normalization,
  file.path(problem_dir, "metric_normalization.csv")
)

knitr::kable(metric_normalization, digits = 6)


# Source chunk: model-tables
allowed_actions <- tidyr::crossing(
  landuse = nature_landuses,
  candidate_table |>
    select(coarse_cell_id, IsProductionForestCell)
) |>
  filter(
    IsProductionForestCell == 0L |
      landuse %in% forest_landuses
  ) |>
  select(landuse, coarse_cell_id)

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
  !anyNA(objective_table),
  all(is.finite(objective_table$BioDiversity)),
  all(is.finite(objective_table$ExistingNatureContact)),
  all(objective_table$ExistingNatureContact >= 0)
)

actions_before_pruning <- nrow(objective_table)
dominance_pruning_applied <-
  isTRUE(prune_dominated_landuse_actions) &&
  spatial_contiguity_bonus == 0

if (dominance_pruning_applied) {
  objective_table <- objective_table |>
    group_by(coarse_cell_id) |>
    arrange(desc(BioDiversity), landuse, .by_group = TRUE) |>
    slice_head(n = 1L) |>
    ungroup() |>
    arrange(landuse, coarse_cell_id)
}

allowed_actions <- objective_table |>
  select(landuse, coarse_cell_id)

stopifnot(
  n_distinct(allowed_actions$coarse_cell_id) == nrow(candidate_table),
  !anyDuplicated(allowed_actions)
)

if (
  spatial_contiguity_bonus > 0 &&
  new_to_new_adjacency_weight > 0
) {
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
} else {
  same_habitat_edges <- tibble(
    landuse = character(),
    i = integer(),
    j = integer()
  )
}

biodiversity_normalizer <- objective_table |>
  group_by(coarse_cell_id) |>
  summarise(best = max(BioDiversity), .groups = "drop") |>
  summarise(total = sum(best)) |>
  pull(total)
biodiversity_normalizer <- max(1, biodiversity_normalizer)
contiguity_normalizer <- max(1, nrow(candidate_edges))
existing_nature_contact_normalizer <- objective_table |>
  group_by(coarse_cell_id) |>
  summarise(
    best = max(ExistingNatureContact),
    .groups = "drop"
  ) |>
  summarise(total = sum(best)) |>
  pull(total)
existing_nature_contact_normalizer <- max(
  1,
  existing_nature_contact_normalizer
)

problem_audit <- tibble(
  polygon_area_ha = sum(subarea_budget$PolygonAreaHa),
  thirty_percent_target_ha =
    sum(subarea_budget$FinalNatureTargetHa),
  current_nature_ha = sum(subarea_budget$CurrentNatureAreaHa),
  required_new_restoration_ha =
    sum(subarea_budget$RestorationRequirementHa),
  candidate_capacity_ha = sum(candidate_table$EligibleAreaHa),
  candidate_cells = nrow(candidate_table),
  candidate_rook_edges = nrow(candidate_edges),
  actions_before_pruning = actions_before_pruning,
  actions_after_pruning = nrow(allowed_actions),
  dominance_pruning_applied = dominance_pruning_applied,
  same_habitat_edge_binaries = nrow(same_habitat_edges),
  existing_nature_actions_with_contact =
    sum(objective_table$ExistingNatureContact > 0),
  existing_nature_contact_edges_available =
    sum(objective_table$ExistingNatureContact),
  existing_nature_contact_km_available =
    sum(objective_table$ExistingNatureContact) *
      fine_edge_length_m / 1000,
  existing_nature_contact_normalizer =
    existing_nature_contact_normalizer,
  spatial_contiguity_bonus = spatial_contiguity_bonus,
  new_to_new_adjacency_weight = new_to_new_adjacency_weight,
  existing_to_new_adjacency_weight =
    existing_to_new_adjacency_weight
)

readr::write_csv(
  problem_audit,
  file.path(problem_dir, "problem_size_and_policy_audit.csv")
)
readr::write_csv(
  tibble(
    pruning_requested = prune_dominated_landuse_actions,
    pruning_applied = dominance_pruning_applied,
    actions_before = actions_before_pruning,
    actions_after = nrow(allowed_actions),
    actions_removed = actions_before_pruning - nrow(allowed_actions)
  ),
  file.path(problem_dir, "dominated_action_audit.csv")
)

knitr::kable(t(problem_audit), digits = 3)
