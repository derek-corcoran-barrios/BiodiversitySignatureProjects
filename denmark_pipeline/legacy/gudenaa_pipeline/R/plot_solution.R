# Maps for the exact validated prefix supplied by the pipeline.
# Uses the existing NaturetypeStyle.sld and current/selected spatial layers.

# ---- Packages -----------------------------------------------------------
required_packages <- c(
  "terra", "dplyr", "purrr", "readr", "tibble", "ggplot2",
  "tidyterra", "xml2", "ragg"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop(
    "Install the missing package(s) first: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
  library(purrr)
  library(readr)
  library(tibble)
  library(ggplot2)
  library(tidyterra)
  library(xml2)
})

# Paths and settings are supplied by gp_plot().
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
# ---- Helpers ------------------------------------------------------------
stop_for_missing_files <- function(paths, description) {
  paths <- unique(paths)
  missing <- paths[is.na(paths) | !file.exists(paths)]
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

landuse_lookup <- function(x) {
  lev <- terra::levels(x)[[1]]
  if (is.null(lev) || !nrow(lev)) {
    stop("The land-use raster has no categorical lookup table.", call. = FALSE)
  }
  numeric_columns <- names(lev)[vapply(lev, is.numeric, logical(1))]
  label_columns <- names(lev)[vapply(
    lev,
    function(z) is.character(z) || is.factor(z),
    logical(1)
  )]
  if (!length(numeric_columns) || !length(label_columns)) {
    stop("Could not identify land-use value and label columns.", call. = FALSE)
  }
  tibble(
    value = as.integer(lev[[numeric_columns[[1]]]]),
    landuse = as.character(lev[[label_columns[[1]]]])
  )
}

codes_for_labels <- function(lookup, requested_labels) {
  missing_labels <- setdiff(requested_labels, lookup$landuse)
  if (length(missing_labels)) {
    stop(
      "Missing land-use classes: ",
      paste(missing_labels, collapse = ", "),
      call. = FALSE
    )
  }
  stats::setNames(
    lookup$value[match(requested_labels, lookup$landuse)],
    requested_labels
  )
}

read_sld_palette <- function(path) {
  stop_for_missing_files(path, "QGIS SLD palette")
  document <- xml2::read_xml(path)
  rules <- xml2::xml_find_all(
    document,
    ".//*[local-name()='Rule']"
  )
  palette <- purrr::map_dfr(rules, function(rule) {
    literal <- xml2::xml_text(xml2::xml_find_first(
      rule,
      ".//*[local-name()='Literal']"
    ))
    colour <- xml2::xml_text(xml2::xml_find_first(
      rule,
      ".//*[local-name()='SvgParameter' and @name='fill']"
    ))
    title <- xml2::xml_text(xml2::xml_find_first(
      rule,
      ".//*[local-name()='Title']"
    ))
    tibble(
      DN = suppressWarnings(as.integer(literal)),
      title = title,
      colour = colour
    )
  }) |>
    filter(!is.na(DN), nzchar(colour)) |>
    distinct(DN, .keep_all = TRUE) |>
    arrange(DN)

  if (!identical(palette$DN, 1:8)) {
    stop("The SLD does not contain exactly DN values 1 through 8.")
  }
  if (
    palette$colour[1] != palette$colour[3] ||
    palette$colour[2] != palette$colour[4] ||
    palette$colour[5] != palette$colour[7] ||
    palette$colour[6] != palette$colour[8]
  ) {
    stop("The SLD rich/poor colour pairs are not identical as expected.")
  }
  palette
}

resolve_solution_prefix <- function() final_solution_prefix

save_plot_pair <- function(plot, stem, width = 10, height = 10) {
  png_path <- file.path(figure_dir, paste0(stem, ".png"))
  pdf_path <- file.path(figure_dir, paste0(stem, ".pdf"))
  ggplot2::ggsave(
    png_path,
    plot = plot,
    device = ragg::agg_png,
    width = width,
    height = height,
    units = "in",
    dpi = png_dpi,
    bg = "white"
  )
  ggplot2::ggsave(
    pdf_path,
    plot = plot,
    device = if (capabilities("cairo")) grDevices::cairo_pdf else "pdf",
    width = width,
    height = height,
    units = "in",
    bg = "white"
  )
  invisible(c(png = png_path, pdf = pdf_path))
}

map_theme <- function() {
  theme_void(base_size = 12) +
    theme(
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      legend.text = element_text(size = 9),
      legend.key.height = grid::unit(4, "mm"),
      legend.key.width = grid::unit(7, "mm"),
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10, colour = "grey30"),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.margin = margin(6, 6, 6, 6)
    )
}

# ---- Resolve and validate the incumbent --------------------------------
solution_prefix <- resolve_solution_prefix()
solution_tag <- basename(solution_prefix)
certificate_cache_path <- file.path(
  figure_dir,
  paste0(solution_tag, "_certificate_geometries_full_patch_extent.gpkg")
)
selected_solution_path <- paste0(solution_prefix, "_selected_cells.tsv")
network_cells_path <- paste0(solution_prefix, "_network_cells.tsv")
network_patches_path <- paste0(solution_prefix, "_network_patches.tsv")
network_summary_path <- paste0(solution_prefix, "_network_summary.tsv")
solution_summary_path <- paste0(solution_prefix, "_solution_summary.tsv")
stage_summary_path <- paste0(solution_prefix, "_stage_summary.tsv")

stop_for_missing_files(
  c(
    fine_landuse_path,
    coarse_landuse_path,
    subareas_path,
    national_patch_path,
    biodiversity_support_path,
    sld_path,
    selected_solution_path,
    network_cells_path,
    network_patches_path,
    network_summary_path,
    solution_summary_path,
    stage_summary_path
  ),
  "map input file(s)"
)

selected_solution <- readr::read_tsv(
  selected_solution_path,
  show_col_types = FALSE
) |>
  transmute(
    coarse_cell_id = as.integer(coarse_cell_id),
    landuse = as.character(landuse)
  )
network_cells <- readr::read_tsv(
  network_cells_path,
  show_col_types = FALSE
) |>
  mutate(coarse_cell_id = as.integer(coarse_cell_id))
network_patches <- readr::read_tsv(
  network_patches_path,
  show_col_types = FALSE
) |>
  mutate(patch_id = as.integer(patch_id))
network_summary <- readr::read_tsv(
  network_summary_path,
  show_col_types = FALSE
)
solution_summary <- readr::read_tsv(
  solution_summary_path,
  show_col_types = FALSE
)
stage_summary <- readr::read_tsv(
  stage_summary_path,
  show_col_types = FALSE
)

if (anyDuplicated(selected_solution$coarse_cell_id)) {
  stop("The pipeline solution assigns more than one land use to a cell.")
}

network_order <- c(
  "Global", "Subarea1", "Subarea2", "Subarea3", "Subarea4"
)
if (!setequal(network_summary$network_name, network_order)) {
  stop("The solution does not contain the expected five certificates.")
}
if (
  any(network_summary$shortfall_ha > 1e-6) ||
  any(network_summary$certificate_built != 1)
) {
  stop("The selected solution is not a zero-shortfall network incumbent.")
}

# ---- Spatial reconstruction --------------------------------------------
subareas <- terra::vect(subareas_path) |>
  terra::makeValid()
subareas$subarea_id <- seq_len(nrow(subareas))
planning_domain <- terra::aggregate(subareas)

fine_landuse_full <- terra::rast(fine_landuse_path)
coarse_landuse_full <- terra::rast(coarse_landuse_path)
if (!same_crs(subareas, fine_landuse_full)) {
  subareas <- terra::project(subareas, terra::crs(fine_landuse_full))
  subareas$subarea_id <- seq_len(nrow(subareas))
  planning_domain <- terra::aggregate(subareas)
}

fine_landuse <- fine_landuse_full |>
  terra::crop(planning_domain, snap = "out") |>
  terra::mask(planning_domain, touches = FALSE)
coarse_landuse <- coarse_landuse_full |>
  terra::crop(planning_domain, snap = "out") |>
  terra::mask(planning_domain, touches = TRUE)

lookup <- landuse_lookup(fine_landuse)
nature_landuses <- c(
  "ForestDryPoor", "ForestDryRich",
  "ForestWetPoor", "ForestWetRich",
  "OpenDryPoor", "OpenDryRich",
  "OpenWetPoor", "OpenWetRich"
)
eligible_landuses <- c("Agriculture", "ProductionForest")
all_landuses <- c(nature_landuses, eligible_landuses, "Urban", "Other")
landuse_codes <- codes_for_labels(lookup, all_landuses)

coarse_cell_id <- terra::init(coarse_landuse, "cell")
coarse_cell_id <- terra::ifel(is.na(coarse_landuse), NA, coarse_cell_id)
names(coarse_cell_id) <- "coarse_cell_id"
coarse_id_on_fine <- align_to_template(
  coarse_cell_id,
  fine_landuse,
  method = "near"
)

subarea_id_on_fine <- terra::rasterize(
  subareas[, "subarea_id"],
  fine_landuse,
  field = "subarea_id",
  background = NA,
  touches = FALSE
) |>
  terra::mask(fine_landuse)

biodiversity_support <- terra::rast(biodiversity_support_path)
biodiversity_support <- align_to_template(
  biodiversity_support,
  fine_landuse,
  method = "bilinear"
)

fine_current_nature <- terra::ifel(
  fine_landuse %in% unname(landuse_codes[nature_landuses]),
  1,
  0
)
fine_eligible <- terra::ifel(
  fine_landuse %in% unname(landuse_codes[eligible_landuses]) &
    !is.na(biodiversity_support),
  1,
  0
)

selected_ids <- selected_solution$coarse_cell_id
fine_restoration <- terra::ifel(
  fine_eligible == 1 & coarse_id_on_fine %in% selected_ids,
  1,
  0
) |>
  terra::mask(fine_landuse)
names(fine_restoration) <- "New restoration"

# Current nature and restoration are mutually exclusive by construction.
nature_status <- terra::ifel(
  fine_restoration == 1,
  2,
  terra::ifel(fine_current_nature == 1, 1, NA)
)
names(nature_status) <- "Nature status"
levels(nature_status) <- data.frame(
  value = c(1L, 2L),
  status = factor(
    c("Current nature", "Newly selected nature"),
    levels = c("Current nature", "Newly selected nature")
  )
)

national_patches <- terra::vect(
  national_patch_path,
  layer = national_patch_layer
) |>
  terra::makeValid()
if (!same_crs(national_patches, fine_landuse)) {
  national_patches <- terra::project(
    national_patches,
    terra::crs(fine_landuse)
  )
}
if (!"patch_id" %in% names(national_patches)) {
  stop("The national patch layer has no patch_id field.")
}
national_patches$patch_id <- as.integer(national_patches$patch_id)

# ---- Exact certificate geometries --------------------------------------
certificate_labels <- c(
  Global = "Global certificate: minimum 5,000 ha",
  Subarea1 = "Subarea 1: minimum 1,000 ha",
  Subarea2 = "Subarea 2: minimum 1,000 ha",
  Subarea3 = "Subarea 3: minimum 1,000 ha",
  Subarea4 = "Subarea 4: minimum 1,000 ha"
)
certificate_palette <- c(
  "Global certificate: minimum 5,000 ha" = "#1F78B4",
  "Subarea 1: minimum 1,000 ha" = "#33A02C",
  "Subarea 2: minimum 1,000 ha" = "#FF7F00",
  "Subarea 3: minimum 1,000 ha" = "#6A3D9A",
  "Subarea 4: minimum 1,000 ha" = "#E31A1C"
)

build_one_certificate <- function(network_name) {
  member_cell_ids <- network_cells |>
    filter(.data$network_name == .env$network_name) |>
    pull(coarse_cell_id) |>
    unique()
  member_patch_ids <- network_patches |>
    filter(.data$network_name == .env$network_name) |>
    pull(patch_id) |>
    unique()

  new_nature_mask <- terra::ifel(
    fine_restoration == 1 & coarse_id_on_fine %in% member_cell_ids,
    1,
    NA
  )
  if (network_name != "Global") {
    subarea_id <- as.integer(sub("Subarea", "", network_name))
    new_nature_mask <- terra::ifel(
      !is.na(new_nature_mask) & subarea_id_on_fine == subarea_id,
      1,
      NA
    )
  }

  new_nature_vector <- terra::as.polygons(
    new_nature_mask,
    aggregate = TRUE,
    values = FALSE,
    na.rm = TRUE
  ) |>
    terra::makeValid()
  new_nature_vector$certificate <- network_name
  new_nature_vector$source <- "New restoration"

  existing_patch_vector <- national_patches[
    national_patches$patch_id %in% member_patch_ids,
  ]
  if (!nrow(existing_patch_vector)) {
    stop("No existing patch geometry found for ", network_name, ".")
  }
  existing_patch_vector$certificate <- network_name
  existing_patch_vector$source <- "Current nature patch"

  pieces <- rbind(
    new_nature_vector[, c("certificate", "source")],
    existing_patch_vector[, c("certificate", "source")]
  )
  dissolved <- terra::aggregate(
    pieces,
    by = "certificate",
    dissolve = TRUE
  )
  dissolved$certificate <- network_name
  dissolved$certificate_label <- unname(certificate_labels[network_name])
  dissolved
}

if (
  reuse_certificate_geometry_cache &&
  file.exists(certificate_cache_path)
) {
  certificate_geometries <- terra::vect(certificate_cache_path) |>
    terra::makeValid()
} else {
  certificate_geometries <- network_order |>
    purrr::map(build_one_certificate) |>
    purrr::reduce(rbind)
  terra::writeVector(
    certificate_geometries,
    certificate_cache_path,
    overwrite = TRUE
  )
}

certificate_geometries$certificate <- as.character(
  certificate_geometries$certificate
)
certificate_geometries$certificate_label <- factor(
  unname(certificate_labels[certificate_geometries$certificate]),
  levels = unname(certificate_labels[network_order])
)
global_certificate <- certificate_geometries[
  certificate_geometries$certificate == "Global",
]
subarea_certificates <- certificate_geometries[
  certificate_geometries$certificate != "Global",
]

# ---- Map 1: filled 5,000/1,000 ha certificates -------------------------
map_1 <- ggplot() +
  geom_spatvector(
    data = global_certificate,
    aes(fill = certificate_label, colour = certificate_label),
    alpha = 0.28,
    linewidth = 0.8
  ) +
  geom_spatvector(
    data = subarea_certificates,
    aes(fill = certificate_label, colour = certificate_label),
    alpha = 0.48,
    linewidth = 0.65
  ) +
  geom_spatvector(
    data = subareas,
    fill = NA,
    colour = "#111111",
    linewidth = 0.45
  ) +
  scale_fill_manual(
    name = "Certified connected areas",
    values = certificate_palette,
    drop = FALSE
  ) +
  scale_colour_manual(
    name = "Certified connected areas",
    values = certificate_palette,
    drop = FALSE
  ) +
  coord_sf(expand = FALSE) +
  guides(
    fill = guide_legend(nrow = 3, byrow = TRUE),
    colour = "none"
  ) +
  labs(
    title = "Connected nature certificates in the incumbent solution",
    subtitle = paste(
      "Final incumbent",
      "quality stage. Linked existing patches retain their full geometry."
    )
  ) +
  map_theme()

save_plot_pair(
  map_1,
  "pipeline_final_01_certificate_areas_full_patch_extent",
  width = 11,
  height = 10
)

# ---- Map 2: certificate contours + current/new nature ------------------
nature_status_palette <- c(
  "Current nature" = "#01846a",
  "Newly selected nature" = "#d58c43"
)

global_certificate_plot <- terra::aggregate(
  global_certificate,
  by = "certificate_label",
  dissolve = TRUE
)

subarea_certificates_plot <- terra::aggregate(
  subarea_certificates,
  by = "certificate_label",
  dissolve = TRUE
)


map_2 <- ggplot() +
  geom_spatraster(
    data = nature_status,
    maxcell = raster_maxcell
  ) +
  geom_spatvector(
    data = global_certificate,
    aes(colour = certificate_label),
    fill = NA,
    linewidth = 1.05
  ) +
  geom_spatvector(
    data = subarea_certificates,
    aes(colour = certificate_label),
    fill = NA,
    linewidth = 0.75
  ) +
  geom_spatvector(
    data = subareas,
    fill = NA,
    colour = "#111111",
    linewidth = 0.4
  ) +
  scale_fill_manual(
    name = "Nature",
    values = nature_status_palette,
    drop = FALSE,
    na.value = "#00000000"
  ) +
  scale_colour_manual(
    name = "Certificate contour",
    values = certificate_palette,
    drop = FALSE
  ) +
  coord_sf(expand = FALSE) +
  guides(
    fill = guide_legend(order = 1, nrow = 1),
    colour = guide_legend(order = 2, nrow = 3, byrow = TRUE)
  ) +
  labs(
    title = "Current and newly selected nature",
    subtitle = paste(
      "The five contours are from the pipeline incumbent; full linked",
      "existing patches remain visible outside the planning polygon."
    )
  ) +
  map_theme()

save_plot_pair(
  map_2,
  "pipeline_final_02_current_and_new_nature_with_certificate_contours",
  width = 11,
  height = 10
)

# ---- Map 3: resulting land use with the QGIS SLD palette ---------------
sld_palette <- read_sld_palette(sld_path)
readr::write_csv(
  sld_palette,
  file.path(figure_dir, "qgis_sld_palette_extracted.csv")
)

selected_destination_codes <- unname(
  landuse_codes[selected_solution$landuse]
)
if (anyNA(selected_destination_codes)) {
  stop("At least one selected destination is absent from FinalLanduse.tif.")
}

destination_coarse <- coarse_cell_id
destination_values <- rep(NA_integer_, terra::ncell(destination_coarse))
destination_values[selected_ids] <- selected_destination_codes
terra::values(destination_coarse) <- destination_values

destination_fine <- align_to_template(
  destination_coarse,
  fine_landuse,
  method = "near"
)
destination_fine <- terra::mask(
  destination_fine,
  terra::ifel(fine_eligible == 1, 1, NA)
)

resulting_landuse <- terra::cover(destination_fine, fine_landuse)

manager_groups <- tibble(
  landuse = c(
    "ForestDryPoor", "ForestDryRich",
    "ForestWetPoor", "ForestWetRich",
    "OpenDryPoor", "OpenDryRich",
    "OpenWetPoor", "OpenWetRich",
    "Agriculture", "Urban", "ProductionForest", "Other"
  ),
  manager_code = c(1L, 1L, 2L, 2L, 3L, 3L, 4L, 4L, 5L, 6L, 7L, 8L),
  manager_label = c(
    "Dry forest nature", "Dry forest nature",
    "Wet forest nature", "Wet forest nature",
    "Dry open nature", "Dry open nature",
    "Wet open nature", "Wet open nature",
    "Agriculture", "Urban", "Production forest", "Other"
  )
) |>
  left_join(lookup, by = "landuse")
if (anyNA(manager_groups$value)) {
  stop("Could not map all manager-facing classes to raster values.")
}

manager_landuse <- terra::subst(
  resulting_landuse,
  from = manager_groups$value,
  to = manager_groups$manager_code, raw = TRUE
)
names(manager_landuse) <- "Resulting land use"
manager_level_table <- manager_groups |>
  distinct(manager_code, manager_label) |>
  arrange(manager_code)
levels(manager_landuse) <- data.frame(
  value = manager_level_table$manager_code,
  landuse = factor(
    manager_level_table$manager_label,
    levels = manager_level_table$manager_label
  )
)

manager_palette <- c(
  "Dry forest nature" = sld_palette$colour[sld_palette$DN == 1],
  "Wet forest nature" = sld_palette$colour[sld_palette$DN == 2],
  "Dry open nature" = sld_palette$colour[sld_palette$DN == 5],
  "Wet open nature" = sld_palette$colour[sld_palette$DN == 6],
  "Agriculture" = "#D9B3B8",
  "Urban" = "#8C8C8C",
  "Production forest" = "#80739B",
  "Other" = "#F2F2F2"
)

map_3 <- ggplot() +
  geom_spatraster(
    data = manager_landuse,
    maxcell = raster_maxcell
  ) +
  geom_spatvector(
    data = subareas,
    fill = NA,
    colour = "#082A38",
    linewidth = 0.42
  ) +
  scale_fill_manual(
    name = "Resulting land use",
    values = manager_palette,
    breaks = names(manager_palette),
    drop = FALSE,
    na.value = "#00000000"
  ) +
  coord_sf(expand = FALSE) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  labs(
    title = "Resulting land use",
    subtitle = paste(
      "Nature colours come directly from",
      "NaturetypeStyle.sld; rich and poor variants are combined."
    )
  ) +
  map_theme()

save_plot_pair(
  map_3,
  "pipeline_final_03_resulting_landuse_qgis_sld_palette",
  width = 9,
  height = 10
)

# ---- Reproducibility audit ---------------------------------------------
readr::write_csv(
  tibble(
    item = c(
      "solution_prefix",
      "selected_cells",
      "global_certified_area_ha",
      "global_required_area_ha",
      "subarea_certificates",
      "certificate_geometry_includes_full_external_patch_extent",
      "qgis_sld",
      "pipeline_completion_status",
      "pipeline_stage_rows"
    ),
    value = c(
      solution_prefix,
      as.character(nrow(selected_solution)),
      as.character(network_summary$certified_area_ha[
        match("Global", network_summary$network_name)
      ]),
      as.character(network_summary$minimum_ha[
        match("Global", network_summary$network_name)
      ]),
      "4",
      "TRUE",
      sld_path,
      "COMPLETE",
      as.character(nrow(stage_summary))
    )
  ),
  file.path(figure_dir, "incumbent_map_source_audit.csv")
)

message(
  "Finished. Maps and certificate geometry are in: ",
  normalizePath(figure_dir, winslash = "/", mustWork = FALSE)
)

# Additional maps for the supplied 09_plot_gudenaa_incumbent.R script.
# Run the original plotting script first in the SAME R session, then:
# source("09b_current_status_and_subarea_closeup.R")
# No AMPL run or certificate reconstruction is needed.
# This add-on does not modify the original script or input files.
# Uses the original script's subarea_id (shapefile row order); verify that
# row 3 is the intended project Subarea 3. It inherits the original spatial
# reconstruction, including its coarse-cell ID convention.
# These maps are presentation products, not independent feasibility audits.

# ---- Settings ----------------------------------------------------------
zoom_subarea_id <- 3L
detail_width_m <- 2400
detail_height_m <- 2400
# Optional fixed extent in the land-use raster CRS: c(xmin, xmax, ymin, ymax).
# NULL chooses a genuine rook-contact cell in the strongest 1 km contact block.
detail_extent_override <- NULL
contact_search_block_m <- 1000

needed <- c("fine_landuse", "fine_current_nature", "fine_restoration",
            "subareas", "manager_groups", "manager_level_table",
            "manager_palette", "nature_status_palette", "figure_dir",
            "map_theme", "raster_maxcell", "png_dpi", "solution_prefix")
missing_objects <- needed[!vapply(needed, exists, logical(1), inherits = TRUE)]
if (length(missing_objects)) {
  stop("Run the original plotting script first in this R session. Missing: ",
       paste(missing_objects, collapse = ", "))
}
if (terra::is.lonlat(fine_landuse)) {
  stop("The detail-window settings require a projected CRS in metres.")
}
stopifnot(terra::compareGeom(fine_landuse, fine_current_nature),
          terra::compareGeom(fine_landuse, fine_restoration))
if (!zoom_subarea_id %in% subareas$subarea_id) stop("Subarea ID not found.")
additional_figure_dir <- file.path(figure_dir, "current_status_and_closeup")
dir.create(additional_figure_dir, recursive = TRUE, showWarnings = FALSE)

save_additional_map <- function(p, stem, width = 9, height = 10) {
  ggplot2::ggsave(file.path(additional_figure_dir, paste0(stem, ".png")),
                  plot = p, device = ragg::agg_png, width = width, height = height,
                  units = "in", dpi = png_dpi, bg = "white")
  ggplot2::ggsave(file.path(additional_figure_dir, paste0(stem, ".pdf")),
                  plot = p,
                  device = if (capabilities("cairo")) grDevices::cairo_pdf else "pdf",
                  width = width, height = height, units = "in", bg = "white")
}

# ---- Current status: same habitat grouping and QGIS colours as Map 3 ----
current_manager_landuse <- terra::subst(fine_landuse,
                                        from = manager_groups$value, to = manager_groups$manager_code, raw = TRUE)
current_levels <- data.frame(value = manager_level_table$manager_code,
                             landuse = factor(manager_level_table$manager_label,
                                              levels = manager_level_table$manager_label))
levels(current_manager_landuse) <- current_levels
names(current_manager_landuse) <- "Current land use"

# Nature-only map. Agriculture/urban/production forest/other are blank.
current_manager_nature <- terra::mask(current_manager_landuse,
                                      terra::ifel(fine_current_nature == 1, 1, NA))
levels(current_manager_nature) <- current_levels[current_levels$value <= 4, ]
names(current_manager_nature) <- "Current nature"

current_map <- function(x, nature_only) {
  palette <- if (nature_only) manager_palette[1:4] else manager_palette
  ggplot2::ggplot() +
    tidyterra::geom_spatraster(data = x, maxcell = raster_maxcell) +
    tidyterra::geom_spatvector(data = subareas, fill = NA,
                               colour = "#082A38", linewidth = 0.5) +
    ggplot2::scale_fill_manual(name = if (nature_only) "Current nature" else
      "Current land use", values = palette, limits = names(palette),
      drop = FALSE, na.value = "transparent") +
    ggplot2::coord_sf(expand = FALSE) +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 2, byrow = TRUE)) +
    ggplot2::labs(title = if (nature_only) "Current nature" else "Current land use",
                  subtitle = "Existing conditions and the four planning subareas",
                  caption = if (nature_only) "Blank areas: other land uses or missing data."
                  else "Same habitat grouping and palette as the resulting-land-use map.") +
    map_theme()
}
map_4_current_nature <- current_map(current_manager_nature, TRUE)
map_4b_current_landuse <- current_map(current_manager_landuse, FALSE)
save_additional_map(map_4_current_nature, "04_current_nature_qgis_palette")
save_additional_map(map_4b_current_landuse, "04b_current_complete_landuse_qgis_palette")

# ---- Subarea overview and a precise contact detail ---------------------
zoom_area <- subareas[subareas$subarea_id == zoom_subarea_id, ]
zoom_label <- paste("Subarea", zoom_subarea_id)
zoom_current <- terra::crop(fine_current_nature, zoom_area, snap = "out")
zoom_new <- terra::crop(fine_restoration, zoom_current)
# Calculate neighbours BEFORE the subarea mask so nature just across an
# internal subarea boundary can also be recognised as an attachment.
pad <- 2 * max(terra::res(fine_landuse))
ze <- as.vector(terra::ext(zoom_current))
padded_current <- terra::crop(fine_current_nature,
                              terra::ext(ze + c(-pad, pad, -pad, pad)), snap = "out")
current_binary <- terra::ifel(is.na(padded_current), 0, padded_current == 1)
rook_window <- matrix(c(NA, 1, NA, 1, NA, 1, NA, 1, NA), nrow = 3)
current_neighbours <- terra::focal(current_binary, w = rook_window,
                                   fun = "sum", na.rm = TRUE, fillvalue = 0)
current_neighbours <- terra::crop(current_neighbours, zoom_current)
contact_cells <- terra::ifel(zoom_new == 1 & current_neighbours > 0, 1, NA)
contact_cells <- terra::mask(contact_cells, zoom_area, touches = FALSE)

make_nature_status <- function(current, restoration) {
  x <- terra::ifel(restoration == 1, 2,
                   terra::ifel(current == 1, 1, NA))
  levels(x) <- data.frame(value = c(1L, 2L), status = factor(
    c("Current nature", "Newly selected nature"),
    levels = c("Current nature", "Newly selected nature")))
  names(x) <- "Nature status"
  x
}
zoom_status <- make_nature_status(zoom_current, zoom_new)
zoom_status <- terra::mask(zoom_status, zoom_area, touches = FALSE)

status_map <- function(x, boundaries, title, subtitle, exact = FALSE) {
  ggplot2::ggplot() +
    tidyterra::geom_spatraster(data = x,
                               maxcell = if (exact) terra::ncell(x) else raster_maxcell) +
    tidyterra::geom_spatvector(data = boundaries, fill = NA,
                               colour = "#082A38", linewidth = 0.6) +
    ggplot2::scale_fill_manual(name = "Nature status",
                               values = nature_status_palette, limits = names(nature_status_palette),
                               drop = FALSE, na.value = "transparent") +
    ggplot2::coord_sf(expand = FALSE) +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 1)) +
    ggplot2::labs(title = title, subtitle = subtitle,
                  caption = paste("Restoration selected at 200 m; displayed on the 10 m land-use grid.",
                                  "Blank areas: other land uses or missing data.", sep = "\n")) +
    map_theme()
}
map_5_subarea_overview <- status_map(zoom_status, zoom_area,
                                     paste(zoom_label, "- current nature and selected restoration"),
                                     "Green: existing nature. Orange: newly selected restoration.")

contact_total <- as.numeric(terra::global(!is.na(contact_cells), "sum",
                                          na.rm = TRUE)[1, 1])
if (is.null(detail_extent_override) && contact_total == 0) {
  save_additional_map(map_5_subarea_overview,
                      paste0("05_subarea", zoom_subarea_id, "_overview"))
  stop("The overview maps were saved, but no rook-contact cell was found in ",
       zoom_label, ". No attachment close-up was invented. Inspect the overview ",
       "and the original reconstruction. You can set detail_extent_override.")
}

if (is.null(detail_extent_override)) {
  fact <- pmax(1L, as.integer(round(contact_search_block_m /
                                      terra::res(contact_cells))))
  blocks <- terra::aggregate(contact_cells, fact = fact, fun = "sum", na.rm = TRUE)
  block_values <- terra::values(blocks, mat = FALSE)
  block_values[!is.finite(block_values)] <- -Inf
  best_block <- which.max(block_values)
  block_xy <- terra::xyFromCell(blocks, best_block)[1, ]
  block_res <- terra::res(blocks)
  block_extent <- terra::ext(block_xy[1] - block_res[1]/2,
                             block_xy[1] + block_res[1]/2, block_xy[2] - block_res[2]/2,
                             block_xy[2] + block_res[2]/2)
  local_contacts <- terra::crop(contact_cells, block_extent, snap = "near")
  ids <- which(terra::values(local_contacts, mat = FALSE) == 1)
  xy <- terra::xyFromCell(local_contacts, ids)
  # Centre on an actual contact, not a block centroid that may lie in a gap.
  centre <- xy[which.min((xy[,1] - block_xy[1])^2 +
                           (xy[,2] - block_xy[2])^2), ]
  chosen_extent <- terra::ext(centre[1] - detail_width_m/2,
                              centre[1] + detail_width_m/2, centre[2] - detail_height_m/2,
                              centre[2] + detail_height_m/2)
} else {
  stopifnot(length(detail_extent_override) == 4L,
            all(is.finite(detail_extent_override)),
            detail_extent_override[1] < detail_extent_override[2],
            detail_extent_override[3] < detail_extent_override[4])
  chosen_extent <- terra::ext(detail_extent_override)
}

# Context may include neighbouring subareas; stays inside the planning domain
# because the upstream fine rasters are already masked to that domain.
detail_current <- terra::crop(fine_current_nature, chosen_extent, snap = "out")
detail_new <- terra::crop(fine_restoration, detail_current)
detail_status <- make_nature_status(detail_current, detail_new)
detail_box <- terra::as.polygons(terra::ext(detail_status),
                                 crs = terra::crs(detail_status))
# Keep the original boundary lines, not crop-created rectangular edges.
detail_boundary_lines <- terra::crop(terra::as.lines(subareas), detail_box)
map_6_attachment_detail <- status_map(detail_status, detail_boundary_lines,
                                      paste(zoom_label, "- restoration adjoining existing nature"),
                                      "Local detail; all raster cells displayed (no plotting subsampling).", exact = TRUE)
map_5_subarea_overview <- map_5_subarea_overview +
  tidyterra::geom_spatvector(data = detail_box, fill = NA, colour = "#222222",
                             linetype = "dashed", linewidth = 0.8) +
  ggplot2::labs(subtitle = "Green: current nature. Orange: restoration. Dashed box: detail view.")
save_additional_map(map_5_subarea_overview,
                    paste0("05_subarea", zoom_subarea_id, "_overview"))
save_additional_map(map_6_attachment_detail,
                    paste0("06_subarea", zoom_subarea_id, "_attachment_detail"), width = 9, height = 8)

e <- as.vector(terra::ext(detail_status))
utils::write.csv(data.frame(solution_prefix = solution_prefix,
                            subarea_id = zoom_subarea_id, subarea_id_definition = "Original script: shapefile row order",
                            detail_selection = if (is.null(detail_extent_override)) "Rook-contact, strongest 1 km block"
                            else "Manual extent",
                            contact_cells_in_subarea = contact_total,
                            xmin = e[1], xmax = e[2], ymin = e[3], ymax = e[4],
                            grid_resolution_x = terra::res(detail_status)[1],
                            grid_resolution_y = terra::res(detail_status)[2]),
                 file.path(additional_figure_dir, "closeup_source_audit.csv"), row.names = FALSE)
message("Finished: ", normalizePath(additional_figure_dir, winslash = "/"))
# In RStudio: print(map_4_current_nature), print(map_5_subarea_overview),
#            print(map_6_attachment_detail).
