required_solution_paths <- c(
  solution_path,
  subarea_solution_path,
  network_cell_solution_path,
  network_patch_solution_path,
  network_root_solution_path,
  network_solution_summary_path,
  summary_path,
  repair_stage_summary_path,
  repair_completion_path
)

solution_is_current <- all(file.exists(required_solution_paths))
if (solution_is_current) {
  completion_status <- trimws(readLines(
    repair_completion_path,
    warn = FALSE,
    n = 1L
  ))
  solution_is_current <- identical(completion_status, "COMPLETE")
}

if (!solution_is_current) {
  message(
    "No current repaired-network solution found. Run the generated .run file, ",
    "then re-knit this notebook."
  )
} else {
  selected_solution <- readr::read_tsv(
    solution_path,
    show_col_types = FALSE
  ) |>
    mutate(coarse_cell_id = as.integer(coarse_cell_id))
  ampl_subarea_solution <- readr::read_tsv(
    subarea_solution_path,
    show_col_types = FALSE
  ) |>
    mutate(subarea_id = as.integer(subarea_id))
  network_cell_solution <- readr::read_tsv(
    network_cell_solution_path,
    show_col_types = FALSE
  ) |>
    mutate(coarse_cell_id = as.integer(coarse_cell_id))
  network_patch_solution <- readr::read_tsv(
    network_patch_solution_path,
    show_col_types = FALSE
  ) |>
    mutate(patch_id = as.integer(patch_id))
  network_root_solution <- readr::read_tsv(
    network_root_solution_path,
    show_col_types = FALSE
  ) |>
    mutate(patch_id = as.integer(patch_id))
  ampl_network_summary <- readr::read_tsv(
    network_solution_summary_path,
    show_col_types = FALSE
  )
  repair_stage_summary <- readr::read_tsv(
    repair_stage_summary_path,
    show_col_types = FALSE
  )
  solution_summary <- readr::read_tsv(
    summary_path,
    show_col_types = FALSE,
    col_types = cols(.default = col_character())
  )

  selected_ids <- selected_solution$coarse_cell_id
  stopifnot(
    nrow(selected_solution) > 0L,
    !anyDuplicated(selected_ids),
    all(selected_ids %in% candidate_cells),
    all(selected_solution$landuse %in% nature_landuses),
    setequal(ampl_network_summary$network_name, required_networks),
    setequal(network_root_solution$network_name, required_networks),
    all(ampl_network_summary$shortfall_ha <=
      network_shortfall_tolerance_ha),
    all(ampl_network_summary$certificate_built == 1L),
    "3_ecological_contiguity_quality" %in%
      repair_stage_summary$stage
  )

  # ---- Exact change from the successful notebook-07 footprint -----------
  footprint_repair_cells <- candidate_table |>
    select(coarse_cell_id, EligibleAreaHa) |>
    left_join(
      warm_start_assignment |>
        select(coarse_cell_id, warm_landuse = landuse),
      by = "coarse_cell_id"
    ) |>
    left_join(
      selected_solution |>
        select(coarse_cell_id, repaired_landuse = landuse),
      by = "coarse_cell_id"
    ) |>
    mutate(
      WarmSelected = as.integer(!is.na(warm_landuse)),
      RepairedSelected = as.integer(!is.na(repaired_landuse)),
      change_type = case_when(
        WarmSelected == 1L & RepairedSelected == 1L ~ "kept",
        WarmSelected == 1L & RepairedSelected == 0L ~ "dropped",
        WarmSelected == 0L & RepairedSelected == 1L ~ "added",
        TRUE ~ "unselected"
      )
    )

  footprint_repair_summary <- footprint_repair_cells |>
    group_by(change_type) |>
    summarise(
      cells = n(),
      area_ha = sum(EligibleAreaHa),
      .groups = "drop"
    ) |>
    arrange(match(
      change_type,
      c("kept", "dropped", "added", "unselected")
    ))

  computed_footprint_change_ha <- footprint_repair_cells |>
    filter(change_type %in% c("dropped", "added")) |>
    summarise(value = sum(EligibleAreaHa)) |>
    pull(value)
  reported_footprint_change_ha <- solution_summary |>
    filter(item == "warm_footprint_change_ha") |>
    pull(value) |>
    as.numeric()
  if (
    length(reported_footprint_change_ha) != 1L ||
      abs(
        computed_footprint_change_ha - reported_footprint_change_ha
      ) > 1e-5
  ) {
    stop("Warm-footprint change audit does not reconcile.", call. = FALSE)
  }

  readr::write_csv(
    footprint_repair_cells |>
      filter(change_type != "unselected"),
    file.path(problem_dir, "warm_footprint_repair_cells.csv")
  )
  readr::write_csv(
    footprint_repair_summary,
    file.path(problem_dir, "warm_footprint_repair_summary.csv")
  )

  # ---- Whole-polygon 30% budget and subarea reference table -------------
  independent_subarea_solution <-
    restorable_area_by_subarea_cell |>
    filter(coarse_cell_id %in% selected_ids) |>
    group_by(subarea_id) |>
    summarise(
      IndependentlyRestoredHa = sum(RestorableAreaHaBySubarea),
      .groups = "drop"
    ) |>
    right_join(subarea_budget, by = "subarea_id") |>
    mutate(
      IndependentlyRestoredHa = coalesce(IndependentlyRestoredHa, 0),
      FinalNatureAreaHa =
        CurrentNatureAreaHa + IndependentlyRestoredHa,
      FinalNaturePercent =
        100 * FinalNatureAreaHa / PolygonAreaHa,
      ReferenceThirtyPercentDifferenceHa =
        FinalNatureAreaHa - FinalNatureTargetHa
    ) |>
    arrange(subarea_id)

  if (max(abs(
    ampl_subarea_solution$restored_ha -
      independent_subarea_solution$IndependentlyRestoredHa
  )) > 1e-5) {
    stop("AMPL and R subarea totals disagree.", call. = FALSE)
  }

  independent_global_solution <-
    independent_subarea_solution |>
    summarise(
      IndependentlyRestoredHa = sum(IndependentlyRestoredHa)
    ) |>
    bind_cols(global_budget) |>
    mutate(
      FinalNatureAreaHa =
        GlobalCurrentNatureAreaHa + IndependentlyRestoredHa,
      FinalNaturePercent =
        100 * FinalNatureAreaHa / GlobalPolygonAreaHa,
      ExcessAboveThirtyPercentHa =
        FinalNatureAreaHa - GlobalFinalNatureTargetHa,
      PassesGlobalThirtyPercentMinimum =
        IndependentlyRestoredHa >=
          GlobalRestorationRequirementHa - 1e-6,
      PassesGlobalThirtyPercentMaximum =
        IndependentlyRestoredHa <=
          GlobalMaximumRestorationHa + 1e-6,
      PassesGlobalThirtyPercentBudget =
        PassesGlobalThirtyPercentMinimum &
          PassesGlobalThirtyPercentMaximum
    )

  if (!independent_global_solution$PassesGlobalThirtyPercentBudget[[1]]) {
    stop("Independent whole-polygon 30% audit failed.", call. = FALSE)
  }

  # ---- Production-forest destination rule -------------------------------
  forest_rule_audit <- selected_solution |>
    select(coarse_cell_id, landuse) |>
    left_join(
      candidate_table |>
        select(
          coarse_cell_id,
          ProductionForestHa,
          AgricultureHa,
          ProductionForestShare,
          IsProductionForestCell
        ),
      by = "coarse_cell_id"
    ) |>
    mutate(
      PassesForestRule =
        IsProductionForestCell == 0L |
        landuse %in% forest_landuses
    )
  if (!all(forest_rule_audit$PassesForestRule)) {
    stop("Production-forest destination audit failed.", call. = FALSE)
  }

  # ---- Recalculate the five AMPL certificates ---------------------------
  independent_certificate_summary <- tibble(
    network_name = required_networks,
    RequiredAreaHa = unname(network_minimums[required_networks])
  ) |>
    left_join(
      bind_rows(
        network_cell_solution |>
          left_join(
            network_cell_area,
            by = c("network_name", "coarse_cell_id")
          ) |>
          transmute(
            network_name,
            area_ha = coalesce(NetworkCellAreaHa, 0),
            member_cells = 1L,
            member_patches = 0L
          ),
        network_patch_solution |>
          left_join(
            network_patch_area,
            by = c("network_name", "patch_id")
          ) |>
          transmute(
            network_name,
            area_ha = coalesce(NetworkPatchAreaHa, 0),
            member_cells = 0L,
            member_patches = 1L
          )
      ) |>
        group_by(network_name) |>
        summarise(
          IndependentlyCertifiedAreaHa = sum(area_ha),
          member_cells = sum(member_cells),
          member_patches = sum(member_patches),
          .groups = "drop"
        ),
      by = "network_name"
    ) |>
    left_join(
      network_root_solution |>
        count(network_name, name = "root_count"),
      by = "network_name"
    ) |>
    mutate(
      PassesArea =
        IndependentlyCertifiedAreaHa >= RequiredAreaHa - 1e-6,
      PassesOneRoot = root_count == 1L,
      PassesNewRestoration = member_cells >= 1L,
      PassesCoarseCertificate =
        PassesArea & PassesOneRoot & PassesNewRestoration
    )

  if (
    anyNA(independent_certificate_summary) ||
    !all(independent_certificate_summary$PassesCoarseCertificate)
  ) {
    stop("Independent AMPL certificate audit failed.", call. = FALSE)
  }
  if (any(!network_cell_solution$coarse_cell_id %in% selected_ids)) {
    stop("A certificate contains an unselected cell.", call. = FALSE)
  }
  if (max(abs(
    ampl_network_summary$certified_area_ha[
      match(
        independent_certificate_summary$network_name,
        ampl_network_summary$network_name
      )
    ] - independent_certificate_summary$IndependentlyCertifiedAreaHa
  )) > 1e-5) {
    stop("AMPL and R certificate totals disagree.", call. = FALSE)
  }

  readr::write_csv(
    independent_subarea_solution,
    file.path(problem_dir, "independent_subarea_solution_audit.csv")
  )
  readr::write_csv(
    independent_global_solution,
    file.path(problem_dir, "independent_global_30pct_solution_audit.csv")
  )
  readr::write_csv(
    forest_rule_audit,
    file.path(problem_dir, "production_forest_solution_audit.csv")
  )
  readr::write_csv(
    independent_certificate_summary,
    file.path(problem_dir, "independent_coarse_certificate_audit.csv")
  )

  # ---- Exact 10 m restoration + full national patch network -------------
  fine_restoration <- terra::ifel(
    fine_eligible == 1 & coarse_id_on_fine %in% selected_ids,
    1,
    0
  )
  fine_restoration <- terra::mask(fine_restoration, fine_landuse)
  names(fine_restoration) <- "RestorationCells"

  final_nature <- terra::ifel(
    fine_nature == 1 | fine_restoration == 1,
    1,
    0
  )
  final_nature <- terra::mask(final_nature, fine_landuse)
  names(final_nature) <- "FinalNatureCells"

  # This clipped raster remains useful for mapping final nature inside the
  # planning polygon. It is not used to truncate national patch area in the
  # 5,000 ha calculation below.
  exact_component_raster <- terra::patches(
    final_nature,
    directions = 4,
    zeroAsNA = TRUE
  )
  names(exact_component_raster) <- "component_id"

  # Restoration is first split into exact 10 m rook components. National
  # existing patches then act as supernodes. This merges restoration pieces
  # that connect through the same current patch, including connections that
  # continue outside the planning boundary.
  restoration_component_raster <- terra::patches(
    fine_restoration,
    directions = 4,
    zeroAsNA = TRUE
  )
  names(restoration_component_raster) <- "restoration_component_id"

  restoration_component_summary <- zonal_sum_table(
    fine_restoration,
    restoration_component_raster,
    "restoration_component_id"
  ) |>
    transmute(
      restoration_component_id = as.integer(round(
        restoration_component_id
      )),
      RestorationAreaHa =
        RestorationCells * fine_cell_area_ha
    ) |>
    filter(RestorationAreaHa > 0) |>
    arrange(desc(RestorationAreaHa))

  exact_pair_multiplier <- expected_number_of_subareas + 1L
  restoration_component_subarea_zone <- terra::ifel(
    is.na(restoration_component_raster) |
      is.na(subarea_id_on_fine),
    NA,
    restoration_component_raster * exact_pair_multiplier +
      subarea_id_on_fine
  )

  restoration_component_by_subarea <- zonal_sum_table(
    fine_restoration,
    restoration_component_subarea_zone,
    "pair_id"
  ) |>
    mutate(
      restoration_component_id = as.integer(floor(
        pair_id / exact_pair_multiplier
      )),
      subarea_id = as.integer(round(
        pair_id -
          restoration_component_id * exact_pair_multiplier
      )),
      RestorationAreaHa =
        RestorationCells * fine_cell_area_ha
    ) |>
    filter(
      subarea_id %in% seq_len(expected_number_of_subareas)
    ) |>
    select(
      restoration_component_id,
      subarea_id,
      RestorationAreaHa
    ) |>
    arrange(subarea_id, desc(RestorationAreaHa))

  restoration_component_polygons <- terra::as.polygons(
    restoration_component_raster,
    aggregate = TRUE,
    values = TRUE,
    na.rm = TRUE
  )
  restoration_component_polygons$restoration_component_id <-
    as.integer(
      restoration_component_polygons$restoration_component_id
    )

  restoration_contact_buffers <- terra::buffer(
    restoration_component_polygons,
    width = contact_buffer_m
  )
  exact_contact_geometry <- terra::intersect(
    restoration_contact_buffers[, "restoration_component_id"],
    patches_near_domain[, "patch_id"]
  )

  if (nrow(exact_contact_geometry)) {
    exact_contact_geometry$contact_proxy_m2 <- terra::expanse(
      exact_contact_geometry,
      unit = "m",
      transform = FALSE
    )
    exact_restoration_patch_links <- exact_contact_geometry |>
      as.data.frame() |>
      transmute(
        restoration_component_id = as.integer(
          restoration_component_id
        ),
        patch_id = as.integer(patch_id),
        contact_proxy_m2 = as.numeric(contact_proxy_m2)
      ) |>
      group_by(restoration_component_id, patch_id) |>
      summarise(
        contact_proxy_m2 = sum(contact_proxy_m2),
        .groups = "drop"
      ) |>
      filter(contact_proxy_m2 >= minimum_contact_proxy_m2) |>
      arrange(restoration_component_id, patch_id)
  } else {
    exact_restoration_patch_links <- tibble(
      restoration_component_id = integer(),
      patch_id = integer(),
      contact_proxy_m2 = numeric()
    )
  }

  exact_linked_patch_table <- exact_restoration_patch_links |>
    distinct(patch_id) |>
    left_join(patch_dictionary_full, by = "patch_id")
  if (anyNA(exact_linked_patch_table$ExistingPatchAreaHa)) {
    stop("Exact audit found a patch missing from the national dictionary.")
  }

  exact_network_vertices <- bind_rows(
    restoration_component_summary |>
      transmute(
        vertex = paste0("r_", restoration_component_id),
        node_type = "restoration",
        node_id = restoration_component_id,
        GlobalNodeAreaHa = RestorationAreaHa
      ),
    exact_linked_patch_table |>
      transmute(
        vertex = paste0("p_", patch_id),
        node_type = "existing_patch",
        node_id = patch_id,
        GlobalNodeAreaHa = ExistingPatchAreaHa
      )
  )

  exact_network_edges <- exact_restoration_patch_links |>
    transmute(
      from = paste0("r_", restoration_component_id),
      to = paste0("p_", patch_id)
    ) |>
    distinct()

  exact_augmented_graph <- igraph::graph_from_data_frame(
    exact_network_edges,
    directed = FALSE,
    vertices = exact_network_vertices
  )
  exact_membership_raw <- igraph::components(
    exact_augmented_graph
  )$membership
  exact_membership <- tibble(
    vertex = names(exact_membership_raw),
    exact_network_component = as.integer(exact_membership_raw)
  )
  exact_network_vertices <- exact_network_vertices |>
    left_join(exact_membership, by = "vertex")

  exact_augmented_component_summary <- exact_network_vertices |>
    group_by(exact_network_component) |>
    summarise(
      restoration_components = sum(node_type == "restoration"),
      existing_national_patches =
        sum(node_type == "existing_patch"),
      RestorationAreaHa = sum(
        GlobalNodeAreaHa[node_type == "restoration"]
      ),
      ExistingPatchFullNationalAreaHa = sum(
        GlobalNodeAreaHa[node_type == "existing_patch"]
      ),
      GlobalNetworkAreaHa = sum(GlobalNodeAreaHa),
      IsAnchoredRestorationNetwork =
        restoration_components > 0 & existing_national_patches > 0,
      .groups = "drop"
    ) |>
    arrange(desc(GlobalNetworkAreaHa))

  exact_restoration_membership <- exact_network_vertices |>
    filter(node_type == "restoration") |>
    transmute(
      restoration_component_id = node_id,
      exact_network_component
    )
  exact_patch_membership <- exact_network_vertices |>
    filter(node_type == "existing_patch") |>
    transmute(
      patch_id = node_id,
      exact_network_component
    )

  exact_augmented_by_subarea <- bind_rows(
    restoration_component_by_subarea |>
      left_join(
        exact_restoration_membership,
        by = "restoration_component_id"
      ) |>
      transmute(
        exact_network_component,
        subarea_id,
        RestorationAreaHa,
        ExistingNatureAreaHa = 0
      ),
    patch_area_by_subarea_all |>
      semi_join(
        exact_linked_patch_table,
        by = "patch_id"
      ) |>
      left_join(exact_patch_membership, by = "patch_id") |>
      transmute(
        exact_network_component,
        subarea_id,
        RestorationAreaHa = 0,
        ExistingNatureAreaHa = ExistingPatchAreaHaBySubarea
      )
  ) |>
    group_by(exact_network_component, subarea_id) |>
    summarise(
      RestorationAreaHa = sum(RestorationAreaHa),
      ExistingNatureAreaHa = sum(ExistingNatureAreaHa),
      NetworkAreaInsideSubareaHa =
        RestorationAreaHa + ExistingNatureAreaHa,
      .groups = "drop"
    ) |>
    left_join(
      exact_augmented_component_summary |>
        select(
          exact_network_component,
          IsAnchoredRestorationNetwork
        ),
      by = "exact_network_component"
    ) |>
    arrange(subarea_id, desc(NetworkAreaInsideSubareaHa))

  exact_anchored_components <- exact_augmented_component_summary |>
    filter(IsAnchoredRestorationNetwork)

  if (nrow(exact_anchored_components)) {
    exact_global_row <- exact_anchored_components |>
      slice_max(GlobalNetworkAreaHa, n = 1, with_ties = FALSE) |>
      transmute(
        scope = "whole_polygon",
        subarea_id = NA_integer_,
        component_id = exact_network_component,
        RequiredAreaHa = global_network_minimum_ha,
        ExactFinalNatureAreaHa = GlobalNetworkAreaHa,
        ExactCurrentNatureAreaHa =
          ExistingPatchFullNationalAreaHa,
        ExactRestorationAreaHa = RestorationAreaHa,
        ExistingAreaInterpretation = "full_national_patch_area",
        PassesExact10m =
          ExactFinalNatureAreaHa >= RequiredAreaHa - 1e-6
      )
  } else {
    exact_global_row <- tibble(
      scope = "whole_polygon",
      subarea_id = NA_integer_,
      component_id = NA_integer_,
      RequiredAreaHa = global_network_minimum_ha,
      ExactFinalNatureAreaHa = 0,
      ExactCurrentNatureAreaHa = 0,
      ExactRestorationAreaHa = 0,
      ExistingAreaInterpretation = "full_national_patch_area",
      PassesExact10m = FALSE
    )
  }

  exact_largest_by_subarea <- exact_augmented_by_subarea |>
    filter(IsAnchoredRestorationNetwork) |>
    group_by(subarea_id) |>
    slice_max(
      NetworkAreaInsideSubareaHa,
      n = 1,
      with_ties = FALSE
    ) |>
    ungroup() |>
    transmute(
      subarea_id,
      component_id = exact_network_component,
      ExactFinalNatureAreaHa = NetworkAreaInsideSubareaHa,
      ExactCurrentNatureAreaHa = ExistingNatureAreaHa,
      ExactRestorationAreaHa = RestorationAreaHa
    )

  exact_network_goal_audit <- bind_rows(
    exact_global_row,
    tibble(
      subarea_id = seq_len(expected_number_of_subareas),
      scope = subarea_labels,
      RequiredAreaHa = subarea_network_minimum_ha
    ) |>
      left_join(exact_largest_by_subarea, by = "subarea_id") |>
      mutate(
        ExactFinalNatureAreaHa = coalesce(
          ExactFinalNatureAreaHa,
          0
        ),
        ExactCurrentNatureAreaHa = coalesce(
          ExactCurrentNatureAreaHa,
          0
        ),
        ExactRestorationAreaHa = coalesce(
          ExactRestorationAreaHa,
          0
        ),
        ExistingAreaInterpretation =
          "area_inside_this_subarea_only",
        PassesExact10m =
          ExactFinalNatureAreaHa >= RequiredAreaHa - 1e-6
      ) |>
      select(
        scope,
        subarea_id,
        component_id,
        RequiredAreaHa,
        ExactFinalNatureAreaHa,
        ExactCurrentNatureAreaHa,
        ExactRestorationAreaHa,
        ExistingAreaInterpretation,
        PassesExact10m
      )
  )

  readr::write_csv(
    restoration_component_summary,
    file.path(problem_dir, "exact_10m_restoration_components.csv")
  )
  readr::write_csv(
    restoration_component_by_subarea,
    file.path(
      problem_dir,
      "exact_10m_restoration_component_by_subarea.csv"
    )
  )
  readr::write_csv(
    exact_restoration_patch_links,
    file.path(
      problem_dir,
      "exact_10m_restoration_to_national_patch_links.csv"
    )
  )
  readr::write_csv(
    exact_augmented_component_summary,
    file.path(
      problem_dir,
      "exact_10m_augmented_network_components.csv"
    )
  )
  readr::write_csv(
    exact_augmented_by_subarea,
    file.path(
      problem_dir,
      "exact_10m_augmented_network_by_subarea.csv"
    )
  )
  readr::write_csv(
    exact_network_goal_audit,
    file.path(problem_dir, "exact_10m_network_goal_audit.csv")
  )

  # ---- Rasters for inspection and the next refinement -------------------
  selected_coarse <- coarse_cell_id
  selected_values <- rep(NA_integer_, terra::ncell(selected_coarse))
  selected_values[candidate_cells] <- 0L
  selected_values[selected_ids] <- 1L
  terra::values(selected_coarse) <- selected_values
  names(selected_coarse) <- "selected_200m"

  landuse_lookup <- tibble(
    warm_start_code = seq_along(nature_landuses),
    landuse = nature_landuses
  )
  destination_coarse <- coarse_cell_id
  destination_values <- rep(
    NA_integer_,
    terra::ncell(destination_coarse)
  )
  destination_values[selected_ids] <- landuse_lookup$warm_start_code[
    match(selected_solution$landuse, landuse_lookup$landuse)
  ]
  terra::values(destination_coarse) <- destination_values
  names(destination_coarse) <- "warm_start_landuse"

  destination_fine <- align_to_template(
    destination_coarse,
    fine_landuse,
    method = "near"
  )
  destination_fine <- terra::mask(
    destination_fine,
    fine_eligible_mask
  )

  terra::writeRaster(
    selected_coarse,
    file.path(problem_dir, "selected_coarse_cells_200m.tif"),
    overwrite = TRUE,
    datatype = "INT1U",
    gdal = c("COMPRESS=DEFLATE", "TILED=YES")
  )
  terra::writeRaster(
    fine_restoration,
    file.path(problem_dir, "repaired_solution_restoration_10m.tif"),
    overwrite = TRUE,
    datatype = "INT1U",
    gdal = c("COMPRESS=DEFLATE", "TILED=YES", "BIGTIFF=YES")
  )
  terra::writeRaster(
    final_nature,
    file.path(problem_dir, "repaired_solution_final_nature_10m.tif"),
    overwrite = TRUE,
    datatype = "INT1U",
    gdal = c("COMPRESS=DEFLATE", "TILED=YES", "BIGTIFF=YES")
  )
  terra::writeRaster(
    exact_component_raster,
    file.path(
      problem_dir,
      "repaired_solution_final_nature_components_10m.tif"
    ),
    overwrite = TRUE,
    datatype = "INT4U",
    gdal = c("COMPRESS=DEFLATE", "TILED=YES", "BIGTIFF=YES")
  )
  terra::writeRaster(
    destination_fine,
    file.path(
      problem_dir,
      "repaired_solution_landuse_warm_start_10m.tif"
    ),
    overwrite = TRUE,
    datatype = "INT1U",
    gdal = c("COMPRESS=DEFLATE", "TILED=YES", "BIGTIFF=YES")
  )
  readr::write_csv(
    landuse_lookup,
    file.path(problem_dir, "warm_start_landuse_lookup.csv")
  )

  print(independent_global_solution)
  print(independent_subarea_solution)
  print(independent_certificate_summary)
  print(exact_network_goal_audit)

  if (!all(exact_network_goal_audit$PassesExact10m)) {
    warning(
      "The AMPL 200 m certificates pass, but one or more exact 10 m ",
      "network goals fail. Use this result as the targeted input to a ",
      "10 m refinement; do not report the failed goal as achieved."
    )
  }
}
