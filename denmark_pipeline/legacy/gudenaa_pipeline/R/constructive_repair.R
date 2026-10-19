stage1_selected <- readr::read_tsv(
  stage1_selected_path,
  show_col_types = FALSE
) |>
  transmute(
    coarse_cell_id = as.integer(coarse_cell_id),
    landuse = as.character(landuse),
    restored_ha = as.numeric(restored_ha)
  ) |>
  arrange(coarse_cell_id)

stage1_network_cells <- readr::read_tsv(
  stage1_network_cells_path,
  show_col_types = FALSE
) |>
  transmute(
    network_name = as.character(network_name),
    coarse_cell_id = as.integer(coarse_cell_id)
  ) |>
  distinct() |>
  arrange(network_name, coarse_cell_id)

stage1_network_patches <- readr::read_tsv(
  stage1_network_patches_path,
  show_col_types = FALSE
) |>
  transmute(
    network_name = as.character(network_name),
    patch_id = as.integer(patch_id),
    is_root = as.integer(is_root)
  ) |>
  distinct() |>
  arrange(network_name, patch_id)

stopifnot(
  !anyDuplicated(stage1_selected$coarse_cell_id),
  setequal(unique(stage1_network_cells$network_name), required_networks),
  setequal(unique(stage1_network_patches$network_name), required_networks),
  all(stage1_network_cells$coarse_cell_id %in%
        stage1_selected$coarse_cell_id),
  all(stage1_network_patches$is_root %in% c(0L, 1L)),
  all(
    stage1_network_patches |>
      count(network_name, wt = is_root, name = "roots") |>
      pull(roots) == 1L
  )
)

invalid_stage1_actions <- stage1_selected |>
  anti_join(allowed_actions, by = c("landuse", "coarse_cell_id"))
if (nrow(invalid_stage1_actions)) {
  stop("The saved Stage-1 incumbent contains a forbidden action.")
}

stage1_selected_ids <- stage1_selected$coarse_cell_id
stage1_cell_sets <- split(
  stage1_network_cells$coarse_cell_id,
  stage1_network_cells$network_name
)
stage1_patch_sets <- split(
  stage1_network_patches$patch_id,
  stage1_network_patches$network_name
)

cell_area_lookup <- stats::setNames(
  candidate_table$EligibleAreaHa,
  candidate_table$coarse_cell_id
)
patch_area_lookup <- stats::setNames(
  patch_table$ExistingPatchAreaHa,
  patch_table$patch_id
)

stage1_selected_area_audit <- stage1_selected |>
  rename(SavedRestoredHa = restored_ha) |>
  left_join(
    candidate_table |>
      select(coarse_cell_id, ModelRestorableAreaHa = EligibleAreaHa),
    by = "coarse_cell_id"
  ) |>
  mutate(
    DifferenceHa = SavedRestoredHa - ModelRestorableAreaHa,
    Passes = !is.na(ModelRestorableAreaHa) & abs(DifferenceHa) <= 1e-6
  )
if (!all(stage1_selected_area_audit$Passes)) {
  stop(
    "Saved Stage-1 restored areas do not match the current candidate ",
    "dictionary. Do not construct a start from mismatched inputs."
  )
}
readr::write_csv(
  stage1_selected_area_audit,
  file.path(problem_dir, "constructive_stage1_selected_area_audit.csv")
)

stage1_global_certificate_ha <-
  sum(unname(cell_area_lookup[
    as.character(stage1_cell_sets[["Global"]])
  ])) +
  sum(unname(patch_area_lookup[
    as.character(stage1_patch_sets[["Global"]])
  ]))

# Directed copies of the undirected graph allow the cost of entering a cell
# to equal its added restoration area. Existing patches and already-selected
# cells have zero entry cost.
repair_vertex_cost <- potential_vertices |>
  transmute(
    vertex,
    entry_cost_ha = if_else(
      node_type == "cell" & !node_id %in% stage1_selected_ids,
      FlowAreaHa,
      0
    )
  )

repair_directed_edges <- bind_rows(
  potential_edges,
  potential_edges |>
    transmute(
      reverse_from = to,
      reverse_to = from
    ) |>
    rename(from = reverse_from, to = reverse_to)
) |>
  distinct() |>
  left_join(
    repair_vertex_cost,
    by = c("to" = "vertex")
  ) |>
  mutate(
    # Candidate areas occur in centi-hectare increments. This negligible
    # positive tie-breaker selects the shorter route among paths with the
    # same added restoration area and avoids gratuitous zero-cost detours.
    weight = entry_cost_ha + 1e-8
  )

repair_graph <- igraph::graph_from_data_frame(
  repair_directed_edges,
  directed = TRUE,
  vertices = potential_vertices
)

global_source_vertices <- c(
  paste0("c_", stage1_cell_sets[["Global"]]),
  paste0("p_", stage1_patch_sets[["Global"]])
)

find_subarea_connection <- function(target_network) {
  target_vertices <- c(
    paste0("c_", stage1_cell_sets[[target_network]]),
    paste0("p_", stage1_patch_sets[[target_network]])
  )

  distance_matrix <- igraph::distances(
    repair_graph,
    v = global_source_vertices,
    to = target_vertices,
    mode = "out",
    weights = igraph::E(repair_graph)$weight
  )

  if (!any(is.finite(distance_matrix))) {
    return(tibble(
      target_network = target_network,
      reachable = FALSE,
      path_added_ha = Inf,
      certified_area_ha = NA_real_,
      path_vertices = list(character()),
      merged_cells = list(integer()),
      merged_patches = list(integer())
    ))
  }

  minimum_distance <- min(distance_matrix, na.rm = TRUE)
  minimum_index <- which(
    distance_matrix == minimum_distance,
    arr.ind = TRUE
  )[1, , drop = FALSE]
  start_vertex <- rownames(distance_matrix)[minimum_index[1, 1]]
  end_vertex <- colnames(distance_matrix)[minimum_index[1, 2]]

  path_vertices <- igraph::shortest_paths(
    repair_graph,
    from = start_vertex,
    to = end_vertex,
    mode = "out",
    weights = igraph::E(repair_graph)$weight,
    output = "vpath"
  )$vpath[[1]] |>
    igraph::as_ids()

  path_cells <- as.integer(stringr::str_remove(
    path_vertices[stringr::str_starts(path_vertices, "c_")],
    "^c_"
  ))
  path_patches <- as.integer(stringr::str_remove(
    path_vertices[stringr::str_starts(path_vertices, "p_")],
    "^p_"
  ))

  merged_cells <- sort(unique(c(
    stage1_cell_sets[["Global"]],
    stage1_cell_sets[[target_network]],
    path_cells
  )))
  merged_patches <- sort(unique(c(
    stage1_patch_sets[["Global"]],
    stage1_patch_sets[[target_network]],
    path_patches
  )))

  # The AMPL formulation deliberately closes patch membership: every
  # existing-nature patch touched by a member cell is itself a member. The
  # connecting path can contain cells that were selected but did not belong
  # to any Stage-1 certificate, so add all of their contacted patches here.
  contacted_patches <- cell_patch_links |>
    filter(coarse_cell_id %in% merged_cells) |>
    pull(patch_id)
  merged_patches <- sort(unique(c(merged_patches, contacted_patches)))

  certified_area_ha <-
    sum(unname(cell_area_lookup[as.character(merged_cells)])) +
    sum(unname(patch_area_lookup[as.character(merged_patches)]))

  tibble(
    target_network = target_network,
    reachable = TRUE,
    path_added_ha = sum(unname(cell_area_lookup[
      as.character(setdiff(path_cells, stage1_selected_ids))
    ])),
    certified_area_ha = certified_area_ha,
    path_vertices = list(path_vertices),
    merged_cells = list(merged_cells),
    merged_patches = list(merged_patches)
  )
}

connection_candidates <- purrr::map_dfr(
  setdiff(required_networks, "Global"),
  find_subarea_connection
)

connection_audit <- connection_candidates |>
  mutate(
    path_nodes = purrr::map_int(path_vertices, length),
    global_member_cells = purrr::map_int(merged_cells, length),
    global_member_patches = purrr::map_int(merged_patches, length),
    chosen = FALSE
  ) |>
  select(
    target_network,
    reachable,
    path_added_ha,
    certified_area_ha,
    path_nodes,
    global_member_cells,
    global_member_patches,
    chosen
  )
readr::write_csv(
  connection_audit,
  file.path(problem_dir, "constructive_connection_candidates.csv")
)

feasible_connections <- connection_candidates |>
  filter(
    reachable,
    certified_area_ha >= global_network_minimum_ha - 1e-6
  ) |>
  arrange(path_added_ha, desc(certified_area_ha), target_network)

if (!nrow(feasible_connections)) {
  stop(
    "No completed Stage-1 subarea certificate can be joined to the ",
    "partial Global certificate to reach 5,000 ha."
  )
}

chosen_connection <- feasible_connections |>
  slice_head(n = 1L)
chosen_target_network <- chosen_connection$target_network[[1]]
constructive_global_cells <- chosen_connection$merged_cells[[1]]
constructive_global_patches <- chosen_connection$merged_patches[[1]]
chosen_path_vertices <- chosen_connection$path_vertices[[1]]

chosen_path_cells <- as.integer(stringr::str_remove(
  chosen_path_vertices[stringr::str_starts(chosen_path_vertices, "c_")],
  "^c_"
))
chosen_path_patches <- as.integer(stringr::str_remove(
  chosen_path_vertices[stringr::str_starts(chosen_path_vertices, "p_")],
  "^p_"
))
addition_ids <- sort(setdiff(chosen_path_cells, stage1_selected_ids))

connection_audit <- connection_audit |>
  mutate(chosen = target_network == chosen_target_network)
readr::write_csv(
  connection_audit,
  file.path(problem_dir, "constructive_connection_candidates.csv")
)

constructive_path_audit <- tibble(vertex = chosen_path_vertices) |>
  mutate(
    node_type = stringr::str_sub(vertex, 1, 1),
    node_id = as.integer(stringr::str_remove(vertex, "^[cp]_")),
    was_selected = node_type == "c" & node_id %in% stage1_selected_ids,
    added_by_08b = node_type == "c" & node_id %in% addition_ids
  ) |>
  left_join(
    potential_vertices |>
      select(vertex, FlowAreaHa),
    by = "vertex"
  )
readr::write_csv(
  constructive_path_audit,
  file.path(problem_dir, "constructive_global_connection_path.csv")
)

# Choose valid destination labels for newly added cells. This is only the
# initial assignment; Stage 3 remains free to improve it.
addition_assignment <- objective_table |>
  filter(coarse_cell_id %in% addition_ids) |>
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
  select(coarse_cell_id, landuse)

if (nrow(addition_assignment) != length(addition_ids)) {
  stop("Could not assign every constructive path addition.")
}

constructive_network_cells <- bind_rows(
  stage1_network_cells |>
    filter(network_name != "Global"),
  tibble(
    network_name = "Global",
    coarse_cell_id = constructive_global_cells
  )
) |>
  distinct() |>
  arrange(network_name, coarse_cell_id)

global_root_patch <- stage1_network_patches |>
  filter(network_name == "Global", is_root == 1L) |>
  pull(patch_id)
stopifnot(length(global_root_patch) == 1L)

constructive_network_patches <- bind_rows(
  stage1_network_patches |>
    filter(network_name != "Global"),
  tibble(
    network_name = "Global",
    patch_id = constructive_global_patches,
    is_root = as.integer(constructive_global_patches == global_root_patch)
  )
) |>
  distinct() |>
  arrange(network_name, patch_id)

invalid_constructive_network_cells <- constructive_network_cells |>
  anti_join(
    eligible_network_cells |>
      select(network_name, coarse_cell_id),
    by = c("network_name", "coarse_cell_id")
  )
invalid_constructive_network_patches <- constructive_network_patches |>
  anti_join(
    eligible_network_patches |>
      select(network_name, patch_id),
    by = c("network_name", "patch_id")
  )
if (
  nrow(invalid_constructive_network_cells) ||
  nrow(invalid_constructive_network_patches)
) {
  stop(
    "The constructive certificate uses a node outside its exact eligible ",
    "potential component."
  )
}

protected_cell_ids <- unique(constructive_network_cells$coarse_cell_id)

stage1_restoration_ha <- sum(unname(cell_area_lookup[
  as.character(stage1_selected_ids)
]))
addition_area_ha <- sum(unname(cell_area_lookup[
  as.character(addition_ids)
]))
pre_removal_restoration_ha <- stage1_restoration_ha + addition_area_ha
minimum_restoration_ha <- global_budget$GlobalRestorationRequirementHa
maximum_restoration_ha <-
  minimum_restoration_ha + global_restoration_tolerance_ha
minimum_removal_ha <- max(
  0,
  pre_removal_restoration_ha - maximum_restoration_ha
)
maximum_removal_ha <-
  pre_removal_restoration_ha - minimum_restoration_ha

removal_candidates <- stage1_selected |>
  filter(!coarse_cell_id %in% protected_cell_ids) |>
  left_join(
    objective_table |>
      select(
        coarse_cell_id,
        landuse,
        BioDiversity,
        ExistingNatureContact
      ),
    by = c("coarse_cell_id", "landuse")
  ) |>
  arrange(
    BioDiversity,
    ExistingNatureContact,
    restored_ha,
    coarse_cell_id
  )

choose_removal_subset <- function(
  candidates,
  minimum_ha,
  maximum_ha
) {
  if (minimum_ha <= 1e-9) return(integer())

  minimum_units <- as.integer(ceiling((minimum_ha - 1e-9) * 100))
  maximum_units <- as.integer(floor((maximum_ha + 1e-9) * 100))
  if (maximum_units < minimum_units) {
    stop("No numerical room exists for the required area removal.")
  }

  item_units <- as.integer(round(candidates$restored_ha * 100))
  reachable <- rep(FALSE, maximum_units + 1L)
  predecessor_sum <- rep(NA_integer_, maximum_units + 1L)
  predecessor_item <- rep(NA_integer_, maximum_units + 1L)
  reachable[[1]] <- TRUE

  for (item in seq_len(nrow(candidates))) {
    units <- item_units[[item]]
    if (units <= 0L || units > maximum_units) next
    for (old_sum in seq.int(maximum_units - units, 0L, by = -1L)) {
      old_index <- old_sum + 1L
      new_sum <- old_sum + units
      new_index <- new_sum + 1L
      if (reachable[[old_index]] && !reachable[[new_index]]) {
        reachable[[new_index]] <- TRUE
        predecessor_sum[[new_index]] <- old_sum
        predecessor_item[[new_index]] <- item
      }
    }
  }

  target_units <- seq.int(minimum_units, maximum_units)
  target_units <- target_units[reachable[target_units + 1L]]
  if (!length(target_units)) {
    stop(
      "Could not find a centi-hectare removal subset that preserves the ",
      "global 30% interval."
    )
  }

  current_sum <- target_units[[1]]
  chosen_items <- integer()
  while (current_sum > 0L) {
    item <- predecessor_item[[current_sum + 1L]]
    if (is.na(item)) stop("Internal subset-sum reconstruction error.")
    chosen_items <- c(chosen_items, item)
    current_sum <- predecessor_sum[[current_sum + 1L]]
  }
  sort(unique(chosen_items))
}

removal_rows <- choose_removal_subset(
  removal_candidates,
  minimum_removal_ha,
  maximum_removal_ha
)
removal_assignment <- removal_candidates[removal_rows, , drop = FALSE]
removal_ids <- removal_assignment$coarse_cell_id

constructive_assignment <- bind_rows(
  stage1_selected |>
    filter(!coarse_cell_id %in% removal_ids) |>
    select(coarse_cell_id, landuse),
  addition_assignment
) |>
  distinct() |>
  arrange(coarse_cell_id)

if (anyDuplicated(constructive_assignment$coarse_cell_id)) {
  stop("Constructive assignment contains duplicate cell decisions.")
}
if (nrow(constructive_assignment |>
  anti_join(allowed_actions, by = c("landuse", "coarse_cell_id")))) {
  stop("Constructive assignment contains a forbidden land-use action.")
}
if (!all(protected_cell_ids %in% constructive_assignment$coarse_cell_id)) {
  stop("Area balancing removed a certificate member cell.")
}

constructive_restoration_ha <- candidate_table |>
  filter(coarse_cell_id %in% constructive_assignment$coarse_cell_id) |>
  summarise(area_ha = sum(EligibleAreaHa)) |>
  pull(area_ha)

if (
  constructive_restoration_ha < minimum_restoration_ha - 1e-6 ||
  constructive_restoration_ha > maximum_restoration_ha + 1e-6
) {
  stop("Constructive repair violates the whole-polygon 30% interval.")
}

readr::write_csv(
  addition_assignment |>
    left_join(
      candidate_table |>
        select(coarse_cell_id, EligibleAreaHa),
      by = "coarse_cell_id"
    ),
  file.path(problem_dir, "constructive_repair_additions.csv")
)
readr::write_csv(
  removal_assignment,
  file.path(problem_dir, "constructive_repair_removals.csv")
)

# Build deterministic spanning-tree flows for the completed certificates.
build_constructive_certificate_start <- function(network_key) {
  member_cells <- constructive_network_cells |>
    filter(.data$network_name == network_key) |>
    pull(coarse_cell_id) |>
    sort()
  member_patches <- constructive_network_patches |>
    filter(.data$network_name == network_key) |>
    pull(patch_id) |>
    sort()
  root_patch <- constructive_network_patches |>
    filter(
      .data$network_name == network_key,
      is_root == 1L
    ) |>
    pull(patch_id)

  if (length(root_patch) != 1L) {
    stop("Certificate ", network_key, " does not have exactly one root.")
  }

  member_vertices <- bind_rows(
    candidate_table |>
      filter(coarse_cell_id %in% member_cells) |>
      transmute(
        vertex = paste0("c_", coarse_cell_id),
        node_type = "cell",
        node_id = coarse_cell_id,
        FlowAreaHa = EligibleAreaHa
      ),
    patch_table |>
      filter(patch_id %in% member_patches) |>
      transmute(
        vertex = paste0("p_", patch_id),
        node_type = "patch",
        node_id = patch_id,
        FlowAreaHa = ExistingPatchAreaHa
      )
  )
  member_vertex_names <- member_vertices$vertex
  member_edges <- potential_edges |>
    filter(
      from %in% member_vertex_names,
      to %in% member_vertex_names
    )
  member_graph <- igraph::graph_from_data_frame(
    member_edges,
    directed = FALSE,
    vertices = member_vertices
  )
  if (igraph::components(member_graph)$no != 1L) {
    stop("Certificate ", network_key, " is disconnected.")
  }

  root_vertex <- paste0("p_", root_patch)
  distances_to_root <- as.numeric(igraph::distances(
    member_graph,
    v = root_vertex,
    to = igraph::V(member_graph),
    mode = "all"
  )[1, ])
  names(distances_to_root) <- igraph::V(member_graph)$name

  parent <- stats::setNames(
    rep(NA_character_, length(member_vertex_names)),
    member_vertex_names
  )
  for (vertex in setdiff(member_vertex_names, root_vertex)) {
    neighbors <- igraph::as_ids(igraph::neighbors(
      member_graph,
      vertex,
      mode = "all"
    ))
    possible_parents <- neighbors[
      distances_to_root[neighbors] ==
        distances_to_root[[vertex]] - 1
    ]
    if (!length(possible_parents)) {
      stop("Could not assign a flow-tree parent for ", vertex, ".")
    }
    parent[[vertex]] <- sort(possible_parents)[[1]]
  }

  subtree_area <- stats::setNames(
    member_vertices$FlowAreaHa[
      match(member_vertex_names, member_vertices$vertex)
    ],
    member_vertex_names
  )
  flow_records <- list()
  flow_index <- 0L
  descending_vertices <- member_vertex_names[
    order(distances_to_root[member_vertex_names], decreasing = TRUE)
  ]
  for (vertex in descending_vertices) {
    if (vertex == root_vertex) next
    parent_vertex <- parent[[vertex]]
    flow_index <- flow_index + 1L
    flow_records[[flow_index]] <- tibble(
      network_name = network_key,
      child = vertex,
      parent = parent_vertex,
      flow_ha = subtree_area[[vertex]]
    )
    subtree_area[[parent_vertex]] <-
      subtree_area[[parent_vertex]] + subtree_area[[vertex]]
  }

  list(
    network_name = network_key,
    cells = member_cells,
    patches = member_patches,
    root_patch = root_patch,
    collected_area_ha = subtree_area[[root_vertex]],
    flows = bind_rows(flow_records) |>
      mutate(
        child_type = stringr::str_sub(child, 1, 1),
        parent_type = stringr::str_sub(parent, 1, 1),
        child_id = as.integer(stringr::str_remove(child, "^[cp]_")),
        parent_id = as.integer(stringr::str_remove(parent, "^[cp]_"))
      )
  )
}

certificate_starts <- purrr::map(
  required_networks,
  build_constructive_certificate_start
)
names(certificate_starts) <- required_networks

start_certificate_cell_members <- constructive_network_cells
start_certificate_patch_members <- constructive_network_patches |>
  select(network_name, patch_id)
start_certificate_roots <- constructive_network_patches |>
  filter(is_root == 1L) |>
  transmute(
    network_name,
    patch_id,
    collected_area_ha = purrr::map_dbl(
      network_name,
      function(k) certificate_starts[[k]]$collected_area_ha
    )
  )
start_certificate_flows <- purrr::map_dfr(
  certificate_starts,
  "flows"
)

invalid_cell_flows <- start_certificate_flows |>
  filter(child_type == "c", parent_type == "c") |>
  transmute(
    network_name,
    i = child_id,
    j = parent_id
  ) |>
  anti_join(
    eligible_network_arcs,
    by = c("network_name", "i", "j")
  )
invalid_cell_patch_flows <- start_certificate_flows |>
  filter(child_type != parent_type) |>
  transmute(
    network_name,
    coarse_cell_id = if_else(
      child_type == "c",
      child_id,
      parent_id
    ),
    patch_id = if_else(
      child_type == "p",
      child_id,
      parent_id
    )
  ) |>
  anti_join(
    eligible_network_cell_patch_edges,
    by = c("network_name", "coarse_cell_id", "patch_id")
  )
if (nrow(invalid_cell_flows) || nrow(invalid_cell_patch_flows)) {
  stop("The constructive flow tree contains an ineligible arc.")
}

start_certificate_area_check <- bind_rows(
  start_certificate_cell_members |>
    left_join(
      network_cell_area,
      by = c("network_name", "coarse_cell_id")
    ) |>
    transmute(
      network_name,
      area_ha = coalesce(NetworkCellAreaHa, 0)
    ),
  start_certificate_patch_members |>
    left_join(
      network_patch_area,
      by = c("network_name", "patch_id")
    ) |>
    transmute(
      network_name,
      area_ha = coalesce(NetworkPatchAreaHa, 0)
    )
) |>
  group_by(network_name) |>
  summarise(CertifiedAreaHa = sum(area_ha), .groups = "drop") |>
  mutate(
    RequiredAreaHa = unname(network_minimums[network_name]),
    Passes = CertifiedAreaHa >= RequiredAreaHa - 1e-6
  )

if (
  nrow(start_certificate_area_check) != length(required_networks) ||
  !all(start_certificate_area_check$Passes)
) {
  stop("The constructive certificate-area audit failed.")
}

constructive_summary <- tibble(
  chosen_connection = chosen_target_network,
  stage1_global_certificate_ha = stage1_global_certificate_ha,
  constructed_global_certificate_ha =
    start_certificate_area_check$CertifiedAreaHa[
      start_certificate_area_check$network_name == "Global"
    ],
  added_cells = length(addition_ids),
  added_restoration_ha = addition_area_ha,
  removed_cells = length(removal_ids),
  removed_restoration_ha = sum(removal_assignment$restored_ha),
  final_restoration_ha = constructive_restoration_ha,
  final_nature_ha =
    global_budget$GlobalCurrentNatureAreaHa + constructive_restoration_ha,
  final_nature_percent = 100 * final_nature_ha /
    global_budget$GlobalPolygonAreaHa,
  all_five_certificates_pass = all(start_certificate_area_check$Passes)
)

readr::write_csv(
  start_certificate_area_check,
  file.path(problem_dir, "constructive_certificate_area_audit.csv")
)
readr::write_csv(
  constructive_assignment |>
    left_join(
      candidate_table |>
        select(coarse_cell_id, RestorableAreaHa = EligibleAreaHa),
      by = "coarse_cell_id"
    ),
  file.path(problem_dir, "constructive_feasible_selected_cells.csv")
)
readr::write_csv(
  start_certificate_cell_members,
  file.path(problem_dir, "constructive_certificate_cells.csv")
)
readr::write_csv(
  constructive_network_patches,
  file.path(problem_dir, "constructive_certificate_patches.csv")
)
readr::write_csv(
  start_certificate_roots,
  file.path(problem_dir, "constructive_certificate_roots.csv")
)
readr::write_csv(
  start_certificate_flows,
  file.path(problem_dir, "constructive_certificate_flow_start.csv")
)
readr::write_csv(
  constructive_summary,
  file.path(problem_dir, "constructive_repair_summary.csv")
)

print(knitr::kable(constructive_summary, digits = 4))
knitr::kable(start_certificate_area_check, digits = 2)
