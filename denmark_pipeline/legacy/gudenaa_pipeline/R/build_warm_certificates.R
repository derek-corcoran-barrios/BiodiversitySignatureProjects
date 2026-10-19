build_certificate_start <- function(network_name, component_id) {
  component_nodes <- warm_start_vertices |>
    filter(warm_component == component_id)

  root_row <- component_nodes |>
    filter(node_type == "patch") |>
    arrange(desc(FlowAreaHa), node_id) |>
    slice_head(n = 1L)
  if (nrow(root_row) != 1L) {
    stop("Certificate has no existing-nature root.", call. = FALSE)
  }

  root_vertex <- root_row$vertex[[1]]
  component_graph <- igraph::induced_subgraph(
    warm_start_graph,
    vids = component_nodes$vertex
  )

  distances_to_root <- as.numeric(igraph::distances(
    component_graph,
    v = root_vertex,
    to = igraph::V(component_graph),
    mode = "all"
  )[1, ])
  names(distances_to_root) <- igraph::V(component_graph)$name
  if (any(!is.finite(distances_to_root))) {
    stop("Internal error: certificate component is disconnected.")
  }

  all_vertices <- component_nodes$vertex
  parent <- stats::setNames(
    rep(NA_character_, length(all_vertices)),
    all_vertices
  )

  for (vertex in setdiff(all_vertices, root_vertex)) {
    neighbors <- igraph::as_ids(igraph::neighbors(
      component_graph,
      vertex,
      mode = "all"
    ))
    possible_parents <- neighbors[
      distances_to_root[neighbors] ==
        distances_to_root[[vertex]] - 1
    ]
    if (!length(possible_parents)) {
      stop("Could not assign a flow-tree parent to ", vertex)
    }
    parent[[vertex]] <- sort(possible_parents)[[1]]
  }

  subtree_area <- stats::setNames(
    component_nodes$FlowAreaHa[match(
      all_vertices,
      component_nodes$vertex
    )],
    all_vertices
  )
  flow_records <- list()
  flow_index <- 0L
  descending_vertices <- all_vertices[
    order(distances_to_root[all_vertices], decreasing = TRUE)
  ]

  for (vertex in descending_vertices) {
    if (vertex == root_vertex) next
    parent_vertex <- parent[[vertex]]
    flow_index <- flow_index + 1L
    flow_records[[flow_index]] <- tibble(
      network_name = network_name,
      child = vertex,
      parent = parent_vertex,
      flow_ha = subtree_area[[vertex]]
    )
    subtree_area[[parent_vertex]] <-
      subtree_area[[parent_vertex]] + subtree_area[[vertex]]
  }

  flow_table <- bind_rows(flow_records) |>
    mutate(
      child_type = str_sub(child, 1, 1),
      parent_type = str_sub(parent, 1, 1),
      child_id = as.integer(str_remove(child, "^[cp]_")),
      parent_id = as.integer(str_remove(parent, "^[cp]_"))
    )

  if (any(
    flow_table$child_type == "p" &
      flow_table$parent_type == "p"
  )) {
    stop("Unexpected patch-to-patch edge in a flow tree.")
  }

  list(
    network_name = network_name,
    component_id = component_id,
    cells = component_nodes |>
      filter(node_type == "cell") |>
      pull(node_id) |>
      as.integer() |>
      sort(),
    patches = component_nodes |>
      filter(node_type == "patch") |>
      pull(node_id) |>
      as.integer() |>
      sort(),
    root_patch = as.integer(root_row$node_id[[1]]),
    collected_area_ha = subtree_area[[root_vertex]],
    flows = flow_table
  )
}

certificate_starts <- map(
  seq_len(nrow(chosen_certificate_components)),
  function(index) {
    build_certificate_start(
      chosen_certificate_components$network_name[[index]],
      chosen_certificate_components$warm_component[[index]]
    )
  }
)
names(certificate_starts) <- map_chr(
  certificate_starts,
  "network_name"
)

start_certificate_cell_members <- map_dfr(
  certificate_starts,
  function(x) tibble(
    network_name = x$network_name,
    coarse_cell_id = x$cells
  )
)
start_certificate_patch_members <- map_dfr(
  certificate_starts,
  function(x) tibble(
    network_name = x$network_name,
    patch_id = x$patches
  )
)
start_certificate_roots <- map_dfr(
  certificate_starts,
  function(x) tibble(
    network_name = x$network_name,
    patch_id = x$root_patch,
    collected_area_ha = x$collected_area_ha
  )
)
start_certificate_flows <- map_dfr(certificate_starts, "flows")

start_flow_structure_audit <- map_dfr(
  certificate_starts,
  function(x) {
    expected_collected <-
      sum(candidate_table$EligibleAreaHa[
        match(x$cells, candidate_table$coarse_cell_id)
      ]) +
      sum(patch_table$ExistingPatchAreaHa[
        match(x$patches, patch_table$patch_id)
      ])

    tibble(
      network_name = x$network_name,
      member_nodes = length(x$cells) + length(x$patches),
      spanning_tree_arcs = nrow(x$flows),
      expected_collected_area_ha = expected_collected,
      start_collected_area_ha = x$collected_area_ha,
      tree_has_n_minus_one_arcs =
        nrow(x$flows) == length(x$cells) + length(x$patches) - 1L,
      collected_area_reconciles =
        abs(x$collected_area_ha - expected_collected) <= 1e-6,
      all_tree_flows_positive = all(x$flows$flow_ha > 0)
    )
  }
)

if (!all(
  start_flow_structure_audit$tree_has_n_minus_one_arcs &
    start_flow_structure_audit$collected_area_reconciles &
    start_flow_structure_audit$all_tree_flows_positive
)) {
  stop("Internal error in certificate flow starts.", call. = FALSE)
}

# Independent membership and area checks before any AMPL files are written.
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
  stop("Internal error in warm-start certificate areas.", call. = FALSE)
}

if (
  nrow(start_certificate_roots) != length(required_networks) ||
  anyDuplicated(start_certificate_roots$network_name)
) {
  stop("Each certificate must have exactly one root.", call. = FALSE)
}

readr::write_csv(
  start_certificate_area_check,
  file.path(problem_dir, "warm_start_certificate_area_check.csv")
)
readr::write_csv(
  start_certificate_roots,
  file.path(problem_dir, "warm_start_certificate_roots.csv")
)
readr::write_csv(
  start_flow_structure_audit,
  file.path(problem_dir, "warm_start_flow_structure_audit.csv")
)

print(knitr::kable(start_certificate_area_check, digits = 2))
knitr::kable(start_flow_structure_audit, digits = 6)
