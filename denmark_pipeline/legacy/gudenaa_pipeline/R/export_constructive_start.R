mip_start_temporary_path <- paste0(mip_start_path, ".building")
# Write the complete zero-shortfall start constructed above. Every binary and
# continuous flow variable is initialized consistently, so this is a genuine
# feasible MIP start rather than a hint.
writeLines(
  c(
    "# Complete feasible start constructed by notebook 08b.",
    "# All five certificate shortfalls are zero.",
    "let {c in Cells} Selected[c] := 0;",
    "let {(l,c) in AllowedActions} LanduseDecision[l,c] := 0;",
    "let {(l,i,j) in SameHabitatEdges} Contiguity[l,i,j] := 0;",
    "let {k in RequiredNetworks} CertificateBuilt[k] := 1;",
    "let {k in RequiredNetworks} NetworkShortfallHa[k] := 0;",
    paste0(
      "let {(k,c) in EligibleNetworkCells} ",
      "NetworkCell[k,c] := 0;"
    ),
    paste0(
      "let {(k,p) in EligibleNetworkPatches} ",
      "NetworkPatch[k,p] := 0;"
    ),
    paste0(
      "let {(k,p) in EligibleNetworkPatches} ",
      "NetworkRoot[k,p] := 0;"
    ),
    paste0(
      "let {(k,i,j) in NetworkCellArcs} ",
      "AreaFlowHa[k,i,j] := 0;"
    ),
    paste0(
      "let {(k,c,p) in NetworkCellPatchEdges} ",
      "CellToPatchFlowHa[k,c,p] := 0;"
    ),
    paste0(
      "let {(k,c,p) in NetworkCellPatchEdges} ",
      "PatchToCellFlowHa[k,c,p] := 0;"
    ),
    paste0(
      "let {(k,p) in EligibleNetworkPatches} ",
      "CollectedAreaHa[k,p] := 0;"
    )
  ),
  mip_start_temporary_path
)

append_start(
  mip_start_temporary_path,
  paste0(
    "let Selected[",
    constructive_assignment$coarse_cell_id,
    "] := 1;"
  )
)
append_start(
  mip_start_temporary_path,
  constructive_assignment |>
    transmute(line = paste0(
      "let LanduseDecision['",
      landuse,
      "',",
      coarse_cell_id,
      "] := 1;"
    )) |>
    pull(line)
)

repair_start_labels_i <- constructive_assignment |>
  transmute(landuse, i = coarse_cell_id)
repair_start_labels_j <- constructive_assignment |>
  transmute(landuse, j = coarse_cell_id)
repair_start_contiguity <- same_habitat_edges |>
  semi_join(repair_start_labels_i, by = c("landuse", "i")) |>
  semi_join(repair_start_labels_j, by = c("landuse", "j"))
append_start(
  mip_start_temporary_path,
  paste0(
    "let Contiguity['",
    repair_start_contiguity$landuse,
    "',",
    repair_start_contiguity$i,
    ",",
    repair_start_contiguity$j,
    "] := 1;"
  )
)

append_start(
  mip_start_temporary_path,
  start_certificate_cell_members |>
    transmute(line = paste0(
      "let NetworkCell['",
      network_name,
      "',",
      coarse_cell_id,
      "] := 1;"
    )) |>
    pull(line)
)
append_start(
  mip_start_temporary_path,
  start_certificate_patch_members |>
    transmute(line = paste0(
      "let NetworkPatch['",
      network_name,
      "',",
      patch_id,
      "] := 1;"
    )) |>
    pull(line)
)
append_start(
  mip_start_temporary_path,
  start_certificate_roots |>
    transmute(line = paste0(
      "let NetworkRoot['",
      network_name,
      "',",
      patch_id,
      "] := 1;"
    )) |>
    pull(line)
)

flow_start_lines <- start_certificate_flows |>
  mutate(line = case_when(
    child_type == "c" & parent_type == "c" ~ paste0(
      "let AreaFlowHa['",
      network_name,
      "',",
      child_id,
      ",",
      parent_id,
      "] := ",
      ampl_number(flow_ha),
      ";"
    ),
    child_type == "c" & parent_type == "p" ~ paste0(
      "let CellToPatchFlowHa['",
      network_name,
      "',",
      child_id,
      ",",
      parent_id,
      "] := ",
      ampl_number(flow_ha),
      ";"
    ),
    child_type == "p" & parent_type == "c" ~ paste0(
      "let PatchToCellFlowHa['",
      network_name,
      "',",
      parent_id,
      ",",
      child_id,
      "] := ",
      ampl_number(flow_ha),
      ";"
    ),
    TRUE ~ NA_character_
  ))
stopifnot(!anyNA(flow_start_lines$line))
append_start(mip_start_temporary_path, flow_start_lines$line)

append_start(
  mip_start_temporary_path,
  start_certificate_roots |>
    transmute(line = paste0(
      "let CollectedAreaHa['",
      network_name,
      "',",
      patch_id,
      "] := ",
      ampl_number(collected_area_ha),
      ";"
    )) |>
    pull(line)
)

commit_generated_file(mip_start_temporary_path, mip_start_path)
