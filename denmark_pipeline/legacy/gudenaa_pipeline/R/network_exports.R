gp_network_export_lines <- function(e) {
  eval(quote(c(
  paste0(
    'printf "coarse_cell_id\\tlanduse\\trestored_ha',
    '\\tbiodiversity\\texisting_nature_contact_edges\\n" > "',
    ampl_path(solution_path), '";'
  ),
  paste0(
    'printf {(l,c) in AllowedActions: ',
    'LanduseDecision[l,c] > 0.5} ',
    '"%d\\t%s\\t%.8f\\t%.8f\\t%.8f\\n", ',
    'c, l, RestorableAreaHa[c], BioDiversity[l,c], ',
    'ExistingNatureContact[l,c] >> "',
    ampl_path(solution_path), '";'
  ),
  "",
  paste0(
    'printf "subarea_id\\tcoarse_cell_id\\trestored_ha\\n" > "',
    ampl_path(selected_subarea_path), '";'
  ),
  paste0(
    'printf {s in Subareas, c in Cells: Selected[c] > 0.5 ',
    'and RestorableAreaHaBySubarea[s,c] > 0} ',
    '"%d\\t%d\\t%.8f\\n", s, c, ',
    'RestorableAreaHaBySubarea[s,c] >> "',
    ampl_path(selected_subarea_path), '";'
  ),
  "",
  paste0(
    'printf "subarea_id\\tpolygon_ha\\treference_30pct_target_ha',
    '\\tcurrent_ha\\treference_30pct_restoration_ha',
    '\\trestored_ha\\tfinal_nature_ha',
    '\\tfinal_nature_percent\\n" > "',
    ampl_path(subarea_solution_path), '";'
  ),
  paste0(
    'printf {s in Subareas} ',
    '"%d\\t%.8f\\t%.8f\\t%.8f\\t%.8f\\t%.8f\\t%.8f',
    '\\t%.8f\\n", ',
    's, SubareaAreaHa[s], FinalNatureTargetHa[s], ',
    'CurrentNatureAreaHa[s], RestorationRequirementHa[s], ',
    'sum {c in Cells} RestorableAreaHaBySubarea[s,c] * Selected[c], ',
    'CurrentNatureAreaHa[s] + sum {c in Cells} ',
    'RestorableAreaHaBySubarea[s,c] * Selected[c], ',
    '100 * (CurrentNatureAreaHa[s] + sum {c in Cells} ',
    'RestorableAreaHaBySubarea[s,c] * Selected[c]) / ',
    'SubareaAreaHa[s] >> "', ampl_path(subarea_solution_path), '";'
  ),
  "",
  paste0(
    'printf "network_name\\tcoarse_cell_id\\tcounted_area_ha',
    '\\ttotal_restorable_ha\\n" > "',
    ampl_path(network_cell_solution_path), '";'
  ),
  paste0(
    'printf {(k,c) in EligibleNetworkCells: ',
    'NetworkCell[k,c] > 0.5} ',
    '"%s\\t%d\\t%.8f\\t%.8f\\n", ',
    'k, c, NetworkCellAreaHa[k,c], RestorableAreaHa[c] >> "',
    ampl_path(network_cell_solution_path), '";'
  ),
  "",
  paste0(
    'printf "network_name\\tpatch_id\\tcounted_area_ha',
    '\\ttotal_patch_ha\\n" > "',
    ampl_path(network_patch_solution_path), '";'
  ),
  paste0(
    'printf {(k,p) in EligibleNetworkPatches: ',
    'NetworkPatch[k,p] > 0.5} ',
    '"%s\\t%d\\t%.8f\\t%.8f\\n", ',
    'k, p, NetworkPatchAreaHa[k,p], ExistingPatchAreaHa[p] >> "',
    ampl_path(network_patch_solution_path), '";'
  ),
  "",
  paste0(
    'printf "network_name\\tpatch_id\\tcollected_area_ha\\n" > "',
    ampl_path(network_root_solution_path), '";'
  ),
  paste0(
    'printf {(k,p) in EligibleNetworkPatches: ',
    'NetworkRoot[k,p] > 0.5} ',
    '"%s\\t%d\\t%.8f\\n", k, p, CollectedAreaHa[k,p] >> "',
    ampl_path(network_root_solution_path), '";'
  ),
  "",
  paste0(
    'printf "network_name\\tminimum_ha\\tcertified_area_ha',
    '\\tmember_cells\\tmember_patches\\troot_count',
    '\\tcollected_total_area_ha\\tshortfall_ha',
    '\\tcertificate_built\\n" > "',
    ampl_path(network_solution_summary_path), '";'
  ),
  paste0(
    'printf {k in RequiredNetworks} ',
    '"%s\\t%.8f\\t%.8f\\t%d\\t%d\\t%d\\t%.8f',
    '\\t%.8f\\t%d\\n", ',
    'k, NetworkMinimumHa[k], ',
    'sum {(q,c) in EligibleNetworkCells: q = k} ',
    'NetworkCellAreaHa[q,c] * NetworkCell[q,c] + ',
    'sum {(q,p) in EligibleNetworkPatches: q = k} ',
    'NetworkPatchAreaHa[q,p] * NetworkPatch[q,p], ',
    'sum {(q,c) in EligibleNetworkCells: q = k} NetworkCell[q,c], ',
    'sum {(q,p) in EligibleNetworkPatches: q = k} NetworkPatch[q,p], ',
    'sum {(q,p) in EligibleNetworkPatches: q = k} NetworkRoot[q,p], ',
    'sum {(q,p) in EligibleNetworkPatches: q = k} ',
    'CollectedAreaHa[q,p], NetworkShortfallHa[k], ',
    'CertificateBuilt[k] >> "',
    ampl_path(network_solution_summary_path), '";'
  ),
  "",
  paste0(
    'printf "item\\tvalue\\n" > "',
    ampl_path(summary_path), '";'
  ),
  paste0(
    'printf "solve_result\\t%s\\n", solve_result >> "',
    ampl_path(summary_path), '";'
  ),
  paste0(
    'printf "solve_result_num\\t%d\\n", solve_result_num >> "',
    ampl_path(summary_path), '";'
  ),
  paste0(
    'printf "repair_status\\tnetwork_incumbent\\n" >> "',
    ampl_path(summary_path), '";'
  ),
  paste0(
    'printf "objective\\t%.10f\\n", RepairQuality >> "',
    ampl_path(summary_path), '";'
  ),
  paste0(
    'printf "warm_footprint_change_ha\\t%.8f\\n", ',
    'sum {c in Cells} RestorableAreaHa[c] * ',
    '(WarmSelected[c] * (1 - Selected[c]) + ',
    '(1 - WarmSelected[c]) * Selected[c]) >> "',
    ampl_path(summary_path), '";'
  ),
  paste0(
    'printf "network_shortfall_total_ha\\t%.8f\\n", ',
    'sum {k in RequiredNetworks} NetworkShortfallHa[k] >> "',
    ampl_path(summary_path), '";'
  ),
  paste0(
    'printf "spatial_contiguity_bonus\\t%.8f\\n", ',
    'SpatialContiguityBonus >> "', ampl_path(summary_path), '";'
  ),
  paste0(
    'printf "normalized_biodiversity\\t%.10f\\n", ',
    '(sum {(l,c) in AllowedActions} BioDiversity[l,c] * ',
    'LanduseDecision[l,c]) / BiodiversityNormalizer >> "',
    ampl_path(summary_path), '";'
  ),
  paste0(
    'printf "normalized_new_to_new_adjacency\\t%.10f\\n", ',
    '(sum {(l,i,j) in SameHabitatEdges} Contiguity[l,i,j]) / ',
    'ContiguityNormalizer >> "', ampl_path(summary_path), '";'
  ),
  paste0(
    'printf "normalized_existing_to_new_adjacency\\t%.10f\\n", ',
    '(sum {(l,c) in AllowedActions} ExistingNatureContact[l,c] * ',
    'LanduseDecision[l,c]) / ExistingNatureContactNormalizer >> "',
    ampl_path(summary_path), '";'
  ),
  paste0(
    'printf "global_existing_patches_joined_beyond_root\\t%.0f\\n", ',
    'sum {(k,p) in EligibleNetworkPatches: k in GlobalNetworks} ',
    'NetworkPatch[k,p] - sum {k in GlobalNetworks} ',
    'CertificateBuilt[k] >> "', ampl_path(summary_path), '";'
  ),
  paste0(
    'printf "normalized_existing_patch_connection\\t%.10f\\n", ',
    '(sum {(k,p) in EligibleNetworkPatches: k in GlobalNetworks} ',
    'NetworkPatch[k,p] - sum {k in GlobalNetworks} ',
    'CertificateBuilt[k]) / PatchConnectionNormalizer >> "',
    ampl_path(summary_path), '";'
  ),
  paste0(
    'printf "selected_cells\\t%.0f\\n", ',
    'sum {c in Cells} Selected[c] >> "',
    ampl_path(summary_path), '";'
  ),
  paste0(
    'printf "restored_ha\\t%.8f\\n", ',
    'sum {c in Cells} RestorableAreaHa[c] * Selected[c] >> "',
    ampl_path(summary_path), '";'
  ),
  paste0(
    'printf "global_polygon_area_ha\\t%.8f\\n", ',
    'GlobalPolygonAreaHa >> "', ampl_path(summary_path), '";'
  ),
  paste0(
    'printf "global_current_nature_ha\\t%.8f\\n", ',
    'GlobalCurrentNatureAreaHa >> "', ampl_path(summary_path), '";'
  ),
  paste0(
    'printf "global_30pct_target_ha\\t%.8f\\n", ',
    'GlobalFinalNatureTargetHa >> "', ampl_path(summary_path), '";'
  ),
  paste0(
    'printf "global_required_restoration_ha\\t%.8f\\n", ',
    'GlobalRestorationRequirementHa >> "', ampl_path(summary_path), '";'
  ),
  paste0(
    'printf "global_restoration_tolerance_ha\\t%.8f\\n", ',
    'GlobalRestorationToleranceHa >> "', ampl_path(summary_path), '";'
  ),
  paste0(
    'printf "global_final_nature_ha\\t%.8f\\n", ',
    'GlobalCurrentNatureAreaHa + sum {c in Cells} ',
    'RestorableAreaHa[c] * Selected[c] >> "',
    ampl_path(summary_path), '";'
  ),
  paste0(
    'printf "global_final_nature_percent\\t%.8f\\n", ',
    '100 * (GlobalCurrentNatureAreaHa + sum {c in Cells} ',
    'RestorableAreaHa[c] * Selected[c]) / GlobalPolygonAreaHa >> "',
    ampl_path(summary_path), '";'
  )
  )), envir = e)
}
