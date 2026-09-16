# Gudenaa_30pct_feasible_network_repair_200m.mod
#
# Lexicographic repair model for a successful 200 m warm-up footprint.
#
# Stage 1 minimizes normalized certificate shortfall. A deliberately complete
# initial solution is always available: the warm footprint, no active network
# certificates, and full shortfall for all five requirements.
#
# Stage 2 is entered only after every shortfall reaches zero. It minimizes the
# area-weighted symmetric difference from the warm footprint while preserving
# the 30%-by-subarea bounds and all five connected certificates.
#
# Stage 3 fixes that minimum repair (within a small numerical tolerance) and
# maximizes the same biodiversity and 1.2-weighted spatial-contiguity objective
# used by notebook 07. Thus connectivity is hard, footprint repair is second,
# and ecological/spatial quality breaks ties among equally small repairs.

# ---- Core sets -----------------------------------------------------------

set Cells ordered;
set Landuses;
set ForestLanduses within Landuses;
set Subareas ordered;
set RequiredNetworks ordered;
set GlobalNetworks within RequiredNetworks;

set AllowedActions within {Landuses, Cells};

# Undirected candidate rook edges, supplied once as (i,j).
set E within {Cells, Cells};
set SameHabitatEdges within {Landuses, Cells, Cells};

set ExistingPatches ordered;
set CellPatchEdges within {Cells, ExistingPatches};

# Network variables are created only inside potential components that can
# possibly meet the corresponding area target. This is an exact reduction:
# a connected certificate cannot use two disconnected potential components.
set EligibleNetworkCells within {RequiredNetworks, Cells};
set EligibleNetworkPatches within {RequiredNetworks, ExistingPatches};
set NetworkCellArcs within {RequiredNetworks, Cells, Cells};
set NetworkCellPatchEdges within {
    RequiredNetworks, Cells, ExistingPatches
};

# ---- Area, warm-footprint, and objective data ----------------------------

param BioDiversity {(l,c) in AllowedActions};
param ExistingNatureContact {(l,c) in AllowedActions} >= 0 default 0;

param RestorableAreaHa {Cells} > 0;
param RestorableAreaHaBySubarea {Subareas, Cells} >= 0 default 0;

param SubareaAreaHa {Subareas} > 0;
param CurrentNatureAreaHa {Subareas} >= 0;
param FinalNatureTargetHa {Subareas} >= 0;
param RestorationRequirementHa {Subareas} >= 0;
param RestorationToleranceHa {Subareas} >= 0 default 4;

param ExistingPatchAreaHa {ExistingPatches} > 0;

# Global uses the complete national area of each reached existing patch.
# Subarea certificates use only the patch fragment physically inside that
# subarea. Restoration always counts only inside the planning polygon.
param NetworkCellAreaHa {RequiredNetworks, Cells} >= 0 default 0;
param NetworkPatchAreaHa {
    RequiredNetworks, ExistingPatches
} >= 0 default 0;
param NetworkMinimumHa {RequiredNetworks} > 0;

param ComponentMaxAreaCellHa {Cells} > 0;
param ComponentMaxAreaPatchHa {ExistingPatches} > 0;

param IsProductionForestCell {Cells} binary default 0;
param WarmSelected {Cells} binary default 0;

param BiodiversityNormalizer > 0 default 1;
param ContiguityNormalizer > 0 default 1;
param ExistingNatureContactNormalizer > 0 default 1;
param PatchConnectionNormalizer > 0 default 1;

param SpatialContiguityBonus >= 0 default 1;
param NewToNewAdjacencyWeight >= 0 default 0.5;
param ExistingToNewAdjacencyWeight >= 0 default 0.5;
param ConnectedExistingPatchBonus >= 0 default 0.25;

# Initially set to the total candidate area, so it does not restrict stages 1
# or 2. The run file replaces it with the stage-2 incumbent plus the stated
# numerical tolerance before stage 3.
param FootprintChangeLimitHa >= 0;
param RepairChangeToleranceHa >= 0 default 0.01;

# ---- Restoration and habitat decisions ----------------------------------

var LanduseDecision {(l,c) in AllowedActions} binary;
var Selected {Cells} binary;
var Contiguity {(l,i,j) in SameHabitatEdges} binary;

# ---- Connected-certificate decisions ------------------------------------

var CertificateBuilt {RequiredNetworks} binary;
var NetworkShortfallHa {k in RequiredNetworks}
    >= 0, <= NetworkMinimumHa[k];

var NetworkCell {(k,c) in EligibleNetworkCells} binary;
var NetworkPatch {(k,p) in EligibleNetworkPatches} binary;
var NetworkRoot {(k,p) in EligibleNetworkPatches} binary;

var AreaFlowHa {(k,i,j) in NetworkCellArcs} >= 0;
var CellToPatchFlowHa {(k,c,p) in NetworkCellPatchEdges} >= 0;
var PatchToCellFlowHa {(k,c,p) in NetworkCellPatchEdges} >= 0;
var CollectedAreaHa {(k,p) in EligibleNetworkPatches} >= 0;

# ---- Lexicographic objectives -------------------------------------------

minimize RelativeNetworkShortfall:
    sum {k in RequiredNetworks}
        NetworkShortfallHa[k] / NetworkMinimumHa[k];

minimize WarmFootprintChange:
    sum {c in Cells} RestorableAreaHa[c] *
      (
        WarmSelected[c] * (1 - Selected[c])
        + (1 - WarmSelected[c]) * Selected[c]
      );

maximize RepairQuality:
    (sum {(l,c) in AllowedActions}
        BioDiversity[l,c] * LanduseDecision[l,c])
        / BiodiversityNormalizer
  + SpatialContiguityBonus *
    (
      NewToNewAdjacencyWeight *
        (sum {(l,i,j) in SameHabitatEdges} Contiguity[l,i,j])
          / ContiguityNormalizer
      + ExistingToNewAdjacencyWeight *
        (sum {(l,c) in AllowedActions}
          ExistingNatureContact[l,c] * LanduseDecision[l,c])
          / ExistingNatureContactNormalizer
    )
  + ConnectedExistingPatchBonus *
      (
        sum {(k,p) in EligibleNetworkPatches: k in GlobalNetworks}
            NetworkPatch[k,p]
        - sum {k in GlobalNetworks} CertificateBuilt[k]
      ) / PatchConnectionNormalizer;

# ---- Restoration and 30%-by-subarea policy ------------------------------

subject to LinkSelectionToLanduse {c in Cells}:
    Selected[c]
    = sum {l in Landuses: (l,c) in AllowedActions}
        LanduseDecision[l,c];

subject to SubareaRestorationLower {s in Subareas}:
    sum {c in Cells}
        RestorableAreaHaBySubarea[s,c] * Selected[c]
    >= RestorationRequirementHa[s];

subject to SubareaRestorationUpper {s in Subareas}:
    sum {c in Cells}
        RestorableAreaHaBySubarea[s,c] * Selected[c]
    <= RestorationRequirementHa[s] + RestorationToleranceHa[s];

subject to LimitWarmFootprintChange:
    sum {c in Cells} RestorableAreaHa[c] *
      (
        WarmSelected[c] * (1 - Selected[c])
        + (1 - WarmSelected[c]) * Selected[c]
      )
    <= FootprintChangeLimitHa;

subject to DefineContiguityFromI {(l,i,j) in SameHabitatEdges}:
    Contiguity[l,i,j] <= LanduseDecision[l,i];
subject to DefineContiguityFromJ {(l,i,j) in SameHabitatEdges}:
    Contiguity[l,i,j] <= LanduseDecision[l,j];
subject to DefineContiguityLower {(l,i,j) in SameHabitatEdges}:
    Contiguity[l,i,j]
    >= LanduseDecision[l,i] + LanduseDecision[l,j] - 1;

# ---- Certificate membership, roots, and area ----------------------------

subject to CertificateUsesSelectedCell {
    (k,c) in EligibleNetworkCells
}:
    NetworkCell[k,c] <= Selected[c];

subject to CellRequiresBuiltCertificate {
    (k,c) in EligibleNetworkCells
}:
    NetworkCell[k,c] <= CertificateBuilt[k];

subject to PatchRequiresBuiltCertificate {
    (k,p) in EligibleNetworkPatches
}:
    NetworkPatch[k,p] <= CertificateBuilt[k];

subject to CertificateContainsNewRestoration {k in RequiredNetworks}:
    sum {(q,c) in EligibleNetworkCells: q = k}
        NetworkCell[q,c]
    >= CertificateBuilt[k];

subject to PatchNeedsMemberContact {
    (k,p) in EligibleNetworkPatches
}:
    NetworkPatch[k,p]
    <= sum {(q,c,r) in NetworkCellPatchEdges:
              q = k and r = p}
           NetworkCell[q,c];

subject to MemberContactIncludesPatch {
    (k,c,p) in NetworkCellPatchEdges
}:
    NetworkPatch[k,p] >= NetworkCell[k,c];

subject to ExactlyOneExistingNatureRoot {k in RequiredNetworks}:
    sum {(q,p) in EligibleNetworkPatches: q = k}
        NetworkRoot[q,p]
    = CertificateBuilt[k];

subject to RootMustBeMemberPatch {
    (k,p) in EligibleNetworkPatches
}:
    NetworkRoot[k,p] <= NetworkPatch[k,p];

subject to RequiredNetworkArea {k in RequiredNetworks}:
    sum {(q,c) in EligibleNetworkCells: q = k}
        NetworkCellAreaHa[q,c] * NetworkCell[q,c]
    + sum {(q,p) in EligibleNetworkPatches: q = k}
        NetworkPatchAreaHa[q,p] * NetworkPatch[q,p]
    + NetworkShortfallHa[k]
    >= NetworkMinimumHa[k];

# ---- Single-commodity flow connectivity ---------------------------------

subject to CellArcNeedsMemberOrigin {
    (k,i,j) in NetworkCellArcs
}:
    AreaFlowHa[k,i,j]
    <= ComponentMaxAreaCellHa[i] * NetworkCell[k,i];

subject to CellArcNeedsMemberDestination {
    (k,i,j) in NetworkCellArcs
}:
    AreaFlowHa[k,i,j]
    <= ComponentMaxAreaCellHa[j] * NetworkCell[k,j];

subject to CellToPatchNeedsMemberCell {
    (k,c,p) in NetworkCellPatchEdges
}:
    CellToPatchFlowHa[k,c,p]
    <= ComponentMaxAreaCellHa[c] * NetworkCell[k,c];

subject to CellToPatchNeedsMemberPatch {
    (k,c,p) in NetworkCellPatchEdges
}:
    CellToPatchFlowHa[k,c,p]
    <= ComponentMaxAreaPatchHa[p] * NetworkPatch[k,p];

subject to PatchToCellNeedsMemberCell {
    (k,c,p) in NetworkCellPatchEdges
}:
    PatchToCellFlowHa[k,c,p]
    <= ComponentMaxAreaCellHa[c] * NetworkCell[k,c];

subject to PatchToCellNeedsMemberPatch {
    (k,c,p) in NetworkCellPatchEdges
}:
    PatchToCellFlowHa[k,c,p]
    <= ComponentMaxAreaPatchHa[p] * NetworkPatch[k,p];

subject to CollectOnlyAtRoot {
    (k,p) in EligibleNetworkPatches
}:
    CollectedAreaHa[k,p]
    <= ComponentMaxAreaPatchHa[p] * NetworkRoot[k,p];

# Flow balances use complete node areas so every positive-area member has a
# path to its one existing-nature root. Certificate thresholds may use a
# smaller subarea-specific coefficient without weakening connectivity.
subject to CandidateCellAreaBalance {
    (k,c) in EligibleNetworkCells
}:
    RestorableAreaHa[c] * NetworkCell[k,c]
    + sum {(q,i,j) in NetworkCellArcs:
              q = k and j = c}
        AreaFlowHa[q,i,j]
    + sum {(q,d,p) in NetworkCellPatchEdges:
              q = k and d = c}
        PatchToCellFlowHa[q,d,p]
    = sum {(q,i,j) in NetworkCellArcs:
              q = k and i = c}
        AreaFlowHa[q,i,j]
    + sum {(q,d,p) in NetworkCellPatchEdges:
              q = k and d = c}
        CellToPatchFlowHa[q,d,p];

subject to ExistingPatchAreaBalance {
    (k,p) in EligibleNetworkPatches
}:
    ExistingPatchAreaHa[p] * NetworkPatch[k,p]
    + sum {(q,c,r) in NetworkCellPatchEdges:
              q = k and r = p}
        CellToPatchFlowHa[q,c,r]
    = sum {(q,c,r) in NetworkCellPatchEdges:
              q = k and r = p}
        PatchToCellFlowHa[q,c,r]
    + CollectedAreaHa[k,p];

# ---- Defensive checks ----------------------------------------------------

check: card(Cells) > 0;
check: card(Landuses) = 8;
check: card(ForestLanduses) = 4;
check: card(Subareas) = 4;
check: card(RequiredNetworks) = 5;
check: card(GlobalNetworks) = 1;
check: card(ExistingPatches) > 0;
check: card(CellPatchEdges) > 0;
check:
    abs(NewToNewAdjacencyWeight + ExistingToNewAdjacencyWeight - 1)
    <= 0.000000001;

check {c in Cells}:
    sum {l in Landuses: (l,c) in AllowedActions} 1 >= 1;

check {c in Cells}:
    abs(
        RestorableAreaHa[c]
        - sum {s in Subareas} RestorableAreaHaBySubarea[s,c]
    ) <= 0.000001;

check {k in RequiredNetworks}:
    sum {(q,c) in EligibleNetworkCells: q = k} 1 >= 1;

check {k in RequiredNetworks}:
    sum {(q,p) in EligibleNetworkPatches: q = k} 1 >= 1;
