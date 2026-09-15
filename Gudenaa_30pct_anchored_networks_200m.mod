# Gudenaa_30pct_anchored_networks_200m.mod
#
# Second-stage 200 m model for the four-part Gudenaa planning polygon.
# It retains the exact 30%-by-subarea accounting from the warm-up and adds
# five connected-network certificates:
#   * Global: at least 5,000 ha of current + restored nature;
#   * Subarea1..Subarea4: at least 1,000 ha inside the corresponding subarea.
#
# Every certificate contains new restoration, is connected through 200 m rook
# edges and existing-nature patch supernodes, and has exactly one root in an
# existing-nature patch. Certificates may overlap: one physical component can
# legitimately satisfy several "at least one area" requirements.

# ---- Sets ----------------------------------------------------------------

set Cells ordered;
set Landuses;
set ForestLanduses within Landuses;
set Subareas ordered;
set RequiredNetworks ordered;
set GlobalNetworks within RequiredNetworks;

set AllowedActions within {Landuses, Cells};

# Undirected candidate-to-candidate rook edges, supplied once as (i,j).
set E within {Cells, Cells};
set NetworkArcs within {Cells, Cells} :=
    E union setof {(i,j) in E} (j,i);

# Optional triples for the same-habitat new-to-new reward.
set SameHabitatEdges within {Landuses, Cells, Cells};

# Each current-nature rook patch is represented by one supernode. A contact
# pair (c,p) is present only when eligible 10 m land in candidate cell c is
# rook-adjacent to current-nature patch p.
set ExistingPatches ordered;
set CellPatchEdges within {Cells, ExistingPatches};

# ---- Area and objective data ---------------------------------------------

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

# These coefficients determine what area counts toward each certificate.
# Global uses total area; a Subarea certificate uses only area inside its own
# subarea, while its connecting path may cross a boundary.
param NetworkCellAreaHa {RequiredNetworks, Cells} >= 0 default 0;
param NetworkPatchAreaHa {RequiredNetworks, ExistingPatches} >= 0 default 0;
param NetworkMinimumHa {RequiredNetworks} > 0;

# Component-specific valid upper bounds for flow variables.
param ComponentMaxAreaCellHa {Cells} > 0;
param ComponentMaxAreaPatchHa {ExistingPatches} > 0;

# Coarse policy flag. The R notebook uses the established 50%-of-eligible-area
# rule to classify production-forest-dominated 200 m units.
param IsProductionForestCell {Cells} binary default 0;

param BiodiversityNormalizer > 0 default 1;
param ContiguityNormalizer > 0 default 1;
param ExistingNatureContactNormalizer > 0 default 1;
param PatchConnectionNormalizer > 0 default 1;

param SpatialContiguityBonus >= 0 default 1;
param NewToNewAdjacencyWeight >= 0 default 0.5;
param ExistingToNewAdjacencyWeight >= 0 default 0.5;
param ConnectedExistingPatchBonus >= 0 default 0.25;

# ---- Core restoration decisions -----------------------------------------

var LanduseDecision {(l,c) in AllowedActions} binary;
var Selected {Cells} binary;
var Contiguity {(l,i,j) in SameHabitatEdges} binary;

# ---- Connected-network certificates -------------------------------------

var NetworkCell {RequiredNetworks, Cells} binary;
var NetworkPatch {RequiredNetworks, ExistingPatches} binary;
var NetworkRoot {RequiredNetworks, ExistingPatches} binary;

var AreaFlowHa {
    k in RequiredNetworks, (i,j) in NetworkArcs
} >= 0;
var CellToPatchFlowHa {
    k in RequiredNetworks, (c,p) in CellPatchEdges
} >= 0;
var PatchToCellFlowHa {
    k in RequiredNetworks, (c,p) in CellPatchEdges
} >= 0;
var CollectedAreaHa {
    k in RequiredNetworks, p in ExistingPatches
} >= 0;

# The first two terms preserve the warm-up objective. The final normalized
# term rewards a Global certificate that joins multiple distinct existing
# nature patches; subtracting its one mandatory root makes a one-patch network
# score zero.
maximize ConservationIndex:
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
      (sum {k in GlobalNetworks, p in ExistingPatches}
          NetworkPatch[k,p] - 1)
        / PatchConnectionNormalizer;

subject to LinkSelectionToLanduse {c in Cells}:
    Selected[c]
    = sum {l in Landuses: (l,c) in AllowedActions}
        LanduseDecision[l,c];

# Current nature plus selected restoration reaches approximately 30% in every
# subarea. This is stronger and clearer than only enforcing 30% in aggregate,
# and automatically gives 30% for the complete polygon.
subject to SubareaRestorationLower {s in Subareas}:
    sum {c in Cells}
        RestorableAreaHaBySubarea[s,c] * Selected[c]
    >= RestorationRequirementHa[s];

subject to SubareaRestorationUpper {s in Subareas}:
    sum {c in Cells}
        RestorableAreaHaBySubarea[s,c] * Selected[c]
    <= RestorationRequirementHa[s] + RestorationToleranceHa[s];

subject to DefineContiguityFromI {(l,i,j) in SameHabitatEdges}:
    Contiguity[l,i,j] <= LanduseDecision[l,i];
subject to DefineContiguityFromJ {(l,i,j) in SameHabitatEdges}:
    Contiguity[l,i,j] <= LanduseDecision[l,j];
subject to DefineContiguityLower {(l,i,j) in SameHabitatEdges}:
    Contiguity[l,i,j]
    >= LanduseDecision[l,i] + LanduseDecision[l,j] - 1;

# ---- Membership, anchoring, and minimum areas ----------------------------

subject to CertificateUsesSelectedCell {k in RequiredNetworks, c in Cells}:
    NetworkCell[k,c] <= Selected[c];

subject to CertificateContainsNewRestoration {k in RequiredNetworks}:
    sum {c in Cells} NetworkCell[k,c] >= 1;

# An included patch must touch a member cell. Conversely, a member cell that
# touches a patch includes that patch in the same certificate.
subject to PatchNeedsMemberContact {
    k in RequiredNetworks, p in ExistingPatches
}:
    NetworkPatch[k,p]
    <= sum {(c,q) in CellPatchEdges: q = p} NetworkCell[k,c];

subject to MemberContactIncludesPatch {
    k in RequiredNetworks, (c,p) in CellPatchEdges
}:
    NetworkPatch[k,p] >= NetworkCell[k,c];

# Exactly one existing-nature root proves that each required area originates
# from existing nature. Other existing patches can be joined to that root.
subject to ExactlyOneExistingNatureRoot {k in RequiredNetworks}:
    sum {p in ExistingPatches} NetworkRoot[k,p] = 1;

subject to RootMustBeMemberPatch {
    k in RequiredNetworks, p in ExistingPatches
}:
    NetworkRoot[k,p] <= NetworkPatch[k,p];

subject to RequiredNetworkArea {k in RequiredNetworks}:
    sum {c in Cells}
        NetworkCellAreaHa[k,c] * NetworkCell[k,c]
    + sum {p in ExistingPatches}
        NetworkPatchAreaHa[k,p] * NetworkPatch[k,p]
    >= NetworkMinimumHa[k];

# ---- Single-commodity flow connectivity ---------------------------------

subject to CellArcNeedsMemberOrigin {
    k in RequiredNetworks, (i,j) in NetworkArcs
}:
    AreaFlowHa[k,i,j]
    <= ComponentMaxAreaCellHa[i] * NetworkCell[k,i];

subject to CellArcNeedsMemberDestination {
    k in RequiredNetworks, (i,j) in NetworkArcs
}:
    AreaFlowHa[k,i,j]
    <= ComponentMaxAreaCellHa[j] * NetworkCell[k,j];

subject to CellToPatchNeedsMemberCell {
    k in RequiredNetworks, (c,p) in CellPatchEdges
}:
    CellToPatchFlowHa[k,c,p]
    <= ComponentMaxAreaCellHa[c] * NetworkCell[k,c];

subject to CellToPatchNeedsMemberPatch {
    k in RequiredNetworks, (c,p) in CellPatchEdges
}:
    CellToPatchFlowHa[k,c,p]
    <= ComponentMaxAreaPatchHa[p] * NetworkPatch[k,p];

subject to PatchToCellNeedsMemberCell {
    k in RequiredNetworks, (c,p) in CellPatchEdges
}:
    PatchToCellFlowHa[k,c,p]
    <= ComponentMaxAreaCellHa[c] * NetworkCell[k,c];

subject to PatchToCellNeedsMemberPatch {
    k in RequiredNetworks, (c,p) in CellPatchEdges
}:
    PatchToCellFlowHa[k,c,p]
    <= ComponentMaxAreaPatchHa[p] * NetworkPatch[k,p];

subject to CollectOnlyAtRoot {
    k in RequiredNetworks, p in ExistingPatches
}:
    CollectedAreaHa[k,p]
    <= ComponentMaxAreaPatchHa[p] * NetworkRoot[k,p];

# Flow uses total node area, not the subarea-specific threshold coefficient.
# Thus every positive-area certificate node must have a path to its root.
subject to CandidateCellAreaBalance {k in RequiredNetworks, c in Cells}:
    RestorableAreaHa[c] * NetworkCell[k,c]
    + sum {(i,j) in NetworkArcs: j = c} AreaFlowHa[k,i,j]
    + sum {(d,p) in CellPatchEdges: d = c}
        PatchToCellFlowHa[k,d,p]
    = sum {(i,j) in NetworkArcs: i = c} AreaFlowHa[k,i,j]
    + sum {(d,p) in CellPatchEdges: d = c}
        CellToPatchFlowHa[k,d,p];

subject to ExistingPatchAreaBalance {
    k in RequiredNetworks, p in ExistingPatches
}:
    ExistingPatchAreaHa[p] * NetworkPatch[k,p]
    + sum {(c,q) in CellPatchEdges: q = p}
        CellToPatchFlowHa[k,c,q]
    = sum {(c,q) in CellPatchEdges: q = p}
        PatchToCellFlowHa[k,c,q]
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

check {s in Subareas}:
    RestorationRequirementHa[s]
    <= sum {c in Cells} RestorableAreaHaBySubarea[s,c]
       + 0.000001;

check {(l,c) in AllowedActions: IsProductionForestCell[c] = 1}:
    l in ForestLanduses;

check {(i,j) in E}: i <> j;
check {(l,i,j) in SameHabitatEdges}: (i,j) in E;
check {(l,i,j) in SameHabitatEdges}: (l,i) in AllowedActions;
check {(l,i,j) in SameHabitatEdges}: (l,j) in AllowedActions;

check {p in ExistingPatches}:
    sum {(c,q) in CellPatchEdges: q = p} 1 >= 1;

check {(i,j) in E}:
    abs(ComponentMaxAreaCellHa[i] - ComponentMaxAreaCellHa[j])
    <= 0.000001;

check {(c,p) in CellPatchEdges}:
    abs(ComponentMaxAreaCellHa[c] - ComponentMaxAreaPatchHa[p])
    <= 0.000001;

check {k in RequiredNetworks}:
    NetworkMinimumHa[k]
    <= sum {c in Cells} NetworkCellAreaHa[k,c]
       + sum {p in ExistingPatches} NetworkPatchAreaHa[k,p]
       + 0.000001;
