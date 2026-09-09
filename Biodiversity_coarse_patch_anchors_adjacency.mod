# Stage B: hard anchored 200 m networks with spatial-coherence rewards.
#
# Every selected restoration component must connect to an activated existing-
# nature anchor and every root must collect at least MinNetworkAreaHa. The
# objective retains the Stage-A biodiversity and adjacency terms. Carbon is
# absent.

set Cells ordered;
set MustIncludeCells within Cells default {};

set Landuses;
set NatureLanduses within Landuses;
set ForestLanduses within NatureLanduses;
set WetLanduses within NatureLanduses;
set AllowedActions within {NatureLanduses, Cells};

set E within {Cells, Cells};
set SameHabitatEdges within {NatureLanduses, Cells, Cells};

set ExistingPatches;
set AnchorPatches within ExistingPatches;
set CellPatchEdges within {Cells, ExistingPatches};
set CellsTouchingExisting within Cells default {};
set ExistingHabitatContacts within {NatureLanduses, Cells} default {};

param BioDiversity {NatureLanduses, Cells} default 0;
param RestorableAreaHa {Cells} > 0;
param AdjacencyWeightHa {(i,j) in E} > 0;
param ExistingPatchAreaHa {ExistingPatches} > 0;
param ComponentMaxAreaCellHa {Cells} > 0;
param ComponentMaxAreaPatchHa {ExistingPatches} > 0;

param IsProductionForestCell {Cells} binary default 0;
param IsHighCarbonCell {Cells} binary default 0;

param CurrentNatureAreaHa >= 0;
param FinalNatureTargetHa > 0 default 25000;
param RestorationToleranceHa >= 0 default 4;
param MinNetworkAreaHa > 0 default 5000;
param MinRestorationAreaHa {NatureLanduses} >= 0 default 0;

param BiodiversityScalePerHa > 0 default 1000;
param AdjacencyNormalizerHa > 0;
param MatchingExistingNormalizerHa > 0 default 1;

param BiodiversityWeight >= 0 default 0.70;
param FootprintAdjacencyWeight >= 0 default 0.10;
param SameHabitatAdjacencyWeight >= 0 default 0.15;
param ExistingNatureContactWeight >= 0 default 0.05;
param MatchingExistingNatureWeight >= 0 default 0;

param RestorationTargetHa :=
    FinalNatureTargetHa - CurrentNatureAreaHa;
param BiodiversityNormalizer :=
    BiodiversityScalePerHa *
    (RestorationTargetHa + RestorationToleranceHa);
param ExistingContactNormalizerHa :=
    RestorationTargetHa + RestorationToleranceHa;

set NetworkArcs within {Cells, Cells} :=
    E union setof {(i,j) in E} (j,i);

check: card(Cells) > 0;
check: card(AllowedActions) > 0;
check: card(E) > 0;
check: card(AnchorPatches) > 0;
check: RestorationTargetHa > 0;
check: RestorationTargetHa <= sum {c in Cells} RestorableAreaHa[c];
check:
    abs(
        BiodiversityWeight +
        FootprintAdjacencyWeight +
        SameHabitatAdjacencyWeight +
        ExistingNatureContactWeight +
        MatchingExistingNatureWeight - 1
    ) <= 1e-8;

check {c in Cells}:
    sum {l in NatureLanduses: (l,c) in AllowedActions} 1 > 0;
check {(i,j) in E}: i <> j;
check {(l,i,j) in SameHabitatEdges}:
    (i,j) in E and
    (l,i) in AllowedActions and
    (l,j) in AllowedActions;
check {c in Cells}: ComponentMaxAreaCellHa[c] >= MinNetworkAreaHa;
check {p in ExistingPatches}:
    ComponentMaxAreaPatchHa[p] >= MinNetworkAreaHa;
check {(i,j) in E}:
    abs(ComponentMaxAreaCellHa[i] - ComponentMaxAreaCellHa[j]) <= 1e-6;
check {(c,p) in CellPatchEdges}:
    abs(ComponentMaxAreaCellHa[c] - ComponentMaxAreaPatchHa[p]) <= 1e-6;
check {p in ExistingPatches}:
    sum {(c,q) in CellPatchEdges: q = p} 1 > 0;
check {(l,c) in AllowedActions: IsProductionForestCell[c] = 1}:
    l in ForestLanduses;
check {(l,c) in AllowedActions: IsHighCarbonCell[c] = 1}:
    l in WetLanduses;

var LanduseDecision {(l,c) in AllowedActions} binary;
var Selected {c in Cells} >= 0, <= 1;
var SelectedAdjacency {(i,j) in E} >= 0, <= 1;
var SameHabitatAdjacency {(l,i,j) in SameHabitatEdges} >= 0, <= 1;

maximize ConservationIndex:
    BiodiversityWeight / BiodiversityNormalizer *
        sum {(l,c) in AllowedActions}
            BioDiversity[l,c] * LanduseDecision[l,c]
  + FootprintAdjacencyWeight / AdjacencyNormalizerHa *
        sum {(i,j) in E}
            AdjacencyWeightHa[i,j] * SelectedAdjacency[i,j]
  + SameHabitatAdjacencyWeight / AdjacencyNormalizerHa *
        sum {(l,i,j) in SameHabitatEdges}
            AdjacencyWeightHa[i,j] * SameHabitatAdjacency[l,i,j]
  + ExistingNatureContactWeight / ExistingContactNormalizerHa *
        sum {c in CellsTouchingExisting}
            RestorableAreaHa[c] * Selected[c]
  + MatchingExistingNatureWeight / MatchingExistingNormalizerHa *
        sum {(l,c) in ExistingHabitatContacts}
            RestorableAreaHa[c] * LanduseDecision[l,c];

subject to LinkSelectedToOneLanduse {c in Cells}:
    Selected[c]
    = sum {l in NatureLanduses: (l,c) in AllowedActions}
        LanduseDecision[l,c];

subject to FinalAreaLower:
    sum {c in Cells} RestorableAreaHa[c] * Selected[c]
    >= RestorationTargetHa;

subject to FinalAreaUpper:
    sum {c in Cells} RestorableAreaHa[c] * Selected[c]
    <= RestorationTargetHa + RestorationToleranceHa;

subject to MinimumRestorationArea {l in NatureLanduses}:
    sum {c in Cells: (l,c) in AllowedActions}
        RestorableAreaHa[c] * LanduseDecision[l,c]
    >= MinRestorationAreaHa[l];

subject to MustIncludeConstraint {c in MustIncludeCells}:
    Selected[c] = 1;

subject to SelectedAdjacencyFromI {(i,j) in E}:
    SelectedAdjacency[i,j] <= Selected[i];
subject to SelectedAdjacencyFromJ {(i,j) in E}:
    SelectedAdjacency[i,j] <= Selected[j];
subject to SameHabitatAdjacencyFromI {(l,i,j) in SameHabitatEdges}:
    SameHabitatAdjacency[l,i,j] <= LanduseDecision[l,i];
subject to SameHabitatAdjacencyFromJ {(l,i,j) in SameHabitatEdges}:
    SameHabitatAdjacency[l,i,j] <= LanduseDecision[l,j];

# ---- Anchored minimum-area networks -------------------------------------

var PatchActive {p in ExistingPatches} binary;
var AnchorRoot {p in AnchorPatches} binary;
var AreaFlowHa {(i,j) in NetworkArcs} >= 0;
var CellToPatchFlowHa {(c,p) in CellPatchEdges} >= 0;
var PatchToCellFlowHa {(c,p) in CellPatchEdges} >= 0;
var CollectedAreaHa {p in AnchorPatches} >= 0;

subject to RootNeedsActiveAnchor {p in AnchorPatches}:
    AnchorRoot[p] <= PatchActive[p];
subject to AbsorbOnlyAtAnchorRoot {p in AnchorPatches}:
    CollectedAreaHa[p]
    <= ComponentMaxAreaPatchHa[p] * AnchorRoot[p];
subject to MinimumAreaAtAnchorRoot {p in AnchorPatches}:
    CollectedAreaHa[p] >= MinNetworkAreaHa * AnchorRoot[p];

subject to CellFlowNeedsSelectedOrigin {(i,j) in NetworkArcs}:
    AreaFlowHa[i,j] <= ComponentMaxAreaCellHa[i] * Selected[i];
subject to CellFlowNeedsSelectedDestination {(i,j) in NetworkArcs}:
    AreaFlowHa[i,j] <= ComponentMaxAreaCellHa[j] * Selected[j];
subject to CellToPatchNeedsSelectedCell {(c,p) in CellPatchEdges}:
    CellToPatchFlowHa[c,p]
    <= ComponentMaxAreaCellHa[c] * Selected[c];
subject to CellToPatchNeedsActivePatch {(c,p) in CellPatchEdges}:
    CellToPatchFlowHa[c,p]
    <= ComponentMaxAreaPatchHa[p] * PatchActive[p];
subject to PatchToCellNeedsSelectedCell {(c,p) in CellPatchEdges}:
    PatchToCellFlowHa[c,p]
    <= ComponentMaxAreaCellHa[c] * Selected[c];
subject to PatchToCellNeedsActivePatch {(c,p) in CellPatchEdges}:
    PatchToCellFlowHa[c,p]
    <= ComponentMaxAreaPatchHa[p] * PatchActive[p];

subject to PatchActivationNeedsSelectedContact {p in ExistingPatches}:
    PatchActive[p]
    <= sum {(c,q) in CellPatchEdges: q = p} Selected[c];

subject to CandidateCellAreaBalance {c in Cells}:
    RestorableAreaHa[c] * Selected[c]
    + sum {(i,j) in NetworkArcs: j = c} AreaFlowHa[i,j]
    + sum {(d,p) in CellPatchEdges: d = c} PatchToCellFlowHa[d,p]
    = sum {(i,j) in NetworkArcs: i = c} AreaFlowHa[i,j]
    + sum {(d,p) in CellPatchEdges: d = c} CellToPatchFlowHa[d,p];

subject to ExistingPatchAreaBalance {p in ExistingPatches}:
    ExistingPatchAreaHa[p] * PatchActive[p]
    + sum {(c,q) in CellPatchEdges: q = p} CellToPatchFlowHa[c,q]
    = sum {(c,q) in CellPatchEdges: q = p} PatchToCellFlowHa[c,q]
    + (if p in AnchorPatches then CollectedAreaHa[p] else 0);
