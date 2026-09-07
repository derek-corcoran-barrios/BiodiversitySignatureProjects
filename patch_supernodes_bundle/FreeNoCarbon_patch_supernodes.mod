# FreeNoCarbon_patch_supernodes.mod
#
# Biodiversity restoration with a hard minimum connected-area rule.
# Existing natural patches are compressed to area-bearing supernodes, while
# only eligible restoration cells remain in the 10 m decision graph.
#
# Data contract
# -------------
# * Cells: eligible restoration cells only (normally Agriculture and
#   ProductionForest). A full raster domain is accepted through CanChange,
#   but a candidate-only domain is much smaller and strongly recommended.
# * E: rook-adjacent candidate-cell pairs. Write each undirected pair once.
# * ExistingPatches: connected current-nature patches that touch at least one
#   candidate cell in Cells.
# * CellPatchEdges: a pair (cell, patch) when the candidate cell shares an
#   edge with that current-nature patch.
# * ExistingPatchAreaHa: the area of the whole connected patch, counted once.
#
# HARD RULE: every selected restoration cell belongs to a connected network
# of selected cells plus current-natural patches with total area at least
# MinPatchAreaHa. Multiple qualifying networks are allowed. Existing-only
# patches are retained outside the decision model and need not reach the
# threshold unless restoration is attached to them.
#
# At 10 x 10 m, CellAreaHa = 0.01 and 5000 ha = 500000 cells. Patch area is
# supplied directly in hectares, so a 5000 ha current patch is represented by
# one supernode rather than 500000 cell nodes.

# ---- Sets and data -------------------------------------------------------

set Cells ordered;
set MustIncludeCells within Cells default {};

set Landuses;
set NatureLanduses within Landuses;

# Candidate-to-candidate rook adjacency. Supply each undirected edge once.
set E within {Cells, Cells};

# Only patches touching the candidate graph need to be supplied.
set ExistingPatches;
set CellPatchEdges within {Cells, ExistingPatches};

# Biodiversity value of assigning nature type l to candidate cell c.
# This can be the user's combined 0--1000 index without redefining it here.
param BioDiversity {Landuses, Cells} default 0;

# Cost of each restoration action. A constant cost of one reproduces a
# cell-count budget.
param TransitionCost {Landuses, Cells} >= 0 default 1;

# Candidate-only .dat files can omit this because the default is one.
param CanChange {Cells} binary default 1;

param CellAreaHa {Cells} > 0 default 0.01;
param ExistingPatchAreaHa {ExistingPatches} > 0;
param MinPatchAreaHa > 0 default 5000;

# Optional contextual objective term. IDWArea_Distance or its normalized
# version may be supplied here. It is a benefit score, not a patch identity
# and not evidence that the hard minimum-area rule is met.
param AreaInfluence {Cells} >= 0 default 0;
param AreaInfluenceWeight >= 0 default 0;

# Optional same-type contacts with current nature, already counted in R.
# For example, value 2 means cell c shares two rook edges with current nature
# of type l. Set SpatialContiguityBonus to zero to disable the whole bonus.
param ExistingNatureContacts {Landuses, Cells} integer >= 0 default 0;
param SpatialContiguityBonus >= 0 default 0;

# Existing and desired final areas by nature type. These optional parameters
# let pre-counted planning-area totals enforce final-area targets without
# adding existing-nature cells as decisions.
param ExistingLanduseAreaHa {Landuses} >= 0 default 0;
param MinFinalLanduseAreaHa {Landuses} >= 0 default 0;
param MinRestorationAreaHa {Landuses} >= 0 default 0;

# Legacy global minimum count per nature type. Prefer area-based targets above
# when cell sizes can vary. It defaults to zero.
param MinLan integer >= 0 default 0;

# b is a cost budget. ExactBudget = 1 preserves the original equality;
# ExactBudget = 0 treats b as an upper bound.
param b >= 0;
param ExactBudget binary default 1;

set ChangeableCells := {c in Cells: CanChange[c] = 1};
set NetworkArcs within {Cells, Cells} :=
    E union setof {(i,j) in E} (j,i);

# Avoid creating |NatureLanduses| x |E| extra binary variables when the
# optional soft contiguity bonus is disabled (the recommended first run).
set BonusEdges within {Cells, Cells} :=
    {(i,j) in E: SpatialContiguityBonus > 0};

# A valid upper bound for any flow or amount collected at one root.
param MaxNetworkAreaHa :=
    sum {c in ChangeableCells} CellAreaHa[c]
    + sum {p in ExistingPatches} ExistingPatchAreaHa[p];

# ---- Defensive checks ---------------------------------------------------

check: card(NatureLanduses) > 0;
check: MaxNetworkAreaHa >= MinPatchAreaHa;
check {c in MustIncludeCells}: CanChange[c] = 1;
check {(i,j) in E}: i <> j;
check {(i,j) in E}: CanChange[i] = 1 and CanChange[j] = 1;
check {(c,p) in CellPatchEdges}: CanChange[c] = 1;

# ---- Restoration decisions and objective -------------------------------

var LanduseDecision {l in NatureLanduses, c in ChangeableCells} binary;
var Contiguity {l in NatureLanduses, (i,j) in BonusEdges} binary;

maximize ConservationIndex:
    sum {l in NatureLanduses, c in ChangeableCells}
        LanduseDecision[l,c] *
        (BioDiversity[l,c] + AreaInfluenceWeight * AreaInfluence[c])
    + SpatialContiguityBonus * (
        sum {l in NatureLanduses, (i,j) in BonusEdges} Contiguity[l,i,j]
        + sum {l in NatureLanduses, c in ChangeableCells}
            ExistingNatureContacts[l,c] * LanduseDecision[l,c]
    );

subject to AtMostOneNatureType {c in ChangeableCells}:
    sum {l in NatureLanduses} LanduseDecision[l,c] <= 1;

subject to MinimumCellsPerNatureType {l in NatureLanduses}:
    sum {c in ChangeableCells} LanduseDecision[l,c] >= MinLan;

subject to MinimumRestorationArea {l in NatureLanduses}:
    sum {c in ChangeableCells}
        CellAreaHa[c] * LanduseDecision[l,c]
    >= MinRestorationAreaHa[l];

subject to MinimumFinalArea {l in NatureLanduses}:
    ExistingLanduseAreaHa[l]
    + sum {c in ChangeableCells}
        CellAreaHa[c] * LanduseDecision[l,c]
    >= MinFinalLanduseAreaHa[l];

subject to BudgetUpper:
    sum {l in NatureLanduses, c in ChangeableCells}
        TransitionCost[l,c] * LanduseDecision[l,c]
    <= b;

subject to BudgetLower:
    sum {l in NatureLanduses, c in ChangeableCells}
        TransitionCost[l,c] * LanduseDecision[l,c]
    >= ExactBudget * b;

subject to MustIncludeConstraint {c in MustIncludeCells}:
    sum {l in NatureLanduses} LanduseDecision[l,c] = 1;

subject to DefineContiguity1 {l in NatureLanduses, (i,j) in BonusEdges}:
    Contiguity[l,i,j] <= LanduseDecision[l,i];

subject to DefineContiguity2 {l in NatureLanduses, (i,j) in BonusEdges}:
    Contiguity[l,i,j] <= LanduseDecision[l,j];

subject to DefineContiguity3 {l in NatureLanduses, (i,j) in BonusEdges}:
    Contiguity[l,i,j] >=
        LanduseDecision[l,i] + LanduseDecision[l,j] - 1;

# ---- Patch-supernode minimum connected area -----------------------------
#
# Every selected cell supplies its own area. An activated current-nature patch
# supplies ExistingPatchAreaHa exactly once, even if many candidate cells touch
# it. Area can move through selected candidate cells and activated patches, and
# is collected at selected-cell roots. Each root must collect at least the
# threshold. Summing the balances over a connected component cancels internal
# flows; therefore every restoration-containing component has enough area.

var PatchActive {p in ExistingPatches} binary;
var NetworkRoot {c in ChangeableCells} binary;
var AreaFlowHa {(i,j) in NetworkArcs} >= 0;
var CellToPatchFlowHa {(c,p) in CellPatchEdges} >= 0;
var PatchToCellFlowHa {(c,p) in CellPatchEdges} >= 0;
var CollectedAreaHa {c in ChangeableCells} >= 0;

subject to RootOnSelectedCell {c in ChangeableCells}:
    NetworkRoot[c]
    <= sum {l in NatureLanduses} LanduseDecision[l,c];

subject to AbsorbOnlyAtRoot {c in ChangeableCells}:
    CollectedAreaHa[c] <= MaxNetworkAreaHa * NetworkRoot[c];

subject to MinimumAreaAtRoot {c in ChangeableCells}:
    CollectedAreaHa[c] >= MinPatchAreaHa * NetworkRoot[c];

subject to CellFlowNeedsSelectedOrigin {(i,j) in NetworkArcs}:
    AreaFlowHa[i,j]
    <= MaxNetworkAreaHa *
        sum {l in NatureLanduses} LanduseDecision[l,i];

subject to CellFlowNeedsSelectedDestination {(i,j) in NetworkArcs}:
    AreaFlowHa[i,j]
    <= MaxNetworkAreaHa *
        sum {l in NatureLanduses} LanduseDecision[l,j];

subject to CellToPatchNeedsSelectedCell {(c,p) in CellPatchEdges}:
    CellToPatchFlowHa[c,p]
    <= MaxNetworkAreaHa *
        sum {l in NatureLanduses} LanduseDecision[l,c];

subject to CellToPatchNeedsActivePatch {(c,p) in CellPatchEdges}:
    CellToPatchFlowHa[c,p]
    <= MaxNetworkAreaHa * PatchActive[p];

subject to PatchToCellNeedsSelectedCell {(c,p) in CellPatchEdges}:
    PatchToCellFlowHa[c,p]
    <= MaxNetworkAreaHa *
        sum {l in NatureLanduses} LanduseDecision[l,c];

subject to PatchToCellNeedsActivePatch {(c,p) in CellPatchEdges}:
    PatchToCellFlowHa[c,p]
    <= MaxNetworkAreaHa * PatchActive[p];

subject to CandidateCellAreaBalance {c in ChangeableCells}:
    CellAreaHa[c] *
        sum {l in NatureLanduses} LanduseDecision[l,c]
    + sum {(i,j) in NetworkArcs: j = c} AreaFlowHa[i,j]
    + sum {(d,p) in CellPatchEdges: d = c} PatchToCellFlowHa[d,p]
    = CollectedAreaHa[c]
    + sum {(i,j) in NetworkArcs: i = c} AreaFlowHa[i,j]
    + sum {(d,p) in CellPatchEdges: d = c} CellToPatchFlowHa[d,p];

subject to ExistingPatchAreaBalance {p in ExistingPatches}:
    ExistingPatchAreaHa[p] * PatchActive[p]
    + sum {(c,q) in CellPatchEdges: q = p} CellToPatchFlowHa[c,q]
    = sum {(c,q) in CellPatchEdges: q = p} PatchToCellFlowHa[c,q];
