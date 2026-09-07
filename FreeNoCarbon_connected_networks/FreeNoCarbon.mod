# FreeNoCarbon.mod -- biodiversity restoration with minimum connected area
#
# Adapted from the supplied FreeNoCarbon.mod. The biodiversity objective,
# optional same-land-use adjacency bonus, exact budget, and minimum number
# of newly selected cells per nature type are retained.
#
# HARD RULE: every newly restored cell must be in a connected component of
# existing nature + new restoration with at least MinPatchAreaHa hectares.
# Different nature types may connect; multiple qualifying networks are allowed.
# Existing nature is retained everywhere. Small existing-only fragments do
# not have to be enlarged. Unselected agriculture/production forest, urban,
# and other non-nature cells cannot connect or contribute area.
#
# At 10 x 10 m, CellAreaHa = 0.01 and 5000 ha = 500000 cells.
# Supply four-neighbour (shared-edge) adjacency in E. The hard connectivity
# rule uses E in both directions, even if each input edge is listed only once.
# Only cells supplied in Cells can contribute area or provide connections.
#
# Run: model FreeNoCarbon.mod; data your_area.dat; option solver highs; solve;
# See FreeNoCarbon_demo.run and FreeNoCarbon_GUIDE.md for examples and tests.

# ---- Original inputs, with explicit nature-action scope ------------------

set Cells;
set MustIncludeCells within Cells default {};
set Landuses;

# Explicitly list the eight restoration nature types here in the data file.
# Do not include Ag, Agriculture, Urban, ProductionForest, or Other.
# Landuses may still contain Ag or other legacy non-nature labels: decisions
# for every label outside NatureLanduses are fixed to zero below.
set NatureLanduses within Landuses;

set E within {Cells, Cells};

# One-hot current class indicators (zero for unlisted current classes).
# For existing nature, set its nature-type indicator to 1 and CanChange to 0.
param Existingnature {Landuses, Cells} binary default 0;

# The supplied combined biodiversity index is used without redefining it.
# As in the original model, these two tables must be supplied for all pairs.
param BioDiversity {Landuses, Cells};
param TransitionCost {Landuses, Cells} >= 0;

# 1 only for eligible agriculture/production-forest restoration cells.
# 0 for current nature, urban, other excluded land, and unavailable cells.
param CanChange {Cells} binary;
param b >= 0;
param MinLan integer >= 0 default 0;
param SpatialContiguityBonus >= 0 default 0;

# ---- New connected-area inputs ------------------------------------------

param CellAreaHa {Cells} > 0 default 0.01;
param MinPatchAreaHa > 0 default 5000;

check: card(NatureLanduses) > 0;
check:
    card(NatureLanduses inter
         {'Ag', 'Agriculture', 'Urban', 'ProductionForest', 'Other'}) = 0;
check {c in Cells}:
    sum {l in Landuses} Existingnature[l,c] <= 1;
check {(i,j) in E}: i <> j;

param IsExistingNature {c in Cells} binary :=
    sum {l in NatureLanduses} Existingnature[l,c];

check {c in Cells}:
    IsExistingNature[c] + CanChange[c] <= 1;

# MustIncludeCells keeps its original meaning: mandatory NEW restoration.
# Existing nature is already retained; do not put it in this set.
check {c in MustIncludeCells}: CanChange[c] = 1;

# Eligible area is a valid upper bound on any flow or root's collected area.
# It is not a target, and separate components cannot exchange it.
param MaxNetworkAreaHa :=
    sum {c in Cells} CellAreaHa[c] *
        (IsExistingNature[c] + CanChange[c]);

set NetworkArcs within {Cells, Cells} :=
    E union (setof {(i,j) in E} (j,i));

# ---- Original restoration decisions and objective -----------------------

var LanduseDecision {l in Landuses, c in Cells} binary;
var Contiguity {l in Landuses, (i,j) in E} binary;

# The original optional bonus is intentionally kept exactly as a sum over E.
# Thus, if E includes both (i,j) and (j,i), new/new same-type adjacency earns
# two bonus terms, just as before. E orientation never changes the HARD rule.
maximize ConservationIndex:
    sum {l in Landuses, c in Cells}
        LanduseDecision[l,c] * BioDiversity[l,c] * CanChange[c]
    + SpatialContiguityBonus * sum {(i,j) in E, l in Landuses} (
        Contiguity[l,i,j] * CanChange[i] * CanChange[j]
        + Existingnature[l,i] * LanduseDecision[l,j] * CanChange[j]
    );

subject to PropotionalUse {c in Cells}:
    sum {l in Landuses} LanduseDecision[l,c] <= 1;

# Unlike multiplication by CanChange in the objective alone, this constraint
# prevents non-restorable cells from being selected as free connecting land.
subject to RestoreOnlyEligible {c in Cells}:
    sum {l in Landuses} LanduseDecision[l,c] <= CanChange[c];

subject to NatureActionsOnly {l in Landuses diff NatureLanduses, c in Cells}:
    LanduseDecision[l,c] = 0;

# This remains a GLOBAL minimum of NEW cells per type, not a per-network
# requirement. Use MinLan = 0 when you do not want mandatory type quotas.
subject to MinimumCellPerLandUse {l in NatureLanduses}:
    sum {c in Cells} LanduseDecision[l,c] >= MinLan;

# Preserved from the original: spend exactly b, not merely at most b.
subject to Budget:
    sum {l in Landuses, c in Cells}
        LanduseDecision[l,c] * TransitionCost[l,c] = b;

subject to MustIncludeConstraint {c in MustIncludeCells}:
    sum {l in Landuses} LanduseDecision[l,c] = 1;

subject to DefineContiguity1 {l in Landuses, (i,j) in E}:
    Contiguity[l,i,j] <= LanduseDecision[l,i];

subject to DefineContiguity2 {l in Landuses, (i,j) in E}:
    Contiguity[l,i,j] <= LanduseDecision[l,j];

subject to DefineContiguity3 {l in Landuses, (i,j) in E}:
    Contiguity[l,i,j] >= LanduseDecision[l,i] + LanduseDecision[l,j] - 1;

# ---- Hard minimum area for every restoration-containing component -------
#
# Single-commodity AREA-COLLECTION certificate (continuous flow in hectares).
# Each participating nature cell supplies its own area exactly once. Area can
# move only between participating neighbouring cells and must be absorbed at
# a root. Every root must absorb at least MinPatchAreaHa, so every participating
# connected component has at least that much area.
#
# Every restored cell MUST participate. Existing nature MAY participate in
# the certificate without requiring isolated existing fragments to qualify.
# NetworkActive = 0 for existing nature does NOT remove or convert that nature.
# Roots can be chosen only on restoration cells; they are bookkeeping points,
# not ecological centres or separate network identifiers.
#
# Necessity: sum AreaBalance over a participating component. Internal flows
# cancel, so its total area equals the area absorbed at its roots. Positive
# area requires a root, and each root absorbs >= MinPatchAreaHa.
# Sufficiency: for any qualifying final nature component containing restoration,
# activate its cells, choose a restored cell as root, and collect each cell's
# area along a spanning tree. All constraints hold. Thus this does not force
# one global network or reject a valid restoration-containing component.

var NetworkActive {c in Cells} binary;
var NetworkRoot {c in Cells} binary;
var AreaFlowHa {(i,j) in NetworkArcs} >= 0;
var CollectedAreaHa {c in Cells} >= 0;

subject to RestorationMustParticipate {c in Cells}:
    NetworkActive[c] >= sum {l in Landuses} LanduseDecision[l,c];

subject to OnlyNatureCanParticipate {c in Cells}:
    NetworkActive[c] <= IsExistingNature[c]
        + sum {l in Landuses} LanduseDecision[l,c];

subject to RootOnRestoration {c in Cells}:
    NetworkRoot[c] <= sum {l in Landuses} LanduseDecision[l,c];

subject to FlowNeedsActiveOrigin {(i,j) in NetworkArcs}:
    AreaFlowHa[i,j] <= MaxNetworkAreaHa * NetworkActive[i];

subject to FlowNeedsActiveDestination {(i,j) in NetworkArcs}:
    AreaFlowHa[i,j] <= MaxNetworkAreaHa * NetworkActive[j];

subject to AbsorbOnlyAtRoot {c in Cells}:
    CollectedAreaHa[c] <= MaxNetworkAreaHa * NetworkRoot[c];

subject to MinimumAreaAtRoot {c in Cells}:
    CollectedAreaHa[c] >= MinPatchAreaHa * NetworkRoot[c];

subject to AreaBalance {c in Cells}:
    CellAreaHa[c] * NetworkActive[c]
    + sum {(i,j) in NetworkArcs: j = c} AreaFlowHa[i,j]
    = CollectedAreaHa[c]
    + sum {(i,j) in NetworkArcs: i = c} AreaFlowHa[i,j];
