# Stage A: fast 200 m biodiversity and spatial-coherence seed model.
#
# This model has no network-flow variables and no hard 5,000 ha requirement.
# It creates an informative seed by balancing biodiversity, compactness,
# same-habitat adjacency, and contact with existing nature. Carbon is absent.

set Cells ordered;
set MustIncludeCells within Cells default {};

set Landuses;
set NatureLanduses within Landuses;
set ForestLanduses within NatureLanduses;
set WetLanduses within NatureLanduses;
set AllowedActions within {NatureLanduses, Cells};

# Undirected rook edges, listed once.
set E within {Cells, Cells};

# Triples (landuse, cell_i, cell_j) for which both endpoints permit the same
# destination. These are presolved in the data generator.
set SameHabitatEdges within {NatureLanduses, Cells, Cells};

# Current-nature contacts available now. The untyped set rewards contact with
# any existing patch. The habitat-specific set is deliberately empty until the
# fine-resolution habitat-contact layer is prepared.
set CellsTouchingExisting within Cells default {};
set ExistingHabitatContacts within {NatureLanduses, Cells} default {};

param BioDiversity {NatureLanduses, Cells} default 0;
param RestorableAreaHa {Cells} > 0;
param AdjacencyWeightHa {(i,j) in E} > 0;

param CurrentNatureAreaHa >= 0;
param FinalNatureTargetHa > 0 default 25000;
param RestorationToleranceHa >= 0 default 4;
param MinRestorationAreaHa {NatureLanduses} >= 0 default 0;

# The biodiversity rasters are normalized to 0--1, combined, multiplied by
# 1,000, and integrated over hectares. This provides a theoretical scaling
# bound for a fixed restoration-area target.
param BiodiversityScalePerHa > 0 default 1000;
param AdjacencyNormalizerHa > 0;
param MatchingExistingNormalizerHa > 0 default 1;

# Initial experimental weights. They must sum to one. MatchingExisting receives
# zero until habitat-specific existing-nature contacts have been exported.
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

check: card(Cells) > 0;
check: card(AllowedActions) > 0;
check: card(E) > 0;
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

var LanduseDecision {(l,c) in AllowedActions} binary;

# Equality to a sum of mutually exclusive binaries makes Selected integral, so
# it need not itself be declared binary.
var Selected {c in Cells} >= 0, <= 1;

# These auxiliaries also need not be binary. Positive objective coefficients
# push them to one exactly when both corresponding decisions equal one.
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
