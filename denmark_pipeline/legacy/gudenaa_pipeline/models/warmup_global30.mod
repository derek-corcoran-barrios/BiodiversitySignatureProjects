# Eight-habitat warm-up. Subarea budgets are diagnostic only.
set Cells ordered;
set Landuses;
set Subareas ordered;
set AllowedActions within {Landuses, Cells};
set E within {Cells, Cells};
set SameHabitatEdges within {Landuses, Cells, Cells};
param BioDiversity {AllowedActions};
param ExistingNatureContact {AllowedActions} >= 0 default 0;
param RestorableAreaHa {Cells} > 0;
param RestorableAreaHaBySubarea {Subareas, Cells} >= 0 default 0;
param SubareaAreaHa {Subareas} > 0;
param CurrentNatureAreaHa {Subareas} >= 0;
param FinalNatureTargetHa {Subareas} >= 0;
param RestorationRequirementHa {Subareas} >= 0;
param RestorationToleranceHa {Subareas} >= 0;
param BiodiversityNormalizer > 0;
param ContiguityNormalizer > 0;
param ExistingNatureContactNormalizer > 0;
param SpatialContiguityBonus >= 0;
param NewToNewAdjacencyWeight >= 0;
param ExistingToNewAdjacencyWeight >= 0;
param GlobalRestorationToleranceHa >= 0 default 16;
param GlobalRestorationRequirementHa :=
    sum {s in Subareas} (FinalNatureTargetHa[s] - CurrentNatureAreaHa[s]);
var Selected {Cells} binary;
var LanduseDecision {AllowedActions} binary;
var Contiguity {SameHabitatEdges} binary;
maximize ConservationIndex:
    sum {(l,c) in AllowedActions} BioDiversity[l,c] * LanduseDecision[l,c]
      / BiodiversityNormalizer
    + SpatialContiguityBonus * (
      NewToNewAdjacencyWeight * sum {(l,i,j) in SameHabitatEdges}
        Contiguity[l,i,j] / ContiguityNormalizer
      + ExistingToNewAdjacencyWeight * sum {(l,c) in AllowedActions}
        ExistingNatureContact[l,c] * LanduseDecision[l,c]
        / ExistingNatureContactNormalizer);
subject to LinkSelectionToLanduse {c in Cells}:
    Selected[c] = sum {l in Landuses: (l,c) in AllowedActions} LanduseDecision[l,c];
subject to GlobalRestorationMinimum:
    sum {c in Cells} RestorableAreaHa[c] * Selected[c] >= GlobalRestorationRequirementHa;
subject to GlobalRestorationMaximum:
    sum {c in Cells} RestorableAreaHa[c] * Selected[c]
      <= GlobalRestorationRequirementHa + GlobalRestorationToleranceHa;
subject to DefineContiguityFromI {(l,i,j) in SameHabitatEdges}:
    Contiguity[l,i,j] <= LanduseDecision[l,i];
subject to DefineContiguityFromJ {(l,i,j) in SameHabitatEdges}:
    Contiguity[l,i,j] <= LanduseDecision[l,j];
subject to DefineContiguityLower {(l,i,j) in SameHabitatEdges}:
    Contiguity[l,i,j] >= LanduseDecision[l,i] + LanduseDecision[l,j] - 1;
check: card(Landuses) = 8;
check: card(Subareas) = 4;
check: GlobalRestorationRequirementHa >= 0;
check {c in Cells}: card({l in Landuses: (l,c) in AllowedActions}) > 0;
check: abs(NewToNewAdjacencyWeight + ExistingToNewAdjacencyWeight - 1) < 1e-9;
