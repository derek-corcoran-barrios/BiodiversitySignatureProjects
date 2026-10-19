# Generic outcomes and zero or more anchored-network certificates.
set Units ordered;
set Outcomes ordered;
set Actions within {Outcomes, Units};
set SameEdges within {Outcomes, Units, Units};
set Networks ordered;
set Nodes {Networks} ordered;
set Arcs dimen 3;
param Area {Units} > 0;
param Cost {Units} >= 0;
param Score {Actions};
param AnchorContact {Units} >= 0 default 0;
param MinimumArea >= 0;
param MaximumArea >= MinimumArea;
param MaximumCost >= 0;
param BioNorm > 0 default 1;
param EdgeNorm > 0 default 1;
param AnchorNorm > 0 default 1;
param SpatialWeight >= 0 default 1.2;
param NewWeight >= 0 default 0.5;
param AnchorWeight >= 0 default 0.5;
param Owner {k in Networks, n in Nodes[k]} integer >= 0;
param RootAllowed {k in Networks, n in Nodes[k]} binary;
param CreditArea {k in Networks, n in Nodes[k]} >= 0;
param FlowCap {k in Networks, n in Nodes[k]} > 0;
param Requirement {Networks} > 0;
param RequireNew binary default 1;
param ReferenceSelected {Units} binary default 0;
param MaximumChange >= 0 default sum {u in Units} Area[u];
param FinalSpatialWeight >= 0 default SpatialWeight;
var X {Actions} binary;
var Selected {Units} binary;
var Adjacent {SameEdges} binary;
var Built {Networks} binary;
var Shortfall {k in Networks} >= 0, <= Requirement[k];
var Member {k in Networks, n in Nodes[k]} binary;
var Root {k in Networks, n in Nodes[k]} binary;
var Supply {k in Networks, n in Nodes[k]} >= 0;
var Flow {Arcs} >= 0;
minimize MinShortfall: sum {k in Networks} Shortfall[k] / Requirement[k];
minimize FootprintChange:
  sum {u in Units} Area[u] * (ReferenceSelected[u] * (1-Selected[u])
    + (1-ReferenceSelected[u]) * Selected[u]);
maximize Quality:
  sum {(o,u) in Actions} Score[o,u]*X[o,u]/BioNorm
  + SpatialWeight * (NewWeight * sum {(o,u,v) in SameEdges} Adjacent[o,u,v]/EdgeNorm
  + AnchorWeight * sum {u in Units} AnchorContact[u]*Selected[u]/AnchorNorm);
subject to ChooseOne {u in Units}: Selected[u] = sum {o in Outcomes: (o,u) in Actions} X[o,u];
subject to AreaLower: sum {u in Units} Area[u]*Selected[u] >= MinimumArea;
subject to AreaUpper: sum {u in Units} Area[u]*Selected[u] <= MaximumArea;
subject to CostUpper: sum {u in Units} Cost[u]*Selected[u] <= MaximumCost;
subject to ChangeLimit:
  sum {u in Units} Area[u] * (ReferenceSelected[u] * (1-Selected[u])
    + (1-ReferenceSelected[u]) * Selected[u]) <= MaximumChange;
subject to EdgeLeft {(o,u,v) in SameEdges}: Adjacent[o,u,v] <= X[o,u];
subject to EdgeRight {(o,u,v) in SameEdges}: Adjacent[o,u,v] <= X[o,v];
subject to EdgeBoth {(o,u,v) in SameEdges}: Adjacent[o,u,v] >= X[o,u]+X[o,v]-1;
subject to Active {k in Networks, n in Nodes[k]: Owner[k,n] > 0}: Member[k,n] <= Selected[Owner[k,n]];
subject to BuiltMember {k in Networks, n in Nodes[k]}: Member[k,n] <= Built[k];
subject to OneRoot {k in Networks}: sum {n in Nodes[k]} Root[k,n] = Built[k];
subject to RootEligible {k in Networks, n in Nodes[k]}: Root[k,n] <= RootAllowed[k,n];
subject to RootMember {k in Networks, n in Nodes[k]}: Root[k,n] <= Member[k,n];
subject to HasRestoration {k in Networks}: sum {n in Nodes[k]: Owner[k,n] > 0} Member[k,n] >= RequireNew*Built[k];
subject to MeetsArea {k in Networks}: sum {n in Nodes[k]} CreditArea[k,n]*Member[k,n] + Shortfall[k] >= Requirement[k];
subject to SourceOnly {k in Networks, n in Nodes[k]}: Supply[k,n] <= FlowCap[k,n]*Root[k,n];
subject to ArcFrom {(k,i,j) in Arcs}: Flow[k,i,j] <= FlowCap[k,i]*Member[k,i];
subject to ArcTo {(k,i,j) in Arcs}: Flow[k,i,j] <= FlowCap[k,j]*Member[k,j];
# Every member consumes one flow unit, including zero-credit lakes/anchors.
subject to Balance {k in Networks, n in Nodes[k]}:
  Supply[k,n] + sum {(q,i,j) in Arcs: q=k and j=n} Flow[q,i,j]
  = Member[k,n] + sum {(q,i,j) in Arcs: q=k and i=n} Flow[q,i,j];
check {u in Units}: card({o in Outcomes: (o,u) in Actions}) > 0;
check {k in Networks, n in Nodes[k]: Owner[k,n] > 0}: Owner[k,n] in Units;
check {(k,i,j) in Arcs}: k in Networks and i in Nodes[k] and j in Nodes[k];
