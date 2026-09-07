# Connected-nature restoration model

`FreeNoCarbon.mod` extends your original AMPL model. Your supplied biodiversity
index, optional same-type contiguity bonus, exact budget, and global minimum
number of restored cells per nature type are retained.

## The hard rule

Every restored cell must belong to a connected component of **existing nature
plus selected restoration** covering at least **5,000 ha**. Several separate
qualifying components are allowed. All eight nature types can connect to each
other; this is not a separate 5,000-ha requirement for each type.

At 10 x 10 m, each cell is 0.01 ha. Therefore 500,000 cells = 5,000 ha = 50 km2.
Existing nature contributes its area without a restoration decision or cost.
An existing-only small fragment is retained but need not reach the threshold.
Unconverted agriculture or production forest, urban land, and other non-nature
cells cannot contribute area or form a bridge.

This model permits single-cell-wide corridors. It does not impose a minimum
corridor width. It counts only the cells supplied in the planning area: include
relevant surrounding existing nature in `Cells` if it should count.

## Files and quick test

Keep these files together:

- `FreeNoCarbon.mod`: production model, default minimum **5,000 ha**.
- `FreeNoCarbon_demo.dat`: synthetic 12-cell dataset, deliberately using a
  **0.05-ha** minimum to test the logic on five-cell networks.
- `FreeNoCarbon_demo.run`: four small tests, using the HiGHS MILP solver.

From a disposable AMPL session in this directory:

```ampl
include FreeNoCarbon_demo.run;
```

The script calls `reset`, so do not run it in an AMPL session holding unsaved
work. An AMPL interpreter and an available MILP solver are required; change the
solver line in the demo if you use Gurobi or CPLEX instead of HiGHS.

The tests check that two disconnected five-cell networks can each qualify;
they cannot be pooled to meet a six-cell minimum; selecting and paying for a
bridge connects them; and area cannot be counted twice. A small, isolated
existing-nature fragment is also retained without being forced to qualify.

**Do not reuse the demo's reduced area threshold for production.**

## Preparing your planning-area data

| Input | Meaning |
| --- | --- |
| `Cells` | IDs of all cells in the supplied planning area, including existing nature. |
| `Landuses` | Action labels; may include the legacy `Ag` label, whose decision is fixed to zero. |
| `NatureLanduses` | Explicit subset of `Landuses` containing the eight restoration nature types. |
| `E` | Neighbour pairs sharing an edge; do not include self-pairs, wrap raster rows, or bridge real gaps. |
| `Existingnature[l,c]` | 1 for a cell's current class and 0 otherwise; at most one 1 per cell. Non-nature current classes may be left all zero. |
| `CanChange[c]` | 1 for eligible agriculture/production-forest cells; 0 for current nature, urban land, and other unavailable cells. |
| `BioDiversity[l,c]` | Your supplied biodiversity score for assigning action `l` to cell `c`. |
| `TransitionCost[l,c]` | Restoration cost for that action and cell, in the same currency/units as `b`. |
| `b` | Budget to spend **exactly**, as in your original model. |
| `MinLan` | Global minimum number of NEW restored cells per nature type; default 0. |
| `MustIncludeCells` | Mandatory restoration cells, not existing nature; default empty. |
| `SpatialContiguityBonus` | Original optional same-type adjacency reward; default 0. |
| `CellAreaHa[c]` | Positive cell area in hectares; default 0.01 for every 10-m cell. |
| `MinPatchAreaHa` | Minimum area of each restoration-containing network; default 5000. |

Supply complete `BioDiversity` and `TransitionCost` tables, including unused
legacy action rows if present in `Landuses`. Do not accidentally default a
missing eligible action's cost to zero. Zero scores and costs for unused
non-nature action rows are harmless because those decisions are fixed to zero.

The additional nature-type declaration for your real data would be:

```ampl
set NatureLanduses :=
    ForestDryPoor ForestWetPoor ForestDryRich ForestWetRich
    OpenDryPoor OpenWetPoor OpenDryRich OpenWetRich;

param MinPatchAreaHa := 5000;
```

These labels must also occur in `Landuses`, and must exactly match the score
and cost table labels. Current nature must use the same labels in
`Existingnature` and have `CanChange = 0`.

For the hard connectivity rule, an edge in either direction is sufficient:
the model constructs both flow directions internally. For the **original soft
bonus**, `E` is used exactly as supplied. Listing both directions counts a
new/new same-type contact twice, as in your original objective; listing only
one direction makes the existing/new bonus orientation-dependent. Keep your
chosen original convention or set the bonus to zero if you only need the hard
network rule. The hard rule never depends on this bonus.

## Run with your own data

```ampl
reset;
model FreeNoCarbon.mod;
data your_area.dat;
option solver highs;
solve;
display solve_result, ConservationIndex;
display {l in NatureLanduses, c in Cells:
    LanduseDecision[l,c] > 0.5} LanduseDecision[l,c];
```

`your_area.dat` must contain your actual AOI, adjacency, existing nature,
eligibility, biodiversity scores, costs, and budget. The supplied example is
only a logic test, not Gudenaa data.

The exact budget equality is deliberately unchanged. If you later choose an
"at most" budget, the deliberate edit is to change `= b` to `<= b` in `Budget`.
Insufficient connected eligible area, an incompatible exact budget, impossible
`MinLan` quotas, or mandatory restoration in a too-small component will correctly
make the problem infeasible. The model never silently lowers the threshold.

## Reading and auditing the result

- `LanduseDecision[l,c] = 1` is a restoration decision to map back to the raster.
- Final nature is `IsExistingNature[c] + sum {l in Landuses} LanduseDecision[l,c]`.
  This retains ALL existing nature, including small untouched fragments.
- `NetworkActive`, `NetworkRoot`, `AreaFlowHa`, and `CollectedAreaHa` are only
  mathematical certificates. A root is not a network ID. More than one root
  can occur in a physical connected component.
- Do **not** use `NetworkActive` as the final nature map: optional existing
  nature can have `NetworkActive = 0` while still being retained physically.

After mapping the solution back to your 10-m grid, make one combined binary
mask of final nature and identify four-neighbour connected components. Sum
their areas and check the 5,000-ha minimum **only for components containing at
least one selected restoration cell**. This independently checks the edge list,
cell IDs, units, and solution tolerances as well as the optimization result.

The area-flow certificate works because each participating cell supplies its
area exactly once, flow cannot cross inactive cells, and each root must absorb
at least the minimum. Within a component all internal flows cancel, so it
cannot invent area or borrow it from a disconnected component. Conversely,
any qualifying final nature component containing restoration can send its
area along a spanning tree to one restored cell chosen as its root.

## Validation performed

Validated with AMPL 20260809 and HiGHS 1.15.1:

- All four bundled synthetic tests passed.
- The actual AMPL model agreed with an independent connected-component
  calculation in 159 additional cases (102 feasible, 57 infeasible), including
  unselected gaps, reversed edges, cycles, untouched existing fragments,
  variable cell areas, and an attempt to select unavailable land.
- A separate eight-nature-type case confirmed that mixed types can form one
  qualifying network while respecting the original global type quotas.
- The biodiversity objective was checked against your uploaded original and
  is identical after removing whitespace. Production defaults were checked.

These are mathematical and software tests on synthetic data, not a solve on
your actual planning area. The final raster component audit remains important
for checking the adjacency and cell-ID data you supply.
