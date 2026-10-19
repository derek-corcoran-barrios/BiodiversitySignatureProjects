# Version 3 — 7 October 2026

## New in this release

- A local Shiny configuration app with seven panels: area, anchors, biodiversity,
  targets, workflow, projection, and review/export.
- A small synthetic landscape for exploring all the controls without national
  datasets or a solver. The app includes a live explanation of the chosen rules.
- Exact 10 m anchor previews and GeoTIFF download, file checks, complete R config
  export, and JSON export/import for editing the same scenario later.
- A general optional repair stage: enforce zero network shortfall and all budgets,
  then minimize hectares added plus removed from the partial solution. This is
  a mixed-integer solve, so a time limit can expire before it finds a repair.
- Separate warm-up and final spatial weights. Set the final weight to zero to
  optimize biodiversity alone while retaining the hard network requirements.
- Free, capped, or fixed final footprints. Capped changes count additions plus
  removals from a reference that stays constant across final continuation stages.
- Tests for repair, footprint policies, phase weights, app reactions, config
  round trips, and map previews. Valid solver thread options are accepted;
  configurations with network requirements and no anchor sources are rejected.

## Updating from version 2

Replace the code folder and use a new run ID. Existing biodiversity layers can
be reused; no species refitting or reprojection is required by this update.
The bundled example configurations already use v3 IDs. Older hand-written
configurations inherit the shared spatial weight if phase weights are omitted,
and default to enabled repair and a free final footprint.

The new repair solves a different optimization objective from the original
Gudenå constructive path heuristic, which remains in `legacy/`. Its speed on
real landscapes has not been benchmarked. Full-scale Gurobi and saved-model
prediction still need a pilot with the actual data. See `VALIDATION.md`.
