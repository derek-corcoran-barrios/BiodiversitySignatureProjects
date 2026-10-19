# Validation — 7 October 2026 (v3)

Executed successfully with a native R interpreter, terra 1.7.65, igraph 1.6.0,
and actual AMPL/HiGHS solves:

- A separate hard-feasibility repair that adds a connector, removes another
  selected unit to maintain the area/cost budget, and minimizes changed hectares.
- A repair with an empty reference selection and an infeasible repair with no
  usable incumbent. Failed attempts are distinguished from accepted solutions.
- Final free, capped and fixed footprints; a one-hectare cap prevents a
  two-hectare swap, while a two-hectare cap permits it. Fixed locations still
  allow better pathway assignments.
- Independent warm-up and final spatial weights, with comparable final-quality
  reporting across phases.
- Configuration defaults and validation, optional subareas, multiple anchor
  sources, single-model destination rules and projection output path generation.
- Solver thread-option validation and rejection of network requirements without
  any enabled anchor source.
- Exact R configuration export/reload and app settings JSON round trips.
- Actual Shiny server reactions to model/footprint changes and invalid budgets.
- A 24 ha synthetic preview using terra, including individual source masks and
  polygon/subarea overlays.
- A live Shiny 1.8.0 app exercised in Chromium: demo loading, input-file checks,
  a separate background R preview process, map rendering, R/JSON/GeoTIFF
  downloads, single/eight model switching, footprint controls, stale-preview
  rejection, settings re-import, invalid-percentage feedback and mobile layout.
  Desktop and mobile screenshots were visually inspected.

- Complete preparation, warm-up, shortfall minimization, ecological optimization,
  export, and independent 10 m audit for a synthetic single-outcome landscape.
- Completed-stage reuse and rejection of changed scientific settings under the
  same run ID.
- Two subareas plus a Global certificate, with no fixed four-subarea assumption.
- A cost-limited case with no network constraints.
- Eight selectable habitat outcomes and the production-forest destination rule.
- Shoreline anchor-layer creation without crediting lake pixels as nature, and
  a complete solve whose root is selected restoration rather than existing nature.
- A counterexample with disconnected fragments inside one planning unit: the
  solver retains positive shortfall and its bound excludes zero. The unit is not
  treated as a bridge across the intervening nonhabitat pixels.
- A spatial footprint comparison with the expected shared hectares and Jaccard.
- Four real terra projection-input tiles, all eight endpoint transformations,
  unchanged nonselectable land and NA domains, and preserved class dictionaries.
- Tile-by-tile climate interpolation agrees with interpolation over the complete
  synthetic area. A different category dictionary is rejected.
- VRT dependency tracking includes the backing tile files.

The R test scripts are `tests/test_spatial_pipeline.R`,
`tests/test_repair_workflow.R`, `tests/test_config_app.R` and
`tests/test_projection_inputs.R`. Their examples are deliberately small enough
to check exact expectations, rather than to estimate national performance.
The older-terra rook-connectivity fallback was the value-patch implementation
exercised here. The newer terra `patches(values = TRUE)` route was checked
against current official documentation but was not executed in this environment.

Not executed:

- Prediction from Derek's actual saved SpeciesPoolR model objects. Those fits,
  thresholds, national occurrence buffers, tree, and national predictor rasters
  were not supplied. The adapter uses projection/metric helper functions extracted
  from the recovered targets file, but its full real-data integration needs a pilot.
- The complete three-metric reduction from all actual species, the targets/crew
  scheduler, or a parallel multiworker national run.
- A v3 full-scale Gurobi solve, real-landscape feasibility, national runtime or
  disk/memory benchmarking, or equivalence to the earlier Gudenå solution.
- The eight-neighbour alternative. Four-neighbour connectivity was tested.
- The app on Windows/macOS or a multi-user server. Background previews use the
  local Rscript executable. The tested deployment was a local Linux session.
- A proof that the new general repair is faster than the original Gudenå
  constructive heuristic. They use different search strategies; a regional
  pilot is needed to measure runtime and tune stage time limits.

The original package and its original validation notes are included under
`legacy/gudenaa_pipeline`. Its prior checks should not be read as validation of
the new formulation or the national projection adapter.

Before committing server time, run a small representative area with the real
saved-model assets. Inspect scenario and metric audits, species coverage,
the anchor masks, the existing-area ledger and the independent connectivity audit.
Then measure a few tiles to choose worker count and tile size. These scientific
and resource checks cannot be inferred from the successful synthetic tests.
