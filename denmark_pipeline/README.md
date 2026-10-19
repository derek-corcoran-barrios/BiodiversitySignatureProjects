# Denmark restoration pipeline — version 3

Use one configuration per area and model. This version separates the national
species-model assets, local or national tiled projections, and regional
restoration optimization. It includes the original Gudenå package under
`legacy/gudenaa_pipeline/` for reproducing the earlier formulation.

## Start here

Extract `denmark_pipeline` beside your project's `Datasets`, `Predictors`,
`Results`, and `Optimization` folders. Open the Shiny configuration editor, or
edit a copy of `configs/area.R`.

```r
install.packages(c("terra", "igraph", "digest", "shiny", "jsonlite")) # once
source("denmark_pipeline/run.R")
launch_config_app()
```

The app explains the active rules, combines anchor sources, checks input files,
previews the current landscape and anchor masks, and downloads a complete R
configuration. Save the app settings JSON to reopen the form later. It includes
a synthetic example that needs neither national data nor a solver. See
`shiny/README.md` for the app workflow and the exact scope of its checks.

```r
source("denmark_pipeline/run.R")

# Build and inspect the anchor layers before undertaking a long solve.
build_anchors("denmark_pipeline/configs/area.R")

# Only needed when the current + eight pathway metric layers do not yet exist.
run_projection("denmark_pipeline/configs/area.R")

# Uses the metric paths selected in the config; no species-model fitting here.
result <- run_area("denmark_pipeline/configs/area.R")
```

To use existing Gudenå metric TIFFs immediately, replace the outcomes assignment:

```r
config$model$outcomes <- dp_eight_outcomes(
  "Results/gudenaa_plant_biodiversity/full/metrics", "tif")
```

Review input paths, network sizes, connectivity context, and solver times before
launching. Set `config$solver$execute <- FALSE` to prepare the first AMPL run
without launching it. Sourcing `run.R` alone launches nothing.

Optimization needs R packages `terra`, `igraph`, `digest`, and AMPL with a Gurobi
or HiGHS MP driver on the executable search path. Projection additionally needs
your tested SpeciesPoolR installation, `maxnet`, `dplyr`, `readr`, `ape`, and
`Matrix`. The optional targets entry point also needs `targets` and `crew`.
Current terra is recommended; a rook-connectivity fallback for older terra is
included but can use more time and memory when building the graph.

## Configuration choices

| Component | Configuration | Meaning |
|---|---|---|
| Planning region | `region$boundary` | One polygon or a set of polygons whose union defines the budget boundary. |
| Subareas | `region$subareas = NULL` or a vector file | Any number; no requirement for four. |
| Subarea IDs | `region$subarea_id` | A unique attribute, or row-based IDs if NULL. |
| Subarea networks | `connectivity$subarea_min_ha` | NULL, one common size, or a named vector indexed by subarea ID. |
| Global network | `connectivity$global_min_ha` | A connected-area requirement, or NULL to omit it. |
| Decision resolution | `grid$planning_unit_m` | Normally 200 m; must be a multiple of 10 m. |
| Anchor sources | `anchors$sources` | Any combination of current nature, vector polygons, or a custom raster. |
| Outcome model | `model$mode` and `model$outcomes` | Eight pathway inputs, or one independently supplied restoration score. |
| Nature target | `budget$target_nature_fraction` | A fraction of the whole planning polygon, or NULL. |
| Restoration budget | `budget$min_restoration_ha`, `max_restoration_ha` | Absolute new-restoration bounds, additional to any nature target. |
| Monetary budget | `budget$cost_cap`, `cost_per_ha` | Optional currency cap and aligned raster in currency per hectare. |
| Phase-specific spatial weights | `objective$warmup_contiguity`, `quality_contiguity` | Separate warm-up and final bonuses; NULL inherits `contiguity`. |
| General repair | `repair$enabled`, `seconds`, `gap` | Attempt hard feasibility while minimizing changed hectares from the partial solution. |
| Final footprint | `quality$footprint`, `max_change_ha` | `free`, `capped`, or `fixed`, relative to the first solution entering the final phase. |
| Projection extent | `projection$boundary` | NULL means the planning area; supply a Denmark boundary for national projection. |
| Parallelism | `projection$workers`, `tile_cells` | One tile per worker; species and scenarios run sequentially within a tile. |

The 30% value is editable. Subareas do not receive separate percentage budgets.
Their connected-network requirements are independent of the whole-polygon
budget. Subareas must not overlap or extend beyond the region; they need not
cover every part of the region. Set both network minima to NULL to optimize
without hard connected-network requirements.

When a nature fraction is supplied, the minimum restoration is
`max(0, fraction * polygon_area - credited_current_nature)`. The default
`target_tolerance_ha = 16` supplies an upper bound on additional restoration.
Set it to NULL for a minimum-only nature target. Existing nature above the
target is retained. Monetary caps do not force the budget to be completely spent.

## Make your own 10 m anchors

The template is the national **modelled categorical land-use raster** used in
the scenario workflow, with its fitted-model class dictionary. The supplied
example points to `Landuse_Denmark_current.tif`. Do not substitute a raster with
different class meanings just because its integer codes happen to match.

```r
config$anchors$sources <- list(
  nature = list(type = "nature"),
  political = list(type = "vector", path = "Datasets/political_anchors.gpkg"),
  lakes = list(type = "vector", path = "Datasets/lakes.gpkg", buffer_m = 10)
)
config$anchors$combine <- "union"
```

`union` means any chosen source can supply a root; `intersection` means a root
must lie in all the source masks. Vector inputs are transformed to the grid CRS
and rasterized by pixel centre. Raster inputs must already match the grid CRS,
resolution, and alignment. Positive raster values are anchors unless an explicit
`values` vector is provided; NA is treated as outside that source.

For a hand-built layer, use:

```r
config$anchors$sources <- list(
  custom = list(type = "raster", path = "Datasets/my_anchors_10m.tif", values = 1)
)
```

Outputs include `anchors_10m.tif`, separate named source bands in
`anchor_sources_10m.tif`, and separate existing-habitat, budget-credit,
network-credit, and eligibility rasters. This keeps source overlap from being
double-counted. A designated polygon supplies eligible roots, not a requirement
to select every pixel or use every designated anchor.

Political farmland remains selectable land. It must actually be selected before
it can act as a root; designating it does not create free existing habitat.
The lake example marks land within one pixel of the shoreline as eligible to
root a network. Water itself is not automatically a terrestrial corridor and
does not automatically contribute to the nature or network area totals.

If a source really represents additional existing, traversable habitat, its
options can explicitly include:

```r
list(type = "raster", path = "Datasets/additional_existing_habitat.tif",
     values = 1, include_as_existing = TRUE,
     credit_budget = TRUE, credit_network = TRUE)
```

These flags assert a modelling interpretation; they do not establish ecological
dispersal through the habitat. Root buffers never enlarge credited area.

## One outcome or eight

`configs/single_outcome.R` is a separate run configuration. It consumes your
independent model's one score raster and makes a restore/do-not-restore decision.
It never derives that score by taking the best or average of the eight pathways.
It uses the same general budget and connectivity machinery, so changing the
underlying scientific benefit model is explicit.

The pathway configuration uses the eight existing destinations: dry/wet ×
poor/rich × forest/open. Optional `allowed_mask` rasters restrict individual
outcomes. An outcome is allowed in a planning unit only if **all of that unit's
selectable pixels** pass its mask. This avoids claiming a whole unit's hectares
when only a fraction can receive the destination. Units with no allowed outcome
are removed, and capacity is recalculated.

The original production-forest rule is available: units with at least 50% of
their selectable area in production forest can only receive forest outcomes.
For an independent single outcome, explicitly disable that destination rule,
as the supplied model has no forest/open choice. You can still restrict its
eligibility with a mask.

Scores are integrated over selectable hectares. `score_type = "endpoint"`
uses destination scores; `"gain"` subtracts supplied `current_scores`.
Automatic metric scales are shared across the outcomes in a run. Supply fixed
`metric_scales` for comparisons requiring identical metric scaling across
areas or different sets of outcomes. The default spatial preference is normalized
same-outcome adjacency plus contact with eligible anchors.

Set `objective$quality_contiguity <- 0` for a biodiversity-only final objective;
hard network requirements still apply. A high warm-up bonus can independently
encourage compact starting patterns. The weights are preferences, not evidence
that any network threshold is satisfied.

Compare validated footprints with:

```r
compare_footprints(
  "PATH_TO_PATHWAY_FINAL/selected_10m.tif",
  "PATH_TO_SINGLE_FINAL/selected_10m.tif",
  output_csv = "footprint_comparison.csv")
```

This reports shared and unique hectares and area-weighted Jaccard overlap.
Shared unselected pixels do not increase agreement. Matching geometry and domain
are checked. Independently optimized alternatives can differ even when their
scores are close, so inspect solver gaps alongside overlap.

## Projection: reuse the national fits

The adapter is grounded in the recovered `targets_plant_biodiversity.R` and
`02_build_gudenaa_scenarios.Rmd`. It expects the saved fitted-model objects,
national status files, thresholds, occurrence buffers, range weights recorded
in those statuses, `phylogeny/plant_tree.rds`, and `audits/tree_coverage.csv`.
The adapter does not download occurrences, fit models, recalibrate thresholds,
recompute national rarity weights, or build a new phylogeny.

The training predictor manifest validates class IDs/names and climate band names.
Projection uses the **original climate source** from
`Predictors/current_predictor_paths.rds`, or the explicit `projection$climate`
path. Each tile is projected directly from that full climate source using the
same bilinear interpolation as the previous scenario builder. Interpolating
climate to 10 m does not add independent 10 m climate information.

For each tile, the code creates current conditions and eight counterfactual
land-use rasters. Every selectable agriculture/production-forest pixel is
changed to the relevant endpoint, while nonselectable pixels and climate remain
unchanged. The optimizer subsequently chooses only the permitted destinations
and locations. These are potential endpoint surfaces, not eight simultaneous
restoration plans or predictions of restoration success/time.

Every scenario uses the same usable national species set, fitted models,
thresholds, occurrence buffers, inverse national range-area weights, and tree.
The masks use the saved national buffers, including sources outside the tile.
Missing predictions inside a species' accessible domain fail instead of becoming
absences. Missing predictor coverage also fails explicitly.

Each scenario produces Richness, PD and RangeRarity: **27 surfaces including the
three current baselines**. Individual species rasters can be removed after a
scenario's summaries pass the inherited arithmetic audits. Tile metrics remain
and are assembled as VRTs; keep the backing tile folders with those VRTs.
Completed tile/scenario metrics have checksummed resume markers.

Start with the next area you actually need. Use a small representative tile set
to measure wall time, RAM, and disk use, then increase workers. `tile_cells = 1000`
means a 10 × 10 km tile. With eight workers, eight tiles are active, with no
nested species-level worker pool. Terra's memory setting is a hint, not a hard
process memory reservation. Allow RAM for models, buffers, PD calculations and
the parent process as well as raster blocks.

For national projection, edit `configs/national_projection.R` with the actual
Denmark boundary and run:

```r
run_projection("denmark_pipeline/configs/national_projection.R")
```

It projects the saved models across Denmark; it does not launch a national
optimization. Local optimizers can read the resulting national metric VRTs.
More cores help until RAM or storage throughput becomes limiting; no national
runtime claim is made without your fitted models and server benchmark.

If you prefer targets:

```r
Sys.setenv(DENMARK_CONFIG = "denmark_pipeline/configs/area.R")
targets::tar_make(script = "denmark_pipeline/_targets.R",
                  store = "_targets_denmark_projection_area1")
```

Use a separate store for each projection configuration. This entry point handles
projection only. Run `run_area()` after the desired metric layers exist. Run all
commands from the project root.

## Connectivity, migration and resume

The generalized model keeps one decision per planning unit, but splits its selectable
pixels into their truly connected 10 m fragments. Selecting one unit therefore
cannot bridge a road or other nonhabitat gap between its disconnected fragments.
Certificates are connected to an eligible root; certificates can overlap.
They do not require all selected restoration to belong to a certificate.
Four-neighbour connectivity is the default. Eight-neighbour connectivity is an
explicit alternative and needs current terra for value-based patches.

The default connectivity context is the planning polygon. Subarea networks must
connect **within their own subarea**. A wider `context_boundary` can supply
external existing habitat to the Global network. `global_area_scope` determines
whether its credited existing hectares are restricted to the planning region or
include that context. A clipped context never represents full national patch
area. To preserve that particular Gudenå policy exactly, continue using the
bundled legacy formulation or provide a context containing the full relevant
patches and choose `within_context`.

The solve sequence is:

1. **Warm-up:** biodiversity + the warm-up spatial bonus, within the area/cost
   budgets. Network requirements are relaxed.
2. **Minimum shortfall:** minimize the sum of network deficits divided by their
   respective requirements. The warm-up is a starting solution.
3. **Repair, if needed and enabled:** fix every shortfall to zero and minimize
   the hectares added plus the hectares removed from the partial solution.
   All network, area, cost, eligibility and destination constraints apply.
   This is a general mixed-integer repair optimization, not the older
   Gudenå-specific path-joining heuristic. Any number of networks is supported.
   Repair is skipped when shortfall is already zero or a solver bound has
   excluded zero shortfall. A time-limited failed repair retains the partial
   solution and proceeds to the next configured shortfall stage, if any.
4. **Final optimization:** maximize biodiversity plus the final spatial bonus
   while enforcing every network requirement. Use the first feasible solution
   as a starting point, with the footprint policy described below.
5. **Independent 10 m audit:** check actual connected components, area/cost
   budgets, legal outcomes and any footprint change limit before acceptance.

The repair objective is
`sum(area of newly selected units) + sum(area of removed units)` relative to the
partial incumbent immediately preceding that attempt. A feasible time-limited
repair is usable; its best bound and incumbent objective are exported. It is
only proven to minimize changes when the solver closes that objective gap.
The repair stage can be expensive and does not guarantee a feasible solution
within its time budget. Full-scale Gurobi performance still needs a pilot.

The final footprint policy refers to **the first feasible solution entering the
final phase**, saved as `quality_reference.csv`. It remains the same across all
final continuation stages:

| `quality$footprint` | Final-stage behavior |
|---|---|
| `free` | Locations and outcomes may change while meeting all hard constraints. Default. |
| `capped` | Additions + removals cannot exceed `quality$max_change_ha`. |
| `fixed` | Selected/unselected locations are fixed; allowed pathway choices remain free. |

Swapping 100 ha for 100 different hectares counts as 200 changed hectares.
With a fixed footprint and one independent outcome, the biodiversity decision
is already determined. Footprint protection is an explicit additional policy,
not an implied consequence of using a warm start.

Example controls:

```r
config$objective$warmup_contiguity <- 1.2
config$objective$quality_contiguity <- 0  # Pure biodiversity in final stage
config$repair <- list(enabled = TRUE, seconds = 1800, gap = 0.005)
config$quality <- list(footprint = "free", max_change_ha = NULL)
# Or: list(footprint = "capped", max_change_ha = 100)
# Or: list(footprint = "fixed", max_change_ha = NULL)
```

The original constructive heuristic and its priority for preserving the
original warm-up footprint remain in `legacy/`. The generalized model also
uses a simpler anchor-contact preference; it does not reproduce the legacy
bonus for joining several existing patches. It is **not a promise of identical
Gudenå results**.

Every solve exports selected actions, network membership, solver statistics,
and all incumbent variable values. A time limit is not automatically called
infeasibility. Positive-shortfall results are retained; a positive lower bound
excluding zero is reported separately. Extend the feasibility or quality time
vectors, or increase the repair time, to add search effort. Completed matching
stages are reused. Checkpoints restart from incumbent values,
not from the old branch-and-bound tree. A changed scientific configuration,
input content, or implementation requires a new run/projection ID.

`latest_validated_solution.txt` is written only after budget checks, legal-action
checks, and the independent 10 m network audit pass. Numerical outputs include
selected cells, selected-pixel raster, network and area audits, stage summaries,
and a restored categorical land-use raster for pathway runs. The original six
Gudenå map layouts remain in the legacy package; the generalized outputs are
ordinary rasters/tables suitable for R or QGIS.

This archive contains code and examples, not the national rasters, fitted species
models, or your independent single-outcome model. See `VALIDATION.md` for the
tests actually executed and the remaining full-data checks.
