# Gudenå eight-habitat pipeline

This package assembles the supplied 07, global-30% 08a, constructive 08b and 09
work into one R entry point. The eight habitat choices remain present in every
optimization stage. The separate yes/no model is not used as an initialization
or substituted for any of these stages.

## Run

1. Extract this folder into the `BiodiversitySignatureProjects` project root,
   beside `Datasets`, `Predictors`, `Results` and `Optimization`.
2. Review `config.R`, particularly input paths, AMPL executable, run ID,
   time limits, and `quality_mode`.
3. Run from the project root:

```r
source("gudenaa_pipeline/run_pipeline.R")
result <- run_gudenaa_pipeline("gudenaa_pipeline/config.R")
```

Or from a terminal:

```bash
Rscript gudenaa_pipeline/run_pipeline.R
```

Required R packages: terra, dplyr, tidyr, purrr, readr, stringr, digest, igraph,
and knitr. Mapping additionally uses ggplot2, tidyterra, xml2 and ragg, plus the
existing `NaturetypeStyle.sld`. A working AMPL installation with the current
Gurobi MP driver and an appropriate license is needed for the real data.
Input rasters, national nature patches, the boundary shapefile, and biodiversity
metric rasters are project inputs and are not bundled here.

Set `execute_solver = FALSE` to prepare the next run file and stop before
launching AMPL. No solver job is launched by extracting or inspecting the files.

## Policy preserved

- Current nature plus new restoration reaches 30% **over the whole polygon**,
  with the existing 16 ha upper tolerance. No subarea-specific 30% constraint
  is imposed, including in the warm-up. Subarea 30% columns are reference values.
- One Global certificate reaches 5,000 ha and each of the four subarea
  certificates reaches 1,000 ha inside its respective subarea.
- Certificates may overlap, include new restoration, and are rooted in existing
  nature. Different habitats may connect. These requirements do not force all
  selected restoration into one component or into the five certificates.
- Global certificate area includes the full national area of each reached
  existing patch, counted once. Subarea certificates and the polygon's 30%
  ledger count only hectares inside the relevant boundary.
- Production-forest-dominated cells (at least 50% of eligible area) receive only
  forest destinations. The decision restores eligible 10 m hectares within a
  200 m cell; it does not assume every selected unit contributes 4 ha.
- The configured area reference is 44,180.77 ha, as in the attached notebooks.
  Actual targets are calculated from the shapefile. Verify the reference if
  using another boundary; the older 73,123.09 ha polygon is not silently assumed.

## Stage progression

1. Prepare the global-budget warm-up, solve it and export selected cells.
2. Prepare the patch graph, exact area ledgers and network model. If the warm-up
   already meets all certificates, construct its full start; otherwise use a
   valid start with explicit shortfalls.
3. Minimize shortfall. Every accepted incumbent is exported, including partial
   network membership and a complete variable-value checkpoint.
4. When needed, try the supplied constructive path-and-area repair. Its current
   heuristic applies when the four subarea certificates are complete and the
   Global certificate remains deficient. If that condition does not hold, or
   repair fails, preserve the incumbent and use the configured general MIP
   continuation. Repair failure is not reported as mathematical infeasibility.
   A positive solver lower bound that excludes the zero-shortfall tolerance is
   reported separately as `ZERO_SHORTFALL_RULED_OUT_BY_BOUND`.
5. Once shortfall is zero, optionally minimize change from the warm-up footprint.
6. Optimize ecological/spatial quality with all network shortfalls fixed to zero.
   `quality_seconds = c(28800, 14400)` means two ecological solves of 8 h and 4 h,
   with a checkpoint between them. These are newly specified continuation stages,
   not a reconstruction of the unattached historical continuation script.
7. Independently reconstruct and audit at 10 m. Only a passing solution updates
   `latest_validated_solution.txt`. A failed fine audit stops with
   `COARSE_FEASIBLE_FINE_REFINEMENT_REQUIRED`; targeted spatial refinement is not
   automated in this version. Maps are generated only after validation.

`quality_mode = "minimum_change"` preserves the original priority: minimize
changes from the warm-up, then limit later footprint change to that incumbent
plus 0.01 ha. A time-limited minimum-change result is not a proof of minimum
change. `quality_mode = "free"` skips that priority and lets the final solver
change locations freely while respecting all ecological constraints. Changing
this scientific choice requires a new run ID. Neither option fixes the initial
selected cells in place.

The final objective remains the original weighted biodiversity plus spatial
contiguity and patch-connection reward. CCI, habitat SDM masks and agricultural
history are not invented or applied without the relevant input layers and
agreed definitions.

## Resume and outputs

Results go under `Optimization/gudenaa_pipeline/<run_id>/`. Existing notebook
output folders are read-only cache sources and are never overwritten.

- `manifest.rds`: configuration plus input and code hashes. Changed data/code or
  scientific settings require a new run ID.
- `warmup/`: shared candidate dictionary, cell/subarea areas, graph and warm-up.
- `network/`: prepared network data and construction diagnostics.
- `stages/<stage_and_hash>/`: run file, AMPL/Gurobi logs, solution tables, solver
  statistics and `checkpoint.run`. A checkpoint reloads incumbent values; it
  does not resume the previous branch-and-bound tree.
- `stage_summary.tsv`: objective, bound, status, shortfall and footprint change
  for every stage used in the current pipeline invocation.
- `final/<selected_stage_and_hash>/`: chosen incumbent, independent audits and
  maps. Distinct accepted stage results have distinct final folders.
- `latest_validated_solution.txt`: the exact prefix to use for downstream work.
- `grid_specification.rds` and `grid_id.txt`: grid, domain and candidate-area
  information for aligning independent models.

Rerun the same command after interruption. Completed stages are reused after
signature checks; an interrupted solve restarts from its previous completed
incumbent. It does not recover an unexported solution from a killed solver.
Time budgets may be extended: append to `feasibility_continuation_seconds` or
`quality_seconds`. Unchanged completed stage hashes are reused. A completed
warm-up is reused; choose a new run ID to intentionally rerun the warm-up.

A `RUNNING.lock` directory prevents concurrent writes. After an abrupt crash,
confirm no R/AMPL job for this run is active before removing a leftover lock.
An unavailable map package or SLD, or a map error, leaves the validated numerical
solution available. Diagnostic files describe map failures separately.

## Compare the separate yes/no model

Solve and validate the binary model independently using its own model and input
score. Both comparisons must refer to the same grid, candidate domain and
restorable hectares. Export **selected cells only**, with `coarse_cell_id` and,
preferably, `restored_ha`. Filter an all-cell decision table to selected rows first.

```r
source("gudenaa_pipeline/R/compare_footprints.R")
root <- "Optimization/gudenaa_pipeline/global30_eight_habitats_v1"
prefix <- readLines(file.path(root, "latest_validated_solution.txt"))
eight_grid_id <- readLines(file.path(root, "grid_id.txt"))
# The independent model should record its verified shared grid identifier.
binary_grid_id <- readLines("PATH_TO_BINARY_RUN/grid_id.txt")

comparison <- compare_gudenaa_footprints(
  eight_habitat_selected = paste0(prefix, "_selected_cells.tsv"),
  binary_selected = "PATH_TO_BINARY_RUN/selected_cells.tsv",
  candidate_dictionary = file.path(root, "warmup/candidate_cell_dictionary.csv"),
  eight_grid_id = eight_grid_id,
  binary_grid_id = binary_grid_id,
  subarea_dictionary = file.path(root, "warmup/candidate_area_by_subarea.csv"),
  output_dir = file.path(root, "comparison_with_binary")
)
comparison$summary
```

Do not copy a grid ID into another run without verifying its spatial definition.
If that model uses a different raster extent or fine-scale restoration mask,
align the spatial footprints first; identical integer cell IDs alone are not
evidence of alignment. When `restored_ha` is supplied, the comparison rejects
incompatible within-cell areas.

The comparison reports hectares selected by both, hectares selected by each
model only, area-weighted Jaccard overlap, and the fraction of each footprint
recovered by the other. It also exports cell-level agreement for mapping and
subarea summaries. Shared unselected cells do not inflate the overlap metric.
High agreement does not show the ecological objectives are equivalent; low
agreement may reflect their different objectives or remaining solver gaps.

## Verification and limitations

The included tests exercise the global budget (including an intentionally
uneven subarea allocation), production-forest habitat restriction, all five
network certificates, positive shortfall when area is insufficient, checkpoint
round-tripping, and area-weighted comparison checks. See `VALIDATION.md` for the
checks actually executed for this package.

The spatial preparation, constructive heuristic and independent audit were
extracted from the supplied scripts rather than rederived. Full-scale performance
and fine-resolution feasibility require a run against the real project data.
The independent binary model was not supplied and is not implemented here.
