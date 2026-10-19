# Start with the next area

1. Extract `denmark_pipeline/` into your existing project root.
2. Open the configuration app (below), or copy `configs/area.R` to a config named for the area. Change the boundary,
   optional subareas, anchor sources, budget, output IDs, and input paths.
3. If the biodiversity layers already exist, point `model$outcomes` at them.
   Otherwise check the saved-fit root and predictor manifests in `projection`.
4. Review the configured solver time limits before starting.

```r
source("denmark_pipeline/run.R")

# Optional visual editor. Install shiny and jsonlite once if needed.
# install.packages(c("shiny", "jsonlite"))
launch_config_app()

# After closing the app, use your downloaded configuration's filename below.

# Inspect this output in QGIS first.
anchor_file <- build_anchors("denmark_pipeline/configs/area.R")

# Skip this line if you supplied existing biodiversity layers.
run_projection("denmark_pipeline/configs/area.R")

result <- run_area("denmark_pipeline/configs/area.R")
```

For a region without subareas, set both `region$subareas` and
`connectivity$subarea_min_ha` to NULL. For the independent yes/no model, copy
`configs/single_outcome.R` and supply its one score layer.

For Denmark-wide prediction, use `configs/national_projection.R` after a small
real-data projection pilot. It uses the same saved species models and the same
eight endpoint definitions. It does not fit models or optimize the whole country.

Start with **Load a synthetic example** in the app to explore the controls.
Save both the R configuration (for running) and the settings JSON (for reopening
the form). See `shiny/README.md`. The app is a local configuration editor; long
solver and projection jobs run through the R commands it displays.

The pipeline now runs warm-up → minimum shortfall → repair when needed → final
optimization → independent audit. Repair seeks zero shortfall while minimizing
changed hectares; it is a general optimization stage and may time out. The
final stage allows free locations by default. Choose `capped` or `fixed` in
`quality$footprint` when footprint preservation is intended. Separate warm-up
and final spatial weights are supported; final weight zero means pure biodiversity.

Use new run/projection IDs when migrating from v2 because the implementation
fingerprint has changed. Existing completed outputs remain available.
Existing biodiversity metric layers can be reused by pointing the configuration
at them; updating the optimizer does not require new species predictions.

The default Global connectivity calculation is confined to the planning region.
Review `context_boundary` and `global_area_scope` if external nature patches
should contribute. The original Gudenå formulation remains in `legacy/`.

See `README.md` for the exact anchor and budget semantics, outputs, resuming,
and differences from the original formulation; see `VALIDATION.md` for what
was and was not run here.
