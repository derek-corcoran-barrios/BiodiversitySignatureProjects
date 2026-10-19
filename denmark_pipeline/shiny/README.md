# Shiny configuration editor

From the project root:

```r
install.packages(c("terra", "igraph", "digest", "shiny", "jsonlite")) # once
source("denmark_pipeline/run.R")
launch_config_app()
```

Alternatively: `shiny::runApp("denmark_pipeline/shiny")`. The default project
root is the parent of `denmark_pipeline`. The app runs locally and reads paths
on the machine running R. It does not upload the national data to a service.

## A first pass

1. Click **Load a synthetic example**. This creates a small artificial landscape
   with two subareas, current nature, production forest, a lake, a political
   polygon, a custom anchor raster, eight sets of metrics and an independent
   single-outcome score. These are invented demonstration data.
2. In **Review & export**, click **Check input files**, then **Build map preview**.
   The map can show current nature, combined anchors, individual anchor sources,
   eligibility and existing-area accounting. The exact anchor raster is downloadable.
3. Explore the anchor switches and the two model choices. The sidebar explains
   the hard requirements; the Workflow tab shows the optimization stages.
4. Replace the demonstration paths and values with the next area's inputs.
   Use a new area/run ID and check the network sizes and area budget.
5. Download the **R configuration** and **editable settings JSON**. Put the R file
   in the project root and use the displayed commands. Reopen the JSON to edit
   the same scenario in the app later.

Demo files are temporary and are not research inputs. To retain them after R
exits, copy the displayed demo folder to a permanent location and update paths.

## What the controls mean

- The percentage target is always for the **whole planning polygon**. Subareas
  only add optional local connected-network requirements.
- Sources mark **eligible roots**, not compulsory selections. Existing nature
  remains existing nature when its anchor switch is off.
- Political designations on farmland remain restoration candidates. Lake
  shores use a 10 m root buffer by default. The area-accounting switches can
  additionally treat the source as existing habitat and credit its hectares;
  this is an explicit modeling choice. The buffer itself receives no credit.
- Union accepts any enabled anchor source; intersection requires their overlap.
  Arbitrarily many additional sources can be supplied as named JSON definitions.
- The single-outcome option uses your independent score raster and disables the
  forest/open destination restriction. It never derives a score from the eight
  pathways. Switching back restores the pathway-specific settings.
- Metric weights are normalized. Optional fixed scales and per-outcome score
  paths or allowed masks are available under advanced biodiversity settings.
- Repair enforces the hard constraints and minimizes additions + removals from
  the preceding partial solution. It runs only when necessary and enabled.
- Final location controls refer to the **feasible solution entering that phase**.
  The reference does not drift during subsequent continuation stages.
- A final spatial weight of zero optimizes biodiversity alone, while keeping
  connectivity and all other hard requirements.
- New projections use the selected saved models and original climate source.
  Their VRT output paths automatically populate the eight outcomes. The app
  does not project the independent single-outcome model.

## Checks and execution

The green settings message means the choices are structurally consistent. It
does **not** certify landscape feasibility. File checks verify presence (and
Shapefile/VRT dependencies), not complete scientific compatibility.

The preview uses the real 10 m mask-building functions. A separate R process
builds the preview, so the interface remains responsive. Its outputs live in
a temporary directory. Large rasters are sampled only for display; the
downloaded anchor raster retains the full 10 m resolution. If settings change,
the app marks the old preview as stale and disables its export until rebuilt.
Preview hectare counts precede minimum-unit and outcome-mask filtering. Full
preparation and final audits remain mandatory.

If the preview takes more than five minutes, the app displays its log path.
Closing a browser session does not forcibly stop a preview already in progress.
The app has no button that starts a long solver or species-projection job.
Use `build_anchors()`, `run_projection()` when required, and `run_area()` from R.

The JSON contains form settings only. Uploaded R scripts are not executed.
The R export is a complete configuration list with round-trip numeric precision.
Existing hand-written R configs can still be run directly; arbitrary R-script
import into the form is not supported.

## Additional anchor example

```json
{
  "reserve": {
    "type": "vector",
    "path": "Datasets/reserves.gpkg",
    "filter_field": "status",
    "filter_values": ["designated"],
    "buffer_m": 0
  },
  "additional_roots": {
    "type": "raster",
    "path": "Datasets/other_anchors_10m.tif",
    "values": [1, 2]
  }
}
```

Use names different from the enabled built-in sources (`nature`, `political`,
`lakes`, `custom`). All rasters must match the land-use grid and CRS.
