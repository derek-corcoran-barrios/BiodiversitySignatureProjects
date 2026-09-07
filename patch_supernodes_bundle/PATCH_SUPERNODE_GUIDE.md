# Using precomputed natural-patch areas in the AMPL restoration model

## Recommendation

Use the already-calculated **connected-patch IDs and hectares** in R, and pass
each current natural patch to AMPL as one supernode. Keep only Agriculture and
ProductionForest cells as 10 m decision nodes.

This is the right division of work:

| Task | R / TroublemakeR | AMPL |
|---|---:|---:|
|Determine which 10 m cells are current nature or restorable|Yes|No|
|Identify rook-connected current natural patches|Yes, once|No|
|Calculate each patch's hectares|Yes, once|No|
|Build candidate-to-candidate and candidate-to-patch edges|Yes|No|
|Choose cells and their future nature types|No|Yes|
|Ensure every selected network totals at least 5000 ha|No|Yes|

At 10 m resolution, one cell is 100 m2 = 0.01 ha, so 5000 ha is indeed
500000 cells. The revised model does **not** instantiate 500000 variables for a
5000 ha existing patch: it supplies that patch's precomputed hectares once.

## Two different meanings of "area"

The uploaded area workflow produces useful but different objects:

1. A stable connected-patch ID plus one exact area in hectares per patch.
   These are valid inputs to the **hard 5000 ha connected-network rule**.
2. `LogAreaRaw1000.tif`, `IDWArea_Distance`, and the smoothed influence layer.
   These summarize nearby habitat context and are valid optional **objective
   scores**. They cannot prove that a selected network is connected or 5000 ha.

Do not use polygon area values themselves as patch IDs. Two disconnected
patches can have the same area. Rasterize a stable ID and retain hectares in a
separate lookup table.

## Why this is much smaller than the old README workflow

The old workflow calls `define_cells()` and `find_connections()` on
`NumNatureTemp`, which contains both current nature and agriculture. Thus every
current-nature 10 m cell and every adjacency enters the MIP.

The new workflow uses:

- `Cells`: only eligible Agriculture and ProductionForest cells;
- `E`: only rook edges between those candidate cells;
- `ExistingPatches`: one element per touching current-natural patch;
- `CellPatchEdges`: only shared-edge contacts between candidates and patches;
- `ExistingPatchAreaHa`: the precomputed hectares of each patch.

The model can still route connectivity through a natural patch. If two selected
cells touch different sides of the same genuinely connected patch, the patch
supernode connects them. Its area is counted once, irrespective of the number
of contacts.

## Files

- `FreeNoCarbon_patch_supernodes.mod`: production AMPL model.
- `patch_supernodes.R`: standalone preparation and writer functions suitable
  for later inclusion in TroublemakeR.
- `FreeNoCarbon_patch_supernodes_demo.dat`: tiny two-network test.
- `FreeNoCarbon_patch_supernodes_demo.run`: AMPL test runner.

## One-time patch catalog

If your existing polygon dataset already has a stable connected-patch ID and
area in hectares, reuse it. Rasterize the **ID**, not the area:

```r
library(terra)
source("patch_supernodes.R")

template_10m <- rast(
  "Predictors/gudenaa_scenarios/Landuse_gudenaa_projection_current.tif"
)

natural_patches <- vect("Datasets/current_nature_connected_patches.gpkg")

# Example expected polygon attributes:
#   patch_id  positive integer, one value per connected patch
#   area_ha   full connected-patch area in hectares
patch_id_10m <- rasterize_precomputed_patch_ids(
  patches = natural_patches,
  template = template_10m,
  patch_id_col = "patch_id",
  filename = "Optimization_inputs/CurrentNature_patch_id_10m.tif",
  overwrite = FALSE
)

patch_areas <- as.data.frame(natural_patches)[, c("patch_id", "area_ha")]
patch_areas <- unique(patch_areas)
stopifnot(!anyDuplicated(patch_areas$patch_id))

readr::write_csv(
  patch_areas,
  "Optimization_inputs/CurrentNature_patch_area_ha.csv"
)
```

The polygon IDs from the area-weighted-distance workflow are reusable only if
each ID is one complete rook-connected patch. Disaggregating multipart features
does not necessarily merge adjacent source polygons. If that condition is
uncertain, generate the catalog once from the final 10 m current-land-use grid:

```r
current <- rast(
  "Predictors/gudenaa_scenarios/Landuse_gudenaa_projection_current.tif"
)

# Classes 1:8 are the eight current nature types. Binarizing first lets
# different but touching nature types belong to the same ecological network.
current_nature <- ifel(as.int(current) %in% 1:8, 1, NA)

patch_id_10m <- patches(
  current_nature,
  directions = 4,
  zeroAsNA = TRUE
)
names(patch_id_10m) <- "patch_id"

writeRaster(
  patch_id_10m,
  "Optimization_inputs/CurrentNature_patch_id_10m.tif",
  overwrite = FALSE,
  datatype = "INT4S",
  gdal = c("COMPRESS=DEFLATE")
)

patch_areas <- freq(patch_id_10m, value = TRUE)
names(patch_areas)[names(patch_areas) == "value"] <- "patch_id"

cell_area_ha <- prod(res(patch_id_10m)) / 10000
stopifnot(isTRUE(all.equal(cell_area_ha, 0.01)))
patch_areas$area_ha <- patch_areas$count * cell_area_ha
patch_areas <- patch_areas[, c("patch_id", "area_ha")]

readr::write_csv(
  patch_areas,
  "Optimization_inputs/CurrentNature_patch_area_ha.csv"
)
```

That connected-component calculation is a one-time checkpoint, not something
to repeat in every regional optimization.

## Important boundary choice

Choose and document one of these meanings before solving:

| Area supplied for a boundary-crossing patch | Interpretation |
|---|---|
|Whole connected patch in the wider landscape|Recommended when ecological connectivity can continue beyond the optimization polygon|
|Only the part inside the optimization polygon|Use when the policy requires all 5000 ha to lie inside the optimization polygon|

For the first interpretation, identify and measure patches on the wider
landscape first, then crop the ID raster to each optimization area while
retaining the full area lookup. For the second, clip first and regenerate patch
IDs and areas inside each area. Do not mix the two interpretations silently.

## Per-area preparation and `.dat` export

The following replaces the old full-domain `define_cells(NumNatureTemp)` and
`find_connections(NumNatureTemp)` calls.

```r
library(terra)
library(TroublemakeR)
source("patch_supernodes.R")

problem_stem <- "RegionProblems/gudenaa_test"
dat_path <- paste0(problem_stem, ".dat")

current <- rast(
  "Predictors/gudenaa_scenarios/Landuse_gudenaa_projection_current.tif"
)
patch_id <- rast("Optimization_inputs/CurrentNature_patch_id_10m.tif")
patch_areas <- readr::read_csv(
  "Optimization_inputs/CurrentNature_patch_area_ha.csv",
  show_col_types = FALSE
)

# If optimizing a subpolygon, crop and mask current and patch_id to the same
# aligned template here. Never bilinearly interpolate class or ID rasters.

candidate <- make_candidate_domain(
  current_landuse = current,
  eligible_values = c(9L, 11L), # Agriculture and ProductionForest
  filename = "Optimization_inputs/Candidate_gudenaa_10m.tif",
  overwrite = FALSE
)

nature_lookup <- data.frame(
  value = 1:8,
  landuse = c(
    "ForestDryPoor", "ForestWetPoor",
    "ForestDryRich", "ForestWetRich",
    "OpenDryPoor", "OpenWetPoor",
    "OpenDryRich", "OpenWetRich"
  )
)

network <- prepare_patch_supernodes(
  candidate_domain = candidate,
  patch_id = patch_id,
  patch_areas = patch_areas,
  current_landuse = as.int(current),
  class_lookup = nature_lookup,
  class_value_col = "value",
  class_name_col = "landuse"
)

network_audit <- audit_patch_supernodes(
  network,
  min_patch_area_ha = 5000
)
print(network_audit)
readr::write_csv(
  network_audit,
  "Optimization_inputs/gudenaa_patch_network_audit.csv"
)

# TroublemakeR writes only eligible cells because candidate is NA elsewhere.
TroublemakeR::define_cells(
  Rasterdomain = candidate,
  name = problem_stem
)

nature_names <- nature_lookup$landuse
TroublemakeR::landuse_names(
  landuses = nature_names,
  name = problem_stem
)
TroublemakeR::write_ampl_lines(
  paste("set NatureLanduses :=", paste(nature_names, collapse = " ")),
  name = problem_stem
)

# This writes E, ExistingPatches, CellPatchEdges, ExistingPatchAreaHa,
# and optional ExistingNatureContacts. Do not also call find_connections().
write_patch_supernodes_dat(
  network,
  dat_path = dat_path,
  append = TRUE,
  write_cells = FALSE
)
```

### Biodiversity input

Prepare eight aligned biodiversity-benefit layers, one per destination nature
type, and mask them to the candidate domain before writing. This preserves your
combined normalized Richness + PD + Rarity definition exactly as requested.

```r
biodiversity_paths <- c(
  ForestDryPoor = "Results/.../Biodiversity_ForestDryPoor.tif",
  ForestWetPoor = "Results/.../Biodiversity_ForestWetPoor.tif",
  ForestDryRich = "Results/.../Biodiversity_ForestDryRich.tif",
  ForestWetRich = "Results/.../Biodiversity_ForestWetRich.tif",
  OpenDryPoor = "Results/.../Biodiversity_OpenDryPoor.tif",
  OpenWetPoor = "Results/.../Biodiversity_OpenWetPoor.tif",
  OpenDryRich = "Results/.../Biodiversity_OpenDryRich.tif",
  OpenWetRich = "Results/.../Biodiversity_OpenWetRich.tif"
)
biodiversity <- rast(unname(biodiversity_paths))
names(biodiversity) <- names(biodiversity_paths)

stopifnot(compareGeom(biodiversity, candidate, lyrs = FALSE))
biodiversity <- mask(biodiversity, candidate)

TroublemakeR::species_suitability(
  Rastercurrent = biodiversity,
  species_names = names(biodiversity),
  parameter = "BioDiversity",
  name = problem_stem
)
```

If one scenario's biodiversity benefit is represented by more than one metric
raster, combine and scale those rasters before this step. AMPL receives the
final score only.

### Optional area-influence score

`IDWArea_Distance` may be included as an objective covariate:

```r
area_influence <- rast("AreaLayers/IDWArea_Distance.tif")
area_influence <- resample(area_influence, candidate, method = "bilinear")
area_influence <- mask(area_influence, candidate)

TroublemakeR::write_cell_param(
  Rasterparam = area_influence,
  parameter = "AreaInfluence",
  name = problem_stem,
  default = 0
)

TroublemakeR::write_ampl_lines(
  "param AreaInfluenceWeight := 100",
  name = problem_stem
)
```

Choose the weight on the same numerical scale as biodiversity. If biodiversity
is 0--1000 and influence is 0--1, a weight of 1 makes influence almost
irrelevant; this is a policy choice, not a connectivity requirement.

### Reusing pre-counted current areas by nature type

Planning-area totals can be passed without retaining current cells:

```r
existing_area <- data.frame(
  landuse = nature_names,
  area_ha = c(...) # totals for THIS planning-area interpretation
)

lines <- c(
  "param ExistingLanduseAreaHa :=",
  sprintf("  %s %.10g", existing_area$landuse, existing_area$area_ha),
  ";"
)
dat_connection <- file(dat_path, open = "a", encoding = "UTF-8")
writeLines(lines, dat_connection, sep = "\n", useBytes = TRUE)
close(dat_connection)
```

The national class totals in the biodiversity-potential document are useful
for a national audit, but must not be inserted as Gudenå or municipal totals.
`MinFinalLanduseAreaHa[l]` means:

```text
pre-counted current area of type l + newly restored area of type l
    >= desired final area of type l
```

This is different from `MinPatchAreaHa`, which applies to each connected final
network regardless of nature type.

### Remaining scalar settings

For a 10 m run with a cell-count budget:

```r
TroublemakeR::write_ampl_lines(
  c(
    "param TransitionCost default 1",
    "param b := 500000",
    "param ExactBudget := 1",
    "param MinPatchAreaHa := 5000",
    "param MinLan := 0",
    "param SpatialContiguityBonus := 0"
  ),
  name = problem_stem
)
```

`b := 500000` means select exactly 500000 cells only when every transition
cost is one and `ExactBudget := 1`. It does not itself impose a 5000 ha minimum
per network; `MinPatchAreaHa` does that.

Start with `SpatialContiguityBonus := 0`. The hard 5000 ha network rule remains
active, while the model avoids creating an additional binary variable for every
nature-type/edge combination. Turn the soft same-type bonus on only after the
hard-rule model is solving acceptably.

## What to remove from the old README loop

For this model, remove or replace these parts:

| Old step | New treatment |
|---|---|
|`define_cells(Rasterdomain = NumNatureTemp)`|Run it on candidate-only Agriculture + ProductionForest raster|
|`find_connections(Rasterdomain = NumNatureTemp)`|Do not call; `write_patch_supernodes_dat()` writes candidate-only `E`|
|`generate_dummy_stack(NumNature)` and `Existingnature` table|Not needed for hard connectivity; patch supernodes replace individual current-nature cells|
|Agriculture as a destination land use (`Ag`)|Do not put it in `NatureLanduses`; restoration decisions have only the eight nature types|
|Full current-nature raster in the MIP|Keep in R only; export patch ID, patch area, and boundary contacts|
|`PotentialCarbon` in `FreeNoCarbon` data|Do not emit it unless using a carbon-enabled `.mod` that declares it|

Keep the regional loop, the candidate masks, the eight biodiversity layers,
transition costs, budget, and any policy-level land-use targets.

## Model behavior to confirm before production

1. **Connectivity rule:** rook/shared-edge, not queen/corner contact.
2. **Cross-type connectivity:** all eight nature types may form one network.
   Type-specific networks would require a different formulation.
3. **Multiple networks:** allowed, but each restoration-containing network must
   reach 5000 ha.
4. **Current patch already over 5000 ha:** a selected cell touching it can
   qualify through that patch. This is usually the desired interpretation.
5. **Small existing-only patch:** it may remain small if no restoration is
   attached. Existing nature is not removed.
6. **Boundary patches:** use either whole-landscape or clipped area consistently.
7. **Patch IDs:** one ID must be one connected component. Separate disconnected
   pieces sharing an ID would create false connectivity.
8. **Geometry:** candidate, patch ID, biodiversity, costs, and optional influence
   rasters must share the exact 10 m grid and CRS. Use nearest-neighbour only for
   categorical class/ID rasters.

## Minimum validation sequence

First run the tiny demo:

```text
ampl FreeNoCarbon_patch_supernodes_demo.run
```

Then run one small real polygon with a deliberately small threshold, verify the
selected components independently in R, restore `MinPatchAreaHa := 5000`, and
only then scale the optimization area. Keep the generated `.dat`, network audit,
patch ID raster, and patch area table as reproducible checkpoints.
