# From the project root: source('denmark_pipeline/tests/test_spatial_pipeline.R')
# Requires terra, igraph, digest, AMPL and a HiGHS or Gurobi driver.
source("denmark_pipeline/run.R")
fixture <- tempfile("denmark_pipeline_test_", tmpdir = getwd()); dir.create(fixture)
fixture <- normalizePath(fixture, winslash = "/")
rat <- data.frame(id = 1:12, class = c(dp_pathways()$class, "Agriculture", "Urban", "ProductionForest", "Other"))
r <- terra::rast(nrows = 6, ncols = 6, xmin = 500000, xmax = 500060, ymin = 6200000, ymax = 6200060, crs = "EPSG:25832")
m <- matrix(rep(c(1,1,9,9,10,9), 6), nrow = 6, byrow = TRUE)
terra::values(r) <- as.vector(t(m)); r <- terra::categories(r, value = rat); names(r) <- "Landuse"
lu <- file.path(fixture, "landuse.tif"); dp_raster(r, lu, "INT2U")
b <- terra::vect(terra::ext(r), crs = terra::crs(r))
boundary <- file.path(fixture, "boundary.gpkg"); terra::writeVector(b, boundary, overwrite = TRUE)
score <- terra::rast(r); terra::values(score) <- as.vector(t(matrix(rep(c(1,1,1,1,0,9), 6), nrow = 6, byrow = TRUE)))
score <- as.numeric(score)
sf <- file.path(fixture, "score.tif"); dp_raster(score, sf)
cfg <- denmark_config(); cfg$project_root <- getwd(); cfg$area_id <- "tiny"
cfg$grid$landuse <- lu; cfg$grid$planning_unit_m <- 20
cfg$region$boundary <- boundary; cfg$region$subareas <- NULL
cfg$eligibility$min_unit_ha <- 0.001; cfg$eligibility$forest_only <- FALSE
cfg$output_root <- file.path(fixture, "runs")
cfg$model$mode <- "single"; cfg$model$outcomes <- list(restore = list(class = NULL, scores = c(Score = sf)))
cfg$model$metric_weights <- c(Score = 1); cfg$objective$contiguity <- 0
cfg$budget$target_nature_fraction <- NULL; cfg$budget$min_restoration_ha <- 0.08; cfg$budget$max_restoration_ha <- 0.08
cfg$connectivity$global_min_ha <- 0.20
cfg$solver$ampl <- Sys.getenv("AMPL_EXECUTABLE", "ampl"); cfg$solver$name <- "highs"
cfg$solver$warmup_seconds <- 20; cfg$solver$feasibility_seconds <- c(20,20); cfg$solver$quality_seconds <- 20
cfg$solver$options <- "tech:outlev=0"; cfg$solver$gap <- 0
one <- run_area(cfg)
stopifnot(abs(sum(one$selected$area_ha)-0.08) < 1e-9)
audit <- read.csv(file.path(one$final, "exact_10m_network_audit.csv"))
stopifnot(nrow(audit)==1L, all(audit$passes))
# Rerunning reuses completed stages; changing scientific settings is rejected.
again <- run_area(cfg); stopifnot(identical(one$final, again$final))
changed <- cfg; changed$budget$max_restoration_ha <- 0.09
stopifnot(inherits(try(run_area(changed), silent = TRUE), "try-error"))

# Arbitrary subarea count: two, with no per-subarea percentage budget.
top <- terra::vect(terra::ext(500000,500060,6200030,6200060), crs=terra::crs(r))
bottom <- terra::vect(terra::ext(500000,500060,6200000,6200030), crs=terra::crs(r))
subareas <- rbind(top,bottom); subareas$label <- c("north","south")
sp <- file.path(fixture,"subareas.gpkg"); terra::writeVector(subareas,sp,overwrite=TRUE)
two <- cfg; two$run_id <- "two_subareas"; two$region$subareas <- sp; two$region$subarea_id <- "label"
two$connectivity$subarea_min_ha <- c(north=0.08,south=0.08)
two_result <- run_area(two)
stopifnot(nrow(read.csv(file.path(two_result$final,"exact_10m_network_audit.csv")))==3L)

# Independent cap-only run with no network constraints.
cap <- cfg; cap$run_id <- "cost_cap"; cap$connectivity$global_min_ha <- NULL
cap$budget$min_restoration_ha <- 0; cap$budget$cost_per_ha <- sf; cap$budget$cost_cap <- 0.08
cap_result <- run_area(cap)
stopifnot(read.csv(file.path(cap_result$final,"solution_summary.csv"))$cost <= 0.08+1e-8)

# All eight outcomes are real choices. A production-forest unit may only get a forest endpoint.
eight <- cfg; eight$run_id <- "eight"; eight$model$mode <- "pathways"; eight$eligibility$forest_only <- TRUE
r8 <- r; raw8 <- m; raw8[1:2,3:4] <- 11
terra::values(r8) <- as.vector(t(raw8)); r8 <- terra::categories(r8, value = rat)
lu8 <- file.path(fixture,"landuse8.tif"); dp_raster(r8,lu8,"INT2U"); eight$grid$landuse <- lu8
p <- dp_pathways(); eight$model$outcomes <- setNames(lapply(seq_len(8),function(i) {
  x <- score + i; x <- as.numeric(x); f<-file.path(fixture,paste0("score",i,".tif"));dp_raster(x,f)
  list(class=p$class[i], scores=c(Score=f))
}),p$outcome)
er <- run_area(eight)
prep <- readRDS(file.path(er$root,"prepared","prepared.rds"))
forest_units <- prep$units$unit[prep$units$forest_share>=0.5]
stopifnot(all(grepl("^forest", prep$actions$outcome[prep$actions$unit %in% forest_units])))
stopifnot(file.exists(file.path(er$final,"restored_landuse_10m.tif")))

# Lakes can root adjacent land without being nature credit or a free corridor.
lake <- terra::vect(terra::ext(500040,500050,6200000,6200060),crs=terra::crs(r))
lp <- file.path(fixture,"lake.gpkg");terra::writeVector(lake,lp,overwrite=TRUE)
acfg <- cfg; acfg$anchors$sources <- list(lakes=list(type="vector",path=lp,buffer_m=10))
af <- file.path(fixture,"anchors");build_anchors(acfg,af)
stopifnot(dp_sum(terra::rast(file.path(af,"anchors_10m.tif")))>0)
stopifnot(abs(dp_sum(terra::rast(file.path(af,"budget_10m.tif")))*0.01-0.12)<1e-8)
stopifnot(dp_sum(terra::rast(file.path(af,"fixed_10m.tif")))==12)
acfg$run_id <- "shoreline_root"
acfg$budget$min_restoration_ha <- 0.06; acfg$budget$max_restoration_ha <- 0.06
acfg$connectivity$global_min_ha <- 0.06
ar <- run_area(acfg)
ap <- readRDS(file.path(ar$root,"prepared","prepared.rds"))
members <- read.delim(file.path(tail(ar$stages,1)[[1]]$folder,"network_members.tsv"))
roots <- members$node[members$root==1]
stopifnot(all(ap$graphs$Global$nodes$unit[match(roots,ap$graphs$Global$nodes$id)]>0))

# The entire area is one planning unit; disconnected fragments must stay disconnected.
broken <- cfg; broken$run_id <- "disconnected"; broken$grid$planning_unit_m <- 100
broken$budget$min_restoration_ha <- 0.18; broken$budget$max_restoration_ha <- 0.18
broken$connectivity$global_min_ha <- 0.29
bad <- try(run_area(broken),silent=TRUE)
stopifnot(inherits(bad,"try-error"))
status <- readLines(file.path(broken$output_root,broken$area_id,broken$run_id,"status.txt"))
stopifnot(status=="ZERO_SHORTFALL_EXCLUDED_BY_SOLVER_BOUND")

# Anchor-only creation and outcome comparisons need no full-scale data.
cmp <- compare_footprints(file.path(one$final,"selected_10m.tif"),file.path(one$final,"selected_10m.tif"))
stopifnot(cmp$area_weighted_jaccard==1,abs(cmp$shared_ha-0.08)<1e-9)
cat("PASS: complete spatial/AMPL workflow, resume, changed-input rejection, one/eight outcomes, optional subareas, budgets, lake roots, disconnected-fragment rejection, comparison.\n")
cat("Fixture retained:",fixture,"\n")
