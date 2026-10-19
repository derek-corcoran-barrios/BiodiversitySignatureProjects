dp_app_demo <- function(folder) {
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  folder <- normalizePath(folder, winslash = "/")
  s <- dp_app_defaults(folder)
  r <- terra::rast(nrows = 40, ncols = 60, xmin = 500000, xmax = 500600, ymin = 6200000, ymax = 6200400, crs = "EPSG:25832")
  xy <- terra::xyFromCell(r, seq_len(terra::ncell(r)))
  x <- (xy[,1]-500000)/600; y <- (xy[,2]-6200000)/400
  values <- rep(9, terra::ncell(r))
  values[((x-.2)/.18)^2+((y-.75)/.24)^2 < 1] <- 1
  values[((x-.8)/.16)^2+((y-.22)/.23)^2 < 1] <- 6
  values[x > .64 & y > .65] <- 11
  values[((x-.52)/.055)^2+((y-.46)/.13)^2 < 1] <- 12
  terra::values(r) <- values
  r <- terra::categories(r, value = data.frame(id = 1:12, class = c(dp_pathways()$class, "Agriculture", "Urban", "ProductionForest", "Other")))
  names(r) <- "Landuse"; s$landuse <- file.path(folder, "landuse.tif"); dp_raster(r, s$landuse, "INT2U")
  b <- terra::vect(terra::ext(r), crs = terra::crs(r)); s$boundary <- file.path(folder, "boundary.gpkg"); terra::writeVector(b, s$boundary)
  a <- terra::vect(terra::ext(500000,500300,6200000,6200400), crs = terra::crs(r))
  b <- terra::vect(terra::ext(500300,500600,6200000,6200400), crs = terra::crs(r)); sub <- rbind(a,b); sub$id <- c("west","east")
  s$subareas <- file.path(folder, "subareas.gpkg"); terra::writeVector(sub, s$subareas)
  lake <- terra::as.polygons(terra::ifel(as.numeric(r) == 12, 1, NA), dissolve = TRUE)
  s$lakes_path <- file.path(folder, "lakes.gpkg"); terra::writeVector(lake, s$lakes_path)
  political <- terra::vect(terra::ext(500230,500290,6200220,6200290), crs = terra::crs(r))
  s$political_path <- file.path(folder, "political.gpkg"); terra::writeVector(political, s$political_path)
  custom <- terra::ifel(as.numeric(r) == 1, 1, 0); s$custom_path <- file.path(folder, "custom_anchors.tif"); dp_raster(custom, s$custom_path, "INT1U")
  s$metrics_dir <- file.path(folder, "metrics"); dir.create(s$metrics_dir); s$baseline_dir <- s$metrics_dir; s$extension <- "tif"
  for (i in seq_len(8)) for (metric in c("Richness", "PD", "RangeRarity")) {
    score <- as.numeric(r)
    v <- 0.2 + 0.8*exp(-((x-(i/10))^2+(y-((9-i)/10))^2)*3)
    terra::values(score) <- v * switch(metric, Richness = 200, PD = 1000, RangeRarity = 5)
    dp_raster(score, file.path(s$metrics_dir, paste0(metric, "_all_to_", dp_pathways()$outcome[i], ".tif")))
    if (i == 1) dp_raster(score*0.6, file.path(s$metrics_dir, paste0(metric, "_current.tif")))
  }
  score <- as.numeric(r); terra::values(score) <- 1 + x + y
  s$single_score <- file.path(folder, "independent_score.tif"); dp_raster(score, s$single_score)
  s$single_baseline <- file.path(folder, "independent_current.tif"); dp_raster(score*0.6, s$single_baseline)
  s$cost_raster <- file.path(folder, "cost_per_ha.tif"); dp_raster(score*50000, s$cost_raster)
  s$area_id <- "synthetic_demo"; s$run_id <- "demo_v3"; s$output_root <- file.path(folder, "runs")
  s$planning_m <- 40; s$min_unit_ha <- 0.01; s$subareas_on <- TRUE; s$subarea_id <- "id"; s$subarea_min <- 3
  s$lakes_on <- TRUE; s$global_min <- 8; s$target_percent <- 40; s$tolerance_ha <- 0.32
  s$warmup_minutes <- 1; s$feasibility_minutes <- "2, 2"; s$quality_minutes <- "2"; s$repair_minutes <- 1
  s
}
