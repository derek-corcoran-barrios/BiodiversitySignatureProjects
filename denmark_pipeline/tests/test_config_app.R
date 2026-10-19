source("denmark_pipeline/run.R")
source("denmark_pipeline/shiny/helpers.R")
source("denmark_pipeline/shiny/demo.R")
source("denmark_pipeline/shiny/app_ui.R")
source("denmark_pipeline/shiny/app_server.R")
s <- dp_app_defaults()
cfg <- dp_app_config(s)
stopifnot(length(cfg$model$outcomes)==8,is.null(cfg$region$subareas),cfg$repair$enabled)
# Letters in valid driver options must not be mistaken for newline escapes.
threads <- s; threads$solver_options <- "tech:threads=8 tech:outlev=1"
stopifnot(identical(dp_app_config(threads)$solver$options, threads$solver_options))
for (invalid in c("tech:outlev=1\nreset;", "tech:outlev=1\rreset;", "'", '\"')) {
  threads$solver_options <- invalid
  stopifnot(inherits(try(dp_app_config(threads), silent=TRUE), "try-error"))
}
# Roots can be omitted only when there are no anchored-network requirements.
no_roots <- s; no_roots$nature_on <- FALSE
stopifnot(inherits(try(dp_app_config(no_roots), silent=TRUE), "try-error"))
no_roots$global_on <- FALSE
stopifnot(length(dp_app_config(no_roots)$anchors$sources)==0)
# Settings JSON round trip and exact executable R export round trip.
js <- tempfile(fileext=".json")
jsonlite::write_json(list(schema_version=1L,settings=s),js,auto_unbox=TRUE)
roundtrip <- dp_app_read_settings(js)
stopifnot(identical(dp_app_config(roundtrip),cfg))
f <- tempfile(fileext=".R");writeLines(dp_app_export(cfg),f)
stopifnot(identical(dp_config(f),cfg))
# Single outcomes never inherit the production-forest pathway restriction.
s$mode <- "single";s$single_score <- "independent.tif";s$footprint <- "fixed"
one <- dp_app_config(s)
stopifnot(length(one$model$outcomes)==1,!one$eligibility$forest_only,one$quality$footprint=="fixed")
# Correct individual subarea mapping, source union and optional area credit.
s$subareas_on <- TRUE;s$subareas <- "subareas.gpkg";s$subarea_id <- "name";s$subarea_goals <- '{"west": 300, "east": 700}'
s$lakes_on <- TRUE;s$lakes_path <- "lakes.gpkg"
s$extra_anchors <- '{"reserve": {"type": "vector", "path": "reserve.gpkg"}}'
many <- dp_app_config(s)
stopifnot(identical(many$connectivity$subarea_min_ha,c(west=300,east=700)),length(many$anchors$sources)==3)
s$lakes_network <- TRUE
stopifnot(inherits(try(dp_app_config(s),silent=TRUE),"try-error"))
s$lakes_existing <- TRUE;stopifnot(dp_app_config(s)$anchors$sources$lakes$credit_network)
s$target_percent <- 110;stopifnot(inherits(try(dp_app_config(s),silent=TRUE),"try-error"))
s <- dp_app_defaults();s$projection_on <- TRUE;s$projection_id <- "test_tiles"
projected <- dp_app_config(s)
stopifnot(all(grepl("test_tiles/metrics/",projected$model$outcomes[[1]]$scores)))
s$quality_weight <- 0;s$warmup_weight <- 3;s$repair_on <- FALSE
stages <- dp_app_config(s)
stopifnot(dp_phase_weight(stages,"warmup")==3,dp_phase_weight(stages,"quality")==0,!stages$repair$enabled)
# Reject imports with unknown keys and non-scalar injected values.
jsonlite::write_json(list(schema_version=1,settings=list(unknown=1)),js,auto_unbox=TRUE)
stopifnot(inherits(try(dp_app_read_settings(js),silent=TRUE),"try-error"))
jsonlite::write_json(list(schema_version=1,settings=list(area_id=list(a=1))),js,auto_unbox=TRUE)
stopifnot(inherits(try(dp_app_read_settings(js),silent=TRUE),"try-error"))
# Real spatial preview with multiple roots, own raster and 2 subareas.
folder <- tempfile("app_demo_",tmpdir=getwd());demo <- dp_app_demo(folder)
demo_cfg <- dp_app_config(demo)
check <- dp_app_file_check(demo_cfg)
stopifnot(all(check$status=="Present"))
preview <- dp_app_make_preview(demo_cfg,file.path(folder,"preview"))
stopifnot(preview$polygon_ha==24,preview$counts[["anchors"]]>0,
  all(c("nature","anchors","source_nature","source_lakes") %in% names(preview$paths)),file.exists(preview$exact_anchors))
custom <- demo;custom$nature_on <- FALSE;custom$lakes_on <- FALSE;custom$custom_on <- TRUE
stopifnot(identical(names(dp_app_config(custom)$anchors$sources),"custom"))
# Render a complete UI document: catches malformed tag nesting and missing inputs.
ui <- dp_app_ui(dp_app_defaults())
html <- htmltools::renderTags(ui)$html
for(id in names(dp_app_defaults())) stopifnot(grepl(paste0('id="',id,'"'),html,fixed=TRUE))
shiny::testServer(dp_app_server, {
  session$setInputs(mode="single",single_score="independent.tif",footprint="fixed",quality_weight=0)
  stopifnot(configuration()$config$model$mode=="single",configuration()$config$quality$footprint=="fixed")
  stopifnot(any(grepl("One independent",as.character(output$summary),fixed=TRUE)))
  session$setInputs(cost_on=TRUE,cost_raster="")
  stopifnot(grepl("cost cap",configuration()$error,fixed=TRUE))
  session$setInputs(cost_raster="cost.tif",subareas_on=FALSE,target_percent=110)
  stopifnot(grepl("Nature fraction",configuration()$error,fixed=TRUE))
  session$setInputs(target_percent=30)
  stopifnot(is.null(configuration()$error))
})
cat("PASS: config validation, independent model, JSON/R round trips, multi-source anchors, optional subareas, projection paths, exact spatial preview, UI controls.\n")
cat("Demo retained:",normalizePath(folder),"\n")
