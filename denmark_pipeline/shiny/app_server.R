dp_app_server <- function(input, output, session) {
  defaults <- getOption("denmark.app.defaults", dp_app_defaults())
  # Each input has a known scalar type. Imported JSON only updates these fields.
  state <- shiny::reactive({
    z <- defaults
    for (nm in names(z)) if (!is.null(input[[nm]])) z[[nm]] <- input[[nm]]
    z
  })
  configuration <- shiny::reactive(tryCatch(list(config = dp_app_config(state()), error = NULL),
    error = function(e) list(config = NULL, error = conditionMessage(e))))
  valid_config <- shiny::reactive({
    x <- configuration(); shiny::validate(shiny::need(is.null(x$error), x$error)); x$config
  })
  import_state <- function(s) {
    selections <- c("directions", "area_scope", "anchor_combine", "mode", "extension", "score_type", "footprint", "solver")
    textareas <- c("subarea_goals", "extra_anchors", "outcome_overrides", "metric_scales")
    for (nm in names(defaults)) {
      if (is.logical(defaults[[nm]])) shiny::updateCheckboxInput(session, nm, value = s[[nm]])
      else if (is.numeric(defaults[[nm]])) shiny::updateNumericInput(session, nm, value = s[[nm]])
      else if (nm %in% selections) shiny::updateSelectInput(session, nm, selected = s[[nm]])
      else if (nm %in% textareas) shiny::updateTextAreaInput(session, nm, value = s[[nm]])
      else shiny::updateTextInput(session, nm, value = s[[nm]])
    }
  }
  output$status <- shiny::renderUI({
    x <- configuration()
    if (!is.null(x$error)) shiny::div(class = "status-error", shiny::strong("Check a setting"), shiny::tags$br(), x$error)
    else shiny::tagList(shiny::div(class = "status-ready", shiny::strong("Settings are consistent")),
      shiny::p(class = "status-sub", "Input coverage and landscape feasibility still need checking."))
  })
  output$summary <- shiny::renderUI({
    cfg <- valid_config()
    shiny::tagList(lapply(dp_app_summary(cfg), function(x) shiny::div(class = "summary-item", x)))
  })
  output$workflow_cards <- shiny::renderUI({
    s <- state()
    steps <- list(c("Spatial warm-up", paste0("Biodiversity and spatial preferences; weight ", s$warmup_weight, ". Network requirements relaxed.")),
      c("Minimum shortfall", "Reduce deficits against Global and subarea network-area targets."),
      c(if (s$repair_on) "Repair if needed" else "Repair disabled", if (s$repair_on)
        "Enforce every network target and minimize changed hectares. Skip when already feasible." else "Continue shortfall stages if deficits remain."),
      c("Final optimization", paste0("Keep all requirements satisfied. Final spatial weight ", s$quality_weight, "; footprint mode: ", s$footprint, ".")),
      c("Independent validation", "Check area, cost, legal outcomes, footprint limits and the actual 10 m connected networks."))
    shiny::div(class = "workflow-list", lapply(seq_along(steps), function(i)
      shiny::div(class = "workflow-step", shiny::span(class = "step-number", i),
        shiny::div(class = "step-body", shiny::strong(steps[[i]][1]), shiny::p(steps[[i]][2])))))
  })
  output$score_files <- shiny::renderTable({
    cfg <- valid_config()
    do.call(rbind, lapply(names(cfg$model$outcomes), function(nm)
      data.frame(Outcome = nm, Metric = names(cfg$model$outcomes[[nm]]$scores), Path = unname(cfg$model$outcomes[[nm]]$scores))))
  }, striped = TRUE, spacing = "s")
  output$tile_summary <- shiny::renderText({
    s <- state()
    sprintf("Each tile is approximately %.2f × %.2f km. Combined terra memory hints: %.1f GB across %s workers.",
      s$tile_cells/100, s$tile_cells/100, s$mem_gb*s$workers, s$workers)
  })
  file_result <- shiny::reactiveVal(NULL); file_signature <- shiny::reactiveVal(NULL)
  shiny::observeEvent(input$check_files, {
    tryCatch({
      cfg <- valid_config(); x <- dp_app_file_check(cfg, isTRUE(state()$projection_on) && cfg$model$mode == "pathways")
      file_result(x); file_signature(digest::digest(cfg))
    }, error = function(e) shiny::showNotification(conditionMessage(e), type = "error", duration = 12))
  })
  output$file_check <- shiny::renderTable({
    shiny::req(file_result())
    shiny::validate(shiny::need(identical(file_signature(), digest::digest(valid_config())), "Settings changed. Run the file check again."))
    file_result()
  }, striped = TRUE, spacing = "s")
  output$file_status <- shiny::renderUI({
    shiny::req(file_result())
    shiny::validate(shiny::need(identical(file_signature(), digest::digest(valid_config())), "Settings changed. Run the file check again."))
    counts <- table(file_result()$status)
    shiny::div(class = "note", paste(paste(names(counts), as.integer(counts), sep = ": "), collapse = " · "))
  })

  preview_result <- shiny::reactiveVal(NULL); preview_signature <- shiny::reactiveVal(NULL)
  preview_message <- shiny::reactiveVal("No preview yet. Load a synthetic example, or supply your own spatial layers.")
  job <- shiny::reactiveVal(NULL)
  shiny::observeEvent(input$preview, {
    if (!is.null(job())) { shiny::showNotification("A preview is already running.", type = "message"); return() }
    tryCatch({
      cfg <- valid_config()
      if (!dir.exists(cfg$project_root)) stop("Project root does not exist on this machine.")
      folder <- tempfile("denmark_preview_"); dir.create(folder)
      request <- list(config = cfg, package = dp_home, folder = folder)
      request_file <- file.path(folder, "request.rds"); saveRDS(request, request_file)
      log <- file.path(folder, "preview.log")
      status <- system2(file.path(R.home("bin"), "Rscript"),
        c("--vanilla", shQuote(file.path(dp_home, "shiny", "preview_worker.R")), shQuote(request_file)),
        stdout = log, stderr = log, wait = FALSE)
      if (status != 0) stop("Could not start preview process. See ", log)
      job(list(folder = folder, signature = digest::digest(cfg), started = Sys.time()))
      preview_message("Building the 10 m masks in the background. You can continue editing.")
    }, error = function(e) preview_message(conditionMessage(e)))
  })
  shiny::observe({
    j <- job(); if (is.null(j)) return()
    shiny::invalidateLater(700, session)
    result_file <- file.path(j$folder, "result.rds")
    if (!file.exists(result_file)) {
      if (as.numeric(difftime(Sys.time(), j$started, units = "mins")) > 5)
        preview_message(paste0("Preview is still running. Progress log: ", file.path(j$folder, "preview.log")))
      return()
    }
    z <- readRDS(result_file); job(NULL)
    if (z$ok) {
      preview_result(z$result); preview_signature(j$signature)
      preview_message("Preview ready. These masks use the 10 m input grid; display may be sampled.")
    } else preview_message(paste("Preview failed:", z$message))
  })
  preview_current <- shiny::reactive({
    shiny::req(preview_result())
    shiny::validate(shiny::need(identical(preview_signature(), digest::digest(valid_config())), "Settings changed. Rebuild the preview to show the current configuration."))
    preview_result()
  })
  output$preview_status <- shiny::renderText(preview_message())
  output$preview_controls <- shiny::renderUI({
    p <- preview_current()
    labels <- c(nature = "Current natural habitat", anchors = "Combined anchors", eligible = "Eligible restoration pixels",
      fixed = "Existing connected habitat", budget = "Existing nature credited to budget", credit = "Existing area credited to networks")
    choices <- setNames(names(p$paths), ifelse(names(p$paths) %in% names(labels), labels[names(p$paths)],
      paste("Anchor source:", sub("^source_", "", names(p$paths)))))
    shiny::selectInput("preview_layer", "Map layer", choices = choices, selected = "anchors")
  })
  output$map <- shiny::renderPlot({
    shiny::validate(shiny::need(!is.null(preview_result()), "The map will appear here after building a preview."))
    p <- preview_current(); layer <- input$preview_layer %or% "anchors"
    shiny::req(layer %in% names(p$paths))
    r <- terra::rast(p$paths[[layer]])
    r <- terra::categories(r, value = data.frame(value = 0:1, presence = c("Outside mask", "Inside mask")))
    old <- setwd(p$project_root); on.exit(setwd(old), add = TRUE)
    par(mar = c(2, 2, 2, 4), bg = "#fbfcf7", fg = "#183a32")
    terra::plot(r, col = c("#e7ecde", "#467957"), axes = FALSE,
      plg = list(title = "Mask"), main = gsub("_", " ", layer))
    terra::plot(dp_vector(p$boundary, terra::crs(r)), add = TRUE, border = "#20392e", lwd = 2)
    if (!is.null(p$subareas)) terra::plot(dp_vector(p$subareas, terra::crs(r)), add = TRUE, border = "#a57529", lwd = 1.5, lty = 2)
  }, res = 100)
  output$preview_ledger <- shiny::renderTable({
    p <- preview_current()
    data.frame(Measure = c("Whole polygon", "Current credited nature", "Eligible pixels before unit/outcome filters", "Area needed for percentage target", "Minimum restoration after explicit lower bound"),
      Hectares = c(p$polygon_ha, p$counts[["budget"]], p$counts[["eligible"]], p$needed_ha, p$minimum_ha))
  }, digits = 3, spacing = "s")
  output$export_controls <- shiny::renderUI({
    valid_config()
    shiny::div(class = "download-row", shiny::downloadButton("config_download", "Download R configuration", class = "btn-primary"),
      if (!is.null(preview_result()) && identical(preview_signature(), digest::digest(valid_config())))
        shiny::downloadButton("anchors_download", "Download anchor raster"))
  })
  output$config_download <- shiny::downloadHandler(
    filename = function() paste0(valid_config()$area_id, "_", valid_config()$run_id, ".R"),
    content = function(file) writeLines(dp_app_export(valid_config()), file, useBytes = TRUE), contentType = "text/plain")
  output$settings_download <- shiny::downloadHandler(filename = "denmark_scenario_settings.json",
    content = function(file) jsonlite::write_json(list(schema_version = 1L, settings = state()), file, auto_unbox = TRUE, pretty = TRUE, digits = NA), contentType = "application/json")
  output$anchors_download <- shiny::downloadHandler(filename = function() paste0(valid_config()$area_id, "_anchors_10m.tif"),
    content = function(file) file.copy(preview_current()$exact_anchors, file, overwrite = TRUE), contentType = "image/tiff")
  output$config_code <- shiny::renderText(paste(dp_app_export(valid_config()), collapse = "\n"))
  output$commands <- shiny::renderText({
    cfg <- valid_config(); filename <- paste0(cfg$area_id, "_", cfg$run_id, ".R")
    q <- function(x) encodeString(x, quote = '"')
    paste(c(paste0("setwd(", q(cfg$project_root), ")"), paste0("source(", q(file.path(dp_home, "run.R")), ")"),
      paste0("build_anchors(", q(filename), ")"),
      if (isTRUE(state()$projection_on) && cfg$model$mode == "pathways") paste0("run_projection(", q(filename), ")"),
      paste0("result <- run_area(", q(filename), ")")), collapse = "\n")
  })
  shiny::observeEvent(input$settings_upload, {
    tryCatch({
      import_state(dp_app_read_settings(input$settings_upload$datapath))
      shiny::showNotification("Settings loaded. File paths are interpreted on this machine.", type = "message")
    }, error = function(e) shiny::showNotification(conditionMessage(e), type = "error", duration = 12))
  })
  demo_message <- shiny::reactiveVal("")
  shiny::observeEvent(input$demo, {
    tryCatch({
      s <- dp_app_demo(tempfile("denmark_demo_")); import_state(s)
      demo_message(paste("Synthetic example loaded. Temporary data:", s$project_root))
      shiny::updateTabsetPanel(session, "panel", selected = "review")
    }, error = function(e) shiny::showNotification(conditionMessage(e), type = "error", duration = 12))
  })
  output$demo_status <- shiny::renderText(demo_message())
}
