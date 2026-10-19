# Orchestration for the existing eight-habitat model. All solves retain habitats.
gp_accepted_result <- function(x) {
  length(x) == 1L && is.finite(x) &&
    ((x >= 0 && x < 100) || (x >= 400 && x < 450))
}
gp_q <- function(x) {
  x <- gsub("\\", "/", x, fixed = TRUE)
  if (any(grepl('["\n\r]', x))) stop("Quotes/newlines are unsupported in AMPL paths.")
  paste0('"', x, '"')
}
gp_write <- function(lines, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (!file.exists(path) || !identical(readLines(path, warn = FALSE), lines)) {
    tmp <- paste0(path, ".building")
    writeLines(lines, tmp, useBytes = TRUE)
    if (!file.copy(tmp, path, overwrite = TRUE)) stop("Cannot write ", path)
    unlink(tmp)
  }
  invisible(path)
}
gp_module <- function(name, e) {
  message("Preparing: ", name)
  sys.source(file.path(gp_package_dir, "R", paste0(name, ".R")), envir = e)
}
gp_read_summary <- function(path) {
  d <- read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  stats::setNames(as.character(d$value), d$item)
}
gp_config <- function(path) {
  e <- new.env(parent = baseenv())
  sys.source(path, envir = e)
  cfg <- e$config
  if (!is.list(cfg)) stop("config.R must define a list called config.")
  if (!cfg$quality_mode %in% c("minimum_change", "free")) stop("Invalid quality_mode.")
  if (!grepl("^[A-Za-z0-9_-]+$", cfg$run_id)) stop("Use letters, digits, _ or - in run_id.")
  if (!identical(as.numeric(cfg$global_nature_target_fraction), 0.30))
    stop("This pipeline implements the agreed whole-polygon 30% target.")
  limits <- unlist(cfg[c("warmup_seconds", "feasibility_seconds",
    "feasibility_continuation_seconds", "repair_verification_seconds",
    "minimum_change_seconds", "quality_seconds")])
  if (any(!is.finite(limits) | limits <= 0)) stop("Time limits must be positive.")
  if (!length(cfg$quality_seconds)) stop("Provide at least one ecological solve time.")
  cfg
}
gp_environment <- function(kind, cfg, root) {
  e <- new.env(parent = globalenv())
  gp_module(paste0("defaults_", kind), e)
  list2env(cfg, envir = e)
  e$expected_number_of_subareas <- 4L
  e$target_nature_proportion <- cfg$global_nature_target_fraction
  e$restoration_tolerance_ha <- cfg$global_restoration_tolerance_ha / 4
  e$problem_dir <- file.path(root, if (kind == "warmup") "warmup" else "network")
  e$problem_name <- if (kind == "warmup") "warmup" else "network"
  e$model_path <- file.path(gp_package_dir, "models",
    if (kind == "warmup") "warmup_global30.mod" else "network_global30.mod")
  e$warmup_dir <- file.path(root, "warmup")
  e$warmup_name <- "warmup"
  e$reuse_matching_checkpoints <- TRUE
  e$solver_name <- cfg$solver
  e$repair_mip_gap <- cfg$mip_gap
  e$repair_feasibility_time_limit_seconds <- cfg$feasibility_seconds
  e$repair_minimum_change_time_limit_seconds <- cfg$minimum_change_seconds
  e$repair_quality_time_limit_seconds <- cfg$quality_seconds[1]
  dir.create(e$problem_dir, recursive = TRUE, showWarnings = FALSE)
  e
}
gp_preflight <- function(cfg) {
  pkgs <- c("terra", "dplyr", "tidyr", "purrr", "readr", "stringr", "digest", "igraph", "knitr")
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install these R packages: ", paste(missing, collapse = ", "))
  shape <- if (tolower(tools::file_ext(cfg$subareas_path)) == "shp")
    paste0(tools::file_path_sans_ext(cfg$subareas_path), c(".shp", ".shx", ".dbf", ".prj")) else cfg$subareas_path
  habitats <- c("forest_dry_poor", "forest_dry_rich", "forest_wet_poor", "forest_wet_rich",
                "open_dry_poor", "open_dry_rich", "open_wet_poor", "open_wet_rich")
  metrics <- as.vector(outer(names(cfg$metric_weights), c("current", paste0("all_to_", habitats)),
                             function(m, s) paste0(m, "_", s, ".tif")))
  files <- c(cfg$fine_landuse_path, cfg$coarse_landuse_path, shape,
             cfg$national_patch_path, cfg$national_patch_dictionary_path,
             file.path(cfg$metrics_dir, metrics))
  if (any(!file.exists(files))) stop("Missing project inputs:\n", paste(files[!file.exists(files)], collapse = "\n"))
  if (cfg$execute_solver && !file.exists(cfg$ampl_executable) && !nzchar(Sys.which(cfg$ampl_executable)))
    stop("AMPL was not found; set ampl_executable in config.R.")
  files
}
gp_lock_manifest <- function(cfg, files, root) {
  # Full content hashes avoid accidental reuse of different data with the same name.
  scientific <- cfg[!names(cfg) %in% c("ampl_executable", "solver_options", "mip_gap",
    "warmup_seconds", "feasibility_seconds", "feasibility_continuation_seconds",
    "repair_verification_seconds", "minimum_change_seconds", "quality_seconds",
    "execute_solver", "make_maps", "sld_path", "raster_maxcell", "png_dpi")]
  sources <- list.files(file.path(gp_package_dir, "R"), full.names = TRUE)
  sources <- c(sources, list.files(file.path(gp_package_dir, "models"), full.names = TRUE))
  signature <- digest::digest(list(config = scientific,
    inputs = tools::md5sum(files), source = tools::md5sum(sources)))
  path <- file.path(root, "manifest.rds")
  if (file.exists(path) && !identical(readRDS(path)$signature, signature))
    stop("Data, code or scientific settings changed. Choose a new run_id; existing results are preserved.")
  if (!file.exists(path)) saveRDS(list(signature = signature, config = cfg,
    input_hashes = tools::md5sum(files), source_hashes = tools::md5sum(sources),
    created = Sys.time()), path)
  signature
}
gp_options <- function(cfg, seconds, logfile, gap = cfg$mip_gap) {
  paste(paste0("lim:time=", as.integer(seconds)), paste0("mip:gap=", gap),
    if (cfg$solver == "gurobi") "mip:start=1" else "", "mip:bestbound=1", cfg$solver_options,
    paste0("tech:logfile=", gp_q(logfile)))
}
gp_solution_paths <- function(e, prefix) {
  fields <- c(solution_path = "_selected_cells.tsv", selected_subarea_path = "_selected_cell_subareas.tsv",
    subarea_solution_path = "_subarea_summary.tsv", network_cell_solution_path = "_network_cells.tsv",
    network_patch_solution_path = "_network_patches.tsv", network_root_solution_path = "_network_roots.tsv",
    network_solution_summary_path = "_network_summary.tsv", summary_path = "_solution_summary.tsv",
    repair_stage_summary_path = "_repair_stage_summary.tsv", repair_completion_path = "_repair_completion_status.txt")
  for (n in names(fields)) assign(n, paste0(prefix, fields[[n]]), envir = e)
}
gp_guard <- function(network = TRUE) {
  binaries <- c("max {c in Cells} abs(Selected[c] - round(Selected[c]))",
    "max {(l,c) in AllowedActions} abs(LanduseDecision[l,c] - round(LanduseDecision[l,c]))")
  if (network) binaries <- c(binaries,
    "max {k in RequiredNetworks} abs(CertificateBuilt[k] - round(CertificateBuilt[k]))",
    "max {(k,c) in EligibleNetworkCells} abs(NetworkCell[k,c] - round(NetworkCell[k,c]))",
    "max {(k,p) in EligibleNetworkPatches} abs(NetworkPatch[k,p] - round(NetworkPatch[k,p]))",
    "max {(k,p) in EligibleNetworkPatches} abs(NetworkRoot[k,p] - round(NetworkRoot[k,p]))")
  c('if not ((solve_result_num >= 0 and solve_result_num < 100) or (solve_result_num >= 400 and solve_result_num < 450)) then { printf "No accepted feasible solver result.\\n"; exit 1; };',
    'if max {i in 1.._ncons} max(0, -_con[i].slack) > 1e-5 then { printf "Constraint residual check failed.\\n"; exit 1; };',
    'if max {j in 1.._nvars} max(0, -_var[j].slack) > 1e-5 then { printf "Variable bound check failed.\\n"; exit 1; };',
    paste0("if max(", paste(binaries, collapse = ", "), ') > 1e-5 then { printf "Integrality check failed.\\n"; exit 1; };'))
}
gp_execute <- function(cfg, runfile, donefile, stdoutfile) {
  if (file.exists(donefile) && identical(readLines(donefile, warn = FALSE), "COMPLETE")) return(TRUE)
  if (!cfg$execute_solver) stop("Prepared AMPL run: ", runfile,
    "\nSet execute_solver = TRUE and rerun, or execute this file with AMPL and rerun.")
  message("Solving: ", runfile, "\nLog: ", stdoutfile)
  status <- system2(cfg$ampl_executable, args = shQuote(runfile), stdout = stdoutfile, stderr = stdoutfile)
  isTRUE(status == 0) && file.exists(donefile) && identical(readLines(donefile, warn = FALSE), "COMPLETE")
}
gp_stage <- function(e, cfg, root, label, seconds, objective, start,
                     hard = FALSE, footprint_limit = NULL) {
  key <- substr(digest::digest(list(label, seconds, objective, hard, footprint_limit,
    cfg$solver, cfg$solver_options, cfg$mip_gap, tools::md5sum(c(e$model_path, e$dat_path, start)))), 1, 12)
  folder <- file.path(root, "stages", paste0(label, "_", key))
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  prefix <- file.path(folder, "solution")
  gp_solution_paths(e, prefix)
  checkpoint <- file.path(folder, "checkpoint.run")
  done <- file.path(folder, "done.txt")
  runfile <- file.path(folder, "solve.run")
  statsfile <- file.path(folder, "stats.tsv")
  other <- setdiff(c("RelativeNetworkShortfall", "WarmFootprintChange", "RepairQuality"), objective)
  code <- c("reset;", paste0("model ", gp_q(e$model_path), ";"),
    paste0("data ", gp_q(e$dat_path), ";"), paste0("include ", gp_q(start), ";"),
    paste0("drop ", other, ";"),
    if (hard) "fix {k in RequiredNetworks} NetworkShortfallHa[k] := 0;",
    if (!is.null(footprint_limit)) paste0("let FootprintChangeLimitHa := ", format(footprint_limit, digits = 17), ";"),
    "suffix bestbound OUT;", paste0("option solver ", cfg$solver, ";"),
    paste0("option ", cfg$solver, "_options '", gp_options(cfg, seconds,
      file.path(folder, "solver.log"), if (objective == "RelativeNetworkShortfall") 0 else cfg$mip_gap), "';"),
    "solve;", "display solve_result_num, solve_result;",
    paste0('printf "item\\tvalue\\nsolve_result_num\\t%d\\n", solve_result_num > ', gp_q(statsfile), ";"),
    gp_guard(TRUE), gp_network_export_lines(e),
    paste0('printf "objective\\t%.17g\\nbest_bound\\t%.17g\\nseconds\\t%.8f\\n", ',
      objective, ", ", objective, ".bestbound, _solve_elapsed_time >> ", gp_q(statsfile), ";"),
    paste0('printf "# Complete incumbent values; not fixed decisions.\\n" > ', gp_q(checkpoint), ";"),
    paste0('printf {j in 1.._nvars} "let %s := %.17g;\\n", _varname[j], _var[j] >> ', gp_q(checkpoint), ";"),
    paste0('printf "COMPLETE\\n" > ', gp_q(done), ";"))
  gp_write(code, runfile)
  if (!gp_execute(cfg, runfile, done, file.path(folder, "ampl.log")))
    stop("Stage did not produce a validated incumbent: ", label, ". Inspect ", folder)
  s <- gp_read_summary(e$summary_path)
  data <- read.delim(statsfile, stringsAsFactors = FALSE)
  list(label = label, folder = folder, prefix = prefix, checkpoint = checkpoint,
       summary = s, stats = stats::setNames(data$value, data$item), objective = objective,
       hard = hard, footprint_limit = footprint_limit)
}
gp_shortfall <- function(stage) as.numeric(stage$summary[["network_shortfall_total_ha"]])
gp_stage_row <- function(s) data.frame(
  stage = s$label, solve_result = s$summary[["solve_result"]],
  solve_result_num = as.numeric(s$summary[["solve_result_num"]]),
  objective_value = as.numeric(s$stats[["objective"]]),
  best_bound = as.numeric(s$stats[["best_bound"]]),
  total_shortfall_ha = gp_shortfall(s),
  footprint_change_ha = as.numeric(s$summary[["warm_footprint_change_ha"]]),
  footprint_limit_ha = if (is.null(s$footprint_limit)) NA_real_ else s$footprint_limit,
  seconds = as.numeric(s$stats[["seconds"]]), stringsAsFactors = FALSE)

run_gudenaa_pipeline <- function(config_path = file.path(gp_package_dir, "config.R")) {
  cfg <- gp_config(config_path)
  old <- setwd(normalizePath(cfg$project_root, mustWork = TRUE)); on.exit(setwd(old), add = TRUE)
  files <- gp_preflight(cfg)
  root <- file.path(cfg$output_root, cfg$run_id)
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  root <- normalizePath(root, winslash = "/")
  # Prevent concurrent writers; an interrupted process may leave a lock to inspect/remove.
  lock <- file.path(root, "RUNNING.lock")
  if (!dir.create(lock, showWarnings = FALSE)) stop("Run is locked: ", lock,
    ". If no process is running, remove this empty lock directory and resume.")
  on.exit(unlink(lock, recursive = TRUE), add = TRUE)
  gp_lock_manifest(cfg, files, root)
  gp_write("RUNNING", file.path(root, "pipeline_status.txt"))
  on.exit({
    status_file <- file.path(root, "pipeline_status.txt")
    if (file.exists(status_file) && identical(readLines(status_file, warn = FALSE), "RUNNING"))
      gp_write(if (cfg$execute_solver) "FAILED_OR_INTERRUPTED_CHECK_LOGS" else
                 "PREPARED_AWAITING_SOLVER", status_file)
  }, add = TRUE)
  gp_write(capture.output(sessionInfo()), file.path(root, "R_session.txt"))
  source(file.path(gp_package_dir, "R", "network_exports.R"), local = environment())
  w <- gp_environment("warmup", cfg, root)
  gp_module("prepare_warmup", w)
  grid_spec <- list(crs = terra::crs(w$coarse_cell_id),
    extent = as.vector(terra::ext(w$coarse_cell_id)), resolution = terra::res(w$coarse_cell_id),
    dimensions = c(nrow(w$coarse_cell_id), ncol(w$coarse_cell_id)),
    candidates = w$candidate_table[c("coarse_cell_id", "EligibleAreaHa")],
    boundary_hash = tools::md5sum(cfg$subareas_path))
  saveRDS(grid_spec, file.path(root, "grid_specification.rds"))
  gp_write(digest::digest(grid_spec), file.path(root, "grid_id.txt"))
  gp_module("export_warmup", w)
  # The warm-up export is retained, but guarded and switched to a global budget.
  warmdone <- file.path(w$problem_dir, "done.txt")
  code <- w$run_lines
  solve_at <- which(code == "solve;")
  code <- append(code, gp_guard(FALSE), after = solve_at)
  code <- append(code, paste0("let GlobalRestorationToleranceHa := ",
                              cfg$global_restoration_tolerance_ha, ";"), after = 4L)
  opt_at <- grep(paste0("^option ", cfg$solver, "_options"), code)
  code[opt_at] <- paste0("option ", cfg$solver, "_options '", gp_options(cfg, cfg$warmup_seconds,
                                                 file.path(w$problem_dir, "solver.log")), "';")
  code <- c(code, paste0('printf "COMPLETE\\n" > ', gp_q(warmdone), ";"))
  gp_write(code, w$run_path)
  if (!gp_execute(cfg, w$run_path, warmdone, file.path(w$problem_dir, "ampl.log")))
    stop("Warm-up did not produce a feasible incumbent. Inspect warmup/ampl.log.")
  n <- gp_environment("network", cfg, root)
  gp_module("prepare_network", n)
  gp_module("export_network_data", n)
  initial_start <- n$mip_start_path
  # If warm-up already meets all certificates, build its complete flow start.
  if (n$complete_feasible_mip_start_available) {
    gp_module("build_warm_certificates", n)
    # The original complete-start exporter expects constructive_assignment.
    n$constructive_assignment <- n$warm_start_assignment
    gp_module("export_constructive_start", n)
  }
  stages <- list()
  record <- function(s) {
    stages[[length(stages) + 1L]] <<- s
    write.table(do.call(rbind, lapply(stages, gp_stage_row)), file.path(root, "stage_summary.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)
    s
  }
  best <- record(gp_stage(n, cfg, root, "01_feasibility",
    if (n$complete_feasible_mip_start_available) cfg$repair_verification_seconds else cfg$feasibility_seconds,
    "RelativeNetworkShortfall", initial_start))
  # Try the established path-and-area repair after each partial solve. When its
  # assumptions do not hold, retain the incumbent and continue the general MIP.
  for (attempt in seq_len(length(cfg$feasibility_continuation_seconds) + 1L)) {
    if (gp_shortfall(best) <= cfg$network_shortfall_tolerance_ha) break
    n$stage1_selected_path <- paste0(best$prefix, "_selected_cells.tsv")
    n$stage1_network_cells_path <- paste0(best$prefix, "_network_cells.tsv")
    patches <- read.delim(paste0(best$prefix, "_network_patches.tsv"))
    roots <- read.delim(paste0(best$prefix, "_network_roots.tsv"))
    patches$is_root <- as.integer(paste(patches$network_name, patches$patch_id) %in%
                                   paste(roots$network_name, roots$patch_id))
    n$stage1_network_patches_path <- file.path(best$folder, "repair_network_patches.tsv")
    write.table(patches, n$stage1_network_patches_path, sep = "\t", quote = FALSE, row.names = FALSE)
    constructed <- tryCatch({
      a <- read.delim(paste0(best$prefix, "_network_summary.tsv"))
      sub <- a[a$network_name != "Global", ]
      if (nrow(sub) != 4L || any(sub$shortfall_ha > cfg$network_shortfall_tolerance_ha) ||
          any(sub$certificate_built != 1)) stop("Subarea certificates still incomplete; use general solver continuation.")
      gp_module("constructive_repair", n)
      n$mip_start_path <- file.path(best$folder, "constructive_start.run")
      gp_module("export_constructive_start", n)
      TRUE
    }, error = function(err) {
      gp_write(conditionMessage(err), file.path(best$folder, "constructive_repair_not_used.txt"))
      message("Constructive repair: ", conditionMessage(err)); FALSE
    })
    if (constructed) {
      repaired <- tryCatch(gp_stage(n, cfg, root, paste0("02_constructive_", attempt),
        cfg$repair_verification_seconds, "RelativeNetworkShortfall", n$mip_start_path),
        error = function(err) { message(conditionMessage(err)); NULL })
      if (!is.null(repaired)) {
        record(repaired)
        if (as.numeric(repaired$stats[["objective"]]) < as.numeric(best$stats[["objective"]])) best <- repaired
      }
    }
    if (gp_shortfall(best) <= cfg$network_shortfall_tolerance_ha) break
    if (attempt <= length(cfg$feasibility_continuation_seconds)) {
      next_stage <- record(gp_stage(n, cfg, root, paste0("03_feasibility_continue_", attempt),
        cfg$feasibility_continuation_seconds[attempt], "RelativeNetworkShortfall", best$checkpoint))
      if (as.numeric(next_stage$stats[["objective"]]) <= as.numeric(best$stats[["objective"]])) best <- next_stage
    }
  }
  if (gp_shortfall(best) > cfg$network_shortfall_tolerance_ha) {
    lower_bound <- max(vapply(stages, function(s)
      if (s$objective == "RelativeNetworkShortfall") as.numeric(s$stats[["best_bound"]]) else -Inf,
      numeric(1)), na.rm = TRUE)
    if (is.finite(lower_bound) && lower_bound >
        cfg$network_shortfall_tolerance_ha / min(n$network_minimums) + 1e-9) {
      gp_write(c("ZERO_SHORTFALL_RULED_OUT_BY_BOUND", best$prefix), file.path(root, "pipeline_status.txt"))
      stop("The solver's positive shortfall lower bound rules out the requested zero-shortfall tolerance under this model.\n",
           "Inspect the area budget, eligibility and network assumptions. Best prefix: ", best$prefix)
    }
    gp_write(c("INCOMPLETE_NETWORK_SHORTFALL", best$prefix), file.path(root, "pipeline_status.txt"))
    stop("Feasibility time budget exhausted; best checkpoint retained. This is not proof of infeasibility.\n",
      "Add a feasibility_continuation_seconds entry and rerun. Best prefix: ", best$prefix)
  }
  limit <- sum(n$candidate_table$EligibleAreaHa)
  if (cfg$quality_mode == "minimum_change") {
    best <- record(gp_stage(n, cfg, root, "04_minimum_change", cfg$minimum_change_seconds,
                            "WarmFootprintChange", best$checkpoint, TRUE, limit))
    limit <- as.numeric(best$summary[["warm_footprint_change_ha"]]) + cfg$repair_change_tolerance_ha
  }
  for (i in seq_along(cfg$quality_seconds)) {
    next_stage <- record(gp_stage(n, cfg, root, paste0("05_quality_", i), cfg$quality_seconds[i],
                                  "RepairQuality", best$checkpoint, TRUE, limit))
    if (as.numeric(next_stage$summary[["objective"]]) + 1e-10 >= as.numeric(best$summary[["objective"]]))
      best <- next_stage
  }
  finaldir <- file.path(root, "final", basename(best$folder))
  dir.create(finaldir, recursive = TRUE, showWarnings = FALSE)
  prefix <- file.path(finaldir, "solution")
  suffix <- c("_selected_cells.tsv", "_selected_cell_subareas.tsv", "_subarea_summary.tsv",
    "_network_cells.tsv", "_network_patches.tsv", "_network_roots.tsv", "_network_summary.tsv", "_solution_summary.tsv")
  if (!all(file.copy(paste0(best$prefix, suffix), paste0(prefix, suffix), overwrite = TRUE))) stop("Final export copy failed.")
  stage_table <- do.call(rbind, lapply(stages, gp_stage_row))
  # Retain the audit's historical quality-stage name in the compatibility table.
  stage_table$stage[grepl("^05_quality", stage_table$stage)] <- "3_ecological_contiguity_quality"
  for (s in c("_stage_summary.tsv", "_repair_stage_summary.tsv"))
    write.table(stage_table, paste0(prefix, s), sep = "\t", quote = FALSE, row.names = FALSE)
  gp_solution_paths(n, prefix)
  gp_write("COMPLETE", n$repair_completion_path)
  n$problem_dir <- finaldir
  gp_module("audit_network", n)
  checks <- read.csv(file.path(finaldir, "exact_10m_network_goal_audit.csv"))
  if (nrow(checks) != 5L || anyNA(checks$PassesExact10m) || !all(checks$PassesExact10m)) {
    gp_write("COARSE_FEASIBLE_FINE_REFINEMENT_REQUIRED", file.path(root, "pipeline_status.txt"))
    stop("The coarse solution failed the 10 m audit. See ",
         file.path(finaldir, "exact_10m_network_goal_audit.csv"),
         "; no validated result was published.")
  }
  gp_write(prefix, file.path(root, "latest_validated_solution.txt"))
  gp_write("VALIDATED", file.path(root, "pipeline_status.txt"))
  if (cfg$make_maps) {
    map_packages <- c("ggplot2", "tidyterra", "xml2", "ragg")
    ready <- file.exists(cfg$sld_path) && all(vapply(map_packages, requireNamespace, logical(1), quietly = TRUE))
    if (ready) {
      maps <- new.env(parent = globalenv()); list2env(cfg, maps)
      maps$problem_dir <- n$problem_dir; maps$problem_name <- "solution"
      maps$final_solution_prefix <- prefix
      maps$biodiversity_support_path <- file.path(cfg$metrics_dir, "Richness_current.tif")
      maps$figure_dir <- file.path(finaldir, "figures")
      maps$reuse_certificate_geometry_cache <- FALSE
      tryCatch(gp_module("plot_solution", maps), error = function(err) {
        gp_write(conditionMessage(err), file.path(finaldir, "map_error.txt"))
        warning("Solution validated, but maps failed: ", conditionMessage(err))
      })
    } else warning("Solution validated. Maps need NaturetypeStyle.sld and ggplot2, tidyterra, xml2, ragg.")
  }
  message("Validated solution: ", prefix)
  invisible(list(root = root, solution_prefix = prefix, stage_summary = stage_table))
}
