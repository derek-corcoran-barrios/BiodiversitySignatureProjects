dp_aq <- function(x) {
  if (!length(x)) return(character())
  x <- gsub("\\", "/", x, fixed = TRUE)
  if (any(grepl('["\n\r]', x))) stop("Unsupported quote/newline in AMPL string.")
  paste0('"', x, '"')
}
dp_num <- function(x) format(as.numeric(x), digits = 17, scientific = FALSE, trim = TRUE)
dp_ampl_data <- function(cfg, p, path) {
  out <- character()
  emit <- function(x) out <<- c(out, x)
  set <- function(n, x) emit(paste0("set ", n, " := ", paste(x, collapse = " "), ";"))
  scalar <- function(n, v) emit(paste0("param ", n, " := ", dp_num(v), ";"))
  param <- function(n, rows) emit(c(paste0("param ", n, " :="), rows, ";"))
  u <- p$units; a <- p$actions; s <- dp_spatial_terms(p)
  set("Units", u$unit); set("Outcomes", dp_aq(names(cfg$model$outcomes)))
  set("Actions", paste0("(", dp_aq(a$outcome), ",", a$unit, ")"))
  se <- character()
  for (o in names(cfg$model$outcomes)) {
    allowed <- a$unit[a$outcome == o]
    e <- s$pairs[s$pairs$from %in% allowed & s$pairs$to %in% allowed, , drop = FALSE]
    if (nrow(e)) se <- c(se, paste0("(", dp_aq(o), ",", e$from, ",", e$to, ")"))
  }
  set("SameEdges", se); set("Networks", dp_aq(names(p$graphs)))
  param("Area", paste(u$unit, dp_num(u$area_ha)))
  param("Cost", paste(u$unit, dp_num(u$cost)))
  param("Score", paste(dp_aq(a$outcome), a$unit, dp_num(a$score)))
  param("AnchorContact", paste(names(s$contact), dp_num(s$contact)))
  scalar("MinimumArea", p$ledger$minimum_restoration_ha); scalar("MaximumArea", p$ledger$maximum_restoration_ha)
  scalar("MaximumCost", p$ledger$cost_cap)
  scalar("BioNorm", max(1e-12, sum(vapply(split(abs(a$score), a$unit), max, numeric(1)))))
  scalar("EdgeNorm", max(1, nrow(s$pairs))); scalar("AnchorNorm", max(1, sum(s$contact)))
  scalar("SpatialWeight", dp_phase_weight(cfg, "quality")); scalar("FinalSpatialWeight", dp_phase_weight(cfg, "quality"))
  scalar("NewWeight", cfg$objective$new_to_new)
  scalar("AnchorWeight", cfg$objective$to_anchors); scalar("RequireNew", as.integer(cfg$connectivity$require_new_restoration))
  arcs <- character(); tables <- list(Owner = character(), RootAllowed = character(), CreditArea = character(), FlowCap = character())
  for (k in names(p$graphs)) {
    g <- p$graphs[[k]]; n <- g$nodes; e <- g$edges
    set(paste0("Nodes[", dp_aq(k), "]"), n$id)
    if (nrow(e)) arcs <- c(arcs, paste0("(", dp_aq(k), ",", e$from, ",", e$to, ")"),
      paste0("(", dp_aq(k), ",", e$to, ",", e$from, ")"))
    if (nrow(n)) {
      tables$Owner <- c(tables$Owner, paste(dp_aq(k), n$id, dp_num(n$unit)))
      tables$RootAllowed <- c(tables$RootAllowed, paste(dp_aq(k), n$id, n$anchor))
      tables$CreditArea <- c(tables$CreditArea, paste(dp_aq(k), n$id, dp_num(n$area_ha)))
      tables$FlowCap <- c(tables$FlowCap, paste(dp_aq(k), n$id, n$flow_cap))
    }
  }
  set("Arcs", arcs)
  if (length(p$graphs)) param("Requirement", vapply(names(p$graphs), function(k)
    paste(dp_aq(k), dp_num(p$graphs[[k]]$minimum)), character(1)))
  for (n in names(tables)) if (length(tables[[n]])) param(n, tables[[n]])
  dp_write(out, path)
}
dp_stage <- function(cfg, root, dat, label, seconds, start = NULL,
                     phase = c("quality", "warmup", "feasibility", "repair"),
                     reference = NULL, allow_no_incumbent = FALSE) {
  phase <- match.arg(phase)
  sig <- digest::digest(list(label, seconds, phase, cfg$solver$name, cfg$solver$options, cfg$solver$gap,
    cfg$objective, if (phase == "repair") cfg$repair$gap,
    if (phase == "quality") cfg$quality, reference,
    tools::md5sum(c(dat, start, file.path(dp_home, "models", "restoration.mod")))))
  folder <- file.path(root, "stages", paste0(label, "_", substr(sig, 1, 12)))
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  done <- file.path(folder, "done.rds"); stats <- file.path(folder, "stats.tsv")
  selected <- file.path(folder, "selected.tsv"); checkpoint <- file.path(folder, "checkpoint.run")
  membership <- file.path(folder, "network_members.tsv"); marker <- file.path(folder, "ampl_complete.txt")
  outputs <- c(stats, selected, checkpoint, membership)
  if (file.exists(done)) {
    old <- readRDS(done)
    if (all(file.exists(outputs)) && identical(unname(tools::md5sum(outputs)), old$hashes)) return(old)
    stop("Completed stage files were modified or removed: ", folder)
  }
  obj <- switch(phase, feasibility = "MinShortfall", repair = "FootprintChange", "Quality")
  if (phase == "repair" && is.null(reference)) stop("Repair needs a reference footprint.")
  ref_cmd <- if (!is.null(reference)) c("let {u in Units} ReferenceSelected[u] := 0;",
    if (length(reference)) paste0("let ReferenceSelected[", reference, "] := 1;"))
  change_cmd <- if (phase == "quality" && cfg$quality$footprint != "free") {
    if (is.null(reference)) stop("Footprint protection needs a reference footprint.")
    if (cfg$quality$footprint == "fixed") "fix {u in Units} Selected[u] := ReferenceSelected[u];"
    else paste0("let MaximumChange := ", dp_num(cfg$quality$max_change_ha), ";")
  }
  weight <- dp_phase_weight(cfg, if (phase == "warmup") "warmup" else "quality")
  gap <- if (phase == "feasibility") 0 else if (phase == "repair") cfg$repair$gap else cfg$solver$gap
  spatial <- "(NewWeight * sum {(o,u,v) in SameEdges} Adjacent[o,u,v]/EdgeNorm + AnchorWeight * sum {u in Units} AnchorContact[u]*Selected[u]/AnchorNorm)"
  run <- c("reset;", paste0("model ", dp_aq(normalizePath(file.path(dp_home, "models", "restoration.mod"), winslash = "/")), ";"),
    paste0("data ", dp_aq(dat), ";"), if (!is.null(start)) paste0("include ", dp_aq(start), ";"),
    paste0("drop ", setdiff(c("Quality", "MinShortfall", "FootprintChange"), obj), ";"),
    ref_cmd, change_cmd, paste0("let SpatialWeight := ", dp_num(weight), ";"),
    if (phase == "warmup") c("fix {k in Networks} Built[k] := 0;", "fix {k in Networks} Shortfall[k] := Requirement[k];"),
    if (phase %in% c("quality", "repair")) "fix {k in Networks} Shortfall[k] := 0;",
    "suffix bestbound OUT;", paste0("option solver ", cfg$solver$name, ";"),
    paste0("option ", cfg$solver$name, "_options '",
      paste0("lim:time=", as.integer(seconds)), " mip:gap=", gap,
      if (cfg$solver$name == "gurobi") " mip:start=1" else "", " mip:bestbound=1 ", cfg$solver$options, "';"),
    "solve;",
    paste0('if not ((solve_result_num >= 0 and solve_result_num < 100) or (solve_result_num >= 400 and solve_result_num < 450)) then { printf "NO_INCUMBENT\\t%d\\n", solve_result_num > ', dp_aq(marker), '; exit 0; };'),
    'if max {i in 1.._ncons} max(0, -_con[i].slack) > 1e-5 then { printf "Constraint residual failure.\\n"; exit 1; };',
    'if max {j in 1.._nvars} max(0, -_var[j].slack) > 1e-5 then { printf "Bound residual failure.\\n"; exit 1; };',
    'if max {j in 1.._nvars: _var[j].int > 0} abs(_var[j] - round(_var[j])) > 1e-5 then { printf "Integrality failure.\\n"; exit 1; };',
    paste0('printf "item\\tvalue\\nsolve_result_num\\t%d\\nobjective\\t%.17g\\nbest_bound\\t%.17g\\nquality\\t%.17g\\nshortfall_ha\\t%.17g\\nseconds\\t%.8f\\nrelative_shortfall\\t%.17g\\nfootprint_change_ha\\t%.17g\\nspatial_weight\\t%.17g\\n", solve_result_num, ',
      obj, ", ", obj, ".bestbound, Quality + (FinalSpatialWeight-SpatialWeight)*", spatial,
      ", sum {k in Networks} Shortfall[k], _solve_elapsed_time, MinShortfall, FootprintChange, SpatialWeight > ", dp_aq(stats), ";"),
    paste0('printf "unit\\toutcome\\tarea_ha\\tscore\\n" > ', dp_aq(selected), ";"),
    paste0('printf {(o,u) in Actions: X[o,u] > 0.5} "%d\\t%s\\t%.12g\\t%.17g\\n", u,o,Area[u],Score[o,u] >> ', dp_aq(selected), ";"),
    paste0('printf "network\\tnode\\troot\\tcredit_ha\\n" > ', dp_aq(membership), ";"),
    paste0('printf {k in Networks, n in Nodes[k]: Member[k,n]>0.5} "%s\\t%d\\t%d\\t%.12g\\n", k,n,Root[k,n],CreditArea[k,n] >> ', dp_aq(membership), ";"),
    paste0('printf "# Incumbent values only: decisions remain free.\\n" > ', dp_aq(checkpoint), ";"),
    paste0('printf {j in 1.._nvars} "let %s := %.17g;\\n", _varname[j], _var[j] >> ', dp_aq(checkpoint), ";"),
    paste0('printf "COMPLETE\\n" > ', dp_aq(marker), ";"))
  runfile <- dp_write(run, file.path(folder, "solve.run"))
  if (!cfg$solver$execute) stop("Prepared: ", runfile, ". Set solver$execute = TRUE to run.")
  # A stale marker from an interrupted R process must not certify a new solve.
  if (file.exists(marker)) unlink(marker)
  message("Solving ", label, "; log: ", file.path(folder, "ampl.log"))
  status <- system2(cfg$solver$ampl, shQuote(runfile), stdout = file.path(folder, "ampl.log"), stderr = file.path(folder, "ampl.log"))
  if (status != 0 || !file.exists(marker)) stop("Solver stage failed: ", folder)
  marker_text <- readLines(marker, warn = FALSE)
  if (length(marker_text) && startsWith(marker_text[1], "NO_INCUMBENT")) {
    code <- as.numeric(strsplit(marker_text[1], "\t", fixed = TRUE)[[1]][2])
    # A genuine infeasibility or time limit can be retained as a repair attempt;
    # driver/license/errors must never be silently interpreted as difficulty.
    if (allow_no_incumbent && is.finite(code) && ((code >= 200 && code < 300) || (code >= 450 && code < 500)))
      return(list(folder = folder, phase = phase, label = label, available = FALSE,
        stats = c(solve_result_num = code, objective = NA, best_bound = NA, quality = NA,
          shortfall_ha = NA, seconds = NA, relative_shortfall = NA, footprint_change_ha = NA, spatial_weight = weight)))
    stop("No usable incumbent (code ", code, "). See ", file.path(folder, "ampl.log"))
  }
  if (!identical(marker_text, "COMPLETE") || !all(file.exists(outputs))) stop("Incomplete solver export: ", folder)
  z <- read.delim(stats, stringsAsFactors = FALSE)
  result <- list(folder = folder, phase = phase, label = label, available = TRUE, stats = setNames(as.numeric(z$value), z$item),
    selected = selected, checkpoint = checkpoint, hashes = unname(tools::md5sum(outputs)))
  dp_save(result, done); result
}
run_area <- function(config) {
  cfg <- dp_config(config); old <- setwd(cfg$project_root); on.exit(setwd(old), add = TRUE)
  inputs <- dp_input_files(cfg)
  scientific <- cfg; scientific$solver <- NULL; scientific$projection <- NULL; scientific$repair <- NULL
  signature <- dp_signature(scientific, inputs)
  root <- file.path(cfg$output_root, cfg$area_id, cfg$run_id)
  dp_manifest(root, signature, cfg); root <- normalizePath(root, winslash = "/")
  lock <- dp_lock(root); on.exit(unlink(lock, recursive = TRUE), add = TRUE)
  statusfile <- file.path(root, "status.txt"); dp_write("RUNNING", statusfile)
  on.exit({ if (file.exists(statusfile) && identical(readLines(statusfile), "RUNNING")) dp_write("INTERRUPTED_OR_FAILED", statusfile) }, add = TRUE)
  dir.create(file.path(root, "prepared"), showWarnings = FALSE)
  prepared_file <- file.path(root, "prepared", "prepared.rds")
  p <- if (file.exists(prepared_file)) readRDS(prepared_file) else prepare_area(cfg, file.path(root, "prepared"))
  dat <- dp_ampl_data(cfg, p, file.path(root, "problem.dat"))
  stages <- list(); record <- function(x) {
    stages[[length(stages) + 1L]] <<- x
    dp_csv(do.call(rbind, lapply(stages, function(s) data.frame(stage = s$label,
      as.list(s$stats), check.names = FALSE))), file.path(root, "stage_summary.csv")); x
  }
  best <- record(dp_stage(cfg, root, dat, "warmup", cfg$solver$warmup_seconds, phase = "warmup"))
  if (length(p$graphs)) for (i in seq_along(cfg$solver$feasibility_seconds)) {
    z <- record(dp_stage(cfg, root, dat, paste0("feasibility_", i), cfg$solver$feasibility_seconds[i], best$checkpoint, "feasibility"))
    if (z$stats[["relative_shortfall"]] <= best$stats[["relative_shortfall"]] || z$stats[["shortfall_ha"]] <= 1e-6) best <- z
    if (best$stats[["shortfall_ha"]] <= 1e-6) break
    tolerance <- sum(vapply(p$graphs, function(g) 1e-6/g$minimum, numeric(1)))
    if (is.finite(z$stats[["best_bound"]]) && z$stats[["best_bound"]] > tolerance) break
    if (isTRUE(cfg$repair$enabled)) {
      ref <- read.delim(best$selected)$unit
      repaired <- record(dp_stage(cfg, root, dat, paste0("repair_", i), cfg$repair$seconds,
        best$checkpoint, "repair", reference = ref, allow_no_incumbent = TRUE))
      if (isTRUE(repaired$available)) { best <- repaired; break }
      if (repaired$stats[["solve_result_num"]] >= 200 && repaired$stats[["solve_result_num"]] < 300) break
    }
  }
  dp_write(if (!length(p$graphs)) "NOT_NEEDED_NO_NETWORKS" else
    if (any(vapply(stages, function(s) s$phase == "repair", logical(1)))) "ATTEMPTED_SEE_STAGE_SUMMARY" else
      if (!isTRUE(cfg$repair$enabled)) "DISABLED" else if (best$stats[["shortfall_ha"]] <= 1e-6)
        "SKIPPED_ALREADY_ZERO_SHORTFALL" else "SKIPPED_ZERO_SHORTFALL_EXCLUDED_BY_BOUND",
    file.path(root, "repair_status.txt"))
  if (best$stats[["shortfall_ha"]] > 1e-6) {
    bounds <- vapply(Filter(function(s) s$phase == "feasibility", stages), function(s) s$stats[["best_bound"]], numeric(1))
    tol <- sum(vapply(p$graphs, function(g) 1e-6/g$minimum, numeric(1)))
    proven <- any(is.finite(bounds) & bounds > tol) || any(vapply(stages, function(s)
      s$phase == "repair" && !isTRUE(s$available) && s$stats[["solve_result_num"]] >= 200 && s$stats[["solve_result_num"]] < 300, logical(1)))
    dp_write(if (proven) "ZERO_SHORTFALL_EXCLUDED_BY_SOLVER_BOUND" else "INCOMPLETE_SHORTFALL_TIME_BUDGET", statusfile)
    dp_write(best$folder, file.path(root, "latest_partial_stage.txt"))
    stop("No zero-shortfall solution. Partial result: ", best$folder,
      if (proven) ". Solver evidence excludes zero shortfall." else ". Extend feasibility_seconds or repair$seconds and rerun.")
  }
  # Hold this reference constant across ALL quality continuation stages.
  reference <- read.delim(best$selected)$unit
  dp_csv(read.delim(best$selected), file.path(root, "quality_reference.csv"))
  for (i in seq_along(cfg$solver$quality_seconds)) {
    z <- record(dp_stage(cfg, root, dat, paste0("quality_", i), cfg$solver$quality_seconds[i], best$checkpoint, "quality", reference = reference))
    if (z$stats[["quality"]] + 1e-10 >= best$stats[["quality"]]) best <- z
  }
  final <- file.path(root, "final", basename(best$folder)); dir.create(final, recursive = TRUE, showWarnings = FALSE)
  selected <- read.delim(best$selected, stringsAsFactors = FALSE)
  audit_solution(cfg, p, selected, final)
  changes <- dp_footprint_change(p$units, reference, selected$unit)
  cap <- if (cfg$quality$footprint == "fixed") 0 else if (cfg$quality$footprint == "capped") cfg$quality$max_change_ha else Inf
  if (changes$total_change_ha > cap + 1e-6) stop("Export fails the final footprint change limit.")
  dp_csv(changes, file.path(final, "footprint_change_audit.csv"))
  file.copy(best$selected, file.path(final, "selected.tsv"), overwrite = TRUE)
  dp_write(final, file.path(root, "latest_validated_solution.txt")); dp_write("VALIDATED", statusfile)
  message("Validated: ", final)
  invisible(list(root = root, final = final, selected = selected, stages = stages))
}
dp_footprint_change <- function(units, reference, selected) {
  area <- function(ids) sum(units$area_ha[match(ids, units$unit)])
  added <- area(setdiff(selected, reference)); removed <- area(setdiff(reference, selected))
  data.frame(added_ha = added, removed_ha = removed, total_change_ha = added + removed)
}
compare_footprints <- function(first_raster, second_raster, output_csv = NULL) {
  a <- terra::rast(first_raster); b <- terra::rast(second_raster); terra::compareGeom(a, b)
  if (dp_any(is.na(a) != is.na(b))) stop("Footprint domains differ; align the planning masks first.")
  area <- prod(terra::res(a))/10000
  both <- dp_sum(a > 0 & b > 0)*area; only_a <- dp_sum(a > 0 & b == 0)*area; only_b <- dp_sum(a == 0 & b > 0)*area
  union <- both + only_a + only_b
  ans <- data.frame(shared_ha = both, first_only_ha = only_a, second_only_ha = only_b,
    area_weighted_jaccard = if (union > 0) both/union else NA_real_)
  if (!is.null(output_csv)) dp_csv(ans, output_csv)
  ans
}
