# Tests the actual AMPL repair and quality formulations on a transparent graph.
source("denmark_pipeline/run.R")
root <- tempfile("repair_workflow_", tmpdir = getwd()); dir.create(root); root <- normalizePath(root, winslash = "/")
cfg <- denmark_config(); cfg$solver$name <- "highs"; cfg$solver$ampl <- Sys.getenv("AMPL_EXECUTABLE", "ampl")
cfg$solver$options <- "tech:outlev=0"; cfg$solver$gap <- 0; cfg$repair$gap <- 0
cfg$objective$warmup_contiguity <- 7; cfg$objective$quality_contiguity <- 0
cfg$model$outcomes <- list(a=list(),b=list())
units <- data.frame(unit=1:4, area_ha=1, cost=1)
actions <- expand.grid(outcome=c("a","b"),unit=1:4,stringsAsFactors=FALSE)
actions$score <- rep(c(1,2,10,8),each=2) + ifelse(actions$outcome=="b",.1,0)
nodes <- data.frame(id=10*(1:5),unit=c(0,1:4),anchor=c(1,rep(0,4)),area_ha=1,flow_cap=5)
edges <- data.frame(from=c(10,20,30,10),to=c(20,30,40,50))
g <- list(nodes=nodes,edges=edges,minimum=3)
p <- list(units=units,actions=actions,graphs=list(Global=g),spatial=g,
  ledger=list(minimum_restoration_ha=2,maximum_restoration_ha=2,cost_cap=2))
dat <- dp_ampl_data(cfg,p,file.path(root,"problem.dat"))
# Units 3 and 4 have the best scores, but 3 is disconnected without 1 and 2.
# Repair must swap 3 for 1. Retaining 4 is the unique minimum-change repair.
repair <- dp_stage(cfg,root,dat,"repair_direct",30,phase="repair",reference=c(3,4))
chosen <- read.delim(repair$selected)
stopifnot(setequal(chosen$unit,c(1,4)),repair$stats[["shortfall_ha"]]==0,
  abs(repair$stats[["footprint_change_ha"]]-2)<1e-9)
stopifnot(all(dp_footprint_change(units,c(3,4),chosen$unit)==data.frame(added_ha=1,removed_ha=1,total_change_ha=2)))
# A different feasible reference (1+2) makes the final-stage choices explicit.
ref <- c(1,2)
for (mode in c("free","capped","fixed")) {
  cfg$quality$footprint <- mode; cfg$quality$max_change_ha <- 1
  z <- dp_stage(cfg,root,dat,paste0("quality_",mode),30,phase="quality",reference=ref)
  x <- read.delim(z$selected)
  stopifnot(setequal(x$unit,if(mode=="free") c(1,4) else ref),all(x$outcome=="b"))
}
cfg$quality$footprint <- "capped"; cfg$quality$max_change_ha <- 2
z <- dp_stage(cfg,root,dat,"quality_cap_two",30,phase="quality",reference=ref)
stopifnot(setequal(read.delim(z$selected)$unit,c(1,4)))
# No reference-selected units is legitimate: repair selects the necessary land.
empty <- dp_stage(cfg,root,dat,"repair_empty",30,phase="repair",reference=integer())
stopifnot(empty$stats[["shortfall_ha"]]==0,empty$stats[["footprint_change_ha"]]==2)
# Repair cannot invent a bridge through disconnected fragments or exceed a budget.
p$ledger$maximum_restoration_ha <- 1; p$ledger$minimum_restoration_ha <- 1
bad_dat <- dp_ampl_data(cfg,p,file.path(root,"impossible.dat"))
bad <- dp_stage(cfg,root,bad_dat,"repair_infeasible",30,phase="repair",reference=4,allow_no_incumbent=TRUE)
stopifnot(!bad$available,bad$stats[["solve_result_num"]]>=200,bad$stats[["solve_result_num"]]<300)
# Both phase weights are used, while the exported 'quality' always uses the
# final-phase weight and is comparable across warm-up and later stages.
cfg$quality$footprint <- "free"
w <- dp_stage(cfg,root,dat,"weight_warmup",30,phase="warmup")
q <- dp_stage(cfg,root,dat,"weight_quality",30,phase="quality",reference=ref)
stopifnot(w$stats[["spatial_weight"]]==7,q$stats[["spatial_weight"]]==0)
for (z in list(w,q)) {
  x <- read.delim(z$selected)
  expected <- sum(x$score)/sum(vapply(split(abs(actions$score),actions$unit),max,numeric(1)))
  stopifnot(abs(z$stats[["quality"]]-expected)<1e-9)
}
cat("PASS: minimum-change repair, budget-preserving swaps, empty reference, failed repair, free/capped/fixed footprints, pathway changes, separate spatial weights.\n")
