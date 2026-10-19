# Compare independent models only after they have been solved and audited.
# Tables contain selected restoration cells only. A full decision table must
# first be filtered to its selected rows. candidate_dictionary is the shared
# spatial dictionary, with coarse_cell_id and EligibleAreaHa.
compare_gudenaa_footprints <- function(eight_habitat_selected, binary_selected,
  candidate_dictionary, eight_grid_id, binary_grid_id,
  subarea_dictionary = NULL, output_dir = NULL) {
  read_table <- function(x) {
    if (is.data.frame(x)) return(x)
    if (tolower(tools::file_ext(x)) == "csv") read.csv(x, stringsAsFactors = FALSE) else
      read.delim(x, stringsAsFactors = FALSE)
  }
  if (length(eight_grid_id) != 1L || is.na(eight_grid_id) || !nzchar(eight_grid_id) ||
      !identical(eight_grid_id, binary_grid_id))
    stop("Both runs must have the same verified grid/domain identifier; regrid before comparing.")
  a <- read_table(eight_habitat_selected); b <- read_table(binary_selected)
  d <- read_table(candidate_dictionary)
  if (!all(c("coarse_cell_id", "EligibleAreaHa") %in% names(d))) stop("Missing candidate columns.")
  if (anyDuplicated(d$coarse_cell_id) || anyNA(d$coarse_cell_id) ||
      any(!is.finite(d$EligibleAreaHa) | d$EligibleAreaHa <= 0)) stop("Invalid candidate dictionary.")
  check_ids <- function(x) {
    if (!"coarse_cell_id" %in% names(x)) stop("Selection table needs coarse_cell_id.")
    if (anyNA(x$coarse_cell_id) || anyDuplicated(x$coarse_cell_id)) stop("Missing or duplicate selected cell IDs.")
    if (!all(x$coarse_cell_id %in% d$coarse_cell_id)) stop("Selected cells are outside the shared candidate dictionary.")
    if ("restored_ha" %in% names(x) && any(!is.finite(x$restored_ha) |
        abs(x$restored_ha - d$EligibleAreaHa[match(x$coarse_cell_id, d$coarse_cell_id)]) > 1e-6))
      stop("The models restore different areas within a cell; compare a common fine-resolution mask instead.")
    x$coarse_cell_id
  }
  ia <- check_ids(a); ib <- check_ids(b)
  d$eight_selected <- d$coarse_cell_id %in% ia
  d$binary_selected <- d$coarse_cell_id %in% ib
  d$agreement <- ifelse(d$eight_selected & d$binary_selected, "both",
    ifelse(d$eight_selected, "eight_only", ifelse(d$binary_selected, "binary_only", "neither")))
  summarise <- function(x, area) {
    both <- sum(area[x$agreement == "both"])
    only_a <- sum(area[x$agreement == "eight_only"])
    only_b <- sum(area[x$agreement == "binary_only"])
    union <- both + only_a + only_b
    safe <- function(num, den) if (den > 0) num/den else NA_real_
    data.frame(eight_selected_ha = both + only_a, binary_selected_ha = both + only_b,
      both_ha = both, eight_only_ha = only_a, binary_only_ha = only_b,
      symmetric_difference_ha = only_a + only_b, union_ha = union,
      area_weighted_jaccard = safe(both, union),
      fraction_of_eight_recovered = safe(both, both + only_a),
      fraction_of_binary_recovered = safe(both, both + only_b))
  }
  overall <- summarise(d, d$EligibleAreaHa)
  by_subarea <- NULL
  if (!is.null(subarea_dictionary)) {
    s <- read_table(subarea_dictionary)
    cols <- c("subarea_id", "coarse_cell_id", "RestorableAreaHaBySubarea")
    if (!all(cols %in% names(s))) stop("Missing subarea dictionary columns.")
    if (anyNA(s[, cols]) || any(!is.finite(s$RestorableAreaHaBySubarea)) ||
        any(s$RestorableAreaHaBySubarea < 0) || anyDuplicated(s[c("subarea_id", "coarse_cell_id")]) ||
        !all(s$coarse_cell_id %in% d$coarse_cell_id)) stop("Invalid subarea dictionary.")
    totals <- tapply(s$RestorableAreaHaBySubarea, s$coarse_cell_id, sum)
    if (anyNA(totals[as.character(d$coarse_cell_id)]) ||
        any(abs(totals[as.character(d$coarse_cell_id)] - d$EligibleAreaHa) > 1e-6))
      stop("Subarea areas do not reconcile with the shared candidate dictionary.")
    s <- merge(s, d[c("coarse_cell_id", "agreement")], by = "coarse_cell_id", all.x = TRUE)
    by_subarea <- do.call(rbind, lapply(split(s, s$subarea_id), function(z)
      cbind(subarea_id = z$subarea_id[1], summarise(z, z$RestorableAreaHaBySubarea))))
    rownames(by_subarea) <- NULL
  }
  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    write.csv(overall, file.path(output_dir, "footprint_comparison.csv"), row.names = FALSE)
    write.csv(d, file.path(output_dir, "cell_agreement.csv"), row.names = FALSE)
    if (!is.null(by_subarea)) write.csv(by_subarea, file.path(output_dir, "subarea_comparison.csv"), row.names = FALSE)
  }
  list(summary = overall, cells = d, subareas = by_subarea)
}
