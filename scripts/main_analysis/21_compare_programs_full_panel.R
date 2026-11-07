# Resolve repository paths before running the recorded analysis blocks.
if (!exists("portfolio_repo_dir", inherits = TRUE)) {
  portfolio_repo_dir <- normalizePath(Sys.getenv("RETINA_PROJECT_DIR", unset = getwd()), mustWork = TRUE)
}

# Display the actual flight-enriched pathways across every terrestrial contrast.
# Reuse the saved NES and FDR values; retain nonsignificant cells for context.

# Label direction separately from whether the terrestrial FDR meets the existing cutoff.
main_panel_evidence <- function(osd_nes, terrestrial_nes, terrestrial_fdr) {
  stopifnot(length(osd_nes) == length(terrestrial_nes),
    length(osd_nes) == length(terrestrial_fdr),
    all(is.finite(osd_nes)), all(is.finite(terrestrial_nes)),
    all(is.finite(terrestrial_fdr) & terrestrial_fdr >= 0 & terrestrial_fdr <= 1),
    all(osd_nes != 0))
  pattern <- ifelse(osd_nes > 0 & terrestrial_nes > 0, "Positive in both",
    ifelse(osd_nes < 0 & terrestrial_nes < 0, "Negative in both",
      ifelse(osd_nes > 0 & terrestrial_nes < 0, "Positive flight / negative terrestrial",
        ifelse(osd_nes < 0 & terrestrial_nes > 0,
          "Negative flight / positive terrestrial", "Zero terrestrial NES"))))
  ifelse(terrestrial_fdr < 0.05, pattern, "Terrestrial FDR >= 0.05")
}

# Draw named pathways with one shared NES color scale and the original FDR stars.
main_panel_plot <- function(nes, fdr, pathway_names, contrast_labels, color_limit, title, file) {
  stopifnot(is.matrix(nes), identical(dim(nes), dim(fdr)), nrow(nes) > 0L,
    nrow(nes) == length(pathway_names), ncol(nes) == length(contrast_labels),
    all(is.finite(nes)), all(is.finite(fdr)), color_limit > 0)
  wrapped_names <- vapply(pathway_names, function(x) paste(strwrap(x, width = 55),
    collapse = "\n"), character(1))
  name_lines <- max(lengths(strsplit(wrapped_names, "\n", fixed = TRUE)))
  colors <- grDevices::colorRampPalette(c("#2166AC", "#F7F7F7", "#B35806"))(101)
  height <- 620 + nrow(nes) * max(76, name_lines * 31)
  grDevices::png(file, width = 3400, height = height, res = 180, pointsize = 14)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(oma = c(5.2, 0, 4.4, 0))
  graphics::layout(matrix(1:3, nrow = 1), widths = c(1.20, 1.30, 0.12))
  y <- rev(seq_len(nrow(nes)))

  # Keep labels and cells in the same original pathway order without clustering.
  graphics::par(mar = c(5, 0.2, 0.5, 0.2))
  graphics::plot.new()
  graphics::plot.window(xlim = c(0, 1), ylim = c(0.5, nrow(nes) + 0.5),
    xaxs = "i", yaxs = "i")
  graphics::text(0.99, y, labels = wrapped_names, adj = c(1, 0.5), cex = 0.86)

  graphics::par(mar = c(5, 0, 0.5, 0))
  graphics::plot.new()
  graphics::plot.window(xlim = c(0.5, ncol(nes) + 0.5),
    ylim = c(0.5, nrow(nes) + 0.5), xaxs = "i", yaxs = "i")
  for (j in seq_len(ncol(nes))) {
    color_index <- pmax(1L, pmin(101L,
      1L + floor((nes[, j] + color_limit) / (2 * color_limit) * 100)))
    graphics::rect(j - 0.5, y - 0.5, j + 0.5, y + 0.5,
      col = colors[color_index], border = "white", lwd = 0.6)
    marked <- fdr[, j] < 0.05
    if (any(marked)) {
      graphics::text(j, y[marked], labels = "*", cex = 1.4,
        col = ifelse(abs(nes[marked, j]) > color_limit * 0.65, "white", "black"))
    }
  }
  graphics::abline(v = 1.5, lwd = 2, col = "#333333")
  graphics::box(col = "#777777")
  label_lines <- strsplit(contrast_labels, "\n", fixed = TRUE)
  for (k in seq_len(max(lengths(label_lines)))) {
    line_labels <- vapply(label_lines, function(x) if (length(x) >= k) x[k] else "", character(1))
    graphics::mtext(line_labels, side = 1, at = seq_len(ncol(nes)),
      line = k - 0.15, cex = 0.72)
  }

  # Use the same symmetric legend in both figures; positive/negative refer to gene ranks.
  graphics::par(mar = c(5, 0.8, 0.5, 2.4))
  graphics::plot.new()
  graphics::plot.window(xlim = c(0, 1), ylim = c(-color_limit, color_limit),
    xaxs = "i", yaxs = "i")
  limits <- seq(-color_limit, color_limit, length.out = 102)
  graphics::rect(0.15, limits[-length(limits)], 0.7, limits[-1],
    col = colors, border = NA)
  legend_at <- pretty(c(-color_limit, color_limit), n = 5)
  legend_at <- legend_at[abs(legend_at) <= color_limit]
  graphics::axis(4, at = legend_at, las = 1, cex.axis = 0.8)
  graphics::mtext("NES", side = 3, line = 0.4, cex = 0.85)
  graphics::mtext(title, side = 3, outer = TRUE, line = 2.5, cex = 1.25, font = 2)
  graphics::mtext(paste(nrow(nes), "pathways | Flight 0G versus HGC reference | All terrestrial contrasts kept separate"),
    side = 3, outer = TRUE, line = 0.9, cex = 0.9)
  graphics::mtext("* FDR < 0.05 in that contrast. Compare each terrestrial color with the flight column.",
    side = 1, outer = TRUE, line = 1, cex = 0.88)
  graphics::mtext("Unstarred cells retain their NES; they do not establish no response. Colors do not measure pathway activity directly.",
    side = 1, outer = TRUE, line = 2.4, cex = 0.83)
  graphics::mtext("Rows include every flight-enriched pathway in this direction. Related pathways can share genes.",
    side = 1, outer = TRUE, line = 3.8, cex = 0.83)
}

# Restore the reviewed comparison checkpoint without rerunning enrichment or similarity.
if (!exists("main_analysis_dir")) {
  main_analysis_dir <- Sys.getenv("MAIN_ANALYSIS_DIR", unset = file.path(portfolio_repo_dir, "results", "main_analysis"))
}
main_panel_input <- file.path(main_analysis_dir,
  "Main_analysis_pathway_comparison_checkpoint_2026-11-06.rds")
stopifnot(file.exists(main_panel_input))
main_panel_context <- readRDS(main_panel_input)
main_panel_reference <- "OSD758__FLT_0G_vs_HGC"
main_panel_manifest <- main_panel_context$contrasts
main_panel_manifest_keys <- paste(main_panel_manifest$Study,
  main_panel_manifest$Comparison, sep = "__")
main_panel_terrestrial <- main_panel_manifest[main_panel_manifest$Role == "Terrestrial", , drop = FALSE]
main_panel_terrestrial_keys <- paste(main_panel_terrestrial$Study,
  main_panel_terrestrial$Comparison, sep = "__")
main_panel_column_keys <- c(main_panel_reference, main_panel_terrestrial_keys)
main_panel_nes <- main_panel_context$NES
main_panel_fdr <- main_panel_context$FDR
main_panel_key <- main_panel_context$pathway_key
stopifnot(length(main_panel_context$checks) > 0L, all(main_panel_context$checks),
  is.matrix(main_panel_nes), is.matrix(main_panel_fdr),
  identical(dimnames(main_panel_nes), dimnames(main_panel_fdr)),
  identical(rownames(main_panel_nes), main_panel_key$Pathway_ID),
  !anyNA(main_panel_key$Pathway_ID), anyDuplicated(main_panel_key$Pathway_ID) == 0L,
  !anyNA(main_panel_key$Pathway_name),
  identical(colnames(main_panel_nes), main_panel_manifest_keys),
  anyDuplicated(main_panel_manifest_keys) == 0L,
  sum(main_panel_manifest$Role == "Reference") == 1L,
  main_panel_manifest$Role[match(main_panel_reference, main_panel_manifest_keys)] == "Reference",
  length(main_panel_terrestrial_keys) == 12L,
  setequal(main_panel_column_keys, main_panel_manifest_keys),
  all(is.finite(main_panel_nes)),
  all(is.finite(main_panel_fdr) & main_panel_fdr >= 0 & main_panel_fdr <= 1))

# Display every pathway meeting the existing flight FDR cutoff, ordered only by flight results.
main_panel_reference_rows <- which(main_panel_fdr[, main_panel_reference] < 0.05)
stopifnot(length(main_panel_reference_rows) > 0L,
  all(main_panel_nes[main_panel_reference_rows, main_panel_reference] != 0))
main_panel_reference_rows <- main_panel_reference_rows[order(
  main_panel_nes[main_panel_reference_rows, main_panel_reference] < 0,
  main_panel_fdr[main_panel_reference_rows, main_panel_reference],
  main_panel_key$Pathway_ID[main_panel_reference_rows])]
main_panel_key <- main_panel_key[main_panel_reference_rows, , drop = FALSE]
rownames(main_panel_key) <- NULL
main_panel_nes <- main_panel_nes[main_panel_reference_rows, main_panel_column_keys, drop = FALSE]
main_panel_fdr <- main_panel_fdr[main_panel_reference_rows, main_panel_column_keys, drop = FALSE]
main_panel_key$Flight_NES <- unname(main_panel_nes[, main_panel_reference])
main_panel_key$Flight_FDR <- unname(main_panel_fdr[, main_panel_reference])

# Keep the matched names, both NES/FDR values and direction labels for every displayed cell.
main_panel_tables <- lapply(seq_along(main_panel_terrestrial_keys), function(i) {
  key <- main_panel_terrestrial_keys[i]
  evidence <- main_panel_evidence(main_panel_key$Flight_NES,
    main_panel_nes[, key], main_panel_fdr[, key])
  data.frame(Study = main_panel_terrestrial$Study[i],
    Comparison = main_panel_terrestrial$Comparison[i], main_panel_key,
    Terrestrial_NES = unname(main_panel_nes[, key]),
    Terrestrial_FDR = unname(main_panel_fdr[, key]),
    NES_direction = ifelse(main_panel_key$Flight_NES * main_panel_nes[, key] > 0,
      "Same", ifelse(main_panel_key$Flight_NES * main_panel_nes[, key] < 0, "Opposite", "Zero_NES")),
    Evidence = evidence, stringsAsFactors = FALSE, row.names = NULL)
})
main_panel_long <- do.call(rbind, main_panel_tables)
rownames(main_panel_long) <- NULL
main_panel_summary <- do.call(rbind, lapply(main_panel_tables, function(x) {
  data.frame(Study = x$Study[1], Comparison = x$Comparison[1],
    Flight_pathways = nrow(x),
    Positive_in_both = sum(x$Evidence == "Positive in both"),
    Negative_in_both = sum(x$Evidence == "Negative in both"),
    Opposite_in_both = sum(x$Terrestrial_FDR < 0.05 & x$NES_direction == "Opposite"),
    Terrestrial_FDR_not_met = sum(x$Terrestrial_FDR >= 0.05),
    Zero_NES_in_both = sum(x$Terrestrial_FDR < 0.05 & x$NES_direction == "Zero_NES"),
    stringsAsFactors = FALSE)
}))
rownames(main_panel_summary) <- NULL

# Verify the display against the original paired values and previously reviewed overlap counts.
main_panel_old_pairs <- main_panel_context$pathway_pairs
main_panel_pair_key <- paste(main_panel_long$Study, main_panel_long$Comparison,
  main_panel_long$Pathway_ID, sep = "__")
main_panel_old_pair_key <- paste(main_panel_old_pairs$Study, main_panel_old_pairs$Comparison,
  main_panel_old_pairs$Pathway_ID, sep = "__")
main_panel_pair_index <- match(main_panel_pair_key, main_panel_old_pair_key)
main_panel_old_summary <- main_panel_context$overlap_summary
main_panel_summary_index <- match(main_panel_terrestrial_keys,
  paste(main_panel_old_summary$Study, main_panel_old_summary$Comparison, sep = "__"))
stopifnot(!anyNA(main_panel_pair_index), !anyNA(main_panel_summary_index))
main_panel_old_summary <- main_panel_old_summary[main_panel_summary_index, , drop = FALSE]
main_panel_checks <- c(
  All_terrestrial_contrasts = identical(paste(main_panel_summary$Study,
    main_panel_summary$Comparison, sep = "__"), main_panel_terrestrial_keys),
  Every_flight_enriched_pathway = nrow(main_panel_long) == nrow(main_panel_key) * 12L &&
    anyDuplicated(main_panel_pair_key) == 0L,
  Original_values = identical(main_panel_long$Flight_NES, main_panel_old_pairs$OSD_NES[main_panel_pair_index]) &&
    identical(main_panel_long$Flight_FDR, main_panel_old_pairs$OSD_FDR[main_panel_pair_index]) &&
    identical(main_panel_long$Terrestrial_NES, main_panel_old_pairs$Terrestrial_NES[main_panel_pair_index]) &&
    identical(main_panel_long$Terrestrial_FDR, main_panel_old_pairs$Terrestrial_FDR[main_panel_pair_index]) &&
    identical(main_panel_long$Pathway_name, main_panel_old_pairs$Pathway_name[main_panel_pair_index]) &&
    identical(main_panel_long$NES_direction, main_panel_old_pairs$NES_direction[main_panel_pair_index]),
  Reviewed_overlap_counts = all(main_panel_summary$Positive_in_both + main_panel_summary$Negative_in_both ==
      main_panel_old_summary$Same_direction_both) &&
    all(main_panel_summary$Opposite_in_both == main_panel_old_summary$Opposite_direction_both) &&
    all(main_panel_summary$Zero_NES_in_both == main_panel_old_summary$Zero_NES_both) &&
    all(main_panel_summary$Terrestrial_FDR_not_met == main_panel_old_summary$OSD_only_FDR005),
  Complete_accounting = all(rowSums(main_panel_summary[, c("Positive_in_both", "Negative_in_both",
    "Opposite_in_both", "Terrestrial_FDR_not_met", "Zero_NES_in_both")]) == nrow(main_panel_key))
)
print(main_panel_summary, row.names = FALSE)
print(main_panel_checks)
stopifnot(all(main_panel_checks))

# Save a readable pathway-by-contrast evidence table; counts refer to contrasts, not independent studies.
main_panel_evidence_matrix <- vapply(main_panel_tables, function(x) x$Evidence,
  character(nrow(main_panel_key)))
dimnames(main_panel_evidence_matrix) <- list(main_panel_key$Pathway_ID, main_panel_terrestrial_keys)
main_panel_wide <- main_panel_key
main_panel_wide$Same_direction_contrasts <- rowSums(main_panel_evidence_matrix == "Positive in both" |
  main_panel_evidence_matrix == "Negative in both")
main_panel_wide$Opposite_direction_contrasts <- rowSums(
  main_panel_evidence_matrix == "Positive flight / negative terrestrial" |
    main_panel_evidence_matrix == "Negative flight / positive terrestrial")
main_panel_wide <- data.frame(main_panel_wide, main_panel_evidence_matrix, check.names = FALSE)
main_panel_files <- file.path(main_analysis_dir, c(
  "Main_analysis_OSD758_programs_full_panel_values.csv",
  "Main_analysis_OSD758_programs_full_panel_evidence.csv",
  "Main_analysis_OSD758_programs_full_panel_summary.csv",
  "Main_analysis_OSD758_positive_programs_full_panel.png",
  "Main_analysis_OSD758_negative_programs_full_panel.png",
  "Main_analysis_programs_full_panel_checkpoint_2026-11-06.rds",
  "Main_analysis_OSD758_programs_full_panel_method_2026-11-06.txt"
))
utils::write.csv(main_panel_long, main_panel_files[1], row.names = FALSE)
utils::write.csv(main_panel_wide, main_panel_files[2], row.names = FALSE)
utils::write.csv(main_panel_summary, main_panel_files[3], row.names = FALSE)

# Preserve stages and genotypes in fixed manifest order, with no disease-level averaging.
main_panel_label_key <- c(OSD758__FLT_0G_vs_HGC = "Flight\n0G / HGC",
  Glaucoma__3D_vs_Naive = "Glaucoma\n3 days", Glaucoma__2W_vs_Naive = "Glaucoma\n2 weeks",
  Glaucoma__6W_vs_Naive = "Glaucoma\n6 weeks", RP__P23H_vs_WT = "RP\nP23H",
  LHON__YARS2_vs_WT = "LHON\nYARS2", LHON__COI_vs_WT = "LHON\nCOI",
  LHON__YARS2_COI_vs_WT = "LHON\nYARS2 +\nCOI", DR__Diabetic_vs_Normal = "DR\nDiabetic",
  AMD__Young_CfhKO_vs_WT = "AMD\nYoung\nCfhKO", AMD__Aged_CfhKO_vs_WT = "AMD\nAged\nCfhKO",
  Autoimmune_uveitis__EAU_vs_CON = "Autoimmune\nuveitis", Light_damage__LD_vs_CTRL = "Light\ndamage")
stopifnot(all(main_panel_column_keys %in% names(main_panel_label_key)))
main_panel_labels <- unname(main_panel_label_key[main_panel_column_keys])
main_panel_color_limit <- max(1, ceiling(max(abs(main_panel_nes)) * 2) / 2)
main_panel_positive <- main_panel_key$Flight_NES > 0
main_panel_negative <- main_panel_key$Flight_NES < 0
main_panel_plot(main_panel_nes[main_panel_positive, , drop = FALSE],
  main_panel_fdr[main_panel_positive, , drop = FALSE], main_panel_key$Pathway_name[main_panel_positive],
  main_panel_labels, main_panel_color_limit, "OSD-758 positive pathways across the terrestrial panel", main_panel_files[4])
main_panel_plot(main_panel_nes[main_panel_negative, , drop = FALSE],
  main_panel_fdr[main_panel_negative, , drop = FALSE], main_panel_key$Pathway_name[main_panel_negative],
  main_panel_labels, main_panel_color_limit, "OSD-758 negative pathways across the terrestrial panel", main_panel_files[5])

# Record display rules and original data provenance without adding a new statistical test.
main_panel_method <- list(
  purpose = "Descriptive inspection of actual shared/opposing pathway programs across all terrestrial contrasts",
  reference = main_panel_reference, selection = "Every pathway with the original flight FDR < 0.05",
  selection_basis = "Flight response only; no selection using a preferred terrestrial disease",
  plot_order = "Positive/negative flight NES split; within direction, increasing flight FDR then pathway ID; no clustering",
  plot_cells = "Original continuous NES for every selected pathway and contrast; * marks that contrast's original FDR < 0.05",
  color_limit = main_panel_color_limit,
  evidence = "Same/opposite nonzero NES signs with FDR < 0.05 separately in both contrasts; zero NES labelled separately",
  nonsignificance = "An unstarred cell is not evidence of no response or a significant between-study difference",
  interpretation = "NES indicates directional enrichment in existing gene rankings, not directly measured pathway activity",
  dependencies = "Overlapping/nested pathways and related contrasts are not independent mechanisms or independent studies",
  scope = "All original pathway correlations still use all 799 pathways; this display does not replace those comparisons",
  inference = "No rerun of enrichment, new FDR adjustment, correlation, overlap significance or rank-uncertainty test",
  documentation = "https://bioconductor.org/packages/release/bioc/vignettes/fgsea/inst/doc/fgsea-tutorial.html"
)
saveRDS(list(pathway_key = main_panel_key, NES = main_panel_nes, FDR = main_panel_fdr,
  values = main_panel_long, evidence = main_panel_wide, summary = main_panel_summary,
  checks = main_panel_checks, method = main_panel_method,
  input_file = main_panel_input, input_md5 = tools::md5sum(main_panel_input)), main_panel_files[6])
utils::capture.output(print(main_panel_method), file = main_panel_files[7])
print(file.exists(main_panel_files))
