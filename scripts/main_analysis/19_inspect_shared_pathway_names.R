# Resolve repository paths before running the recorded analysis blocks.
if (!exists("portfolio_repo_dir", inherits = TRUE)) {
  portfolio_repo_dir <- normalizePath(Sys.getenv("RETINA_PROJECT_DIR", unset = getwd()), mustWork = TRUE)
}

# Inspect the actual pathway names underlying the reviewed overlap counts.
# Select pathways with the original FDR < 0.05 in both responses for this display.

# Restore the saved pathway comparison without repeating enrichment or correlation.
if (!exists("main_analysis_dir")) {
  main_analysis_dir <- Sys.getenv("MAIN_ANALYSIS_DIR", unset = file.path(portfolio_repo_dir, "results", "main_analysis"))
}
main_shared_names_input_file <- file.path(
  main_analysis_dir, "Main_analysis_pathway_comparison_checkpoint_2026-11-06.rds"
)
stopifnot(file.exists(main_shared_names_input_file))
main_shared_names_context <- readRDS(main_shared_names_input_file)
main_shared_pair_rows <- main_shared_names_context$pathway_pairs
stopifnot(
  length(main_shared_names_context$checks) > 0L, all(main_shared_names_context$checks),
  all(c("Study", "Comparison", "Pathway_ID", "Pathway_name", "OSD_NES", "OSD_FDR",
    "Terrestrial_NES", "Terrestrial_FDR", "NES_direction", "Both_enriched") %in% names(main_shared_pair_rows)),
  identical(main_shared_pair_rows$Both_enriched,
    main_shared_pair_rows$OSD_FDR < 0.05 & main_shared_pair_rows$Terrestrial_FDR < 0.05)
)

# Keep every both-enriched pathway record across all terrestrial contrasts.
main_both_enriched_pathway_names <- main_shared_pair_rows[
  main_shared_pair_rows$Both_enriched, , drop = FALSE
]
main_both_enriched_pathway_names$Direction_pattern <- with(main_both_enriched_pathway_names,
  ifelse(OSD_NES > 0 & Terrestrial_NES > 0, "Positive in both",
    ifelse(OSD_NES < 0 & Terrestrial_NES < 0, "Negative in both",
      ifelse(OSD_NES > 0 & Terrestrial_NES < 0, "Positive OSD / negative terrestrial",
        ifelse(OSD_NES < 0 & Terrestrial_NES > 0, "Negative OSD / positive terrestrial", "Zero NES")))))
main_shared_direction_levels <- c(
  "Positive in both", "Negative in both", "Positive OSD / negative terrestrial",
  "Negative OSD / positive terrestrial", "Zero NES"
)
main_shared_contrast_info <- main_shared_names_context$contrasts[
  main_shared_names_context$contrasts$Role == "Terrestrial", , drop = FALSE
]
main_shared_contrast_keys <- paste(
  main_shared_contrast_info$Study, main_shared_contrast_info$Comparison, sep = "__"
)
main_both_enriched_pathway_names <- main_both_enriched_pathway_names[order(
  match(paste(main_both_enriched_pathway_names$Study,
    main_both_enriched_pathway_names$Comparison, sep = "__"), main_shared_contrast_keys),
  match(main_both_enriched_pathway_names$Direction_pattern, main_shared_direction_levels),
  main_both_enriched_pathway_names$OSD_FDR,
  main_both_enriched_pathway_names$Terrestrial_FDR,
  main_both_enriched_pathway_names$Pathway_ID
), , drop = FALSE]
rownames(main_both_enriched_pathway_names) <- NULL

# Separate matching and opposing directions while retaining the original pathway IDs and values.
main_same_direction_pathway_names <- main_both_enriched_pathway_names[
  main_both_enriched_pathway_names$NES_direction == "Same", , drop = FALSE
]
main_opposite_direction_pathway_names <- main_both_enriched_pathway_names[
  main_both_enriched_pathway_names$NES_direction == "Opposite", , drop = FALSE
]

# Begin manual interpretation with the leading observed 2-week glaucoma contrast.
main_glaucoma_2W_shared_pathways <- main_both_enriched_pathway_names[
  main_both_enriched_pathway_names$Study == "Glaucoma" &
    main_both_enriched_pathway_names$Comparison == "2W_vs_Naive", , drop = FALSE
]
main_glaucoma_2W_overlap_record <- main_shared_names_context$overlap_summary[
  main_shared_names_context$overlap_summary$Study == "Glaucoma" &
    main_shared_names_context$overlap_summary$Comparison == "2W_vs_Naive", , drop = FALSE
]
stopifnot(
  nrow(main_glaucoma_2W_overlap_record) == 1L,
  nrow(main_glaucoma_2W_shared_pathways) == main_glaucoma_2W_overlap_record$Both_FDR005,
  sum(main_glaucoma_2W_shared_pathways$NES_direction == "Same") ==
    main_glaucoma_2W_overlap_record$Same_direction_both,
  sum(main_glaucoma_2W_shared_pathways$NES_direction == "Opposite") ==
    main_glaucoma_2W_overlap_record$Opposite_direction_both
)
cat("OSD-758 flight 0G versus HGC and glaucoma 2W versus naive\n")
print(c(
  Enriched_in_both = nrow(main_glaucoma_2W_shared_pathways),
  Same_direction = sum(main_glaucoma_2W_shared_pathways$NES_direction == "Same"),
  Opposite_direction = sum(main_glaucoma_2W_shared_pathways$NES_direction == "Opposite")
))

# Print all focused pathway names grouped by NES sign; only display precision is rounded.
for (direction in main_shared_direction_levels) {
  displayed <- main_glaucoma_2W_shared_pathways[
    main_glaucoma_2W_shared_pathways$Direction_pattern == direction,
    c("Pathway_name", "OSD_NES", "Terrestrial_NES", "OSD_FDR", "Terrestrial_FDR"), drop = FALSE
  ]
  if (nrow(displayed) > 0L) {
    cat("\n", direction, "\n", sep = "")
    names(displayed)[names(displayed) == "Terrestrial_NES"] <- "Glaucoma_NES"
    names(displayed)[names(displayed) == "Terrestrial_FDR"] <- "Glaucoma_FDR"
    print(displayed, digits = 4, row.names = FALSE)
  }
}

# Save complete-precision name lists for the focused contrast and every terrestrial comparison.
main_shared_names_files <- file.path(main_analysis_dir, c(
  "Main_analysis_OSD758_Glaucoma_2W_shared_pathway_names.csv",
  "Main_analysis_OSD758_both_enriched_pathway_names_all_contrasts.csv",
  "Main_analysis_OSD758_same_direction_pathway_names_all_contrasts.csv",
  "Main_analysis_OSD758_opposite_direction_pathway_names_all_contrasts.csv"
))
utils::write.csv(main_glaucoma_2W_shared_pathways, main_shared_names_files[1], row.names = FALSE)
utils::write.csv(main_both_enriched_pathway_names, main_shared_names_files[2], row.names = FALSE)
utils::write.csv(main_same_direction_pathway_names, main_shared_names_files[3], row.names = FALSE)
utils::write.csv(main_opposite_direction_pathway_names, main_shared_names_files[4], row.names = FALSE)
print(file.exists(main_shared_names_files))
