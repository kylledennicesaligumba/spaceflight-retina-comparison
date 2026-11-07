# Review how the reference mapping issues overlap with the existing flight DE results.
# Keep the original results and every mapping-audit row unchanged.

# Confirm that the mapping audit remains aligned with the complete ranked reference.
stopifnot(
  exists("main_osd_mapping_audit"), exists("main_ranked_tables"),
  dir.exists(main_analysis_dir),
  identical(
    main_osd_mapping_audit$Source_ID,
    as.character(main_ranked_tables[["OSD758__FLT_0G_vs_HGC"]]$gene_id)
  )
)

# Attach the existing statistics to a separate review table without recalculating them.
main_osd_mapping_review <- cbind(
  main_osd_mapping_audit,
  main_ranked_tables[["OSD758__FLT_0G_vs_HGC"]][,
    c("logFC", "PValue", "FDR", "rank_score"), drop = FALSE
  ]
)

# Separate repeated one-match targets from the other one-match rows for this review.
# These groups describe mapping status and do not define the final eligible genes.
main_osd_mapping_review$Review_group <- main_osd_mapping_review$Mapping_status
main_osd_mapping_review$Review_group[
  main_osd_mapping_review$Mapping_status == "One_match" &
    !main_osd_mapping_review$Target_repeated
] <- "One_match_no_repeat"
main_osd_mapping_review$Review_group[
  main_osd_mapping_review$Target_repeated
] <- "One_match_shared_target"
main_osd_review_groups <- c(
  "One_match_no_repeat", "One_match_shared_target", "Multiple_matches", "No_match"
)

# Count all rows and previously significant genes within each disjoint mapping group.
main_osd_mapping_impact <- data.frame(
  Mapping_group = main_osd_review_groups,
  Rows = vapply(main_osd_review_groups, function(g) {
    sum(main_osd_mapping_review$Review_group == g)
  }, numeric(1)),
  DEGs_FDR005 = vapply(main_osd_review_groups, function(g) {
    sum(main_osd_mapping_review$Review_group == g & main_osd_mapping_review$FDR < 0.05)
  }, numeric(1)),
  row.names = NULL
)
stopifnot(
  sum(main_osd_mapping_impact$Rows) == nrow(main_osd_mapping_review),
  sum(main_osd_mapping_impact$DEGs_FDR005) == sum(main_osd_mapping_review$FDR < 0.05)
)
print(main_osd_mapping_impact, row.names = FALSE)

# Review unresolved rows with the lowest existing P values first.
# This display ordering does not change the original gene ranking or DEG calls.
main_osd_unresolved_review <- subset(
  main_osd_mapping_review, Review_group != "One_match_no_repeat"
)
main_osd_unresolved_review <- main_osd_unresolved_review[
  order(main_osd_unresolved_review$PValue), , drop = FALSE
]
print(head(
  main_osd_unresolved_review[,
    c("Source_ID", "Review_group", "Candidate_Entrez", "logFC", "FDR"),
    drop = FALSE
  ], n = 10
), row.names = FALSE)

# Save the impact summary and full unresolved review for documented follow-up.
main_osd_mapping_review_files <- file.path(
  main_analysis_dir,
  c("OSD758_Entrez_mapping_impact_summary.csv", "OSD758_Entrez_unresolved_mapping_review.csv")
)
write.csv(main_osd_mapping_impact, main_osd_mapping_review_files[1], row.names = FALSE)
write.csv(main_osd_unresolved_review, main_osd_mapping_review_files[2], row.names = FALSE)
print(file.exists(main_osd_mapping_review_files))
