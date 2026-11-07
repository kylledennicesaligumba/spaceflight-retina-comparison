# Import the complete saved response tables and review their gene-identifier fields.
# Run after confirming the 13 input paths in the main-analysis setup.

# Confirm that the reviewed input manifest and output folder are available.
stopifnot(
  exists("main_input_files"), dir.exists(main_analysis_dir),
  all(file.exists(main_input_files$File))
)

# Read each full ranking without changing its original column headings.
main_ranked_tables <- lapply(
  main_input_files$File, utils::read.csv,
  check.names = FALSE, stringsAsFactors = FALSE
)
names(main_ranked_tables) <- paste(
  main_input_files$Study, main_input_files$Comparison, sep = "__"
)

# Show the imported row and column counts for every contrast.
main_import_summary <- data.frame(
  Study = main_input_files$Study,
  Comparison = main_input_files$Comparison,
  Rows = vapply(main_ranked_tables, nrow, integer(1)),
  Columns = vapply(main_ranked_tables, ncol, integer(1)),
  row.names = NULL
)
print(main_import_summary, row.names = FALSE)

# Save all column names alongside the import counts for later mapping decisions.
main_column_inventory <- main_import_summary
main_column_inventory$Column_names <- vapply(
  main_ranked_tables, function(x) paste(names(x), collapse = " | "), character(1)
)
main_inventory_file <- file.path(main_analysis_dir, "Main_analysis_table_inventory.csv")
write.csv(main_column_inventory, main_inventory_file, row.names = FALSE)

# Preview columns and three identifier rows from one contrast per study.
# The review identifies source-label formats; no gene mapping is performed here.
for (i in which(!duplicated(main_input_files$Study))) {
  cat("\n", names(main_ranked_tables)[i], "\n", sep = "")
  print(names(main_ranked_tables[[i]]))
  main_id_preview_columns <- grep(
    "gene|ensembl|entrez|symbol|accid|^id$|^x$",
    names(main_ranked_tables[[i]]), ignore.case = TRUE, value = TRUE
  )
  print(head(main_ranked_tables[[i]][, main_id_preview_columns, drop = FALSE], n = 3))
}

# Confirm that the import inventory was saved in the main-analysis folder.
print(file.exists(main_inventory_file))
