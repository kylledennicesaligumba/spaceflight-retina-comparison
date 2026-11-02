# Restore the reviewed stage after running the 2 November save script.
precheck_project_dir <- normalizePath(getwd(), mustWork = TRUE)
uveitis_checkpoint_file <- file.path(
  precheck_project_dir, "results", "autoimmune_uveitis",
  "Autoimmune_Uveitis_GSE284505_precheck_checkpoint_2026-11-02.RData"
)
stopifnot(file.exists(uveitis_checkpoint_file))
load(uveitis_checkpoint_file, envir = .GlobalEnv)
setwd(uveitis_output_dir)

# Reopen XML reader objects from the preserved metadata file.
uveitis_metadata_xml <- xml2::read_xml(uveitis_metadata_file)
uveitis_sample_nodes <- xml2::xml_find_all(
  uveitis_metadata_xml, "/geo:MINiML/geo:Sample", ns = uveitis_ns
)
uveitis_characteristic_nodes <- xml2::xml_find_all(
  uveitis_sample_nodes, ".//geo:Characteristics", ns = uveitis_ns
)
uveitis_protocol_nodes <- xml2::xml_find_all(
  uveitis_sample_nodes,
  ".//geo:Growth-Protocol | .//geo:Treatment-Protocol | .//geo:Extract-Protocol",
  ns = uveitis_ns
)

# Check restored counts, annotations and the filtered/TMM exploration matrix.
stopifnot(
  identical(dim(uveitis_counts_selected), c(56859L, 6L)),
  identical(dim(uveitis_metadata_selected), c(6L, 6L)),
  identical(dim(uveitis_gene_annotation), c(56859L, 9L)),
  identical(dim(uveitis_dge_qc), c(18187L, 6L)),
  identical(dim(uveitis_log_cpm_qc), c(18187L, 6L)),
  identical(colnames(uveitis_counts_selected), uveitis_metadata_selected$sample_id),
  identical(rownames(uveitis_counts_selected), uveitis_raw$AccID),
  identical(uveitis_gene_annotation$AccID, rownames(uveitis_counts_selected)),
  identical(uveitis_dge_qc$samples$sample_id, colnames(uveitis_dge_qc)),
  identical(uveitis_dge_qc$genes$AccID, rownames(uveitis_dge_qc)),
  identical(colnames(uveitis_log_cpm_qc), colnames(uveitis_dge_qc))
)

# Reopen the formatting helper for the saved projection plots when needed.
source(file.path(precheck_project_dir, "scripts", "14_Uveitis_projection_style.R"))
cat("Checkpoint restored through annotation, filtering, TMM, MDS and PCA.\n")
print(getwd())
print(dim(uveitis_counts_selected))
print(dim(uveitis_gene_annotation))
print(dim(uveitis_log_cpm_qc))
cat("All six samples remain retained; see docs/autoimmune-uveitis-precheck.md for the completed record.\n")
