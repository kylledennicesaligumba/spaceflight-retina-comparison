# Portable repository copy; run from the repository root.
# Original working files are preserved. GitHub sync did not execute this code.
# Reload the completed LHON precheck after running 08_LHON_GSE304707_save_final.R.
# Make the final local checkpoint and its saved XML path available first.
# This restores existing results without refitting; see docs/running-the-lhon-precheck.md.
# The partial 29 October checkpoint and its old reload script remain archived.

# Load the packages used by the saved count and plotting objects.
stopifnot(all(vapply(
  c("edgeR", "ggplot2", "ggrepel", "xml2"),
  requireNamespace, logical(1), quietly = TRUE
)))

# Target the completed-stage checkpoint explicitly.
lhon_reload_file <- file.path("results", "lhon", "LHON_GSE304707_precheck_checkpoint_2026-10-30.RData")
if (!file.exists(lhon_reload_file)) {
  stop("The final checkpoint is not saved yet. Run 08_LHON_GSE304707_save_final.R in the completed R session first. The older 29 October file is only a partial backup.")
}
load(lhon_reload_file)

# Reopen XML document/node objects, which were excluded from serialization.
# This restores access to existing metadata rather than rerunning expression QC.
lhon_metadata_xml <- xml2::read_xml(lhon_metadata_file)
lhon_sample_nodes <- xml2::xml_find_all(
  lhon_metadata_xml,
  "/geo:MINiML/geo:Sample",
  ns = lhon_ns
)

# Check the saved dimensions and the calculated projection modes.
stopifnot(
  identical(dim(lhon_counts_all), c(49714L, 12L)),
  identical(dim(lhon_dge_qc), c(17860L, 12L)),
  identical(dim(lhon_log_cpm_qc), c(17860L, 12L)),
  identical(lhon_mds$gene.selection, "pairwise"),
  identical(lhon_pca$gene.selection, "common")
)

# Print the restored matrix sizes and alignment checks explicitly.
cat("Original count matrix:\n")
print(dim(lhon_counts_all))
cat("Filtered QC object:\n")
print(dim(lhon_dge_qc))
print(c(
  sample_order = identical(lhon_dge_qc$samples$sample_id, colnames(lhon_dge_qc)),
  gene_order = identical(lhon_dge_qc$genes$gene_name, rownames(lhon_dge_qc))
))

# These commands redisplay existing plot objects when needed:
# lhon_draw_boxplot()
# lhon_draw_density()
# lhon_draw_mds()
# lhon_draw_pca()
# lhon_tmm_summary
# The LHON processed-data precheck is closed with all 12 libraries retained.
