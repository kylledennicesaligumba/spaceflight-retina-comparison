# Reload the saved Light Damage checkpoint without touching AMD objects.
ld_output_dir <- file.path(normalizePath(getwd(), mustWork = TRUE), "results", "light_damage")
load(file.path(ld_output_dir, "Light_Damage_GSE314043_precheck_checkpoint_2026-11-03.RData"))

# Confirm the restored sample and gene orders before using the QC objects.
stopifnot(
  identical(colnames(ld_counts_selected), ld_metadata_selected$sample_id),
  identical(rownames(ld_counts_selected), ld_gene_annotation$Geneid),
  identical(colnames(ld_log_cpm_qc), ld_metadata_selected$sample_id)
)

# Redraw a saved QC plot when needed; uncomment the desired command.
# ld_draw_boxplot()
# ld_draw_density()
# ld_draw_mds()
# ld_draw_pca()
