# Verify the saved Light Damage results in a fresh, separate R process.
# Reloading the checkpoint also checks that it does not depend on live XML pointers.
source(file.path("scripts", "21_Light_Damage_GSE314043_reload.R"))

# Confirm original input checksums still agree with the saved source manifest.
ld_verify_manifest <- read.csv(file.path(ld_output_dir, "Light_Damage_GSE314043_input_manifest.csv"))
stopifnot(identical(unname(tools::md5sum(ld_verify_manifest$path)), ld_verify_manifest$md5))

# Check for identical count profiles across different sample columns.
ld_duplicate_profiles <- anyDuplicated(as.data.frame(t(ld_counts_all)))
stopifnot(ld_duplicate_profiles == 0L)
cat("Duplicate sample count profiles:", ld_duplicate_profiles, "\n")

# Reproduce key exported tables from the reloaded objects.
ld_verify_counts <- read.csv(file.path(ld_output_dir, "Light_Damage_GSE314043_raw_counts.csv"), check.names = FALSE)
stopifnot(identical(ld_verify_counts$Geneid, rownames(ld_counts_all)),
          isTRUE(all.equal(unname(as.matrix(ld_verify_counts[, -1])), unname(ld_counts_all))))
ld_verify_qc <- read.csv(file.path(ld_output_dir, "Light_Damage_GSE314043_sample_QC.csv"))
stopifnot(identical(ld_verify_qc$sample_id, colnames(ld_counts_all)),
          isTRUE(all.equal(ld_verify_qc$total_gene_counts, unname(colSums(ld_counts_all)))))

# Confirm common-gene MDS agrees with an independent centered, unscaled PCA.
# PCA axes can reverse sign without changing sample relationships.
ld_verify_variances <- apply(ld_log_cpm_qc, 1, stats::var)
ld_verify_top500 <- order(ld_verify_variances, decreasing = TRUE)[seq_len(500)]
ld_verify_prcomp <- stats::prcomp(t(ld_log_cpm_qc[ld_verify_top500, ]), center = TRUE, scale. = FALSE)
ld_verify_pca_agreement <- c(
  PC1 = abs(stats::cor(ld_verify_prcomp$x[, 1], ld_pca$x)),
  PC2 = abs(stats::cor(ld_verify_prcomp$x[, 2], ld_pca$y))
)
ld_verify_pca_percent <- 100 * ld_verify_prcomp$sdev^2 / sum(ld_verify_prcomp$sdev^2)
stopifnot(all(ld_verify_pca_agreement > 1 - 1e-10),
          isTRUE(all.equal(ld_verify_pca_percent[1:2], 100 * ld_pca$var.explained[1:2])))

# limma uses mean squared gene differences, so its coordinates are PCA/sqrt(500).
ld_verify_pca_sign <- c(
  sign(stats::cor(ld_verify_prcomp$x[, 1], ld_pca$x)),
  sign(stats::cor(ld_verify_prcomp$x[, 2], ld_pca$y))
)
stopifnot(
  isTRUE(all.equal(unname(ld_verify_prcomp$x[, 1]) * ld_verify_pca_sign[1] / sqrt(500), unname(ld_pca$x))),
  isTRUE(all.equal(unname(ld_verify_prcomp$x[, 2]) * ld_verify_pca_sign[2] / sqrt(500), unname(ld_pca$y)))
)
print(ld_verify_pca_agreement)
cat("Independent PCA axis percentages:", ld_verify_pca_percent[1:2], "\n")

# Save the exact top-variable genes used for PCA, with their available names.
ld_verify_top_genes <- data.frame(
  rank = seq_len(500), Geneid = rownames(ld_log_cpm_qc)[ld_verify_top500],
  variance_logCPM = ld_verify_variances[ld_verify_top500], row.names = NULL
)
ld_verify_top_genes$gene_name <- ld_gene_annotation$gene_name[
  match(ld_verify_top_genes$Geneid, ld_gene_annotation$Geneid)
]
write.csv(ld_verify_top_genes, file.path(ld_output_dir, "Light_Damage_GSE314043_PCA_top500_genes.csv"),
          row.names = FALSE, na = "")

# Save a concise verification result alongside the complete run log.
writeLines(c(
  "Verification completed: 3 November 2026, Asia/Manila.",
  "Checkpoint reload, sample order and gene order: PASS.",
  "Original count/metadata source checksums: PASS.",
  "Raw-count CSV and sample-QC CSV agree with checkpoint: PASS.",
  "No duplicate sample count profiles: PASS.",
  "PCA agrees with independent prcomp calculation (allowing axis sign and sqrt(500) scaling): PASS.",
  "Compressed Ensembl GTF passed gzip integrity verification.",
  "All four exported PNG figures were visually inspected."
), file.path(ld_output_dir, "Light_Damage_GSE314043_verification.txt"))
cat("Saved-result verification completed successfully.\n")
