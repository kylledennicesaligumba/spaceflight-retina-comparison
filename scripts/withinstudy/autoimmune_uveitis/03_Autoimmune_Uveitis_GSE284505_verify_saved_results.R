# Repository paths: begin in the repository root, or set RETINA_PROJECT_DIR.
# Existing model calculations are preserved; downloaded checkpoints stay local.
if (!exists("portfolio_repo_dir", inherits = TRUE)) {
  portfolio_repo_dir <- normalizePath(Sys.getenv("RETINA_PROJECT_DIR", unset = getwd()), mustWork = TRUE)
}

# Verify saved within-study objects against the precheck and exported results.
# This checks saved results without fitting another model or changing the analysis.
verify_folder <- file.path(portfolio_repo_dir, "results", "withinstudy", "autoimmune_uveitis")
verify_stem <- "Autoimmune_Uveitis_GSE284505"
verify_prefix <- "uveitis"
verify_env <- new.env()
verify_names <- load(file.path(verify_folder, paste0(verify_stem, "_withinstudy_final_checkpoint_2026-11-05.RData")), envir = verify_env)
verify_get <- function(name) get(paste0(verify_prefix, "_", name), envir = verify_env)
verify_precheck <- new.env()
load(verify_get("precheck_file"), envir = verify_precheck)

# Confirm that no counts, saved normalization factors or sample assignments changed.
verify_dge <- verify_get("within_dge")
verify_original <- get(paste0(verify_prefix, "_dge_qc"), envir = verify_precheck)
verify_results <- verify_get("results")
verify_ranked <- verify_get("ranked")
verify_dots <- verify_get("dot_data")
stopifnot(
  identical(verify_dge$counts, verify_original$counts),
  identical(verify_dge$samples$norm.factors, verify_original$samples$norm.factors),
  identical(colnames(verify_dge), verify_get("within_metadata")$sample_id),
  setequal(verify_results$gene_id, rownames(verify_dge)),
  !any(vapply(mget(verify_names, envir = verify_env), is.environment, logical(1)))
)

# Check the original QL statistics, BH calculation, DEG membership and full ranking.
stopifnot(
  isTRUE(all.equal(verify_results$F, verify_get("test")$table[verify_results$gene_id, "F"])),
  isTRUE(all.equal(verify_results$FDR, p.adjust(verify_results$PValue, "BH"))),
  identical(verify_get("degs")$gene_id, verify_results$gene_id[verify_results$FDR < 0.05]),
  isTRUE(all.equal(verify_ranked$rank_score, with(verify_ranked, sign(logFC) * sqrt(pmax(F, 0))))),
  all(diff(verify_ranked$rank_score) <= 0),
  nrow(verify_ranked) == nrow(verify_results)
)

# Check heatmap scaling, significant-gene selections and one dot per gene-library pair.
stopifnot(
  identical(dim(verify_get("heat_z")), c(30L, ncol(verify_dge))),
  isTRUE(all.equal(rowMeans(verify_get("heat_z")), rep(0, 30), check.attributes = FALSE, tolerance = 1e-10)),
  isTRUE(all.equal(apply(verify_get("heat_z"), 1, sd), rep(1, 30), check.attributes = FALSE)),
  all(verify_get("heat_genes")$FDR < 0.05),
  all(verify_get("dot_genes")$FDR < 0.05),
  identical(verify_get("dot_genes")$gene_id, verify_get("md_labels")$gene_id),
  nrow(verify_dots) == 10L * ncol(verify_dge),
  anyDuplicated(verify_dots[c("gene_id", "sample_id")]) == 0L
)
verify_cpm <- edgeR::cpm(verify_dge, log = FALSE)
verify_dot_values <- as.numeric(verify_cpm[cbind(match(verify_dots$gene_id, rownames(verify_cpm)), match(verify_dots$sample_id, colnames(verify_cpm)))])
stopifnot(isTRUE(all.equal(verify_dots$CPM, verify_dot_values)))

# Check the source checkpoint checksum, exported statistics and four saved PNG files.
verify_manifest <- read.csv(file.path(verify_folder, paste0(verify_stem, "_input_manifest.csv")))
verify_csv <- read.csv(file.path(verify_folder, paste0(verify_stem, "_EAU_vs_CON_all_genes.csv")))
verify_plots <- file.path(verify_folder, paste0(verify_stem, c("_QL_dispersion.png", "_MD_plot.png", "_DEG_expression_heatmap_top30.png", "_DEG_expression_dotplots_top10.png")))
stopifnot(
  identical(verify_manifest$md5, unname(tools::md5sum(verify_manifest$path))),
  identical(verify_csv$gene_id, verify_results$gene_id),
  isTRUE(all.equal(verify_csv$logFC, verify_results$logFC)),
  isTRUE(all.equal(verify_csv$FDR, verify_results$FDR)),
  all(file.exists(verify_plots)), all(file.info(verify_plots)$size > 0)
)

# Exercise the actual reload helper in an independent environment.
source(file.path(portfolio_repo_dir, "scripts", "withinstudy", "autoimmune_uveitis", paste0("02_", verify_stem, "_reload_final.R")), local = new.env())
writeLines(c(
  "Fresh-process checkpoint restoration and reload helper: PASS.",
  "Counts and TMM factors identical to precheck: PASS.",
  "Sample/gene alignment, QL statistics, BH FDR and DEG membership: PASS.",
  "Signed full-gene ranking: PASS.",
  "Heatmap row z-scores and DEG selection: PASS.",
  "Dot-plot keys, gene selection and normalized CPM values: PASS.",
  "Result CSV agrees with saved objects: PASS.",
  "Input checkpoint checksum unchanged: PASS.",
  "All four PNG files present and nonempty: PASS.",
  "No saved precheck environment or unrelated study objects: PASS."
), file.path(verify_folder, paste0(verify_stem, "_verification_2026-11-05.txt")))
cat(verify_stem, ": all saved-result checks passed.\n")
