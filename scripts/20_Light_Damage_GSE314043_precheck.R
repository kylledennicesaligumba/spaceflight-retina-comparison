# Light Damage GSE314043: reproducible count-matrix prequality check.
# Run in a separate R process with Rscript --vanilla to keep AMD work separate.
# All created objects use ld_ names; downloaded source files stay unchanged.
# Completed on 3 November 2026, Asia/Manila; no differential testing is performed.

# Set input and output locations without changing R's working directory.
ld_dir <- Sys.getenv(
  "LIGHT_DAMAGE_DATA_DIR",
  unset = file.path(normalizePath(getwd(), mustWork = TRUE), "data", "raw", "Light_Damage_GSE314043")
)
ld_output_dir <- file.path(normalizePath(getwd(), mustWork = TRUE), "results", "light_damage")
ld_metadata_dir <- file.path(ld_output_dir, "metadata")
ld_annotation_dir <- file.path(ld_output_dir, "annotation")
dir.create(ld_metadata_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(ld_annotation_dir, recursive = TRUE, showWarnings = FALSE)
ld_count_file <- file.path(ld_dir, "GSE314043_genecount.txt")
ld_metadata_archive <- file.path(ld_dir, "GSE314043_family")
stopifnot(file.exists(ld_count_file), file.exists(ld_metadata_archive))

# Check required packages without installing or changing the user's R setup.
ld_packages <- c("xml2", "edgeR", "limma", "ggplot2", "ggrepel")
stopifnot(all(vapply(ld_packages, requireNamespace, logical(1), quietly = TRUE)))

# Import the already-uncompressed tab-separated count file, preserving labels.
ld_raw <- read.delim(ld_count_file, check.names = FALSE, stringsAsFactors = FALSE)
print(dim(ld_raw))
print(colnames(ld_raw))
print(head(ld_raw, 5))
stopifnot(colnames(ld_raw)[1] == "Geneid")

# Check identifiers and count storage before constructing the count matrix.
stopifnot(
  !anyNA(ld_raw),
  all(nzchar(trimws(ld_raw$Geneid))),
  anyDuplicated(ld_raw$Geneid) == 0L,
  anyDuplicated(colnames(ld_raw)[-1]) == 0L,
  all(vapply(ld_raw[, -1, drop = FALSE], is.numeric, logical(1)))
)
ld_counts_all <- as.matrix(ld_raw[, -1, drop = FALSE])
rownames(ld_counts_all) <- ld_raw$Geneid
stopifnot(
  all(is.finite(ld_counts_all)), all(ld_counts_all >= 0),
  all(ld_counts_all == round(ld_counts_all)), all(colSums(ld_counts_all) > 0)
)
ld_integrity <- data.frame(
  check = c("gene_rows", "sample_columns", "missing_entries", "blank_gene_ids",
            "duplicate_gene_ids", "finite_counts", "nonnegative_counts",
            "whole_counts", "minimum_count", "maximum_count", "all_zero_genes"),
  value = c(nrow(ld_counts_all), ncol(ld_counts_all), sum(is.na(ld_raw)),
            sum(!nzchar(trimws(ld_raw$Geneid))), sum(duplicated(ld_raw$Geneid)),
            all(is.finite(ld_counts_all)), all(ld_counts_all >= 0),
            all(ld_counts_all == round(ld_counts_all)), range(ld_counts_all),
            sum(rowSums(ld_counts_all) == 0))
)
print(ld_integrity)

# List the archive and extract only the known XML member to the output folder.
ld_archive_members <- utils::untar(ld_metadata_archive, list = TRUE)
print(ld_archive_members)
stopifnot("GSE314043_family.xml" %in% ld_archive_members)
utils::untar(ld_metadata_archive, files = "GSE314043_family.xml", exdir = ld_metadata_dir)
ld_metadata_file <- file.path(ld_metadata_dir, "GSE314043_family.xml")
ld_metadata_xml <- xml2::read_xml(ld_metadata_file)
ld_ns <- c(geo = "http://www.ncbi.nlm.nih.gov/geo/info/MINiML")
ld_sample_nodes <- xml2::xml_find_all(ld_metadata_xml, "/geo:MINiML/geo:Sample", ld_ns)

# Read the recorded characteristics without inventing unreported design fields.
ld_get_characteristic <- function(tag) {
  xml2::xml_text(xml2::xml_find_first(
    ld_sample_nodes, paste0(".//geo:Characteristics[@tag='", tag, "']"), ld_ns
  ), trim = TRUE)
}
ld_metadata_all <- data.frame(
  geo_accession = xml2::xml_attr(ld_sample_nodes, "iid"),
  title = xml2::xml_text(xml2::xml_find_first(ld_sample_nodes, "geo:Title", ld_ns), trim = TRUE),
  organism = xml2::xml_text(xml2::xml_find_first(ld_sample_nodes, ".//geo:Organism", ld_ns), trim = TRUE),
  tissue = ld_get_characteristic("tissue"),
  cell_type = ld_get_characteristic("cell type"),
  genotype_recorded = ld_get_characteristic("genotype"),
  treatment_recorded = ld_get_characteristic("treatment"),
  batch_recorded = ld_get_characteristic("batch"),
  stringsAsFactors = FALSE
)

# The deposited genotype field contains library names, not biological genotypes.
# Use these exact names as the key and validate them against count columns.
ld_metadata_all$sample_id <- ld_metadata_all$genotype_recorded
stopifnot(
  !anyNA(ld_metadata_all$sample_id), anyDuplicated(ld_metadata_all$sample_id) == 0L,
  setequal(ld_metadata_all$sample_id, colnames(ld_counts_all))
)
ld_match <- match(colnames(ld_counts_all), ld_metadata_all$sample_id)
ld_metadata_matched <- ld_metadata_all[ld_match, , drop = FALSE]
rownames(ld_metadata_matched) <- NULL
stopifnot(identical(ld_metadata_matched$sample_id, colnames(ld_counts_all)))

# Retain all deposited libraries for LD-vehicle versus CTRL-vehicle.
# Replicate numbers in different groups do not establish paired animals.
ld_metadata_selected <- ld_metadata_matched
ld_metadata_selected$condition <- ifelse(
  grepl("^CTRL-vehicle-[1-5]$", ld_metadata_selected$sample_id), "CTRL",
  ifelse(grepl("^LD-vehicle-[1-5]$", ld_metadata_selected$sample_id), "LD", NA_character_)
)
ld_metadata_selected$replicate_label <- sub(".*-", "", ld_metadata_selected$sample_id)
ld_counts_selected <- ld_counts_all[, ld_metadata_selected$sample_id, drop = FALSE]
ld_qc_group <- factor(ld_metadata_selected$condition, levels = c("CTRL", "LD"))
stopifnot(
  !anyNA(ld_qc_group), all(table(ld_qc_group) == 5),
  all(ld_metadata_selected$tissue == "retina"),
  all(ld_metadata_selected$organism == "Mus musculus"),
  identical(colnames(ld_counts_selected), ld_metadata_selected$sample_id)
)
print(ld_metadata_selected)
print(table(ld_qc_group))

# Save the source protocol text so design and processing claims remain auditable.
ld_protocol_text <- unique(xml2::xml_text(xml2::xml_find_all(
  ld_metadata_xml,
  "/geo:MINiML/geo:Series/geo:Summary | /geo:MINiML/geo:Series/geo:Overall-Design | /geo:MINiML/geo:Sample//geo:Growth-Protocol | /geo:MINiML/geo:Sample//geo:Treatment-Protocol | /geo:MINiML/geo:Sample//geo:Extract-Protocol | /geo:MINiML/geo:Sample//geo:Data-Processing",
  ld_ns
), trim = TRUE))
writeLines(ld_protocol_text, file.path(ld_metadata_dir, "Light_Damage_GSE314043_deposited_protocols.txt"))

# Summarise assigned gene counts and positive-count genes for every library.
ld_total_counts <- colSums(ld_counts_selected)
ld_detected_genes <- colSums(ld_counts_selected > 0)
ld_sample_qc <- data.frame(
  sample_id = colnames(ld_counts_selected), geo_accession = ld_metadata_selected$geo_accession,
  condition = as.character(ld_qc_group), total_gene_counts = ld_total_counts,
  detected_genes = ld_detected_genes, zero_genes = colSums(ld_counts_selected == 0),
  row.names = NULL
)
print(ld_sample_qc)

# Calculate initial all-gene log-CPM before expression filtering and TMM.
# CPM accounts for library totals; prior.count = 2 keeps zero counts finite.
ld_log_cpm <- edgeR::cpm(ld_counts_selected, log = TRUE, prior.count = 2)
ld_palette <- c(CTRL = "#0072B2", LD = "#D55E00")
ld_plot_colors <- unname(ld_palette[as.character(ld_qc_group)])
ld_draw_boxplot <- function() {
  par(mar = c(9, 4.5, 4, 1))
  boxplot(ld_log_cpm, las = 2, col = ld_plot_colors, outline = FALSE,
          main = "Light Damage: initial expression distributions",
          ylab = "Log2 counts per million (CPM)", cex.axis = 0.8)
}
ld_boxplot_file <- file.path(ld_output_dir, "Light_Damage_GSE314043_logCPM_boxplot.png")
png(ld_boxplot_file, width = 2400, height = 1800, res = 300)
ld_draw_boxplot()
dev.off()

# Estimate and draw one density curve per library using all supplied gene rows.
ld_density_curves <- lapply(seq_len(ncol(ld_log_cpm)), function(i) density(ld_log_cpm[, i]))
names(ld_density_curves) <- colnames(ld_log_cpm)
ld_density_colors <- hcl.colors(ncol(ld_log_cpm), "Dark 3")
ld_draw_density <- function() {
  par(mar = c(4.5, 4.5, 4, 1))
  plot(ld_density_curves[[1]],
       xlim = range(vapply(ld_density_curves, function(d) range(d$x), numeric(2))),
       ylim = c(0, max(vapply(ld_density_curves, function(d) max(d$y), numeric(1)))),
       col = ld_density_colors[1], lwd = 2,
       main = "Light Damage: initial expression density",
       xlab = "Log2 counts per million (CPM)", ylab = "Density")
  for (i in 2:length(ld_density_curves)) {
    lines(ld_density_curves[[i]], col = ld_density_colors[i], lwd = 2)
  }
  legend("topright", names(ld_density_curves), col = ld_density_colors,
         lty = 1, lwd = 2, bty = "n", cex = 0.7)
}
ld_density_file <- file.path(ld_output_dir, "Light_Damage_GSE314043_logCPM_density.png")
png(ld_density_file, width = 2400, height = 1800, res = 300)
ld_draw_density()
dev.off()

# Validate mouse Ensembl gene identifiers while preserving the original IDs.
ld_gene_labels <- rownames(ld_counts_selected)
ld_ensembl_prefix <- grepl("^ENSMUSG[0-9]{11}$", ld_gene_labels)
print(table(Mouse_Ensembl_gene_ID = ld_ensembl_prefix))
stopifnot(all(ld_ensembl_prefix))

# Use the fixed Ensembl release-116 GRCm39 GTF after REST requests failed.
# The study reports GRCm39 but does not report its original annotation release.
# This adds gene names/biotypes without renaming or discarding count rows.
ld_ensembl_release <- 116L
ld_gtf_url <- "https://ftp.ensembl.org/pub/release-116/gtf/mus_musculus/Mus_musculus.GRCm39.116.gtf.gz"
ld_gtf_file <- file.path(ld_annotation_dir, "Mus_musculus.GRCm39.116.gtf.gz")
if (!file.exists(ld_gtf_file)) {
  ld_download_timeout <- getOption("timeout")
  options(timeout = 300)
  download.file(ld_gtf_url, ld_gtf_file, mode = "wb", method = "libcurl")
  options(timeout = ld_download_timeout)
}

# Stream the compressed GTF in small chunks, keeping only gene-feature records.
# This avoids loading all exon/transcript records into the R workspace.
ld_gtf_connection <- gzfile(ld_gtf_file, open = "rt")
ld_gtf_gene_lines <- character()
repeat {
  ld_gtf_chunk <- readLines(ld_gtf_connection, n = 100000L, warn = FALSE)
  if (!length(ld_gtf_chunk)) break
  ld_gtf_gene_lines <- c(ld_gtf_gene_lines, ld_gtf_chunk[grepl("\tgene\t", ld_gtf_chunk, fixed = TRUE)])
}
close(ld_gtf_connection)
ld_gtf_genes <- read.delim(text = paste(ld_gtf_gene_lines, collapse = "\n"),
                          header = FALSE, sep = "\t", quote = "", comment.char = "",
                          stringsAsFactors = FALSE)
stopifnot(ncol(ld_gtf_genes) == 9L)
colnames(ld_gtf_genes) <- c("chromosome", "source", "feature", "start", "end",
                           "score", "strand", "frame", "attributes")

# Extract named GTF attributes, preserving absent attributes as NA.
ld_gtf_attribute <- function(tag) {
  ld_pattern <- paste0('(?:^|; )', tag, ' "([^"]*)"')
  ld_matches <- regexec(ld_pattern, ld_gtf_genes$attributes, perl = TRUE)
  vapply(regmatches(ld_gtf_genes$attributes, ld_matches), function(x) {
    if (length(x) >= 2L) x[2] else NA_character_
  }, character(1))
}
ld_gtf_genes$ensembl_id <- ld_gtf_attribute("gene_id")
ld_gtf_genes$gene_name <- ld_gtf_attribute("gene_name")
ld_gtf_genes$biotype <- ld_gtf_attribute("gene_biotype")
ld_gtf_genes$gene_version <- ld_gtf_attribute("gene_version")
stopifnot(!anyNA(ld_gtf_genes$ensembl_id), anyDuplicated(ld_gtf_genes$ensembl_id) == 0L)

# Match the annotation to every original gene row, preserving unreturned IDs.
ld_annotation_match <- match(ld_gene_labels, ld_gtf_genes$ensembl_id)
ld_gene_annotation <- data.frame(
  Geneid = ld_gene_labels,
  ld_gtf_genes[ld_annotation_match,
               c("ensembl_id", "gene_name", "biotype", "gene_version", "chromosome", "start", "end", "strand"),
               drop = FALSE],
  annotation_species = "Mus musculus", annotation_assembly = "GRCm39",
  annotation_release = ld_ensembl_release,
  annotation_source = "Ensembl release-116 GRCm39 GTF",
  row.names = NULL, stringsAsFactors = FALSE
)
ld_gene_annotation$lookup_status <- ifelse(
  is.na(ld_gene_annotation$ensembl_id), "not_found",
  ifelse(is.na(ld_gene_annotation$gene_name) | !nzchar(ld_gene_annotation$gene_name),
         "found_without_name", "found_with_name")
)
stopifnot(identical(ld_gene_annotation$Geneid, rownames(ld_counts_selected)))
print(table(ld_gene_annotation$lookup_status))
print(head(ld_gene_annotation[, c("Geneid", "gene_name", "biotype", "annotation_assembly")], 10))
write.csv(ld_gene_annotation[ld_gene_annotation$lookup_status != "found_with_name", ],
          file.path(ld_annotation_dir, "unmapped_or_unnamed_gene_ids.csv"), row.names = FALSE, na = "")
write.csv(data.frame(source_url = ld_gtf_url, release = ld_ensembl_release,
                      assembly = "GRCm39", md5 = unname(tools::md5sum(ld_gtf_file))),
          file.path(ld_annotation_dir, "annotation_provenance.csv"), row.names = FALSE)

# Assemble raw counts, aligned metadata and annotations into an edgeR object.
ld_dge <- edgeR::DGEList(counts = ld_counts_selected, samples = ld_metadata_selected,
                         group = ld_qc_group, genes = ld_gene_annotation)
stopifnot(identical(ld_dge$samples$sample_id, colnames(ld_dge)),
          identical(ld_dge$genes$Geneid, rownames(ld_dge)))

# Apply the same expression-filter settings used for RP, LHON and uveitis.
# This creates a QC copy; every original count row remains in ld_counts_all.
ld_keep_genes <- edgeR::filterByExpr(ld_dge, group = ld_qc_group,
                                   min.count = 10, min.total.count = 15)
print(table(Retain_for_QC = ld_keep_genes))
ld_dge_qc <- ld_dge[ld_keep_genes, , keep.lib.sizes = FALSE]

# Estimate TMM scaling factors and calculate filtered, normalized log-CPM.
ld_dge_qc <- edgeR::normLibSizes(ld_dge_qc, method = "TMM")
stopifnot(all(is.finite(ld_dge_qc$samples$norm.factors)),
          all(ld_dge_qc$samples$norm.factors > 0))
ld_log_cpm_qc <- edgeR::cpm(ld_dge_qc, normalized.lib.sizes = TRUE,
                           log = TRUE, prior.count = 2)
stopifnot(all(is.finite(ld_log_cpm_qc)),
          identical(colnames(ld_log_cpm_qc), ld_metadata_selected$sample_id),
          identical(rownames(ld_log_cpm_qc), ld_dge_qc$genes$Geneid))
print(ld_dge_qc$samples[, c("sample_id", "group", "lib.size", "norm.factors")])
ld_sample_qc$filtered_library_size <- ld_dge_qc$samples$lib.size
ld_sample_qc$tmm_factor <- ld_dge_qc$samples$norm.factors
ld_sample_qc$effective_library_size <- with(ld_sample_qc, filtered_library_size * tmm_factor)

# Calculate pairwise top-500 MDS and common top-500 PCA, as in the uveitis script.
# Common-gene MDS equals PCA on centered, unscaled top-variable gene log-CPM.
ld_mds <- limma::plotMDS(ld_log_cpm_qc, top = 500, gene.selection = "pairwise", plot = FALSE)
ld_pca <- limma::plotMDS(ld_log_cpm_qc, top = 500, gene.selection = "common", plot = FALSE)

# Use the familiar labelled points, condition colours and outside legend.
ld_projection_plot <- function(projection, title, axis_names) {
  ld_plot_data <- data.frame(axis1 = projection$x, axis2 = projection$y,
                            group = ld_qc_group, label = ld_metadata_selected$sample_id)
  ld_percentages <- 100 * projection$var.explained[projection$dim.plot]
  ggplot2::ggplot(ld_plot_data, ggplot2::aes(axis1, axis2, colour = group)) +
    ggplot2::geom_point(size = 3) +
    ggrepel::geom_text_repel(ggplot2::aes(label = label), size = 3.4, seed = 1,
                            box.padding = 0.6, point.padding = 0.4,
                            min.segment.length = 0, max.overlaps = Inf,
                            max.time = 2, show.legend = FALSE) +
    ggplot2::scale_colour_manual(values = ld_palette, labels = c("Control", "Light Damage")) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.25)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = 0.25)) +
    ggplot2::labs(title = title,
                  x = sprintf("%s (%.1f%%)", axis_names[1], ld_percentages[1]),
                  y = sprintf("%s (%.1f%%)", axis_names[2], ld_percentages[2]),
                  colour = "Condition") +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(legend.position = "right", legend.title = ggplot2::element_text(face = "bold"),
                   plot.title = ggplot2::element_text(hjust = 0.5))
}
ld_mds_plot <- ld_projection_plot(ld_mds, "Light Damage: MDS after TMM",
                                  c("Leading logFC dim 1", "Leading logFC dim 2"))
ld_pca_plot <- ld_projection_plot(ld_pca, "Light Damage: PCA after TMM", c("PC1", "PC2"))
ld_draw_mds <- function() print(ld_mds_plot)
ld_draw_pca <- function() print(ld_pca_plot)
ld_mds_file <- file.path(ld_output_dir, "Light_Damage_GSE314043_TMM_MDS.png")
ld_pca_file <- file.path(ld_output_dir, "Light_Damage_GSE314043_TMM_PCA_top500.png")
ggplot2::ggsave(ld_mds_file, ld_mds_plot, width = 9, height = 6, dpi = 300)
ggplot2::ggsave(ld_pca_file, ld_pca_plot, width = 9, height = 6, dpi = 300)

# Save projection coordinates, gene filtering decisions and sample summaries.
ld_coordinates <- data.frame(sample_id = ld_metadata_selected$sample_id,
                              condition = as.character(ld_qc_group),
                              MDS1 = ld_mds$x, MDS2 = ld_mds$y,
                              PC1 = ld_pca$x, PC2 = ld_pca$y)
ld_variance <- data.frame(axis = seq_along(ld_pca$var.explained),
                           MDS_percent = 100 * ld_mds$var.explained,
                           PCA_percent = 100 * ld_pca$var.explained)
ld_filter_table <- data.frame(Geneid = rownames(ld_counts_selected), retain_for_qc = ld_keep_genes)
ld_export_csv <- function(object, name) {
  write.csv(object, file.path(ld_output_dir, paste0("Light_Damage_GSE314043_", name, ".csv")),
            row.names = FALSE, na = "")
}
ld_export_csv(ld_raw, "raw_counts")
ld_export_csv(ld_integrity, "integrity_checks")
ld_export_csv(ld_metadata_all, "metadata_all")
ld_export_csv(ld_metadata_matched, "metadata_matched")
ld_export_csv(ld_metadata_selected, "metadata_selected")
ld_export_csv(ld_sample_qc, "sample_QC")
ld_export_csv(ld_gene_annotation, "gene_annotation")
ld_export_csv(ld_filter_table, "expression_filter")
ld_export_csv(ld_coordinates, "MDS_PCA_coordinates")
ld_export_csv(ld_variance, "MDS_PCA_variance")
saveRDS(ld_counts_all, file.path(ld_output_dir, "Light_Damage_GSE314043_raw_counts.rds"))
saveRDS(ld_log_cpm_qc, file.path(ld_output_dir, "Light_Damage_GSE314043_TMM_logCPM.rds"))
capture.output(sessionInfo(), file = file.path(ld_output_dir, "Light_Damage_GSE314043_sessionInfo.txt"))

# Record input checksums so the downloaded sources can be verified unchanged.
ld_input_manifest <- data.frame(path = c(ld_count_file, ld_metadata_archive),
                                 md5 = unname(tools::md5sum(c(ld_count_file, ld_metadata_archive))))
ld_export_csv(ld_input_manifest, "input_manifest")

# Save only Light Damage objects, excluding XML pointers and transient responses.
# XML documents can be reread from ld_metadata_file when needed after reloading.
ld_checkpoint_objects <- ls(pattern = "^ld_")
ld_checkpoint_objects <- setdiff(ld_checkpoint_objects, c(
  "ld_metadata_xml", "ld_sample_nodes", "ld_release_response", "ld_lookup_response",
  "ld_all_lookup_result", "ld_batch_cache", "ld_checkpoint_objects",
  "ld_gtf_connection", "ld_gtf_gene_lines", "ld_gtf_chunk", "ld_gtf_genes",
  "ld_gtf_attribute", "ld_get_characteristic"
))
save(list = ld_checkpoint_objects,
     file = file.path(ld_output_dir, "Light_Damage_GSE314043_precheck_checkpoint_2026-11-03.RData"))
stopifnot(all(file.exists(c(ld_boxplot_file, ld_density_file, ld_mds_file, ld_pca_file))))
print(ld_coordinates)
print(head(ld_variance, 2))
cat("Light Damage prequality check completed successfully.\n")
