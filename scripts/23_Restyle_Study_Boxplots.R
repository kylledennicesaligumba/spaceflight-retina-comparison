# Restyle the six existing initial-expression boxplots to match Light Damage.
# Run this file with Source in RStudio, or with Rscript. Light Damage is excluded.
# Saved checkpoints are loaded into private environments; the working R session
# and original analysis scripts are not changed. No filtering, TMM, or DE is run.
# The original sample selections, expression scales, and boxplot statistics stay
# unchanged. Only sample order, displayed labels, colours, and layout are updated.

local({
  # Start in the repository root; override PREQUALITYCHECKS_DIR for existing outputs.
  base_dir <- Sys.getenv(
    "PREQUALITYCHECKS_DIR",
    unset = file.path(normalizePath(getwd(), mustWork = TRUE), "results")
  )
  stopifnot(dir.exists(base_dir))
  restyle_date <- "2026-11-03"
  colours <- c(blue = "#0072B2", orange = "#D55E00",
               green = "#009E73", purple = "#CC79A7")

  # Use the reviewed checkpoint for each study, rather than sourcing its analysis.
  studies <- list(
    glaucoma = list(
      folder = "glaucoma",
      checkpoint = "Glaucoma_GSE241782_precheck_checkpoint_2026-10-28.RData",
      matrix = "glaucoma_log_cpm", metadata = "glaucoma_metadata_selected",
      id = "sample_id", dimensions = c(55487L, 8L),
      output = "Glaucoma_GSE241782_logCPM_boxplot.png",
      title = "Glaucoma: initial expression distributions",
      subtitle = "Retina | naive baseline and glaucoma time points | F: female; M: male",
      ylab = "Log2 counts per million (CPM)"
    ),
    rp = list(
      folder = "retinitis_pigmentosa",
      checkpoint = "RP_GSE156533_precheck_checkpoint_2026-10-29.RData",
      matrix = "rp_log_cpm", metadata = "rp_metadata_matched",
      id = "sample_id", dimensions = c(54532L, 7L),
      output = "RP_GSE156533_logCPM_boxplot.png",
      title = "Retinitis Pigmentosa: initial expression distributions",
      subtitle = "Whole retina | postnatal day 90 | F: female; M: male",
      ylab = "Log2 counts per million (CPM)"
    ),
    lhon = list(
      folder = "lhon",
      checkpoint = "LHON_GSE304707_precheck_checkpoint_2026-10-30.RData",
      matrix = "lhon_log_cpm", metadata = "lhon_metadata_matched",
      id = "sample_id", dimensions = c(49714L, 12L),
      output = "LHON_GSE304707_logCPM_boxplot.png",
      title = "LHON: initial expression distributions",
      subtitle = "Retina | 12 months | all samples recorded as male",
      ylab = "Log2 counts per million (CPM)"
    ),
    dr = list(
      folder = "diabetic_retinopathy",
      checkpoint = "DR_GSE87433_precheck_checkpoint_2026-10-31.RData",
      matrix = "dr_expr_all", metadata = "dr_metadata_matched",
      id = "geo_accession", dimensions = c(24973L, 17L),
      output = "DR_GSE87433_expression_boxplot.png",
      title = "Diabetic Retinopathy: initial expression distributions",
      subtitle = "Whole retina | normal control, diabetic control, and diabetic + transplantation",
      ylab = "Processed expression value"
    ),
    uveitis = list(
      folder = "autoimmune_uveitis",
      checkpoint = "Autoimmune_Uveitis_GSE284505_precheck_checkpoint_2026-11-02.RData",
      matrix = "uveitis_log_cpm", metadata = "uveitis_metadata_selected",
      id = "sample_id", dimensions = c(56859L, 6L),
      output = "Autoimmune_Uveitis_GSE284505_logCPM_boxplot.png",
      title = "Autoimmune Uveitis: initial expression distributions",
      subtitle = "Retina | CON: healthy control; EAU: experimental autoimmune uveitis",
      ylab = "Log2 counts per million (CPM)"
    ),
    amd = list(
      folder = "amd",
      checkpoint = "AMD_GSE38671_precheck_checkpoint_2026-11-02.RData",
      matrix = "amd_expr_log2", metadata = "amd_metadata_all",
      id = "geo_accession", dimensions = c(35556L, 16L),
      output = "AMD_GSE38671_log2_expression_boxplot.png",
      title = "AMD model: initial expression distributions",
      subtitle = "CFH-deficiency model | neuroretina | young and aged WT and Cfh-/- groups",
      ylab = "Log2 expression value"
    )
  )

  # Match sample descriptions to matrix columns before using any group labels.
  align_metadata <- function(x, metadata, id_column) {
    ids <- as.character(metadata[[id_column]])
    stopifnot(!anyNA(ids), !anyDuplicated(ids),
              !is.null(colnames(x)), !anyDuplicated(colnames(x)))
    index <- match(colnames(x), ids)
    stopifnot(!anyNA(index))
    metadata <- metadata[index, , drop = FALSE]
    stopifnot(identical(as.character(metadata[[id_column]]), colnames(x)))
    metadata
  }

  # Export the same base-R boxplot style used for Light Damage.
  # Individual outlier dots are hidden, as before; all input rows are retained.
  save_plot <- function(stats, labels, fills, cfg, path) {
    png(path, width = 2400, height = 1800, res = 300)
    on.exit(dev.off(), add = TRUE)
    par(mar = c(9, 4.5, 4.5, 1))
    stats$names <- labels
    bxp(stats, las = 2, pars = list(boxfill = fills), outline = FALSE,
        ylab = cfg$ylab, cex.axis = 0.8)
    title(main = cfg$title, line = 2.5, cex.main = 1.15)
    mtext(cfg$subtitle, side = 3, line = 0.8, cex = 0.64)
  }

  # Process one study at a time so checkpoints never overwrite one another.
  reports <- list()
  for (study in names(studies)) {
    cfg <- studies[[study]]
    folder <- file.path(base_dir, cfg$folder)
    checkpoint <- file.path(folder, cfg$checkpoint)
    stopifnot(file.exists(checkpoint))
    saved <- new.env(parent = baseenv())
    load(checkpoint, envir = saved)
    stopifnot(exists(cfg$metadata, envir = saved, inherits = FALSE))
    metadata <- get(cfg$metadata, envir = saved, inherits = FALSE)

    # AMD's checkpoint predates its plot. Reproduce only the already reviewed
    # neuroretina selection and log2 conversion from its existing precheck script.
    # These microarray values are not CPM and are not normalised again here.
    if (study == "amd") {
      metadata <- metadata[metadata$tissue == "neuroretina", , drop = FALSE]
      x <- saved$amd_expr_all[, metadata$geo_accession, drop = FALSE]
      stopifnot(all(x > 0))
      x <- log2(x)
    } else {
      stopifnot(exists(cfg$matrix, envir = saved, inherits = FALSE))
      x <- get(cfg$matrix, envir = saved, inherits = FALSE)
    }
    stopifnot(is.matrix(x), is.numeric(x), all(is.finite(x)),
              identical(dim(x), cfg$dimensions))
    metadata <- align_metadata(x, metadata, cfg$id)

    # Use recorded treatment, time, genotype, age, and sex fields for labels.
    # Control groups appear first; the original order is retained within groups.
    if (study == "glaucoma") {
      treatment <- as.character(metadata$treatment)
      stopifnot(all(treatment %in% c("Na\u00efve", "Glaucoma")),
                all(metadata$sex %in% c("Female", "Male")))
      treatment[treatment == "Na\u00efve"] <- "Naive"
      group <- paste(treatment, metadata$time, sep = "-")
      group_levels <- c("Naive-0D", "Glaucoma-3D", "Glaucoma-2W", "Glaucoma-6W")
      labels <- paste0(group, "-", substr(metadata$sex, 1, 1),
                       " (", metadata$sample_id, ")")
      palette <- setNames(unname(colours), group_levels)
    } else if (study == "rp") {
      genotype_names <- c("RhoWT (wild-type)" = "WT",
                          "RhoP23H/WT (mutant)" = "P23H")
      group <- unname(genotype_names[as.character(metadata$genotype)])
      group_levels <- c("WT", "P23H")
      stopifnot(all(metadata$sex %in% c("Female", "Male")))
      labels <- paste0(metadata$sample_id, " (", substr(metadata$sex, 1, 1), ")")
      palette <- setNames(unname(colours[1:2]), group_levels)
    } else if (study == "lhon") {
      group <- as.character(metadata$genotype)
      group_levels <- c("WT", "YARS2", "COI", "YARS2/COI")
      labels <- ifelse(group == "WT", metadata$sample_id,
                       paste0(group, " (", metadata$sample_id, ")"))
      palette <- setNames(unname(colours), group_levels)
    } else if (study == "dr") {
      group_names <- c("normoglycemic_no" = "Normal control",
                       "diabetic_no" = "Diabetic control",
                       "diabetic_yes" = "Diabetic + transplantation")
      group <- unname(group_names[paste(metadata$disease,
                                        metadata$transplantation, sep = "_")])
      group_levels <- unname(group_names)
      replicate <- sub("^.*biological rep ([0-9]+)$", "\\1", metadata$title)
      stopifnot(all(grepl("^[0-9]+$", replicate)))
      label_prefix <- c("Normal control" = "Normal-control",
                        "Diabetic control" = "Diabetic-control",
                        "Diabetic + transplantation" = "Diabetic+transplant")
      labels <- paste0(unname(label_prefix[group]), "-", replicate)
      palette <- setNames(unname(colours[1:3]), group_levels)
    } else if (study == "uveitis") {
      group <- as.character(metadata$genotype)
      group_levels <- c("CON", "EAU")
      labels <- as.character(metadata$sample_id)
      palette <- setNames(unname(colours[1:2]), group_levels)
    } else if (study == "amd") {
      group <- paste(metadata$age_group, metadata$genotype, sep = "-")
      group_levels <- c("young-WT", "young-Cfh-/-", "aged-WT", "aged-Cfh-/-")
      labels <- sub("^NR ", "", metadata$title)
      palette <- setNames(unname(colours), group_levels)
    }

    # Check the labels, then arrange whole columns without changing their values.
    group <- factor(group, levels = group_levels)
    stopifnot(!anyNA(group), !anyNA(labels), !anyDuplicated(labels),
              length(labels) == ncol(x), all(table(group) > 0))
    plot_order <- order(group, seq_along(group))
    original_stats <- boxplot(x, plot = FALSE)
    plot_stats <- boxplot(x[, plot_order, drop = FALSE], plot = FALSE)
    stopifnot(
      identical(unname(plot_stats$stats),
                unname(original_stats$stats[, plot_order, drop = FALSE])),
      identical(unname(plot_stats$n), unname(original_stats$n[plot_order]))
    )

    # Save the exact original PNG once before replacing its usual filename.
    output <- file.path(folder, cfg$output)
    backup_dir <- file.path(folder, "archive",
                            paste0("boxplots_before_restyle_", restyle_date))
    backup <- file.path(backup_dir, cfg$output)
    if (file.exists(output) && !file.exists(backup)) {
      dir.create(backup_dir, recursive = TRUE, showWarnings = FALSE)
      stopifnot(file.copy(output, backup, overwrite = FALSE))
    }

    # Draw to a temporary file first so a failed export cannot replace a PNG.
    temporary_png <- tempfile(pattern = "boxplot_", tmpdir = folder, fileext = ".png")
    save_plot(plot_stats, labels[plot_order],
              unname(palette[as.character(group[plot_order])]), cfg, temporary_png)
    stopifnot(file.exists(temporary_png), file.info(temporary_png)$size > 0,
              file.copy(temporary_png, output, overwrite = TRUE))
    unlink(temporary_png)

    # Preserve a readable link from every displayed label to its original ID.
    # The additional metadata columns document the recorded experimental setup.
    key <- data.frame(
      plot_position = seq_along(plot_order),
      original_column = plot_order,
      original_sample_id = colnames(x)[plot_order],
      plot_label = labels[plot_order],
      plot_group = as.character(group[plot_order]),
      fill_colour = unname(palette[as.character(group[plot_order])]),
      metadata[plot_order, , drop = FALSE],
      check.names = FALSE, row.names = NULL
    )
    key_file <- file.path(folder, sub("\\.png$", "_sample_key.csv", cfg$output))
    write.csv(key, key_file, row.names = FALSE, na = "", fileEncoding = "UTF-8")

    # Report the inputs and confirm that every box retains its original statistics.
    reports[[study]] <- data.frame(
      study = study, checkpoint = checkpoint, expression = cfg$matrix,
      features = nrow(x), samples = ncol(x),
      groups = paste(names(table(group)), as.integer(table(group)),
                     sep = ": ", collapse = "; "),
      box_statistics_unchanged = TRUE,
      output = output, original_backup = backup, sample_key = key_file
    )
    message(cfg$folder, ": saved ", ncol(x), " boxes; statistics unchanged.")
    rm(saved, x)
  }

  # Keep the script and a short export record in the main prequalitychecks folder.
  report <- do.call(rbind, reports)
  report_file <- file.path(base_dir, paste0("Boxplot_restyle_record_", restyle_date, ".csv"))
  write.csv(report, report_file, row.names = FALSE, fileEncoding = "UTF-8")
  message("Finished: six updated plots. Light Damage was not replotted or changed.")
})
