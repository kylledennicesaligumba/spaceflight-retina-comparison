# Portable repository copy. Start R in the repository root.
# This file records the interactive workflow; later result/closure notes
# supersede earlier pending notes. See the execution guides under docs/.
# Original working files are preserved. Repository sync did not run analysis.
# Restyle RP MDS and glaucoma MDS/PCA using the approved RP PCA appearance.
# Prepared 29 October 2026 for manual source() in the current RP R session.
# Reuses existing coordinates. Does not filter, normalize or fit projections.
# Overwrites the three existing PNG exports after successful construction.
# If glaucoma objects are absent, loads the dated checkpoint in an isolated
# environment; it does not load that workspace into the current RP workspace.

stopifnot(
  requireNamespace("ggplot2", quietly = TRUE),
  requireNamespace("ggrepel", quietly = TRUE)
)

retina_projection_plot <- function(
  projection, samples, palette, title, legend_title, axis_names
) {
  # Match sample annotations to the order stored with the projection.
  projection_ids <- colnames(projection$distance.matrix.squared)
  if (!is.null(projection_ids)) {
    sample_index <- match(projection_ids, as.character(samples$sample_id))
    stopifnot(!anyNA(sample_index))
    samples <- samples[sample_index, , drop = FALSE]
  }
  stopifnot(
    nrow(samples) == length(projection$x),
    length(projection$x) == length(projection$y),
    all(as.character(samples$group) %in% names(palette))
  )

  plot_data <- data.frame(
    axis1 = projection$x,
    axis2 = projection$y,
    group = samples$group,
    label = paste(samples$sample_id, substr(samples$sex, 1, 1))
  )
  percentages <- 100 * projection$var.explained[projection$dim.plot]

  ggplot2::ggplot(
    plot_data, ggplot2::aes(x = axis1, y = axis2, colour = group)
  ) +
    ggplot2::geom_point(size = 3) +
    ggrepel::geom_text_repel(
      ggplot2::aes(label = label),
      size = 3.5, seed = 1,
      box.padding = 0.6, point.padding = 0.4,
      min.segment.length = 0,
      max.overlaps = Inf, max.time = 2,
      show.legend = FALSE
    ) +
    ggplot2::scale_colour_manual(values = palette) +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = 0.2)
    ) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = 0.2)
    ) +
    ggplot2::labs(
      title = title,
      x = sprintf("%s (%.0f%%)", axis_names[1], percentages[1]),
      y = sprintf("%s (%.0f%%)", axis_names[2], percentages[2]),
      colour = legend_title
    ) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      legend.position = "right",
      legend.title = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(hjust = 0.5)
    )
}

retina_save_projection <- function(plot, path) {
  png(path, width = 2400, height = 1800, res = 300)
  on.exit(dev.off(), add = TRUE)
  print(plot)
  invisible(path)
}

stopifnot(
  exists("rp_mds", envir = .GlobalEnv, inherits = FALSE),
  exists("rp_dge_qc", envir = .GlobalEnv, inherits = FALSE)
)

# Restore glaucoma objects only if they are not already available.
retina_glaucoma_objects <- c(
  "glaucoma_mds", "glaucoma_pca", "glaucoma_dge_qc"
)
if (all(vapply(
  retina_glaucoma_objects,
  exists, logical(1), envir = .GlobalEnv, inherits = FALSE
))) {
  retina_glaucoma_env <- .GlobalEnv
} else {
  retina_glaucoma_env <- new.env(parent = baseenv())
  invisible(load(
    file.path("results", "glaucoma", "Glaucoma_GSE241782_precheck_checkpoint_2026-10-28.RData"),
    envir = retina_glaucoma_env
  ))
}
stopifnot(all(vapply(
  retina_glaucoma_objects,
  exists, logical(1), envir = retina_glaucoma_env, inherits = FALSE
)))

retina_glaucoma_samples <- get(
  "glaucoma_dge_qc", envir = retina_glaucoma_env
)$samples
retina_glaucoma_palette <- setNames(
  c("#0072B2", "#D55E00", "#009E73", "#CC79A7"),
  levels(factor(retina_glaucoma_samples$group))
)
stopifnot(
  length(levels(factor(retina_glaucoma_samples$group))) == 4L,
  nrow(retina_glaucoma_samples) == 8L
)

# RP MDS: retain the WT/P23H colours.
rp_mds_plot <- retina_projection_plot(
  rp_mds, rp_dge_qc$samples,
  palette = c(WT = "#0072B2", P23H = "#D55E00"),
  title = "RP retinal libraries: MDS after TMM",
  legend_title = "Genotype",
  axis_names = c("Leading logFC dim 1", "Leading logFC dim 2")
)
rp_draw_mds <- function() print(rp_mds_plot)

# Glaucoma MDS: retain the four treatment/time groups and original colours.
glaucoma_mds_plot <- retina_projection_plot(
  get("glaucoma_mds", envir = retina_glaucoma_env),
  retina_glaucoma_samples,
  palette = retina_glaucoma_palette,
  title = "Glaucoma retinal libraries: MDS after TMM",
  legend_title = "Treatment / time",
  axis_names = c("Leading logFC dim 1", "Leading logFC dim 2")
)
glaucoma_draw_mds <- function() print(glaucoma_mds_plot)

# Glaucoma PCA: reuse the existing common-top-500 PCA coordinates.
glaucoma_pca_plot <- retina_projection_plot(
  get("glaucoma_pca", envir = retina_glaucoma_env),
  retina_glaucoma_samples,
  palette = retina_glaucoma_palette,
  title = "Glaucoma retinal libraries: PCA after TMM",
  legend_title = "Treatment / time",
  axis_names = c("PC1", "PC2")
)
glaucoma_draw_pca <- function() print(glaucoma_pca_plot)

retina_projection_plots <- list(
  RP_MDS = rp_mds_plot,
  Glaucoma_MDS = glaucoma_mds_plot,
  Glaucoma_PCA = glaucoma_pca_plot
)
retina_projection_files <- c(
  RP_MDS = file.path("results", "retinitis_pigmentosa", "RP_GSE156533_TMM_MDS.png"),
  Glaucoma_MDS = file.path("results", "glaucoma", "Glaucoma_GSE241782_TMM_MDS.png"),
  Glaucoma_PCA = file.path("results", "glaucoma", "Glaucoma_GSE241782_TMM_PCA_top500.png")
)
stopifnot(all(dir.exists(dirname(retina_projection_files))))

# Display and export each plot; the save helper closes its own PNG device.
for (plot_name in names(retina_projection_plots)) {
  print(retina_projection_plots[[plot_name]])
  retina_save_projection(
    retina_projection_plots[[plot_name]],
    retina_projection_files[[plot_name]]
  )
}

print(vapply(retina_projection_files, file.exists, logical(1)))
# Expected: TRUE for RP_MDS, Glaucoma_MDS and Glaucoma_PCA.
# To redisplay separately:
# rp_draw_mds()
# glaucoma_draw_mds()
# glaucoma_draw_pca()
# These exports do not save a new RData checkpoint.
