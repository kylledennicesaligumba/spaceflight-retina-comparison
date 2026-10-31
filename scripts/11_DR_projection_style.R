# Portable repository copy; start R in the repository root.
# This is a processed MICROARRAY precheck: no additional log2 transformation,
# CPM calculation or TMM normalization is applied to the supplied values.
# DR row labels are Entrez IDs; supplied platform annotation links them to gene symbols.
# Original working files are preserved; GitHub sync did not execute this code.
# Reusable DR projection formatting, following the approved RP/LHON style.
# Loading this file defines functions and colours; it does not analyse data.
stopifnot(
  requireNamespace("limma", quietly = TRUE),
  requireNamespace("ggplot2", quietly = TRUE),
  requireNamespace("ggrepel", quietly = TRUE)
)

# Keep the same group colours in both projection plots.
dr_projection_palette <- c(
  "Normal control" = "#0072B2",
  "Diabetic control" = "#D55E00",
  "Diabetic + transplantation" = "#009E73"
)

dr_projection_plot <- function(projection, samples, title, axis_names) {
  # Match the plotted sample descriptions to the calculated coordinate order.
  projection_ids <- colnames(projection$distance.matrix.squared)
  if (!is.null(projection_ids)) {
    sample_index <- match(projection_ids, as.character(samples$geo_accession))
    stopifnot(!anyNA(sample_index))
    samples <- samples[sample_index, , drop = FALSE]
  }
  stopifnot(
    nrow(samples) == length(projection$x),
    length(projection$x) == length(projection$y),
    all(as.character(samples$group) %in% names(dr_projection_palette))
  )

  # Use calculated coordinates for points; move only the GSM text labels.
  plot_data <- data.frame(
    axis1 = projection$x,
    axis2 = projection$y,
    group = factor(as.character(samples$group), levels = names(dr_projection_palette)),
    label = as.character(samples$geo_accession)
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
    ggplot2::scale_colour_manual(values = dr_projection_palette) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.2)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = 0.2)) +
    ggplot2::labs(
      title = title,
      x = sprintf("%s (%.0f%%)", axis_names[1], percentages[1]),
      y = sprintf("%s (%.0f%%)", axis_names[2], percentages[2]),
      colour = "Group"
    ) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      legend.position = "right",
      legend.title = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(hjust = 0.5)
    )
}

# Export a supplied plot at the same resolution as the previous figures.
dr_save_projection <- function(plot, path) {
  png(path, width = 2400, height = 1800, res = 300)
  on.exit(dev.off(), add = TRUE)
  print(plot)
  invisible(path)
}
