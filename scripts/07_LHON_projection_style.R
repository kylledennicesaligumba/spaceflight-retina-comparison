# Portable repository copy; run from the repository root.
# Original working files are preserved. GitHub sync did not execute this code.
# Reusable LHON projection formatting, adapted from the approved RP style.
# Loading this file defines a function; it does not load or analyse datasets.
# The main precheck script calculates projections and explicitly saves PNGs.

stopifnot(
  requireNamespace("ggplot2", quietly = TRUE),
  requireNamespace("ggrepel", quietly = TRUE)
)

lhon_projection_plot <- function(projection, samples, title, axis_names) {
  # Use consistent genotype colours for MDS and PCA.
  palette <- c(
    WT = "#0072B2", YARS2 = "#D55E00",
    COI = "#009E73", "YARS2/COI" = "#CC79A7"
  )

  # Align descriptive sample rows to the projection's recorded sample order.
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

  # Keep points at their calculated coordinates; move only the labels.
  plot_data <- data.frame(
    axis1 = projection$x,
    axis2 = projection$y,
    group = factor(as.character(samples$group), levels = names(palette)),
    label = paste(samples$sample_id, toupper(substr(samples$sex, 1, 1)))
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
      colour = "Genotype"
    ) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      legend.position = "right",
      legend.title = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(hjust = 0.5)
    )
}
