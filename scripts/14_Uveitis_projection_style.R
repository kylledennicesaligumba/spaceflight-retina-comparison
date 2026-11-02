# Uveitis projection formatting, following the approved RP/LHON style.
# Loading this file defines a function; the main script calculates and saves plots.
stopifnot(
  requireNamespace("limma", quietly = TRUE),
  requireNamespace("ggplot2", quietly = TRUE),
  requireNamespace("ggrepel", quietly = TRUE)
)

uveitis_projection_plot <- function(projection, samples, title, axis_names) {
  # Use the same condition colours for MDS and PCA.
  palette <- c(CON = "#0072B2", EAU = "#D55E00")

  # Match sample descriptions to the projection's recorded sample order.
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

  # Keep points at their calculated coordinates and use the recorded sample IDs.
  plot_data <- data.frame(
    axis1 = projection$x,
    axis2 = projection$y,
    group = factor(as.character(samples$group), levels = names(palette)),
    label = as.character(samples$sample_id)
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
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.2)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = 0.2)) +
    ggplot2::labs(
      title = title,
      x = sprintf("%s (%.0f%%)", axis_names[1], percentages[1]),
      y = sprintf("%s (%.0f%%)", axis_names[2], percentages[2]),
      colour = "Condition"
    ) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      legend.position = "right",
      legend.title = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(hjust = 0.5)
    )
}
