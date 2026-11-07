# OSD-758 projection formatting, adapted from the approved RP/LHON/uveitis style.
# Loading this file defines a plotting function without loading expression data.
stopifnot(
  requireNamespace("limma", quietly = TRUE),
  requireNamespace("ggplot2", quietly = TRUE),
  requireNamespace("ggrepel", quietly = TRUE)
)

osd758_projection_plot <- function(projection, samples, title, axis_names) {
  # Keep the same condition colours in MDS and PCA.
  palette <- c(HGC = "#0072B2", FLT_0G = "#D55E00")

  # Match sample descriptions to the projection's recorded sample order.
  projection_ids <- colnames(projection$distance.matrix.squared)
  if (!is.null(projection_ids)) {
    sample_index <- match(projection_ids, as.character(samples$sample))
    stopifnot(!anyNA(sample_index))
    samples <- samples[sample_index, , drop = FALSE]
  }
  stopifnot(
    nrow(samples) == length(projection$x),
    length(projection$x) == length(projection$y),
    all(is.finite(projection$x)),
    all(is.finite(projection$y)),
    all(as.character(samples$condition) %in% names(palette))
  )

  # Preserve calculated point coordinates and shorten only the displayed IDs.
  plot_data <- data.frame(
    axis1 = projection$x,
    axis2 = projection$y,
    condition = factor(as.character(samples$condition), levels = names(palette)),
    label = sub("^MHU-8_RTN_(FLT|HGC)_", "", as.character(samples$sample))
  )
  percentages <- 100 * projection$var.explained[projection$dim.plot]
  stopifnot(length(percentages) == 2, all(is.finite(percentages)))

  ggplot2::ggplot(
    plot_data, ggplot2::aes(x = axis1, y = axis2, colour = condition)
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
    ggplot2::scale_colour_manual(
      values = palette, labels = c(HGC = "HGC", FLT_0G = "Flight 0G")
    ) +
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
