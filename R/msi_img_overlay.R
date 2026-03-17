#' Rebuild a Cardinal image as a ggplot object
#'
#' This function extracts pixel information from \code{Cardinal::image()} and
#' reconstructs the image using \pkg{ggplot2}. It supports both ion images
#' defined by m/z values and variable images such as \code{"TIC"},
#' \code{"ROI"}, or \code{"Stage"}.
#'
#' @param msi_data A Cardinal MSI object.
#' @param value A single image target. If numeric, it is treated as an m/z
#'   value and passed to \code{Cardinal::image(..., mz = value)}. If character,
#'   it is treated as an image variable name and passed to
#'   \code{Cardinal::image(msi_data, value)}.
#' @param color A color vector passed to \code{Cardinal::image()}. If the vector
#'   has names, it is treated as a discrete color mapping and
#'   \code{scale_fill_manual()} will be used. Otherwise,
#'   \code{scale_fill_gradientn()} will be used.
#' @param xlim Optional numeric vector of length 2 specifying the x-axis range.
#'   If \code{NULL}, it will be inferred from \code{Cardinal::coord(msi_data)}.
#' @param ylim Optional numeric vector of length 2 specifying the y-axis range.
#'   If \code{NULL}, it will be inferred from \code{Cardinal::coord(msi_data)}.
#' @param show_legend Logical; whether to display the legend. Default is
#'   \code{FALSE}.
#' @param show_axes Logical; whether to display axis lines, ticks, and axis text.
#'   Default is \code{TRUE}. If \code{FALSE}, all axis elements are hidden.
#'
#' @return A \code{ggplot} object.
#'
#' @importFrom rlang .data
#'
#' @export
image2ggplot <- function(
    msi_data,
    value,
    color,
    xlim = NULL,
    ylim = NULL,
    show_legend = FALSE,
    show_axes = TRUE
) {
  if (is.null(xlim) || is.null(ylim)) {
    coord_df <- Cardinal::coord(msi_data)
    if (is.null(xlim)) xlim <- range(coord_df[, 1], na.rm = TRUE)
    if (is.null(ylim)) ylim <- range(coord_df[, 2], na.rm = TRUE)
  }

  if (is.numeric(value)) {
    img <- Cardinal::image(
      msi_data,
      mz = value,
      xlim = xlim,
      ylim = ylim,
      smooth = "gaussian",
      enhance = "histogram",
      col = color
    )
  } else if (is.character(value) && length(value) == 1) {
    img <- Cardinal::image(
      msi_data,
      value,
      xlim = xlim,
      ylim = ylim,
      smooth = "gaussian",
      enhance = "histogram",
      col = color
    )
  } else {
    stop("`value` must be either a numeric m/z or a single character variable name.")
  }

  df <- data.frame(
    x = img$plots[[1]]$marks$pixels$encoding$x,
    y = img$plots[[1]]$marks$pixels$encoding$y,
    fill = img$plots[[1]]$marks$pixels$encoding$color
  )

  full_grid <- expand.grid(
    x = seq(xlim[1], xlim[2], by = 1),
    y = seq(ylim[1], ylim[2], by = 1)
  )

  df <- merge(full_grid, df, by = c("x", "y"), all.x = TRUE)

  is_discrete <- !is.null(names(color))

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$x, y = .data$y, fill = .data$fill)) +
    ggplot2::geom_raster() +
    ggplot2::scale_x_continuous(
      limits = xlim,
      expand = c(0, 0)
    ) +
    ggplot2::scale_y_reverse(
      limits = rev(ylim),
      expand = c(0, 0)
    ) +
    ggplot2::coord_fixed() +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = if (show_legend) "right" else "none",
      axis.title = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank()
    )

  if (is_discrete) {
    p <- p + ggplot2::scale_fill_manual(
      values = color,
      na.translate = FALSE
    )
  } else {
    p <- p + ggplot2::scale_fill_gradientn(
      colors = color,
      na.value = "transparent"
    )
  }

  if (show_axes) {
    p <- p + ggplot2::theme(
      axis.line = ggplot2::element_line(color = "black"),
      axis.ticks = ggplot2::element_line(color = "black"),
      axis.text = ggplot2::element_text(color = "black")
    )
  } else {
    p <- p + ggplot2::theme(
      axis.line = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank()
    )
  }

  p
}


#' Overlay two MSI images reconstructed by ggplot2
#'
#' This function rebuilds two MSI images using \code{image2ggplot()} and
#' overlays them using \pkg{cowplot}. The second layer legend can be displayed
#' automatically without requiring the user to manually generate a separate
#' legend plot.
#'
#' @param msi_data01 A Cardinal MSI object for the background layer.
#' @param value01 A single image target for the first layer.
#' @param color01 A color vector for the first layer.
#' @param msi_data02 A Cardinal MSI object for the overlay layer.
#' @param value02 A single image target for the second layer.
#' @param color02 A color vector for the second layer.
#' @param xlim Optional numeric vector of length 2 specifying the x-axis range.
#' @param ylim Optional numeric vector of length 2 specifying the y-axis range.
#' @param show_legend Logical; whether to display the legend for the second
#'   layer. Default is \code{TRUE}.
#' @param show_axes Logical; whether to display axis lines, ticks, and axis text.
#'   Default is \code{TRUE}. If \code{FALSE}, all axis elements are hidden.
#'
#' @return A combined plot object.
#'
#' @export
msi_img_overlay <- function(
    msi_data01,
    value01,
    color01,
    msi_data02,
    value02,
    color02,
    xlim = NULL,
    ylim = NULL,
    show_legend = TRUE,
    show_axes = TRUE
) {
  if (is.null(xlim) || is.null(ylim)) {
    coord1 <- Cardinal::coord(msi_data01)
    coord2 <- Cardinal::coord(msi_data02)

    if (is.null(xlim)) {
      xlim <- range(c(coord1[, 1], coord2[, 1]), na.rm = TRUE)
    }
    if (is.null(ylim)) {
      ylim <- range(c(coord1[, 2], coord2[, 2]), na.rm = TRUE)
    }
  }

  p1 <- image2ggplot(
    msi_data = msi_data01,
    value = value01,
    color = color01,
    xlim = xlim,
    ylim = ylim,
    show_legend = FALSE,
    show_axes = show_axes
  )

  p2 <- image2ggplot(
    msi_data = msi_data02,
    value = value02,
    color = color02,
    xlim = xlim,
    ylim = ylim,
    show_legend = FALSE,
    show_axes = show_axes
  )

  combined <- cowplot::ggdraw(p1) +
    cowplot::draw_plot(p2, x = 0, y = 0, width = 1, height = 1)

  if (show_legend) {
    p2_leg <- image2ggplot(
      msi_data = msi_data02,
      value = value02,
      color = color02,
      xlim = xlim,
      ylim = ylim,
      show_legend = TRUE,
      show_axes = show_axes
    )

    legend <- cowplot::get_legend(
      p2_leg + ggplot2::theme(
        legend.position = "right",
        legend.justification = c(1, 0.5)
      )
    )

    combined <- cowplot::plot_grid(
      combined,
      legend,
      ncol = 2,
      rel_widths = c(1, 0.06)
    )
  }

  combined
}

