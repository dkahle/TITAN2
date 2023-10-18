#' Plot community level change
#'
#' Plot community level change
#'
#' @inheritParams plot-sumz
#' @inheritParams plot_taxa_ridges
#' @param sumz Plot sum z values? (Default = \code{TRUE})
#' @param points Plot points on sum z values? (Default = \code{FALSE})
#' @param ribbon Fill in axes under sum z values?  (Default = \code{TRUE})
#' @param density Plot densities of change points?  (Default = \code{TRUE})
#' @param change_points Plot ranges of change points? These are taken from the
#'   titan object. (Default = \code{TRUE})
#' @param alpha1,alpha2 Transparency of Z- and Z+ values, respectively.
#' @param legend.position ggplot2 legend position in relative coordinates of sum
#'   z plot. (Default = \code{c(.9, .9)})
#' @return A plot
#' @references Baker, ME and RS King.  2010. A new method for detecting and
#'   interpreting biodiversity and ecological community thresholds.  Methods in
#'   Ecology and Evolution 1(1): 25:37.
#' @references King, RS and ME Baker  2010. Considerations for identifying and
#'   interpreting ecological community thresholds. Journal of the North American
#'   Benthological Association 29(3):998-1008.
#' @note Should not be used with output objects from TITAN v1.0.
#' @author M. Baker and R. King
#' @seealso [plot_sumz()], [plot_cps()]
#' @export
#' @examples
#'
#' data(glades.titan)
#'
#' plot_sumz_density(glades.titan)
#' plot_sumz_density(glades.titan, trans = "log10")
#'
#'
#'
plot_sumz_density <- function(
  titan.out,
  filter = TRUE,
  sumz = TRUE,
  points = FALSE,
  ribbon = TRUE, density = TRUE,
  change_points = TRUE,
  sumz1 = TRUE, sumz2 = TRUE,
  xlim = c(min(titan.out$env), max(titan.out$envcls)),
  xlabel = "Environmental Gradient",
  y1label = NULL,
  y2label = "Density",
  alpha1 = 0.65,
  alpha2 = 0.5,
  col1 = "steelblue4", col2 = "red",
  trans = "identity",
  legend.position = c(.9, .9),
  ...,
  axis.text.x, axis.text.y,
  axis.title.x, axis.title.y
) {

  # cran guard
  X1 <- NULL; rm(X1)
  X2 <- NULL; rm(X2)


  # deal with axis font defaulting
  if (missing(axis.text.x)) axis.text.x <- if (is.null(theme_get()$axis.text.x$size)) theme_get()$axis.text.x$size else theme_get()$text$size
  if (missing(axis.text.y)) axis.text.y <- if (is.null(theme_get()$axis.text.y$size)) theme_get()$axis.text.y$size else theme_get()$text$size
  if (missing(axis.title.x)) axis.title.x <- if (is.null(theme_get()$axis.title.x$size)) theme_get()$axis.title.x$size else theme_get()$titletext$size
  if (missing(axis.title.y)) axis.title.y <- if (is.null(theme_get()$axis.title.y$size)) theme_get()$axis.title.y$size else theme_get()$titletext$size


  if (filter) {
    ivz <- titan.out$ivz.f
    sumz1max <- titan.out$sumz.cp[3, 1]
    sumz1quant <- titan.out$sumz.cp[3, c(2, 6)]
    sumz2max <- titan.out$sumz.cp[4, 1]
    sumz2quant <- titan.out$sumz.cp[4, c(2, 6)]
    sumz1lab <- "fsumz-"
    sumz2lab <- "fsumz+"
    maxPsumz <- titan.out$maxFsumz
    if (is.null(y1label)) y1label <- "Filtered Sum(z)"
  } else {
    ivz <- titan.out$ivz
    sumz1max <- titan.out$sumz.cp[1, 1]
    sumz1quant <- titan.out$sumz.cp[1, c(2, 6)]
    sumz2max <- titan.out$sumz.cp[2, 1]
    sumz2quant <- titan.out$sumz.cp[2, c(2, 6)]
    sumz1lab <- "sumz-"
    sumz2lab <- "sumz+"
    maxPsumz <- titan.out$maxSumz
    if (is.null(y1label)) y1label <- "Unfiltered Sum(z)"
  }


  if (sumz) {
    if (sumz1 & sumz2) {
      ggplot() +
        (if (points) geom_point(aes(x = titan.out$envcls, y = ivz[, 1], color = "Z-")) else geom_blank()) +
        (if (ribbon) geom_ribbon(aes(x = titan.out$envcls, ymin = 0, ymax = ivz[, 1], fill = "Z-"), alpha = alpha1) else geom_blank()) +
        geom_line(aes(x = titan.out$envcls, y = ivz[, 1], color = "Z-")) +
        (if (points) geom_point(aes(x = titan.out$envcls, y = ivz[, 2], color = "Z+")) else geom_blank()) +
        (if (ribbon) geom_ribbon(aes(x = titan.out$envcls, ymin = 0, ymax = ivz[, 2], fill = "Z+"), alpha = alpha2) else geom_blank()) +
        geom_line(aes(x = titan.out$envcls, y = ivz[, 2], color = "Z+")) +
        scale_x_continuous(xlabel, trans = trans, limits = xlim) +
        scale_color_manual("",
          breaks = c("Z-", "Z+"),
          values = c("Z-" = col1, "Z+" = col2),
          guide = "none"
        ) +
        scale_fill_manual("",
          breaks = c("Z-", "Z+"),
          values = c("Z-" = col1, "Z+" = col2)
        ) +
        labs(y = y1label, legend = "") +
        theme_classic() +
        theme(
          axis.text.x = element_text(size = axis.text.x),
          axis.text.y = element_text(size = axis.text.y),
          axis.title.x = element_text(size = axis.title.x),
          axis.title.y = element_text(size = axis.title.y),
          plot.margin = structure(
            c(0, 5.5, 5.5, 5.5),
            class = c("margin", "unit"),
            valid.unit = 8L,
            unit = "pt"
          ),
          legend.position = legend.position
        ) -> pbot
    } else {
      if (sumz1) {
        ggplot() +
          (if (points) geom_point(aes(x = titan.out$envcls, y = ivz[, 1]), color = col1) else geom_blank()) +
          (if (ribbon) geom_ribbon(aes(x = titan.out$envcls, ymin = 0, ymax = ivz[, 1]), fill = col1, alpha = alpha1) else geom_blank()) +
          geom_line(aes(x = titan.out$envcls, y = ivz[, 1]), color = col1) +
          scale_x_continuous(xlabel, trans = trans, limits = xlim) +
          labs(y = y1label) +
          theme_classic() +
          theme(
            axis.text.x = element_text(size = axis.text.x),
            axis.text.y = element_text(size = axis.text.y),
            axis.title.x = element_text(size = axis.title.x),
            axis.title.y = element_text(size = axis.title.y),
            plot.margin = structure(
              c(0, 5.5, 5.5, 5.5),
              class = c("margin", "unit"),
              valid.unit = 8L,
              unit = "pt"
            ),
            legend.position = legend.position
          ) -> pbot
      } else {
        if (sumz2) {
          ggplot() +
            (if (points) geom_point(aes(x = titan.out$envcls, y = ivz[, 2]), color = col2) else geom_blank()) +
            (if (ribbon) geom_ribbon(aes(x = titan.out$envcls, ymin = 0, ymax = ivz[, 2]), fill = col2, alpha = alpha2) else geom_blank()) +
            geom_line(aes(x = titan.out$envcls, y = ivz[, 2]), color = col2) +
            scale_x_continuous(xlabel, trans = trans, limits = xlim) +
            labs(y = y1label) +
            theme_classic() +
            theme(
              axis.text.x = element_text(size = axis.text.x),
              axis.text.y = element_text(size = axis.text.y),
              axis.title.x = element_text(size = axis.title.x),
              axis.title.y = element_text(size = axis.title.y),
              plot.margin = structure(
                c(0, 5.5, 5.5, 5.5),
                class = c("margin", "unit"),
                valid.unit = 8L,
                unit = "pt"
              ),
              legend.position = legend.position
            ) -> pbot
        }
      }
    }
  }
  maxFsumz <- data.frame(titan.out$maxFsumz)

  if (density) {
    ggplot(maxFsumz) +
      (if (sumz1) geom_density(aes(x = X1), color = col1, fill = col1, alpha = alpha1) else geom_blank()) +
      (if (sumz2) geom_density(aes(x = X2), color = col2, fill = col2, alpha = alpha2) else geom_blank()) +
      scale_x_continuous(if (sumz) "" else xlabel, trans = trans, limits = xlim) +
      labs(y = y2label) +
      theme_classic() +
      theme(
        axis.text.x = (if (sumz) element_blank() else element_text(size = axis.text.x)),
        axis.text.y = element_text(size = axis.text.y),
        axis.title.x = element_text(size = axis.title.x),
        axis.title.y = element_text(size = axis.title.y),
        plot.margin = structure(
          c(if (change_points) 0 else 5.5, 5.5, 0, 5.5),
          class = c("margin", "unit"),
          valid.unit = 8L,
          unit = "pt"
        )
      ) -> pmid
  }

  if (change_points) {
    ggplot() +
      (if (sumz1) geom_segment(aes(x = sumz1quant[1], y = sumz1lab, xend = sumz1quant[2], yend = sumz1lab), col = col1) else geom_blank()) +
      (if (sumz1) geom_point(aes(x = sumz1max, y = sumz1lab), col = col1, alpha = alpha1, size = 5) else geom_blank()) +
      (if (sumz2) geom_segment(aes(x = sumz2quant[1], y = sumz2lab, xend = sumz2quant[2], yend = sumz2lab), col = col2) else geom_blank()) +
      (if (sumz2) geom_point(aes(x = sumz2max, y = sumz2lab), col = col2, alpha = alpha2, size = 5) else geom_blank()) +
      scale_x_continuous(if (sumz) "" else xlabel, trans = trans, limits = xlim) +
      scale_y_discrete("") +
      theme_classic() +
      theme(
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = axis.text.y),
        axis.title.x = element_text(size = axis.title.x),
        axis.title.y = element_text(size = axis.title.y),
        plot.margin = structure(
          c(5.5, 5.5, 0, 5.5),
          class = c("margin", "unit"),
          valid.unit = 8L,
          unit = "pt"
        )
      ) -> ptop
  }


  if (density) {
    if (change_points) {
      if (sumz) {
        plot_grid(ptop, pmid, pbot, ncol = 1, rel_heights = c(2, 3, 6), align = "v")
      } else {
        plot_grid(ptop, pmid, ncol = 1, rel_heights = c(2, 3), align = "v")
      }
    } else {
      if (sumz) {
        plot_grid(pmid, pbot, ncol = 1, rel_heights = c(3, 6), align = "v")
      } else {
        plot_grid(pmid, ncol = 1)
      }
    }
  } else {
    if (sumz) {
      if (change_points) {
        plot_grid(ptop, pbot, ncol = 1, rel_heights = c(2, 6), align = "v")
      } else {
        plot_grid(pbot)
      }
    }
  }
}
