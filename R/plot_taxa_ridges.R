#' Taxon change point ridge plots
#'
#' @param titan.out A TITAN output object.
#' @param z1 A logical specifying whether decreasing taxa should be plotted.
#' @param z2 A logical specifying whether increasing taxa should be plotted.
#' @param xlabel A character string for the x axis label.
#' @param n_ytaxa The maximum number of taxa to be plotted.
#' @param printspp A logical specifying whether the sppmax table should be
#'   printed.
#' @param z1_fill_low,z1_fill_high,z2_fill_low,z2_fill_high Respective fill
#'   colors passed to [scale_fill_gradient()]
#' @param grid The `grid` argument of [theme_ridges()]. Setting this to
#'   `FALSE` removes horizontal lines from the plot, see examples.
#' @param trans a scale transformation to be applied to the x-axis through
#'   ggplot2. e.g. `"log10"` or `"sqrt"`.
#' @param xaxis Logical; should the x-axis be plotted?
#' @param xlim x axis limits, e.g. \code{xlim = c(0,10)}.
#' @param axis.text.x,axis.text.y,axis.title.x,axis.title.y Font sizes of
#'   respective axis text. Passed as the size argument to the ggplot2 thematic
#'   elements of the same name. Defaults to current default ggplot2 theme.
#' @param bw The bandwidth of used in the kernel density estimate; see
#'   [density()].
#' @param rel_heights Argument to pass to [cowplot::plot_grid()] as an override
#'   for the defaults
#' @param breaks Argument to pass to [geom_density_ridges()].
#' Determines values and labels of breaks (ticks) on the x-axis.
#' @param d_lines Argument to pass to [geom_density_ridges()].
#' Short for density lines, a logical determining whether the median change
#' point values from estimated from each taxon's density function, is plotted.
#' If `TRUE`, a dashed, black vertical line will be plotted.
#' @param ... Arguments to pass to [geom_density_ridges()]
#' @return A plot of decreasing and/or increasing taxon-specific change points
#'   along the environmental gradient.
#' @references Baker, ME and RS King.  2010. A new method for detecting and
#'   interpreting biodiversity and ecological community thresholds.  Methods in
#'   Ecology and Evolution 1(1): 25:37.
#' @references King, RS and ME Baker  2010. Considerations for identifying and
#'   interpreting ecological community thresholds. Journal of the North American
#'   Benthological Association 29(3):998-1008.
#' @note Should not be used with output objects from TITAN v1.0.
#' @author M. Baker, R. King, D. Kahle
#' @seealso [plot_sumz()], [plot_cps()]
#' @export
#' @examples
#'
#' data(glades.titan)
#'
#'
#' # basic usage
#' plot_taxa_ridges(glades.titan)
#'
#'
#' \dontrun{ reduce R CMD check time
#'
#' # manipulating the x axis thematic components
#' plot_taxa_ridges(glades.titan, grid = FALSE)
#' plot_taxa_ridges(glades.titan, xaxis = TRUE)
#'
#'
#' # applying scale transformations
#' plot_taxa_ridges(glades.titan, trans = "log10", xlim = c(1,150))
#'
#'
#' # removing z2
#' plot_taxa_ridges(glades.titan, z2 = FALSE, trans = "sqrt", xlim = c(0, 150))
#'
#'
#' # more styling
#' plot_taxa_ridges(glades.titan,
#'   axis.text.x = 12,
#'   axis.text.y = 4,
#'   axis.title.x = 14,
#'   xlabel = expression(paste("Surface water total phosphorus ("*mu*"g/l)"))
#' )
#'
#' }
#'
#'
#'
plot_taxa_ridges <- function(
    titan.out,
    z1 = TRUE, z2 = TRUE,
    z1_fill_low="lightblue", z1_fill_high="black",
    z2_fill_low="pink", z2_fill_high="red",
    # pur.cut = titan.out$arguments[[7]],
    # rel.cut = titan.out$arguments[[8]],
    xlabel = "Environmental Gradient",
    n_ytaxa = 90,
    printspp = FALSE,
    grid = TRUE,
    d_lines = FALSE,
    trans = "identity",
    xaxis = !grid,
    xlim,
    breaks = ggplot2::waiver(),
    bw,
    rel_heights,
    ...,
    axis.text.x, axis.text.y,
    axis.title.x, axis.title.y
) {

  # cran guard
  filter <- NULL; rm(filter)
  maxgrp <- NULL; rm(maxgrp)
  obsiv.prob <- NULL; rm(obsiv.prob)
  p05 <- NULL; rm(p05)
  p95 <- NULL; rm(p95)
  y <- NULL; rm(y)
  id <- NULL; rm(id)
  zscore <- NULL; rm(zscore)
  zenv.cp <- NULL; rm(zenv.cp)
  chk_pts <- NULL; rm(chk_pts)
  desc <- NULL; rm(desc)
  purity <- NULL; rm(purity)
  reliability <- NULL; rm(reliability)
  z.median <- NULL; rm(z.median)
  `5%` <- NULL; rm(`5%`)
  `95%` <- NULL; rm(`95%`)
  Q5 <- NULL; rm(Q5)
  Q50 <- NULL; rm(Q50)
  Q95 <- NULL; rm(Q95)

  # deal with axis font defaulting
  if (missing(axis.text.x)) axis.text.x <- if (is.null(theme_get()$axis.text.x$size)) theme_get()$axis.text.x$size else theme_get()$text$size
  if (missing(axis.text.y)) axis.text.y <- if (is.null(theme_get()$axis.text.y$size)) theme_get()$axis.text.y$size else theme_get()$text$size
  if (missing(axis.title.x)) axis.title.x <- if (is.null(theme_get()$axis.title.x$size)) theme_get()$axis.title.x$size else theme_get()$titletext$size
  if (missing(axis.title.y)) axis.title.y <- if (is.null(theme_get()$axis.title.y$size)) theme_get()$axis.title.y$size else theme_get()$titletext$size

  imax <- titan.out$arguments[["imax"]]
  boot <- titan.out$arguments[["boot"]] > 0.5

  sppmax <- titan.out %>%
    pluck("sppmax") %>%
    tibble::as_tibble() %>%
    mutate(id = row.names(titan.out$sppmax)) %>%
    rename("Q5" = "5%","Q10" = "10%", "Q50" = "50%", "Q90"="90%", "Q95" = "95%")

  # next 3 lines are unnecessary
  # sppmax$filter <- 0
  # reliable_taxa_ndcs <- with(sppmax, which(purity>=pur.cut & reliability>=rel.cut))
  # sppmax$filter[reliable_taxa_ndcs] <- sppmax$maxgrp[reliable_taxa_ndcs]

  num.dcr <- sum(sppmax$filter == 1L)
  num.ncr <- sum(sppmax$filter == 2L)

  cli_alert("There are {sum(num.dcr, num.ncr)} indicator taxa of {n_ytaxa} possible for plotting")
  cat(paste(
    col_red("\u2198"), col_black(glue("Number of Decreasers = {num.dcr}")), "\n"
  ))
  cat(paste(
    cli::col_blue("\u2197"), col_black(glue("Number of Increasers = {num.ncr}"))
  ))

  sppmax_filter_index <- which(sppmax$filter > 0) #changed name of index for pure/reliable taxa (more intuitive)

  sppmax_ny <- sppmax %>%
    arrange(desc(purity), desc(reliability), desc(z.median)) %>%
    slice(1:min(n(), n_ytaxa))

  n_filter_decr <- sum(sppmax_ny$filter == 1L)
  n_filter_incr <- sum(sppmax_ny$filter == 2L)

  gdf <- sppmax_filter_index %>% #use of new pure/reliable filter index
    map_dfr(~
              cbind(
                tibble::tibble(
                  "chk_pts" = sort(titan.out$metricArray[.x,2,])
                ),
                sppmax %>% slice(.x)
              )
    ) %>%
    mutate(filter = case_when(
      (filter == 1) ~ -1L,
      (filter == 2) ~ +1L
    ))

  gdf <- gdf %>%
    arrange(desc(purity), desc(reliability), desc(z.median)) %>%
    slice(1:min(n(), n_ytaxa*length(titan.out$metricArray[1,1,]))) %>%
    filter(Q5 < chk_pts, chk_pts < Q95) #changed to Q5 and Q95 to reflect renaming

  # compute pooled statistics
  if (missing(xlim)) xlim <- grDevices::extendrange(range(gdf$chk_pts), f = 0.05)
  xmin <- xlim[1]
  xmax <- xlim[2]
  if (missing(bw)) bw <- stats::bw.nrd0(do.call(trans, list(gdf$chk_pts)))


  if (z2) {

    ptop <- ggplot(
      filter(gdf, filter == -1L),
      aes(x = chk_pts, y = reorder(id, -Q50, median), fill = zscore)
    ) +
      ggridges::geom_density_ridges(
        quantile_lines = d_lines, #made an argument
        vline_size = 0.25,
        quantiles = 2,
        vline_color = "black", #suggest using black
        vline_linetype = 2, #distinguish from 50% by linetype = 2
        color = NA,
        alpha = 0.6,
        from = do.call(trans, list(xmin)),
        to = do.call(trans, list(xmax)),
        bandwidth = bw
      ) +
      geom_segment(aes(x=Q50, xend=Q50, y=as.numeric(reorder(id, -Q50, median)), yend=as.numeric(reorder(id, -Q50, median))+.9), linewidth=0.25, color="black"
                   #changed to black, removed "linetype"
      ) +
      scale_x_continuous(limits = xlim, breaks=breaks, expand = c(0,0), trans = trans) +
      #added breaks=breaks as argument
      scale_y_discrete("") +
      scale_fill_gradient("Z-Score", low = z1_fill_low, high = z1_fill_high) +
      ggridges::theme_ridges(center_axis_labels = TRUE, grid = grid) +
      theme(
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = axis.text.y),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = axis.title.y),
        axis.line = (if (xaxis) { #added if (xaxis) {...} to all 3 places needed to plot x if xaxis=TRUE
          do.call(element_line, theme_get()$axis.line)
        } else {
          element_blank()
        }),
        axis.line.y = element_blank(),
        plot.margin = structure(
          c(5.5, 5.5, 0, 5.5),
          class = c("margin", "unit"),
          valid.unit = 8L,
          unit = "pt"
        )
      )

    pbottom <- ggplot(
      filter(gdf, filter == +1L),
      aes(x = chk_pts, y = reorder(id, Q50, median), fill = zscore)
    ) +
      ggridges::geom_density_ridges(
        quantile_lines = d_lines, #made an argument
        vline_size = 0.25,
        quantiles = 2,
        vline_color = "black", #suggest using black
        vline_linetype = 2, #distinguish from 50% by linetype = 2
        color = NA,
        alpha = 0.6,
        from = do.call(trans, list(xmin)),
        to = do.call(trans, list(xmax)),
        bandwidth = bw
      ) +
      geom_segment(aes(x=Q50, xend=Q50, y=as.numeric(reorder(id, Q50, median)), yend=as.numeric(reorder(id, Q50, median))+.9), linewidth=0.25, color="black" #changed to black, removed linetype
      ) +
      ggridges::theme_ridges(center_axis_labels = TRUE, grid = grid) +
      scale_x_continuous(xlabel, limits = xlim, breaks=breaks, expand = c(0,0), trans = trans) +
      scale_y_discrete("") +
      scale_fill_gradient("Z-Score", low = z2_fill_low, high = z2_fill_high) +
      theme(
        axis.text.x = element_text(size = axis.text.x),
        axis.text.y = element_text(size = axis.text.y),
        axis.title.x = element_text(size = axis.title.x),
        axis.title.y = element_text(size = axis.title.y),
        axis.line = (if (xaxis) {#added if (xaxis) {...} to all 3 places needed to plot x if xaxis=TRUE
          do.call(element_line, theme_get()$axis.line)
        } else {
          element_blank()
        }),
        axis.line.y = element_blank()
      )

  } else {

    if (z1) {

      ptop <- ggplot(
        filter(gdf, filter == -1),
        aes(x = chk_pts, y = reorder(id, -Q50, median), fill = zscore)
      ) +
        ggridges::theme_ridges(center_axis_labels = TRUE, grid = grid) +
        ggridges::geom_density_ridges(
          quantile_lines = d_lines,  #made an argument
          vline_size = 0.25,
          quantiles = 2,
          vline_color = "black", #suggest using black
          vline_linetype = 2, #distinguish from 50% by linetype = 2
          color = NA,
          alpha = 0.6,
          from = do.call(trans, list(xmin)),
          to = do.call(trans, list(xmax)),
          bandwidth = bw
        )+
        geom_segment(aes(x=Q50, xend=Q50, y=as.numeric(reorder(id, -Q50, median)), yend=as.numeric(reorder(id, -Q50, median))+.9), linewidth=0.25, color="black" #changed to black, removed linetype
        ) +
        scale_x_continuous(xlabel, limits = xlim, breaks=breaks, expand = c(0,0), trans = trans) +
        scale_y_discrete("") +
        scale_fill_gradient("Z-Score", low = z1_fill_low, high = z1_fill_high) +
        theme(
          axis.text.x = element_text(size = axis.text.x),
          axis.text.y = element_text(size = axis.text.y),
          axis.title.x = element_text(size = axis.title.x),
          axis.title.y = element_text(size = axis.title.y),
          axis.line = (if (xaxis) {#added if (xaxis) {...} to all 3 places needed to plot x if xaxis=TRUE
            do.call(element_line, theme_get()$axis.line)
          } else {
            element_blank()
          }),
          axis.line.y = element_blank()
        )


    }

  }

  if (printspp) print(as.data.frame(sppmax_ny))

  if (z1) {
    if (z2) {
      if (missing(rel_heights)) {
        rel_heights <- prop.table(c(n_filter_decr, n_filter_incr))
        if (min(rel_heights) < .2) {
          rel_heights <- if (n_filter_decr > n_filter_incr) c(.8, .2) else c(.2, .8)
        }
      }
      plot_grid(ptop, pbottom, ncol = 1, rel_heights = rel_heights, align = "v")
    } else {
      plot_grid(ptop)
    }
  } else {
    if (z2) plot_grid(pbottom)
  }
}
