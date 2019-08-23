#' Taxon change point ridge plots
#'
#' Taxon change point ridge plots
#'
#' @param titan.out A TITAN output object.
#' @param z1 A logical specifying whether decreasing taxa should be plotted.
#' @param z2 A logical specifying whether increasing taxa should be plotted.
#' @param xlabel A character string for the x axis label.
#' @param ytxt.sz The relative size of the taxa label along the y axis.
#' @param n_ytaxa The maximum number of taxa to be plotted.
#' @param printspp A logical specifying whether the sppmax table should be
#'   printed.
#' @param z1_fill_low,z1_fill_high,z2_fill_low,z2_fill_high Respective fill
#'   colors passed to [scale_fill_gradient()]
#' @param pur.cut pur.cut
#' @param rel.cut rel.cut
#' @param grid The \code{grid} argument of [theme_ridges()].
#' @param xaxis Logical; should the x-axis be plotted?
#' @param bw The bandwidth of used in the kernel density estimate; see
#'   [density()].
#' @param xlim x axis limits, e.g. \code{xlim = c(0,10)}
#' @param trans a scale transformation to be applied to the x-axis through
#'   ggplot2. e.g. \code{"log10"} or \code{"sqrt"}.
#' @param ... ...
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
#' plot_taxa_ridges(glades.titan, ytxt.sz = 8)
#' plot_taxa_ridges(glades.titan, ytxt.sz = 8, grid = FALSE)
#' plot_taxa_ridges(glades.titan, ytxt.sz = 8, xaxis = TRUE)
#' plot_taxa_ridges(glades.titan, ytxt.sz = 8, trans = "log10", xlim = c(1,150))
#' plot_taxa_ridges(glades.titan, ytxt.sz = 8, trans = "sqrt", xlim = c(0,150))
#'
#' # example with x-axis styling
#' # plot_taxa_ridges(glades.titan,
#' #   xlabel = expression(paste("Surface water total phosphorus ("*mu*"g/l)")),
#' #   z2 = FALSE
#' # )
#'
#'
plot_taxa_ridges <- function(
  titan.out,
  z1 = TRUE, z2 = TRUE,
  z1_fill_low="light blue", z1_fill_high="black",
  z2_fill_low="pink", z2_fill_high="red",
  pur.cut=titan.out$arguments[[7]],
  rel.cut=titan.out$arguments[[8]],
  xlabel = "Environmental Gradient",
  n_ytaxa = 90,
  ytxt.sz = 10,
  printspp = FALSE,
  grid = TRUE,
  xaxis = !grid,
  trans = "identity",
  xlim, bw,
  ...
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

  imax <- titan.out$arguments[["imax"]]
  boot <- titan.out$arguments[["boot"]] > 0.5

  sppmax <- titan.out %>%
    pluck("sppmax") %>%
    tibble::as_tibble() %>%
    mutate(id = row.names(titan.out$sppmax)) #%>%
    #dplyr::rename("Q5" = "5%", "Q95" = "95%")

  sppmax$filter <- 0
  pure_and_reliable_ndcs <- with(sppmax, which(purity>=pur.cut & reliability>=rel.cut))
  sppmax$filter[pure_and_reliable_ndcs] <- sppmax$maxgrp[pure_and_reliable_ndcs]

  num.dcr <- sum(sppmax$filter == 1L)
  num.ncr <- sum(sppmax$filter == 2L)

  message(glue("There are {sum(num.dcr, num.ncr)} indicator taxa of {n_ytaxa} possible for plotting"))
  message(glue("Number of Decreasers = {num.dcr}"))
  message(glue("Number of Increasers = {num.ncr}"))

  reliable_taxa_ndcs <- which(sppmax$filter > 0)

  sppmax_ny <- sppmax %>%
   dplyr::arrange(desc(purity),desc(reliability),desc(z.median)) %>%
   dplyr::slice(1:min(dplyr::n(), n_ytaxa))

  n_filter_decr <- with(sppmax_ny, sum(filter == 1L))
  n_filter_incr <- with(sppmax_ny, sum(filter == 2L))

  gdf <- reliable_taxa_ndcs %>%
    map(~
      cbind(
        tibble::tibble(
          "chk_pts" = round(sort(titan.out$metricArray[.x,2,]), digits = 2)
        ),
        sppmax %>% slice(.x)
      )
    ) %>%
    bind_rows() %>%
    mutate(filter = case_when(
      (filter == 1) ~ -1L,
      (filter == 2) ~ +1L
    ))

  gdf <- gdf %>%
    dplyr::arrange(desc(purity), desc(reliability), desc(z.median)) %>%
    dplyr::slice(1:min(dplyr::n(), n_ytaxa*length(titan.out$metricArray[1,1,]))) %>%
    dplyr::filter(`5%` < chk_pts, chk_pts < `95%`)

  # compute pooled statistics
  if (missing(xlim)) xlim <- grDevices::extendrange(range(gdf$chk_pts), f = 0.05)
  xmin <- xlim[1]
  xmax <- xlim[2]
  if (missing(bw)) bw <- stats::bw.nrd0(do.call(trans, list(gdf$chk_pts)))


  if (z2) {

    ptop <- ggplot(
      filter(gdf, filter == -1L),
      aes(x = chk_pts, y = reorder(id, -chk_pts, median), fill = zscore)
    ) +
      ggridges::geom_density_ridges(
        quantile_lines = TRUE,
        vline_size = 0.25,
        quantiles = 2,
        vline_color = "black",
        color = NA,
        alpha = 0.6,
        from = do.call(trans, list(xmin)),
        to = do.call(trans, list(xmax)),
        bandwidth = bw,
        ...
      ) +
      scale_x_continuous(limits = xlim, expand = c(0,0), trans = trans) +
      scale_y_discrete("") +
      scale_fill_gradient("Z-Score", low = z1_fill_low, high = z1_fill_high) +
      ggridges::theme_ridges(center_axis_labels = TRUE, grid = grid) +
      theme(
        axis.text.y = element_text(size = ytxt.sz),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        plot.margin = structure(
          c(5.5, 5.5, 0, 5.5),
          class = c("margin", "unit"),
          valid.unit = 8L,
          unit = "pt"
        )
      )

    pbottom <- ggplot(
      filter(gdf, filter == +1L),
      aes(x = chk_pts, y = reorder(id, chk_pts, median), fill = zscore)
    ) +
    ggridges::geom_density_ridges(
      quantile_lines = TRUE,
      vline_size = 0.25,
      quantiles = 2,
      vline_color = "black",
      color = NA,
      alpha = 0.6,
      from = do.call(trans, list(xmin)),
      to = do.call(trans, list(xmax)),
      bandwidth = bw,
      ...
    ) +
    ggridges::theme_ridges(center_axis_labels = TRUE, grid = grid) +
    theme(
      axis.text.y = element_text(size = ytxt.sz),
      axis.line = (if (xaxis) {
        do.call(element_line, theme_get()$axis.line)
      } else {
        element_blank()
      }),
      axis.line.y = element_blank()
    ) +
    scale_x_continuous(xlabel, limits = xlim, expand = c(0,0), trans = trans) +
    scale_y_discrete("") +
    labs(x = xlabel, y = "") +
    scale_fill_gradient("Z-Score", low = z2_fill_low, high = z2_fill_high)

  } else {

    if (z1) {

      ptop <- ggplot(
        filter(gdf, filter == -1),
        aes(x = chk_pts, y = reorder(id, -chk_pts, median), fill = zscore)
      ) +
        ggridges::theme_ridges(center_axis_labels = TRUE, grid = grid) +
        ggridges::geom_density_ridges(
          quantile_lines = TRUE,
          vline_size = 0.25,
          quantiles = 2,
          vline_color = "black",
          color = NA,
          alpha = 0.6,
          from = do.call(trans, list(xmin)),
          to = do.call(trans, list(xmax)),
          bandwidth = bw,
          ...
        ) +
        scale_x_continuous(xlabel, limits = xlim, expand = c(0,0), trans = trans) +
        scale_y_discrete("") +
        scale_fill_gradient("Z-Score", low = z1_fill_low, high = z1_fill_high) +
        theme(axis.text.y = element_text(size = ytxt.sz))

      }

    }

  if (printspp) print(as.data.frame(sppmax_ny))

  if (z1) {
    if (z2) {
      plot_grid(ptop, pbottom, ncol = 1, rel_heights = c(n_filter_decr, n_filter_incr), align = "v")
    } else {
      plot_grid(ptop)
    }
  } else {
    if (z2) plot_grid(pbottom)
  }
}
