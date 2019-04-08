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
#' @param xlim x axis limits, e.g. \code{xlim =c(0,10)}
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
#' @seealso [plot_sumz()], [plotCPs()]
#' @keywords TITAN kwd2
#' @export
#' @examples
#'
#' data(glades.titan)
#'
#' plot_taxa_ridges(glades.titan, ytxt.sz=8)
#' plot_taxa_ridges(glades.titan, ytxt.sz=8, grid = FALSE)
#' plot_taxa_ridges(glades.titan, ytxt.sz=8, xaxis = TRUE)
#'
#' plot_taxa_ridges(glades.titan,
#'   xlabel = expression(paste("Surface water total phosphorus ("*mu*"g/l)")),
#'   z2 = FALSE
#' )
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



  imax <- titan.out$arguments[["imax"]]
  boot <- titan.out$arguments[["boot"]] > 0.5

  sppmax <- titan.out %>%
    pluck("sppmax") %>%
    tibble::as_tibble() %>%
    mutate(id = row.names(titan.out$sppmax))

  sppmax$filter <- 0
  sppmax$filter[which(sppmax$purity>=pur.cut & sppmax$reliability>=rel.cut)] <- sppmax$maxgrp[which(sppmax$purity>=pur.cut & sppmax$reliability>=rel.cut)]


  num.dcr <- sum(sppmax$filter==1)
  num.ncr <- sum(sppmax$filter==2)

  message(paste("There are", sum(num.dcr,num.ncr), "indicator taxa of", n_ytaxa, "possible for plotting",sep=" "))
  message(paste("Number of Decreasers=", num.dcr,sep=""))
  message(paste("Number of Increasers=", num.ncr,sep=""))

  sppmax <- sppmax %>%
   dplyr::arrange(desc(purity),desc(reliability),desc(z.median)) %>%
   dplyr::slice(1:min(dplyr::n(),n_ytaxa))

  reliable_taxa_ndcs <- which(sppmax$filter > 0)

  gdf <- reliable_taxa_ndcs %>%
    map(~
      cbind(
        data_frame(chk_pts = round(sort(titan.out$metricArray[.x,2,]), digits = 2)),
        sppmax %>% slice(.x)
      )
    ) %>%
    bind_rows() %>%
    mutate(filter = case_when(
      (filter == 1) ~ -1L,
      (filter == 2) ~ +1L
    ))

  # compute pooled statistics
  if (missing(xlim)) xlim <- grDevices::extendrange(range(gdf$chk_pts), f = 0.05)
  xmin <- xlim[1]
  xmax <- xlim[2]
  if(missing(bw)) bw <- stats::bw.nrd0(gdf$chk_pts)


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
        from = xmin, to = xmax, bandwidth = bw,
        ...
      ) +
      scale_x_continuous(limits = xlim, expand = c(0,0)) +
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
      from = xmin, to = xmax, bandwidth = bw,
      ...
    ) +
    ggridges::theme_ridges(center_axis_labels = TRUE, grid = grid) +
    theme(
      axis.text.y = element_text(size=ytxt.sz),
      axis.line = (if (xaxis) {
        do.call(element_line, theme_get()$axis.line)
      } else {
        element_blank()
      }),
      axis.line.y = element_blank()
    ) +
    scale_x_continuous(xlabel, limits = xlim, expand = c(0,0)) +
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
          from = xmin, to = xmax, bandwidth = bw,
          ...
        ) +
        scale_x_continuous(xlabel, limits = xlim, expand = c(0,0)) +
        scale_y_discrete("") +
        scale_fill_gradient("Z-Score", low = z1_fill_low, high = z1_fill_high) +
        theme(axis.text.y = element_text(size = ytxt.sz))

      }

    }

  if (printspp) print(as.data.frame(sppmax))

  if (z1) {
    if (z2) {
      plot_grid(ptop, pbottom, ncol = 1, rel_heights = c(num.dcr, num.ncr), align = "v")
    } else {
      plot_grid(ptop)
    }
  } else {
    if (z2) plot_grid(pbottom)
  }
}
