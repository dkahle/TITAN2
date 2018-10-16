#' Taxon change point ridge plots
#'
#' Taxon change point ridge plots
#'
#' @param titan.out A TITAN output object.
#' @param z1 A logical specifying whether decreasing taxa should be plotted.
#' @param z2 A logical specifying whether decreasing taxa should be plotted.
#' @param prob95 A logical specifying whether change-point locations should be
#'   plotted on the basis of their 5th (for increasers) and 95th (for
#'   decreasers) quantile versus their observed values.
#' @param z.med A logical specifying whether (1) change point magnitudes should
#'   be obtained from the median z score across bootstrap replicates and (2)
#'   whether the locations should be plotted on the basis of the 50th quantile
#'   of change-point locations (i.e., if prob95=FALSE).
#' @param xlabel A character string for the x axis label.
#' @param all A logical specifying whether all taxa with p<0.05 should be
#'   plotted.
#' @param xlim X axis limits.
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
#' @author M. Baker and R. King
#' @seealso \code{\link{plotSumz}}, \code{\link{plotCPs}}
#' @keywords TITAN kwd2
#' @export
#' @examples
#'
#' data(glades.titan)
#'
#' plotTaxa(glades.titan, xlabel = "Surface Water TP (ug/l)")
#' ggplotTaxa(glades.titan, xlabel = "Surface Water TP (ug/l)")
#' taxa_ridges(glades.titan, xlabel = "Surface Water TP (ug/l)")
#' taxa_ridges(glades.titan, xlabel = "Surface Water TP (ug/l)", xlim = c(0, 125))
#'
taxa_ridges <- function(titan.out, z1 = TRUE, z2 = TRUE, prob95 = FALSE,
  z.med = FALSE, xlabel = "Environmental Gradient", all = FALSE, xlim, ...) {

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

  imax <- titan.out$arguments[["imax"]]
  boot <- titan.out$arguments[["boot"]] > 0.5
  if (all) boot <- FALSE

  sppmax <- titan.out %>%
    pluck("sppmax") %>%
    as_data_frame() %>%
    mutate(id = row.names(titan.out$sppmax))

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
      (filter == 1) ~ -1,
      (filter == 2) ~ +1
    ))

  # nbd = rev(dim(titan.out$metricArray))[1],
  # taxa = sppmax$id[.x],
  # z.Scores =

  if (missing(xlim)) xlim <- range(gdf$chk_pts)

  gdf %>%
    filter(filter == -1) %>%
    ggplot(aes(x = chk_pts, y = reorder(id, -chk_pts, median), fill = zscore)) +
      geom_density_ridges(
        scale = 2, alpha = 0.6, size = 0.05, color="gray50",
        quantile_lines = TRUE, quantiles = 2, vline_color = "black", vline_size = .25
      ) +
      scale_x_continuous("", limits = xlim) +
      scale_y_discrete("Taxa") +
      scale_fill_gradient("Z-Score", low = "light blue", high = "blue") +
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        plot.margin = structure(c(5.5, 5.5, 0, 5.5), class = c("margin", "unit"), valid.unit = 8L, unit = "pt")
      ) ->
    ptop

  gdf %>%
    filter(filter == +1) %>%
    ggplot(aes(x = chk_pts, y = reorder(id, chk_pts, median), fill = zscore)) +
      geom_density_ridges(
        scale = 2, alpha = 0.6, size = 0.05, color="gray50",
        quantile_lines = TRUE, quantiles = 2, vline_color = "black", vline_size = .25
      ) +
      scale_x_continuous(xlabel,  limits = xlim) +
      scale_y_discrete("Taxa") +
      scale_fill_gradient("Z-Score", low = "pink", high = "red") ->
    pbottom


  plot_grid(ptop, pbottom, ncol = 1, align = "v")

}



