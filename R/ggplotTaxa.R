#' Plots taxon-specific change points
#'
#' Creates a plot of taxon-specific change points with optional quantiles
#' conveying uncertainty resulting from bootstrapped samples and optional
#' filtering by pure and reliable taxa.
#'
#' The fuction assumes that TITAN objects contain bootstrap summaries and
#' filtering information and automatically determines whether this is the case.
#' Without bootstrap summaries, only observed change-point locations and z-score
#' magnitudes will be plotted.  The plotting function automatically interprets
#' whether observed change-point values were obtained using IndVal or z-score
#' maxima.  The interval option is for turning off the intervals for TITAN
#' objects that contain bootstrap information. The prob95 is recommended for
#' communicating uncertainty involving management or policy action, whereas the
#' z.med option is recommended for increasingly robust estimates (by
#' incorporating uncertainty associated with the sample) of taxon-specific
#' change-point locations beyond those provided by the default (i.e., observed
#' values).
#'
#' @param titan.out A TITAN output object.
#' @param z1 A logical specifying whether decreasing taxa should be plotted.
#' @param z2 A logical specifying whether decreasing taxa should be plotted.
#' @param interval A logical specifying whether quantiles of bootstrapped change
#'   points should be plotted.
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
#'
#' plotTaxa(glades.titan, xlabel = "Surface Water TP (ug/l)", log = "x")
#' ggplotTaxa(glades.titan, xlabel = "Surface Water TP (ug/l)") + ggplot2::scale_x_log10()
#'
ggplotTaxa <- function(titan.out, z1 = TRUE, z2 = TRUE, interval = TRUE, prob95 = FALSE,
  z.med = FALSE, xlabel = "Environmental Gradient", all = FALSE, ...) {

  filter <- NULL; rm(filter)
  maxgrp <- NULL; rm(maxgrp)
  obsiv.prob <- NULL; rm(obsiv.prob)
  p05 <- NULL; rm(p05)
  p95 <- NULL; rm(p95)
  y <- NULL; rm(y)
  id <- NULL; rm(id)
  zscore <- NULL; rm(zscore)
  zenv.cp <- NULL; rm(zenv.cp)

  imax <- titan.out$arguments[[5]]
  boot <- titan.out$arguments[[3]] > 0.5
  if (all) boot <- FALSE

  ## SUBSET TAXA INTO 2 MATRICES BY GROUP ID (1 OR 2)
  sppmax <- as_data_frame(titan.out$sppmax) %>%
    mutate(id = row.names(titan.out$sppmax))

  # sppmax %>%
  #   filter(filter != 0) %>%
  #   ggplot(aes(id, ienv.cp)) +
  #     geom_linerange(aes(x = reorder(id, `50%`), ymin = `5%`, ymax = `95%`, color = factor(filter))) +
  #     geom_point(aes(x = reorder(id, `50%`), y = `50%`, color = factor(filter), size = z.median)) +
  #     coord_flip() + theme_bw() +
  #     xlab("") +
  #     ylab("Surface Water TP (ug/L)")

  if (boot) {
    if (z1) spp1 <- sppmax %>% filter(filter == 1)
    if (z2) spp2 <- sppmax %>% filter(filter == 2)
  } else {
    if (z1) spp1 <- sppmax %>% filter(maxgrp == 1, obsiv.prob <= 0.05)
    if (z2) spp2 <- sppmax %>% filter(maxgrp == 2, obsiv.prob <= 0.05)
  }

  ## Check length of spp1 and spp2
  if (z1 && nrow(spp1) < 1) stop("z1 is empty, set z1=FALSE, change significance criteria, or\nset boot=FALSE if bootstrapping was not used to generate the titan object")
  if (z2 && nrow(spp2) < 1) stop("z2 is empty, set z2=FALSE, change significance criteria, or\nset boot=FALSE if bootstrapping was not used to generate the titan object")

  ## ADD DUMMY GRAPH TO SET PLOT LIMITS
  sppsub.gt <- if (z1 && z2) {
    if (nrow(spp1) >= nrow(spp2)) spp1 else spp2
  } else {
    if (z1) spp1 else spp2
  }

  rank_first <- function(x) rank(x, ties.method = "first")
  descending_rank <- function(x) max(rank_first(x)) - rank_first(x) + 1

  ## DETERMINE RANK ORDER OF SYMBOLS ON Y AXIS
  if (boot) {
    if (prob95) {
      if (z1) yvals1 <- descending_rank(spp1$"95%")
      if (z2) yvals2 <- rank_first(spp2$"5%") + 0.5
    } else {
      if (z.med) {
        if (z1) yvals1 <- descending_rank(spp1$"50%")
        if (z2) yvals2 <- rank_first(spp2$"50%") + 0.5
      } else {
        if (imax) {
          if (z1) yvals1 <- descending_rank(spp1$ienv.cp)
          if (z2) yvals2 <- rank_first(spp2$ienv.cp) + 0.5
        } else {
          if (z1) yvals1 <- descending_rank(spp1$zenv.cp)
          if (z2) yvals2 <- rank_first(spp2$zenv.cp) + 0.5
        }
      }
    }
  } else {
    if (imax) {
      if (z1) yvals1 <- descending_rank(spp1$ienv.cp)
      if (z2) yvals2 <- rank_first(spp2$ienv.cp) + 0.5
    } else {
      if (z1) yvals1 <- descending_rank(spp1$zenv.cp)
      if (z2) yvals2 <- rank_first(spp2$zenv.cp) + 0.5
    }
  }


  spp1$y <- yvals1
  spp2$y <- yvals2
  spp <- rbind(spp1, spp2) # bind_rows(spp1, spp2)

  names(spp) <- c("ienv.cp", "zenv.cp", "freq", "maxgrp", "IndVal", "obsiv.prob", "zscore", "p05", "p10", "p50", "p90", "p95", "purity", "reliability", "z.median", "filter", "id", "y")

  # ybreaks <- spp$y
  # ylabels <- spp$id

  ggplot(spp) +
    geom_segment(aes(
        x = p05, xend = p95,
        y = y, yend = y,
        color = factor(filter), size = zscore, alpha = abs(zscore)#, fill = factor(filter)
    ), lineend = "round") +
    # geom_point(aes(
    #     x = zenv.cp, y = y,
    #     color = factor(filter), size = zscore, alpha = .5
    # )) +
    # geom_segment(aes(x = zenv.cp, xend = zenv.cp, y = y-.3, yend = y+.3), size = .25) +
    geom_point(aes(x = zenv.cp, y = y), size = .5) +
    geom_text(
      aes(x = (p05+p95)/2, y = y, label = paste0(id," (", round(zscore), ")")),
      # aes(x = (p05+p95)/2, y = y, label = glue("{id}({round(zscore)})")),
      size = 2
    ) +
    theme_bw() +
    scale_alpha("Magnitude\nof z-score", guide = FALSE) +
    scale_size(guide = FALSE) +
    # scale_y_continuous(
    #   breaks = ybreaks,
    #   labels = ylabels
    # ) +
    scale_y_continuous("", breaks = NULL) +
    xlab("Surface Water TP (ug/L)") +
    scale_color_discrete("", labels = c("1" = "z-", "2" = "z+")) +
    guides(colour = guide_legend(override.aes = list(size = 3)))

}








