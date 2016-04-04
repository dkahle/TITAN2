#' Plots taxon-specific change points
#'
#' Creates a plot of taxon-specific change points with optional
#' quantiles conveying uncertainty resulting from bootstrapped
#' samples and optional filtering by pure and reliable taxa.
#'
#' The fuction assumes that TITAN objects contain bootstrap
#' summaries and filtering information and automatically determines
#' whether this is the case.  Without bootstrap summaries, only
#' observed change-point locations and z-score magnitudes will be
#' plotted.  The plotting function automatically interprets whether
#' observed change-point values were obtained using IndVal or
#' z-score maxima.  The interval option is for turning off the
#' intervals for TITAN objects that contain bootstrap information.
#' The prob95 is recommended for communicating uncertainty involving
#' management or policy action, whereas the z.med option is
#' recommended for increasingly robust estimates (by incorporating
#' uncertainty associated with the sample) of taxon-specific
#' change-point locations beyond those provided by the default
#' (i.e., observed values).
#'
#' @param titan.out A TITAN output object.
#' @param z1 A logical specifying whether decreasing taxa should be
#'   plotted.
#' @param z2 A logical specifying whether decreasing taxa should be
#'   plotted.
#' @param interval A logical specifying whether quantiles of
#'   bootstrapped change points should be plotted.
#' @param prob95 A logical specifying whether change-point locations
#'   should be plotted on the basis of their 5th (for increasers)
#'   and 95th (for decreasers) quantile versus their observed
#'   values.
#' @param z.med A logical specifying whether (1) change point
#'   magnitudes should be obtained from the median z score across
#'   bootstrap replicates and (2) whether the locations should be
#'   plotted on the basis of the 50th quantile of change-point
#'   locations (i.e., if prob95=FALSE).
#' @param xlabel A character string for the x axis label.
#' @param log A graphical argument specifying whether an axis should
#'   be log scaled.
#' @param at A graphical argument controlling the locatino of the x
#'   axis label.
#' @param xmin A graphical argument specifying the value of the x
#'   axis minimum.
#' @param xmax A graphical argument specifying the value of the x
#'   axis maximum.
#' @param tck A graphical argument specifying the scale of axis tick
#'   marks.
#' @param bty A graphical argument specyfying the box type around
#'   the plot.
#' @param ntick A graphical argument specifying the default number
#'   of axis tick marks.
#' @param prtty A logical specifying whether pretty() should be used
#'   for axis plotting.
#' @param dig A numeric argument controlling the number of decimal
#'   digits in axis labels.
#' @param leg.x A graphical argument specifying the x coordinate of
#'   the legend.
#' @param leg.y A graphical argument specifying the y coordinate of
#'   the legend.
#' @param cex.taxa A graphical argument specifying the scaling of
#'   the taxa names.
#' @param cex A graphical argument specifying the scaling of the
#'   figure.
#' @param cex.axis A graphical argument specifying the scaling of
#'   the axes.
#' @param cex.leg A graphical argument specifying the scaling of the
#'   legend.
#' @param cex.lab A graphical argument specifying the scaling of the
#'   lables.
#' @param legend A logical specifying whether or not to plot the
#'   legend.
#' @param col1 A graphical argument specifying the color of group 1
#'   symbols.
#' @param fil1 A graphical argument specifying the color of group 1
#'   fills.
#' @param col2 A graphical argument specifying the color of group 2
#'   symbols.
#' @param fil2 A graphical argument specifying the color of group 2
#'   fills.
#' @param write A logical specifying whether summary tables are
#'   written to screen.
#' @param all A logical specifying whether all taxa with p<0.05
#'   should be plotted.
#' @param ... An argument for passing generic plotting function
#'   parameters.
#' @return A plot of decreasing and/or increasing taxon-specific
#'   change points along the environmental gradient.
#' @references Baker, ME and RS King.  2010. A new method for
#'   detecting and interpreting biodiversity and ecological
#'   community thresholds.  Methods in Ecology and Evolution 1(1):
#'   25:37.
#' @references King, RS and ME Baker  2010. Considerations for
#'   identifying and interpreting ecological community thresholds.
#'   Journal of the North American Benthological Association
#'   29(3):998-1008.
#' @note Should not be used with output objects from TITAN v1.0.
#' @author M. Baker and R. King
#' @seealso \code{\link{plotSumz}}, \code{\link{plotCPs}}
#' @keywords TITAN kwd2
#' @export
#' @examples
#'
#' data(glades.titan)
#' plotTaxa(glades.titan, xlabel = "Surface Water TP (ug/l)")
#' plotTaxa(glades.titan, xlabel = "Surface Water TP (ug/l)", log = "x")
#' ggplotTaxa(glades.titan, xlabel = "Surface Water TP (ug/l)")
#' ggplotTaxa(glades.titan, xlabel = "Surface Water TP (ug/l)") + scale_x_log10()
#'
ggplotTaxa <- function(titan.out, z1 = TRUE, z2 = TRUE, interval = TRUE, prob95 = FALSE,
  z.med = FALSE, xlabel = "Environmental Gradient", log = "", at = NULL,
  xmin = min(titan.out$sppmax[,"5%"]),
  xmax = 1.05*max(titan.out$sppmax[,"95%"]),
  tck = 0.025, bty = "u", ntick = 6, prtty = TRUE,
  dig = 1, leg.x = 0.8, leg.y = 10, cex.taxa = 0.75, cex = 1.25,
  cex.axis = 1.25, cex.leg = 1.25, cex.lab = 1.25, legend = TRUE,
  col1 = "black", fil1 = "black", col2 = "black", fil2 = "white",
  write = FALSE, all = FALSE, ...) {

  imax <- titan.out$arguments[[5]]
  boot <- titan.out$arguments[[3]] > 0.5
  if (all) boot <- FALSE

  ## SUBSET TAXA INTO 2 MATRICES BY GROUP ID (1 OR 2)
  sppmax <- as.data.frame(titan.out$sppmax)
  sppmax$id <- row.names(titan.out$sppmax)

  # names(sppmax) <- c("ienv.cp", "zenv.cp", "freq", "maxgrp",
  #   "IndVal", "obsiv.prob", "zscore",
  #   "p05", "p10", "p50", "p90", "p95",
  #   "purity", "reliability", "z.median", "filter", "id"
  # )
  #
  # sppmax %>%
  #   filter(filter != 0) %>%
  #   ggplot(aes(id, ienv.cp)) +
  #     geom_linerange(aes(x = reorder(id, p50), ymin = p05, ymax = p95, color = factor(filter))) +
  #     geom_point(aes(x = reorder(id, p50), y = p50, color = factor(filter), size = z.median)) +
  #     coord_flip() + theme_bw() +
  #     xlab("") +
  #     ylab("Surface Water TP (ug/L)")


  if (boot) {
    if (z1) spp1 <- subset(sppmax, filter == 1)
    if (z2) spp2 <- subset(sppmax, filter == 2)
  } else {
    if (z1) spp1 <- subset(sppmax, maxgrp == 1 & obsiv.prob <= 0.05)
    if (z2) spp2 <- subset(sppmax, maxgrp == 2 & obsiv.prob <= 0.05)
  }

  ## Check length of spp1 and spp2
  if (z1 && nrow(spp1) < 1) stop(
    "z1 is empty, set z1=FALSE, ",
    "change significance criteria, ",
    "or set boot=FALSE if bootstrapping was not used ",
    "to generate the titan object"
  )
  if (z2 && nrow(spp2) < 1) stop(
    "z2 is empty, set z2=FALSE, ",
    "change significance criteria, ",
    "or set boot=FALSE if bootstrapping was not used ",
    "to generate the titan object"
  )


  ## SET GRAPH WINDOW, LEAVE MARGIN FOR TAXA LABELS
  # par(mar = c(8, 8, 1, 8), oma = c(0, 3, 0, 3))

  ## ADD DUMMY GRAPH TO SET PLOT LIMITS
  if (z1 && z2) {
    if (nrow(spp1) >= nrow(spp2)) {
      sppsub.gt <- spp1
    } else {
      sppsub.gt <- spp2
    }
  } else {
    if (z1) {
      sppsub.gt <- spp1
    } else {
      sppsub.gt <- spp2
    }
  }

  rank_first <- function(x) rank(x, ties.method = "first")

  descending_rank <- function(x){
    ranks <- rank_first(x)
    max(ranks) - ranks + 1
  }

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

  # if (boot) {
  #   if (prob95) {
  #     if (z.med) {
  #       if (z1) spp1$point_size <- spp1$z.median
  #       if (z2) spp2$point_size <- spp2$z.median
  #     } else {
  #       if (z1) symbols(spp1$"95%", yvals1, circles = spp1$zscore, inches = 0.1, add = TRUE, xlim = c(0,5), fg = col1, bg =  grpcol, lwd = 2)
  #       if (z2) symbols(spp2$"5%", yvals2, circles = spp2$zscore, inches = 0.1, add = TRUE, xlim = c(0,5), fg = col2, bg = grpcol2, lwd = 2)
  #     }
  #   } else {
  #     if (z.med) {
  #       if (z1) symbols(spp1$"50%", yvals1, circles = spp1$z.median, inches = 0.1, add = TRUE, xlim = c(0,5), fg = col1, bg =  grpcol, lwd = 2)
  #       if (z2) symbols(spp2$"50%", yvals2, circles = spp2$z.median, inches = 0.1, add = TRUE, xlim = c(0,5), fg = col2, bg = grpcol2, lwd = 2)
  #     } else {
  #       if (imax) {
  #         if (z1) symbols(spp1$ienv.cp, yvals1, circles = spp1$zscore, inches = 0.1, add = TRUE, xlim = c(0,5), fg = col1, bg =  grpcol, lwd = 2)
  #         if (z2) symbols(spp2$ienv.cp, yvals2, circles = spp2$zscore, inches = 0.1, add = TRUE, xlim = c(0,5), fg = col2, bg = grpcol2, lwd = 2)
  #       } else {
  #         if (z1) symbols(spp1$zenv.cp, yvals1, circles = spp1$zscore, inches = 0.1, add = TRUE, xlim = c(0,5), fg = col1, bg = grpcol, lwd = 2)
  #         if (z2) symbols(spp2$zenv.cp, yvals2, circles = spp2$zscore, inches = 0.1, add = TRUE, xlim = c(0,5), fg = col2, bg = grpcol2, lwd = 2)
  #       }
  #     }
  #   }
  # } else {
  #   if (imax) {
  #     if (z1) symbols(spp1$ienv.cp, yvals1, circles = spp1$zscore, inches = 0.1, add = TRUE, xlim = c(0, 5), fg = col1, bg = grpcol, lwd = 2)
  #     if (z2) symbols(spp2$ienv.cp, yvals2, circles = spp2$zscore, inches = 0.1, add = TRUE, xlim = c(0, 5), fg = col2, bg = grpcol2, lwd = 2)
  #   } else {
  #     if (z1) symbols(spp1$zenv.cp, yvals1, circles = spp1$zscore, inches = 0.1, add = TRUE, xlim = c(0, 5), fg = col1, bg = grpcol, lwd = 2)
  #     if (z2) symbols(spp2$zenv.cp, yvals2, circles = spp2$zscore, inches = 0.1, add = TRUE, xlim = c(0, 5), fg = col2, bg = grpcol2, lwd = 2)
  #   }
  # }

  names(spp) <- c("ienv.cp", "zenv.cp", "freq", "maxgrp",
    "IndVal", "obsiv.prob", "zscore",
    "p05", "p10", "p50", "p90", "p95",
    "purity", "reliability", "z.median", "filter", "id", "y"
  )

  # ybreaks <- spp$y
  # ylabels <- spp$id

  p <- ggplot(spp) +
    geom_segment(aes(
        x = p05, xend = p95,
        y = y, yend = y,
        color = factor(filter), size = zscore, alpha = abs(zscore)
    ), lineend = "round") +
    # geom_point(aes(
    #     x = zenv.cp, y = y,
    #     color = factor(filter), size = zscore, alpha = .5
    # )) +
    geom_segment(aes(x = zenv.cp, xend = zenv.cp, y = y-.3, yend = y+.3), size = .25) +
    geom_text(aes(x = (p05+p95)/2, y = y, label = id), size = 2) +
    theme_bw() +
    scale_alpha("Magnitude\nof z-score") +
    scale_size(guide = FALSE) +
    # scale_y_continuous(
    #   breaks = ybreaks,
    #   labels = ylabels
    # ) +
    scale_y_continuous("", breaks = NULL) +
    xlab("Surface Water TP (ug/L)") +
    scale_color_discrete("", labels = c("1" = "z-", "2" = "z+"))

  return(p)











  ## ADD INTERVALS AS LINE SEGMENTS
  if (boot && interval) {
    if (z1) segments(spp1$"5%", yvals1, spp1$"95%", yvals1, col = col1, lwd = 2)
    if (z2) segments(spp2$"5%", yvals2, spp2$"95%", yvals2, col = col2, lwd = 2, lty = 3)
  }


  ## CREATE VECTOR FOR COLOR CODING MAX GROUP ASSIGNMENTS

  if (z1) {
    grpcol <- rep(fil1, nrow(spp1))
    grpcol[spp1$maxgrp > 1.5] <- col2
  }

  if (z2) {
    grpcol2 <- rep(fil1, nrow(spp2))
    grpcol2[spp2$maxgrp > 1.5] <- fil2
  }



  ## ADD TAXA SYMBOLS IN PROPORTION TO Z SCORE OR MEDIAN BOOT.Z,
  ## COLORED BY MAX GROUP ASSIGMENT
  if (boot) {
    if (prob95) {
      if (z.med) {
        if (z1) symbols(spp1$"95%", yvals1, circles = spp1$z.median, inches = 0.1, add = TRUE, xlim = c(0,5), fg = col1, bg =  grpcol, lwd = 2)
        if (z2) symbols(spp2$"5%", yvals2, circles = spp2$z.median, inches = 0.1, add = TRUE, xlim = c(0,5), fg = col2, bg = grpcol2, lwd = 2)
      } else {
        if (z1) symbols(spp1$"95%", yvals1, circles = spp1$zscore, inches = 0.1, add = TRUE, xlim = c(0,5), fg = col1, bg =  grpcol, lwd = 2)
        if (z2) symbols(spp2$"5%", yvals2, circles = spp2$zscore, inches = 0.1, add = TRUE, xlim = c(0,5), fg = col2, bg = grpcol2, lwd = 2)
      }
    } else {
      if (z.med) {
        if (z1) symbols(spp1$"50%", yvals1, circles = spp1$z.median, inches = 0.1, add = TRUE, xlim = c(0,5), fg = col1, bg =  grpcol, lwd = 2)
        if (z2) symbols(spp2$"50%", yvals2, circles = spp2$z.median, inches = 0.1, add = TRUE, xlim = c(0,5), fg = col2, bg = grpcol2, lwd = 2)
      } else {
        if (imax) {
          if (z1) symbols(spp1$ienv.cp, yvals1, circles = spp1$zscore, inches = 0.1, add = TRUE, xlim = c(0,5), fg = col1, bg =  grpcol, lwd = 2)
          if (z2) symbols(spp2$ienv.cp, yvals2, circles = spp2$zscore, inches = 0.1, add = TRUE, xlim = c(0,5), fg = col2, bg = grpcol2, lwd = 2)
        } else {
          if (z1) symbols(spp1$zenv.cp, yvals1, circles = spp1$zscore, inches = 0.1, add = TRUE, xlim = c(0,5), fg = col1, bg = grpcol, lwd = 2)
          if (z2) symbols(spp2$zenv.cp, yvals2, circles = spp2$zscore, inches = 0.1, add = TRUE, xlim = c(0,5), fg = col2, bg = grpcol2, lwd = 2)
        }
      }
    }
  } else {
    if (imax) {
      if (z1) symbols(spp1$ienv.cp, yvals1, circles = spp1$zscore, inches = 0.1, add = TRUE, xlim = c(0, 5), fg = col1, bg = grpcol, lwd = 2)
      if (z2) symbols(spp2$ienv.cp, yvals2, circles = spp2$zscore, inches = 0.1, add = TRUE, xlim = c(0, 5), fg = col2, bg = grpcol2, lwd = 2)
    } else {
      if (z1) symbols(spp1$zenv.cp, yvals1, circles = spp1$zscore, inches = 0.1, add = TRUE, xlim = c(0, 5), fg = col1, bg = grpcol, lwd = 2)
      if (z2) symbols(spp2$zenv.cp, yvals2, circles = spp2$zscore, inches = 0.1, add = TRUE, xlim = c(0, 5), fg = col2, bg = grpcol2, lwd = 2)
    }
  }

  ## ADD TAXA NAMES TO Y-AXES
  # if (z1) axis(2, at = yvals1, labels = rownames(spp1), las = 1, mgp = c(1, 0.5, 0), cex.axis = cex.taxa, tck = tck)
  # if (z2) axis(4, at = yvals2, labels = rownames(spp2), las = 1, mgp = c(1, 0.5, 0), cex.axis = cex.taxa, tck = tck)


  ## ADD X-AXIS WITH SAME RANGE AS DUMMY PLOT
  if (log == "x") {
    axis(1, at = at, mgp = c(2, 0.5, 0), cex.axis = cex.axis, tck = tck)
  } else {
    if (prtty) {
      axis(1, pretty(xmin:xmax, ntick), mgp = c(2, 0.5, 0), cex.axis = cex.axis, tck = tck)
    } else {
      axis(1, at = seq(from = round(xmin, digits = dig), to = round(xmax, digits = dig), by = round((xmax - xmin)/4, digits = dig)), mgp = c(2, 0.5, 0), cex.axis = cex.axis, tck = tck)
    }
  }

  ## ADD X-AXIS LABEL
  mtext(xlabel, side = 1, line = 2, cex = cex)

  ## ADD LEGEND TO FIGURE
  if (z1 && z2) {
    leg <- c("z-", "z+")
    fill.leg <- c(fil1, fil2)
    legend(titan.out$envcls[length(titan.out$envcls) * leg.x], leg.y, leg, fill = fill.leg, ncol = 1, bty = "n", plot = TRUE, cex = cex.leg)
  }
  box(which = "plot", bty = bty)

  ## WRITE spp1 AND/OR 2 TO FILE AND PRINT TO CONSOLE if(z1)
  ## {write.table(spp1,'spp1.txt')} if(z2)
  ## {write.table(spp2,'spp2.txt')}
  if (z1 & z2 & write) {
    sppsub <- list(spp1, spp2)
    names(sppsub) <- c("spp1", "spp2")
    return(sppsub)
  }
  if (z1 & write) return(spp1)
  if (z2 & write) return(spp2)



}








