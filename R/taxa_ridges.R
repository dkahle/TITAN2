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
#' @param ytxt.sz The relative size of the taxa label along the y axis.
#' @param n_ytaxa The maximum number of taxa to be plotted.
#' @param printspp A logical specifying whether the sppmax table should be
#'   printed.
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


###REQUIRES tidyverse, cowplot, ggridges

##z1, z2 not functional
##xlim problem
##scale x log10, sqrt, or continuous

taxa_ridges <- function(titan.out, z1 = TRUE, z2 = TRUE, z1_fill_low="light blue", z1_fill_high="black", z2_fill_low="pink", z2_fill_high="red", pur.cut=titan.out$arguments[[7]], rel.cut=titan.out$arguments[[8]], xlabel = "Environmental Gradient", n_ytaxa=90, ytxt.sz=10, all = FALSE, printspp=T, xmin, xmax, ...) {

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

  sppmax$filter <- 0
  sppmax$filter[which(sppmax$purity>=pur.cut & sppmax$reliability>=rel.cut)] <- sppmax$maxgrp[which(sppmax$purity>=pur.cut & sppmax$reliability>=rel.cut)]


  num.dcr <- sum(sppmax$filter==1)
  num.ncr <- sum(sppmax$filter==2)

  message(paste("There are", sum(num.dcr,num.ncr), "indicator taxa of", n_ytaxa, "possible for plotting",sep=" "))
  message(paste("Number of Decreasers=", num.dcr,sep=""))
  message(paste("Number of Increasers=", num.ncr,sep=""))

  sppmax <- sppmax %>%
   arrange(desc(purity),desc(reliability),desc(z.median)) %>%
   slice(1:min(n(),n_ytaxa))

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

  ####if (missing(xlim)) xlim <- range(gdf$chk_pts)

  if(z2){
    gdf %>%
      filter(filter == -1) %>%
      ggplot(aes(x = chk_pts, y = reorder(id, -chk_pts, median), fill = zscore)) +
      stat_density_ridges(quantile_lines=T, vline_size=0.25, geom="density_ridges",quantiles=2, vline_color="black", size=0.05, color="gray50",alpha=0.6) +
      xlim(xmin, xmax) +
      scale_y_discrete("") +
      scale_fill_gradient("Z-Score", low = z1_fill_low, high = z1_fill_high) +
      theme(
        axis.text.y = element_text(size=ytxt.sz),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        plot.margin = structure(c(5.5, 5.5, 0, 5.5), class = c("margin", "unit"), valid.unit = 8L, unit = "pt")
      ) ->
      ptop

    gdf %>%
      filter(filter == +1) %>%
      ggplot(aes(x = chk_pts, y = reorder(id, chk_pts, median), fill = zscore)) +
      stat_density_ridges(quantile_lines=T, vline_size=0.25, geom="density_ridges",quantiles=2, vline_color="black", size=0.05, color="gray50",alpha=0.6) +
      theme(axis.text.y = element_text(size=ytxt.sz)) +
      xlim(xmin, xmax) +
      labs(x=xlabel)+
      scale_y_discrete("") +
      scale_fill_gradient("Z-Score", low = z2_fill_low, high = z2_fill_high) ->
      pbottom
  }else{
    if(z1){
      gdf %>%
        filter(filter == -1) %>%
        ggplot(aes(x = chk_pts, y = reorder(id, -chk_pts, median), fill = zscore)) +
        stat_density_ridges(quantile_lines=T, vline_size=0.25, geom="density_ridges",quantiles=2, vline_color="black", size=0.05, color="gray50",alpha=0.6) +
        xlim(xmin, xmax) +
        scale_y_discrete("") +
        scale_fill_gradient("Z-Score", low = z1_fill_low, high = z1_fill_high) +
        theme(axis.text.y = element_text(size=ytxt.sz)) +
        labs(x=xlabel) ->
        ptop
    }
    }

  if(printspp) {print(as.data.frame(sppmax))}

  if(z1){
    if(z2){
      plot_grid(ptop, pbottom, ncol = 1, rel_heights=c(num.dcr,num.ncr), align = "v")
    }else{
      plot_grid(ptop)
    } }else{
      if(z2){plot_grid(pbottom)}
    }
  }
