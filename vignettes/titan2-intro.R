## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo=TRUE, collapse=TRUE, error=TRUE, comment = "#")

## -----------------------------------------------------------------------------
library("TITAN2")

## -----------------------------------------------------------------------------
data(glades.taxa)
str(glades.taxa, list.len = 5)

## -----------------------------------------------------------------------------
data(glades.env)
str(glades.env)

## ---- eval = FALSE------------------------------------------------------------
#  glades.titan <- titan(glades.env, glades.taxa)

## ---- eval = FALSE------------------------------------------------------------
#  glades.titan <- titan(glades.env, glades.taxa,
#    minSplt = 5, numPerm = 250, boot = TRUE, nBoot = 500, imax = FALSE,
#    ivTot = FALSE, pur.cut = 0.95, rel.cut = 0.95, ncpus = 8, memory = FALSE
#  )

## -----------------------------------------------------------------------------
data(glades.titan)
str(glades.titan, 1)

## ---- echo = FALSE------------------------------------------------------------
message("100% occurrence detected 1 times (0.8% of taxa), use of TITAN less than ideal for this data type")
message("Taxa frequency screen complete")

## ---- echo = FALSE------------------------------------------------------------
message("Determining partitions along gradient")
message("Calculating observed IndVal maxima and class values")
message("Calculating IndVals using mean relative abundance")
message("Permuting IndVal scores")
message("IndVal $z$ score calculation complete")
message("Summarizing Observed Results")
message("Estimating taxa change points using z-score maxima")

## ---- echo = FALSE------------------------------------------------------------
message("Bootstrap resampling in sequence...")
message(1*1)
message(2*1)
message(3*1)

## ---- echo = FALSE------------------------------------------------------------
message("Bootstrap resampling in parallel using 2 CPUs...no index will be printed to screen")

## -----------------------------------------------------------------------------
glades.titan$sumz.cp

## -----------------------------------------------------------------------------
head(glades.titan$sppmax)

## -----------------------------------------------------------------------------
str(glades.titan, max.level = 1, give.attr = FALSE)

## ---- fig.height = 6,fig.width = 8--------------------------------------------
plot_sumz_density(glades.titan)

## ---- fig.height = 6,fig.width = 8--------------------------------------------
plot_sumz_density(glades.titan, ribbon = FALSE, points = TRUE)

## ---- fig.height = 6,fig.width = 8--------------------------------------------
plot_sumz_density(glades.titan, 
  ribbon = TRUE, points = FALSE, sumz1 = FALSE, change_points = FALSE, 
  xlabel = expression(paste("Surface Water Total Phosphorus ("*mu*"g/l)"))
)

## ---- fig.height = 6,fig.width = 8--------------------------------------------
plot_sumz(glades.titan, filter = TRUE)

## ---- fig.height = 10,fig.width = 8-------------------------------------------
plot_taxa_ridges(glades.titan, axis.text.y = 8)

## ---- fig.height = 10,fig.width = 8-------------------------------------------
plot_taxa_ridges(glades.titan, 
  xlabel = expression(paste("Surface water total phosphorus ("*mu*"g/l)")), 
  n_ytaxa = 50
)

## ---- fig.height = 6,fig.width = 8--------------------------------------------
plot_taxa_ridges(glades.titan, 
  xlabel = expression(paste("Surface water total phosphorus ("*mu*"g/l)")), 
  z2 = FALSE
)

## ---- fig.height = 6,fig.width = 8--------------------------------------------
plot_taxa_ridges(glades.titan, 
  xlabel = expression(paste("Surface water total phosphorus ("*mu*"g/l)")), 
  z2 = FALSE, grid = FALSE
)

## ---- fig.height = 10,fig.width = 8-------------------------------------------
plot_taxa_ridges(
  glades.titan, 
  axis.text.x = 12, axis.text.y = 8, axis.title.x = 14,  
  rel_heights = c(0.45, 0.55), xaxis = TRUE, d_lines = TRUE, trans = "log10", 
  xlim = c(4, 200), breaks = c( 10, 20, 40, 80, 160), 
  xlabel = expression(paste("Surface water total phosphorus ("*mu*"g/l)"))
)

## ---- fig.height = 8,fig.width = 9--------------------------------------------
plot_taxa(glades.titan, xlabel = "Surface Water TP (ug/l)")

## ---- fig.height = 8,fig.width = 9--------------------------------------------
plot_taxa(
  glades.titan, z.med = FALSE,
  xlabel = expression(paste("Surface water total phosphorus ("*mu*"g/l)"))
)

## ---- fig.height = 8,fig.width = 9--------------------------------------------
plot_taxa(
  glades.titan, 
  xlabel = expression(paste("Surface water total phosphorus ("*mu*"g/l)")), 
  cex.taxa = 0.5, cex = 1.25, cex.axis = 1.1, legend = FALSE, 
  col1 = "black", col2 = "red", fil2 = "red"
)

## ---- fig.height = 8,fig.width = 9--------------------------------------------
plot_taxa(
  glades.titan, 
  xlabel = expression(paste("Surface water total phosphorus ("*mu*"g/l)")), 
  cex.taxa = 0.5, cex = 1.25, cex.axis = 1.1, legend = FALSE,  
  prob95 = TRUE, col1 = "black", col2 = "red", fil2 = "red"
)

## ---- fig.height = 10, fig.width = 10-----------------------------------------
plot_cps(glades.titan)

## ---- fig.height = 5,fig.width = 8--------------------------------------------
plot_cps(glades.titan, taxaID = "ENALCIVI", xlabel = "Surface Water TP (ug/l)")

## ---- fig.height = 5,fig.width = 8--------------------------------------------
plot_cps(glades.titan, taxaID = "ENALCIVI", cp.trace = TRUE, xlabel = "Surface Water TP (ug/l)")

## ---- fig.height = 5,fig.width = 8--------------------------------------------
plot_cps(glades.titan, taxaID = "OSTRASP5", cp.trace = TRUE, xlabel = "Surface Water TP (ug/l)")

## ---- fig.height = 6,fig.width = 8--------------------------------------------
plot_cps(glades.titan, taxa.dist = FALSE, xlabel = "Surface Water TP (ug/l)")

## ---- fig.height = 6,fig.width = 8--------------------------------------------
plot_cps(glades.titan, taxa.dist = FALSE, xlabel = "Surface Water TP (ug/l)", stacked = TRUE)

