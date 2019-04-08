## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo=TRUE, collapse=TRUE, error=TRUE, comment = "#")

## ---- echo=FALSE, warning=FALSE, message=FALSE, comment=FALSE------------
library(TITAN2)
library(cowplot)
library(ggridges)
library(tidyverse)

## ------------------------------------------------------------------------
library("TITAN2")

## ------------------------------------------------------------------------
data(glades.taxa)
str(glades.taxa, list.len = 5)

## ------------------------------------------------------------------------
data(glades.env)
str(glades.env)

## ---- eval = F-----------------------------------------------------------
#  glades.titan <- titan(glades.env, glades.taxa)

## ---- eval = F-----------------------------------------------------------
#  glades.titan <- titan(glades.env, glades.taxa,
#    minSplt = 5, numPerm = 250, boot = TRUE, nBoot = 500, imax = FALSE,
#    ivTot = FALSE, pur.cut = 0.95, rel.cut = 0.95, ncpus = 1, memory = FALSE
#  )

## ------------------------------------------------------------------------
data(glades.titan)

## ---- echo = FALSE-------------------------------------------------------
message("100% occurrence detected 1 times (0.8% of taxa), use of TITAN less than ideal for this data type")
message("Taxa frequency screen complete")

## ---- echo = FALSE-------------------------------------------------------
message("Determining partitions along gradient")
message("Calculating observed IndVal maxima and class values")
message("Calculating IndVals using mean relative abundance")
message("Permuting IndVal scores")
message("IndVal $z$ score calculation complete")
message("Summarizing Observed Results")
message("Estimating taxa change points using z-score maxima")

## ---- echo = FALSE-------------------------------------------------------
message("Bootstrap resampling in sequence...")
message(1*1)
message(2*1)
message(3*1)

## ---- echo = FALSE-------------------------------------------------------
message("Bootstrap resampling in parallel using 2 CPUs...no index will be printed to screen")

## ------------------------------------------------------------------------
glades.titan$sumz.cp

## ------------------------------------------------------------------------
head(glades.titan$sppmax)

## ------------------------------------------------------------------------
str(glades.titan, max.level = 1, give.attr = FALSE)

## ---- echo = FALSE-------------------------------------------------------
graphics.off()

## ---- fig.height = 6,fig.width = 8---------------------------------------
plot_sumz_density(glades.titan)

## ---- echo = FALSE-------------------------------------------------------
graphics.off()

## ---- fig.height = 6,fig.width = 8---------------------------------------
plot_sumz_density(glades.titan, ribbon=FALSE, points=TRUE)

## ---- echo = FALSE-------------------------------------------------------
graphics.off()

## ---- fig.height = 6,fig.width = 8---------------------------------------
plot_sumz_density(glades.titan, ribbon=TRUE, points=FALSE, sumz1=FALSE, change_points = FALSE, xlabel=expression(paste("Surface Water Total Phosphorus ("*mu*"g/l)")))

## ---- echo = FALSE-------------------------------------------------------
graphics.off()

## ---- fig.height = 6,fig.width = 8---------------------------------------
plotSumz(glades.titan, filter = TRUE)

## ---- echo = FALSE-------------------------------------------------------
graphics.off()

## ---- fig.height = 10,fig.width = 8--------------------------------------
plot_taxa_ridges(glades.titan, ytxt.sz=8)

## ---- echo = FALSE-------------------------------------------------------
graphics.off()

## ---- fig.height = 10,fig.width = 8--------------------------------------
plot_taxa_ridges(glades.titan, xlabel = expression(paste("Surface water total phosphorus ("*mu*"g/l)")), n_ytaxa=50)

## ---- echo = FALSE-------------------------------------------------------
graphics.off()

## ---- fig.height = 6,fig.width = 8---------------------------------------
plot_taxa_ridges(glades.titan, xlabel = expression(paste("Surface water total phosphorus ("*mu*"g/l)")), z2=FALSE)

## ---- echo = FALSE-------------------------------------------------------
graphics.off()

## ---- fig.height = 6,fig.width = 8---------------------------------------
plot_taxa_ridges(glades.titan, xlabel = expression(paste("Surface water total phosphorus ("*mu*"g/l)")), z2=FALSE, grid=FALSE)

## ---- echo = FALSE-------------------------------------------------------
graphics.off()

## ---- fig.height = 8,fig.width = 8---------------------------------------
plotTaxa(glades.titan, xlabel = "Surface Water TP (ug/l)")

## ---- echo = FALSE-------------------------------------------------------
graphics.off()

## ---- fig.height = 8,fig.width = 8---------------------------------------
plotTaxa(glades.titan,xlabel = "Surface Water TP (ug/l)",z.med = T)

## ---- echo = FALSE-------------------------------------------------------
graphics.off()

## ---- fig.height = 8,fig.width = 8---------------------------------------
plotTaxa(glades.titan, xlabel = "Surface Water TP (ug/l)", z.med = F, prob95 = T)

## ---- fig.height = 10,fig.width = 10-------------------------------------
plotCPs(glades.titan)

## ---- echo = FALSE-------------------------------------------------------
graphics.off()

## ---- fig.height = 5,fig.width = 8---------------------------------------
plotCPs(glades.titan, taxaID = "ENALCIVI", xlabel = "Surface Water TP (ug/l)")

## ---- echo = FALSE-------------------------------------------------------
graphics.off()

## ---- fig.height = 5,fig.width = 8---------------------------------------
plotCPs(glades.titan, taxaID = "ENALCIVI", cp.trace = TRUE, xlabel = "Surface Water TP (ug/l)")

## ---- echo = FALSE-------------------------------------------------------
graphics.off()

## ---- fig.height = 5,fig.width = 8---------------------------------------
plotCPs(glades.titan, taxaID = "OSTRASP5", cp.trace = TRUE, xlabel = "Surface Water TP (ug/l)")

## ---- echo = FALSE-------------------------------------------------------
graphics.off()

## ---- fig.height = 6,fig.width = 8---------------------------------------
plotCPs(glades.titan, taxa.dist = FALSE, xlabel = "Surface Water TP (ug/l)")

## ---- echo = FALSE-------------------------------------------------------
graphics.off()

## ---- fig.height = 6,fig.width = 8---------------------------------------
plotCPs(glades.titan, taxa.dist = FALSE, xlabel = "Surface Water TP (ug/l)", stacked = TRUE)

