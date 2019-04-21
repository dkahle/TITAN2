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

