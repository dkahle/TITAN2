<!-- README.md is generated from README.Rmd. Please edit that file -->

# TITAN2

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/TITAN2)](https://cran.r-project.org/package=TITAN2)
[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/dkahle/TITAN2?branch=master&svg=true)](https://ci.appveyor.com/project/dkahle/TITAN2)
<!-- badges: end -->

**TITAN2** is the second R implementation of Threshold Indicator Taxa
ANalysis. It is an R package source controlled with Git on GitHub and
[distributed on CRAN](https://cran.r-project.org/package=TITAN2).

To learn more about **TITAN2**, [check out the vignette
here](https://github.com/dkahle/TITAN2/blob/master/vignettes/titan2-intro.pdf)
(you can click Download to view it in a separate window).

*Note: a previous version of this readme stated that you could read the
vignette; however, the vignette is not built when the package is
downloaded from GitHub, so just access it as above.*

## Installation

You can install **TITAN2** in either of two ways. At the present time,
we recommend installing **TITAN2** from GitHub, as it has several new
features, e.g. `plot_taxa_ridges()`.

-   From Github (dev version):

``` r
if (!requireNamespace("devtools")) install.packages("devtools")
devtools::install_github("dkahle/TITAN2")
```

-   From CRAN:

``` r
install.packages("TITAN2")
```

## Acknowledgements

This work continues to be supported by the [Department of Geography and
Environmental Systems (UMBC)](https://ges.umbc.edu/), [Department of
Biology (Baylor)](https://biology.artsandsciences.baylor.edu/), and
[Department of Statistical Science
(Baylor)](https://statistics.artsandsciences.baylor.edu/).
