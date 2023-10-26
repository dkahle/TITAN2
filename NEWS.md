# TITAN2 2.4.3

MINOR CHANGES

* `plot_taxa_ridges()`: 
 
  - Changed the rank order of taxa plotted on y-axis to match the rank order of
the 2nd quartile (median/50th percentile) of the bootstrap change point
estimates. The previous version plotted order of taxa based on the 2nd quartile
(median) of the estimated values from the density function, which did not always
match the median of the n bootstrap change-point estimates.

  - Changed the location of the vertical line on x for each taxa to match the
2nd quartile (median/50th percentile) of the bootstrap change point estimates
The previous version plotted the location of the vertical line on x based on the
2nd quartile (median) of the estimated values from the density function, which
did not always match the median of the n bootstrap change-point estimates.

  - Added new argument `d_lines`. Short for "density lines", a logical
determining whether the median change point values estimated from each taxon's
density function, is plotted as a vertical dashed line. When `d_lines = TRUE`,
line is added to the plot but it does not replace the vertical solid line
corresponding to the median of the bootstrap change point estimates.

  - Added new argument `breaks`. Determines values and labels of breaks (ticks)
on the x-axis, allowing the user greater control of appearance of axis ticks,
labels, and gridlines.

  - Removed `pur.cut` and `rel.cut` arguments. Plotting function code relies on
output from `titan()`, which is where users specify pur.cut and rel.cut.


* `plot_taxa()`: 

  - Changed default plotting of taxa along y to match the median of the
bootstrap change point estimates. This aligns with the default plotting in
`plot_taxa_ridges()`.

  - Changed default plotting of locations of symbols (points sized by z-scores)
to match the median of the bootstrap change point estimates.This aligns with the
default plotting in `plot_taxa_ridges()`.

  - Changed default cex parameters to smaller values to allow for all taxa names
in glades.titan to be plotted in default runs.

* `plot_sumz_density()`: `guide = FALSE` has been changed to `guide = "none"` to
  conform to **ggplot2**'s changes.


# TITAN2 2.4.2

MINOR CHANGES

* The snow dependency was removed. Now cluster types default to those of 
  `parallel::makeCluster()`.


MINOR BUG FIXES

* A bug was removed from `plot_taxa_ridges()` that rounded the x axis values to 
  the nearest hundredth (`round(..., digits = 2)`), which had the effect of 
  heaping small-scale data, making it impossible to differentiate different 
  taxa and meaningfully impacting the graphic. This rounding has been removed.


# TITAN2 2.4.1

MINOR BUG FIXES

* For some reason the color `"light blue"` began erroring in 
  `plot_taxa_ridges()`. This has been fixed to `"lightblue"` to provide clean 
  CRAN checks.
  
* The relative heights of the graphics in `plot_taxa_ridges()` now have better
  defaults. These can be edited with the new `rel_heights` argument to that 
  function.
