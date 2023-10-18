# TITAN2 2.4.3

MINOR CHANGES

* `plot_taxa_ridges()`: 

* `plot_taxa()`: 


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
