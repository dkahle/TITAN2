# TITAN2 2.4.1

MINOR BUG FIXES

* For some reason the color `"light blue"` began erroring in 
  `plot_taxa_ridges()`. This has been fixed to `"lightblue"` to provide clean 
  CRAN checks.
  
* The relative heights of the graphics in `plot_taxa_ridges()` now have better
  defaults. These can be edited with the new `rel_heights` argument to that 
  function.
