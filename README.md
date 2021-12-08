# inters

  <!-- badges: start -->
  [![R-CMD-check](https://github.com/mattblackwell/inters/workflows/R-CMD-check/badge.svg)](https://github.com/mattblackwell/inters/actions)
  <!-- badges: end -->

`inters` is an R package with tools to estimate interactions. Currently, it implements the post-double-selection approach to estimating interactions described in this [forthcoming paper by Matthew Blackwell and Michael Olson][lasso-paper]. To install the development version of `inters`, run the following code in R:
```R
require(devtools)
install_github("mattblackwell/inters", build_vignettes = TRUE)
```

You can find an example of how to use the package in [our vignette](articles/post-double-selection.html).

[lasso-paper]: http://www.mattblackwell.org/files/papers/lasso-inters.pdf
