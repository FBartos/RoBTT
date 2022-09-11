README
================

<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/RoBTT)](https://CRAN.R-project.org/package=RoBTT)
<!-- badges: end -->

<!-- 
#[![R-CRAN-check](https://github.com/FBartos/RoBTT/workflows/R-CMD-check/badge.svg)](https://github.com/FBartos/RoBTT/actions)
[![R-tests](https://github.com/FBartos/RoBTT/workflows/R-CMD-tests/badge.svg)](https://github.com/FBartos/RoBTT/actions)
#[![Codecov test coverage](https://codecov.io/gh/FBartos/RoBTT/branch/master/graph/badge.svg)](https://app.codecov.io/gh/FBartos/RoBTT?branch=master)
-->

# Robust Bayesian T-Test (RoBTT)

This package provides an implementation of Bayesian model-averaged
t-test that allows users to draw inference about the presence vs absence
of the effect, heterogeneity of variances, and outliers. The RoBTT
packages estimates model ensembles of models created as a combination of
the competing hypotheses and uses Bayesian model-averaging to combine
the models using posterior model probabilities. Users can obtain the
model-averaged posterior distributions and inclusion Bayes factors which
account for the uncertainty in the data generating process. User can
define a wide range of informative priors for all parameters of
interest. The package provides convenient functions for summary,
visualizations, and fit diagnostics.

See our pre-print Maier et al. (forthcoming) introducing the method.

We also prepared a vignette that illustrate functionality of the
package:

-   [Introduction to
    RoBTT](https://fbartos.github.io/RoBTT/articles/Introduction_to_RoBTT.html)

## Installation

The release version can be installed from CRAN:

``` r
install.packages("RoBTT")
```

and the development version of the package can be installed from GitHub:

``` r
devtools::install_github("fbartos/RoBTT")
```

### References

<div id="refs" class="references csl-bib-body hanging-indent"
line-spacing="2">

<div id="ref-maier2022bayesian" class="csl-entry">

Maier, M., Bartoš, F., Quintana, D. S., Bergh, D. van den, Marsman, M.,
Ly, A., & Wagenmakers, E.-J. (forthcoming). Bayesian model-averaged
t-tests. In *PsyArxiv*.

</div>

</div>
