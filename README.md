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
t-tests that allows users to draw inference about the presence vs
absence of the effect, heterogeneity of variances, and outliers. The
RoBTT packages estimates model ensembles of models created as a
combination of the competing hypotheses and uses Bayesian
model-averaging to combine the models using posterior model
probabilities. Users can obtain the model-averaged posterior
distributions and inclusion Bayes factors which account for the
uncertainty in the data generating process. User can define a wide range
of informative priors for all parameters of interest. The package
provides convenient functions for summary, visualizations, and fit
diagnostics.

See our manuscripts for more information about the methodology:

- Maier et al. (2022) introduces a robust Bayesian t-test that
  model-averages over normal and t-distributions to account for the
  uncertainty about potential outliers,
- Godmann et al. (2024) introduces a truncated Bayesian t-test that
  accounts for outlier exclusion when estimating the models.

We also prepared vignettes that illustrate functionality of the package:

- [Introduction to
  RoBTT](https://fbartos.github.io/RoBTT/articles/Introduction_to_RoBTT.html)
- [Truncated
  T-Tests](https://fbartos.github.io/RoBTT/articles/Truncated_t_test.html)

## Installation

The release version can be installed from CRAN:

``` r
install.packages("RoBTT")
```

and the development version of the package can be installed from GitHub:

``` r
devtools::install_github("FBartos/RoBTT")
```

### References

<div id="refs" class="references csl-bib-body hanging-indent"
line-spacing="2">

<div id="ref-godmann2024how" class="csl-entry">

Godmann, H. R., Bartoš, F., & Wagenmakers, E.-J. (2024). *A truncated
t-test: Excluding outliers without biasing the Bayes factor*.

</div>

<div id="ref-maier2022bayesian" class="csl-entry">

Maier, M., Bartoš, F., Quintana, D. S., Bergh, D. van den, Marsman, M.,
Ly, A., & Wagenmakers, E.-J. (2022). *Model-averaged Bayesian t-tests*.
<https://doi.org/10.31234/osf.io/d5zwc>

</div>

</div>
