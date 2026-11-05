## development version
### Fixes
- spurious BayesTools version warning with development versions of BayesTools

## version 1.4.0
### Changes
- update to BayesTools 0.4.0 (seeded results differ from earlier versions and models fitted with earlier versions need to be refitted)
- `interpret()` reports the conditional estimate of the degrees of freedom

### Fixes
- prior densities in `plot()` with `transform_rho = TRUE`

## version 1.3.1
- vignette and references updates 

## version 1.3.0
- adds `truncation` argument to `RoBTT()` function to allow for truncation of the normal/t-distribution models

## version 1.2.2
### Features
- `update.RoBTT()` function to refit models with convergence warnings or change prior model probabilities

## version 1.2.1
Compatibility update for the rstan 2.26.

## version 1.2.0
### Changes
- when specifying prior distribution for the nu parameter, coerce Spike(Inf) to NULL
- better handling of NULL in prior distribution settings and better coercing to normal likelihood

### Fixes
- conditional posterior distributions for degrees of freedom

## version 1.1.0
### Features
- `check_setup()` function
- `diagnostics()` function
- `rho2logsdr` transformations

## version 1.0.3
Compatibility update for the rstan 2.31.
(correct CRAN commit)

## version 1.0.1
Compatibility update for the rstan 2.26.
