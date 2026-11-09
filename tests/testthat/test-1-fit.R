context("(1) Fitting and updating functions")
skip_on_cran()

# test objects
saved_files <- paste0("fit_", 1:4, ".RDS")
saved_fits  <- list()
for(i in seq_along(saved_files)){
  saved_fits[[i]] <- readRDS(file = file.path(testthat::test_path("../results/fits"), saved_files[i]))
}

# A refit matches a saved fit only up to Monte Carlo error: compilers, math libraries
# and even the state of the R process change the draws. Each bound is about twice the
# largest difference from seed 1 across 48 seeds of these four fits.
saved_fit_bounds <- c(log_marglik = 0.08, log_BF = 0.04, post_probs = 0.011, mean = 0.04, nu_mean = 0.11)
saved_fit_values <- function(fit){
  fit_summary <- summary(fit, conditional = TRUE)
  estimates   <- fit_summary$estimates
  conditional <- fit_summary$estimates_conditional
  is_nu       <- rownames(estimates) == "nu"
  values      <- list(
    log_marglik = vapply(fit$models, function(model) model$marglik$logml, numeric(1)),
    log_BF      = vapply(fit$RoBTT$inference, attr, numeric(1), "log_BF"),
    post_probs  = fit$RoBTT$inference[[1]]$post_probs,
    mean        = c(estimates$Mean[!is_nu], conditional$Mean[!is_nu]),
    nu_mean     = conditional$Mean[is_nu]
  )
  return(values[lengths(values) > 0])
}
expect_saved_fit <- function(saved_fit, fit){
  saved_values <- saved_fit_values(saved_fit)
  fit_values   <- saved_fit_values(fit)
  expect_identical(lengths(fit_values), lengths(saved_values))
  for(value in names(saved_values)){
    expect_lt(max(abs(fit_values[[value]] - saved_values[[value]])), saved_fit_bounds[[value]], label = value)
  }
}
try_parallel <- function(x, rep = 3){
  temp_fit <- NULL
  i        <- 0
  while(is.null(temp_fit) & i < rep){
    temp_fit <- tryCatch(eval(x), error = function(e) NULL)
    i        <- i + 1
  }
  return(temp_fit)
}


# create mock data
mean1 <- 3
mean2 <- 2.5
sd1   <- 1.7
sd2   <- 2.2
N1    <- 10
N2    <- 15
x1    <- c(1, 0.5, -3, 2, 5)
x2    <- c(0, 0.6,  9, 4, 2) 


test_that("Default model works", {
  
  set.seed(1)
  
  fit1 <- suppressWarnings(RoBTT(x1 = x1, x2 = x2, seed = 1, parallel = FALSE))
  expect_saved_fit(saved_fits[[1]], fit1)
  
  
  fit2 <- suppressWarnings(RoBTT(mean1 = mean1, mean2 = mean2, sd1 = sd1, sd2 = sd2, N1 = N1, N2 = N2, prior_nu = NULL, seed = 1, parallel = FALSE))
  expect_saved_fit(saved_fits[[2]], fit2)
  
  fit3 <- suppressWarnings(RoBTT(x1 = x1, x2 = x2, seed = 1, parallel = FALSE, 
                                 prior_delta      = prior("cauchy", list(0, 1/sqrt(2)), list(0, Inf)),
                                 prior_rho        = prior("beta",   list(3, 3)),
                                 prior_nu         = prior("exp",    list(1)),
                                 prior_delta_null = prior("normal", list(0, 0.15), list(0, Inf))))
  expect_saved_fit(saved_fits[[3]], fit3)
  
  # update prior model probs
  fit2u <- update(fit2, prior_weights = c(2, 1, 2, 2))
  
  fit1u <- update(fit1)
  
  # fit a truncated t-test
  fit4 <- suppressWarnings(RoBTT(x1 = x1, x2 = x2, seed = 1, truncation = list(x1 = c(-4, 6), x2 = c(-1, 10)),parallel = FALSE))
  expect_saved_fit(saved_fits[[4]], fit4)
})


#### creating / updating the test settings ####
if(FALSE){
  saved_fits <- list(fit1, fit2, fit3, fit4)

  for(i in 1:length(saved_fits)){
    saveRDS(saved_fits[[i]], file = file.path("tests/results/fits", paste0("fit_",i,".RDS")),   compress  = "xz")
  }
}
