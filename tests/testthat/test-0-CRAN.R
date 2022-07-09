context("(0) Basic tests for CRAN")
# These are just a very rudimentary tests that don't require time or saved files.
# The full range of tests is run locally.


test_that("Basic functionality works", {

  # fit a default model
  fit_default <- suppressWarnings(RoBTT(x1 = c(1, 2, 3), x2 = c(0, 1, 2), chains = 1, warmup = 50, iter = 100, seed = 1))

  expect_equal(TRUE, is.RoBTT(fit_default))

  expect_equal(
    capture_output_lines(fit_default, print = TRUE, width = 150),
    c("Call:"                                                           ,
      "RoBTT(x1 = c(1, 2, 3), x2 = c(0, 1, 2), chains = 1, iter = 100, ",
      "    warmup = 50, seed = 1)"                                      ,
      ""                                                                ,
      "Estimates:"                                                      ,
      "     delta        rho "                                          ,
      "-0.1988558  0.5066339 " 
    )
  )

  expect_equal(
    capture_output_lines(summary(fit_default), print = TRUE, width = 150),
    c("Call:"                                                                                       ,
      "RoBTT(x1 = c(1, 2, 3), x2 = c(0, 1, 2), chains = 1, iter = 100, "                            ,
      "    warmup = 50, seed = 1)"                                                                  ,
      ""                                                                                            ,
      "Robust Bayesian t-test"                                                                      ,
      "Components summary:"                                                                         ,
      "              Models Prior prob. Post. prob. Inclusion BF"                                   ,
      "Effect           4/8       0.500       0.436        0.775"                                   ,
      "Heterogeneity    4/8       0.500       0.417        0.717"                                   ,
      "Outliers         4/8       0.500       0.401        0.668"                                   ,
      ""                                                                                            ,
      "Model-averaged estimates:"                                                                   ,
      "        Mean Median  0.025 0.975"                                                            ,
      "delta -0.199  0.000 -1.543 0.366"                                                            ,
      "rho    0.507  0.500  0.191 0.862"                                                            ,
      "nu       Inf    Inf  2.085   Inf"                                                            ,
      "\033[0;31mModel (1): Minimum effective sample size was low (23).\033[0m"                     ,
      "\033[0;31mModel (1): Maximum R-hat was large (1.18).\033[0m"                                 ,
      "\033[0;31mModel (2): Minimum effective sample size was low (9).\033[0m"                      ,
      "\033[0;31mModel (2): Maximum R-hat was large (1.13).\033[0m"                                 ,
      "\033[0;31mModel (3): Minimum effective sample size was low (10).\033[0m"                     ,
      "\033[0;31mThere were another 8 warnings. To see all warnings call 'check_RoBTT(fit)'.\033[0m"
  )
  )
})
