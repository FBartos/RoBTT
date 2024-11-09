context("(3) Plot functions")
skip_on_cran()

# the plotting functions are imported from BayesTools and tested henceforth
# test objects - assuming that the fit function worked properly
saved_files <- paste0("fit_", 1:4, ".RDS")
saved_fits  <- list()
for(i in seq_along(saved_files)){
  saved_fits[[i]] <- readRDS(file = file.path(testthat::test_path("../results/fits"), saved_files[i]))
}


test_that("Parameter plots work", {

  i <- 1
  # ggplot
  expect_doppelganger(paste0("ggplot_delta1_",i), plot(saved_fits[[i]], "delta", plot_type = "ggplot"))
  expect_doppelganger(paste0("ggplot_delta2_",i), plot(saved_fits[[i]], "delta", prior = TRUE, plot_type = "ggplot"))

  # default base plot
  expect_doppelganger(paste0("plot_delta1_",i), function()plot(saved_fits[[i]], "delta"))
  expect_doppelganger(paste0("plot_delta2_",i), function()plot(saved_fits[[i]], "delta", prior = TRUE))
  expect_doppelganger(paste0("plot_delta3_",i), function()plot(saved_fits[[i]], "delta", conditional = TRUE))
  expect_doppelganger(paste0("plot_delta4_",i), function()plot(saved_fits[[i]], "delta", conditional = TRUE, prior = TRUE))

  # additional settings
  expect_doppelganger(paste0("plot_delta5_",i), function()plot(saved_fits[[i]], "delta", prior = TRUE, dots_prior = list(col = "blue", lty = 2), col = "red", lty = 2, xlim = c(0, 1), main = "Title"))

  ### heterogeneity
  i <- 1
  # default ggplot2
  expect_doppelganger(paste0("ggplot_rho1_",i), plot(saved_fits[[i]], "rho", plot_type = "ggplot"))
  expect_doppelganger(paste0("ggplot_rho2_",i), plot(saved_fits[[i]], "rho", prior = TRUE, plot_type = "ggplot"))
  
  # default base plot
  expect_doppelganger(paste0("plot_rho1_",i), function()plot(saved_fits[[i]], "rho"))
  expect_doppelganger(paste0("plot_rho2_",i), function()plot(saved_fits[[i]], "rho", prior = TRUE))
  expect_doppelganger(paste0("plot_rho3_",i), function()plot(saved_fits[[i]], "rho", conditional = TRUE))
  expect_doppelganger(paste0("plot_rho4_",i), function()plot(saved_fits[[i]], "rho", conditional = TRUE, prior = TRUE))
  
  ### heterogeneity (standard deviation ratio)
  i <- 1
  # default ggplot2
  expect_doppelganger(paste0("ggplot_lsdr1_",i), plot(saved_fits[[i]], "rho", plot_type = "ggplot", transform_rho = TRUE))
  expect_doppelganger(paste0("ggplot_lsdr2_",i), plot(saved_fits[[i]], "rho", prior = TRUE, plot_type = "ggplot", transform_rho = TRUE))
  
  # default base plot
  expect_doppelganger(paste0("plot_lsdr1_",i), function()plot(saved_fits[[i]], "rho", transform_rho = TRUE))
  expect_doppelganger(paste0("plot_lsdr2_",i), function()plot(saved_fits[[i]], "rho", prior = TRUE, transform_rho = TRUE))
  expect_doppelganger(paste0("plot_lsdr3_",i), function()plot(saved_fits[[i]], "rho", conditional = TRUE, transform_rho = TRUE))
  expect_doppelganger(paste0("plot_lsdr4_",i), function()plot(saved_fits[[i]], "rho", conditional = TRUE, prior = TRUE, transform_rho = TRUE))
  

  ### degrees of freedom
  i <- 1
  # default ggplot2
  expect_doppelganger(paste0("ggplot_nu1_",i), plot(saved_fits[[i]], "nu", plot_type = "ggplot"))
  expect_doppelganger(paste0("ggplot_nu2_",i), plot(saved_fits[[i]], "nu", prior = TRUE, plot_type = "ggplot"))
  
  # default base plot
  expect_doppelganger(paste0("plot_nu1_",i), function()plot(saved_fits[[i]], "nu"))
  expect_doppelganger(paste0("plot_nu2_",i), function()plot(saved_fits[[i]], "nu", prior = TRUE))
  expect_doppelganger(paste0("plot_nu3_",i), function()plot(saved_fits[[i]], "nu", conditional = TRUE))
  expect_doppelganger(paste0("plot_nu4_",i), function()plot(saved_fits[[i]], "nu", conditional = TRUE, prior = TRUE))

})
