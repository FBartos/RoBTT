context("(4) Diagnostics functions")
skip_on_cran()

# the plotting functions are imported from BayesTools and tested henceforth
# test objects - assuming that the fit function worked properly
saved_files <- paste0("fit_", 1:4, ".RDS")
saved_fits  <- list()
for(i in seq_along(saved_files)){
  saved_fits[[i]] <- readRDS(file = file.path(testthat::test_path("results/fits"), saved_files[i]))
}


test_that("Parameter plots work", {

  p1 <- diagnostics(saved_fits[[1]], parameter = "mu",    type = "trace")
  p2 <- diagnostics(saved_fits[[1]], parameter = "sigma", type = "autocorelation")
  p3 <- diagnostics(saved_fits[[1]], parameter = "delta", type = "density")
  p4 <- diagnostics(saved_fits[[1]], parameter = "rho",   type = "trace")
  p5 <- diagnostics(saved_fits[[1]], parameter = "nu",    type = "autocorelation")
  
  # ggplot
  for(i in 1:8){
    
    expect_doppelganger(paste0("diagnostics_mu_",i),    p1[[i]])
    expect_doppelganger(paste0("diagnostics_sigma_",i), p2[[i]])
    
    if(i %in% 1:4){
      expect_null(p3[[i]])
    }else{
      expect_doppelganger(paste0("diagnostics_delta_",i), p3[[i]])
    }
    
    if(i %in% c(1,2,5,6)){
      expect_null(p4[[i]])
    }else{
      expect_doppelganger(paste0("diagnostics_rho_",i), p4[[i]])
    }
    
    if(i %in% c(1,3,5,7)){
      expect_null(p5[[i]])
    }else{
      expect_doppelganger(paste0("diagnostics_nu_",i), p5[[i]])
    }
    
  }

})
