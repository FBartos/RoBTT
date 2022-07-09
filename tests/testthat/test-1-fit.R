context("(1) Fitting and updating functions")
skip_on_cran()
skip_on_covr()

# test objects
saved_files <- paste0("fit_", 1:13, ".RDS")
saved_fits  <- list()
for(i in seq_along(saved_files)){
  saved_fits[[i]] <- readRDS(file = file.path("../results/fits", saved_files[i]))
}

# functions simplifying the comparison
clean_all  <- function(fit, only_samples = TRUE){
  if(only_samples){
    fit$data     <- NULL
    fit$add_info <- NULL
    fit$control  <- NULL
    fit$models   <- NULL
  }
  return(fit)
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

  fit1 <- suppressWarnings(RoBTT(x1 = x1, x2 = x2, seed = 1, parallel = FALSE))
  expect_equal(clean_all(saved_fits[[1]]), clean_all(fit1))
  
  
  fit2 <- suppressWarnings(RoBTT(mean1 = mean1, mean2 = mean2, sd1 = sd1, sd2 = sd2, N1 = N1, N2 = N2, likelihood = "normal", seed = 1, parallel = FALSE))
  expect_equal(clean_all(saved_fits[[2]]), clean_all(fit2))
  
  fit3 <- suppressWarnings(RoBTT(x1 = x1, x2 = x2, seed = 1, parallel = FALSE, 
                                 prior_delta      = prior("cauchy", list(0, 1/sqrt(2)), list(0, Inf)),
                                 prior_rho        = prior("beta",   list(3, 3)),
                                 prior_nu         = prior("exp",    list(1)),
                                 prior_delta_null = prior("normal", list(0, 0.15), list(0, Inf)),
                                 likelihood       = c("normal", "t"),))
  expect_equal(clean_all(saved_fits[[3]]), clean_all(fit3))

})


#### creating / updating the test settings ####
if(FALSE){
  saved_fits <- list(fit1, fit2, fit3)

  for(i in 1:length(saved_fits)){
    saveRDS(saved_fits[[i]], file = file.path("tests/results/fits", paste0("fit_",i,".RDS")),   compress  = "xz")
  }
}
