context("(2) Print and summary functions")
skip_on_cran()

# the summary tables and print functions are imported from BayesTools and tested henceforth
# test objects - assuming that the fit function worked properly
saved_files <- paste0("fit_", 1:3, ".RDS")
saved_fits  <- list()
for(i in seq_along(saved_files)){
  saved_fits[[i]] <- readRDS(file = file.path("../results/fits", saved_files[i]))
}

test_that("Print functions work", {

  # testing consistency across all model specifications
  for(i in 1:length(saved_fits)){
    expect_equal(
      capture_output_lines(saved_fits[[i]], print = TRUE, width = 150),
      read.table(file = file.path("../results/print", paste0(i, ".txt")), header = FALSE, blank.lines.skip = FALSE)[,1])
  }
  
})

test_that("Summary functions work", {

  # testing consistency across all model specifications
  for(i in 1:length(saved_fits)){
    expect_equal(
      capture_output_lines(summary(saved_fits[[i]]), print = TRUE, width = 150),
      read.table(file = file.path("../results/summary", paste0(i, ".txt")), header = FALSE, blank.lines.skip = FALSE)[,1])
  }

  # all options
  expect_equal(
    capture_output_lines(summary(saved_fits[[1]], conditional = TRUE, logBF = TRUE, BF01 = TRUE, group_estimates = TRUE, probs = c(0.10, 0.50, .90)), print = TRUE, width = 150),
    c("Call:"                                                          ,
      "RoBTT(x1 = x1, x2 = x2, parallel = FALSE, seed = 1)"            ,
      ""                                                               ,
      "Robust Bayesian t-test"                                         ,
      "Components summary:"                                            ,
      "              Models Prior prob. Post. prob. log(Exclusion BF)" ,
      "Effect           4/8       0.500       0.330             0.707" ,
      "Heterogeneity    4/8       0.500       0.373             0.520" ,
      "Outliers         4/8       0.500       0.547            -0.190" ,
      ""                                                               ,
      "Model-averaged estimates:"                                      ,
      "       Mean Median   0.1   0.5   0.9"                           ,
      "delta 0.098  0.000 0.000 0.000 0.484"                           ,
      "rho   0.517  0.500 0.414 0.500 0.689"                           ,
      "nu      Inf  4.433 2.228 4.433   Inf"                           ,
      "\033[0;31mModel (7): There were 1 divergent transitions.\033[0m",
      ""                                                               ,
      "Model-averaged group parameter estimates:"                      ,
      "          Mean Median   0.1   0.5   0.9"                        ,
      "mu[1]    1.693  1.670 0.279 1.670 3.177"                        ,
      "mu[2]    2.045  1.954 0.588 1.954 3.659"                        ,
      "sigma[1] 5.032  3.817 2.474 3.817 8.003"                        ,
      "sigma[2] 5.220  3.953 2.606 3.953 8.147"                        ,
      "\033[0;31mModel (7): There were 1 divergent transitions.\033[0m",
      ""                                                               ,
      "Conditional estimates:"                                         ,
      "       Mean Median    0.1   0.5   0.9"                          ,
      "delta 0.274  0.233 -0.271 0.233 0.875"                          ,
      "rho   0.548  0.559  0.270 0.559 0.806"                          ,
      "nu    3.033  2.740  2.128 2.740 4.349"                          ,
      "\033[0;31mModel (7): There were 1 divergent transitions.\033[0m"
  ))
})

test_that("Models summary functions work", {

  # testing consistency across all model specifications
  for(i in 1:length(saved_fits)){
    expect_equal(
      capture_output_lines(summary(saved_fits[[i]], type = "models"), print = TRUE, width = 150),
      read.table(file = file.path("../results/summary.models", paste0(i, ".txt")), header = FALSE, blank.lines.skip = FALSE)[,1])
  }

  # test short names & no spikes
  expect_equal(
    capture_output_lines(summary(saved_fits[[1]], type = "models", short_name = TRUE, remove_spike_0 = TRUE), print = TRUE, width = 150),
    c("Call:"                                                                                               ,
      "RoBTT(x1 = x1, x2 = x2, parallel = FALSE, seed = 1)"                                                 ,
      ""                                                                                                    ,
      "Robust Bayesian t-test"                                                                              ,
      "Models overview:"                                                                                    ,
      " Model Distribution Prior delta Prior rho Prior nu Prior prob. log(marglik) Post. prob. Inclusion BF",
      "     1       normal                S(0.5)                0.125       -24.40       0.175        1.484",
      "     2            t                S(0.5)     E(1)       0.125       -24.07       0.243        2.250",
      "     3       normal               B(1, 1)                0.125       -24.96       0.100        0.774",
      "     4            t               B(1, 1)     E(1)       0.125       -24.54       0.152        1.255",
      "     5       normal  C(0, 0.71)    S(0.5)                0.125       -24.82       0.115        0.910",
      "     6            t  C(0, 0.71)    S(0.5)     E(1)       0.125       -25.02       0.094        0.726",
      "     7       normal  C(0, 0.71)   B(1, 1)                0.125       -25.42       0.063        0.471",
      "     8            t  C(0, 0.71)   B(1, 1)     E(1)       0.125       -25.50       0.058        0.433",
      "\033[0;31mModel (7): There were 1 divergent transitions.\033[0m"            
    ))

})

test_that("Diagnostics summary functions work", {

  # testing consistency across all model specifications
  for(i in 1:length(saved_fits)){
    expect_equal(
      capture_output_lines(summary(saved_fits[[i]], type = "diagnostics"), print = TRUE, width = 200),
      read.table(file = file.path("../results/summary.diagnostics", paste0(i, ".txt")), header = FALSE, blank.lines.skip = FALSE)[,1])
  }
  
})

test_that("Individual summary functions work", {

  # testing consistency across all model specifications
  for(i in 1:length(saved_fits)){
    expect_equal(
      capture_output_lines(summary(saved_fits[[i]], type = "individual"), print = TRUE, width = 150),
      read.table(file = file.path("../results/summary.individual", paste0(i, ".txt")), header = FALSE, blank.lines.skip = FALSE)[,1])
  }

})

test_that("Interpret functions work", {

  # testing consistency across all model specifications
  for(i in 1:length(saved_fits)){
    expect_equal(
      gsub("(.{80})", "\\1\\\n", interpret(saved_fits[[i]])),
      read.table(file = file.path("../results/interpret", paste0(i, ".txt")), header = FALSE, blank.lines.skip = FALSE)[,1])
  }

})

#### creating / updating the test settings ####
if(FALSE){

  saved_files <- paste0("fit_", 1:3, ".RDS")
  saved_fits  <- list()
  for(i in seq_along(saved_files)){
    saved_fits[[i]] <- readRDS(file = file.path("tests/results/fits", saved_files[i]))
  }

  # generate print files
  for(i in seq_along(saved_fits)){
    write.table(capture_output_lines(saved_fits[[i]], print = TRUE, width = 150), file = file.path("tests/results/print", paste0(i, ".txt")), row.names = FALSE, col.names = FALSE)
  }

  # generate summary files
  for(i in seq_along(saved_fits)){
    write.table(capture_output_lines(summary(saved_fits[[i]]), print = TRUE, width = 150), file = file.path("tests/results/summary", paste0(i, ".txt")), row.names = FALSE, col.names = FALSE)
  }

  # generate summary.models files
  for(i in seq_along(saved_fits)){
    write.table(capture_output_lines(summary(saved_fits[[i]], type = "models"), print = TRUE, width = 150), file = file.path("tests/results/summary.models", paste0(i, ".txt")), row.names = FALSE, col.names = FALSE)
  }

  # generate summary.diagnostics files
  for(i in seq_along(saved_fits)){
    write.table(capture_output_lines(summary(saved_fits[[i]], type = "diagnostics"), print = TRUE, width = 200), file = file.path("tests/results/summary.diagnostics", paste0(i, ".txt")), row.names = FALSE, col.names = FALSE)
  }

  # generate summary.individual files
  for(i in seq_along(saved_fits)){
    write.table(capture_output_lines(summary(saved_fits[[i]], type = "individual"), print = TRUE, width = 150), file = file.path("tests/results/summary.individual", paste0(i, ".txt")), row.names = FALSE, col.names = FALSE)
  }

  # generate summary.individual files
  for(i in seq_along(saved_fits)){
    write.table(gsub("(.{80})", "\\1\\\n", interpret(saved_fits[[i]])), file = file.path("tests/results/interpret", paste0(i, ".txt")), row.names = FALSE, col.names = FALSE)
  }

}
