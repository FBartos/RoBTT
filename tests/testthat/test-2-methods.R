context("(2) Print and summary functions")
skip_on_cran()

# the summary tables and print functions are imported from BayesTools and tested henceforth
# test objects - assuming that the fit function worked properly
saved_files <- paste0("fit_", 1:4, ".RDS")
saved_fits  <- list()
for(i in seq_along(saved_files)){
  saved_fits[[i]] <- readRDS(file = file.path(testthat::test_path("results/fits"), saved_files[i]))
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
      "Heterogeneity    4/8       0.500       0.371             0.526" ,
      "Outliers         4/8       0.500       0.548            -0.192" ,
      ""                                                               ,
      "Model-averaged estimates:"                                      ,
      "       Mean Median   0.1   0.5   0.9"                           ,
      "delta 0.092  0.000 0.000 0.000 0.467"                           ,
      "rho   0.517  0.500 0.414 0.500 0.686"                           ,
      "nu      Inf  4.425 2.242 4.425   Inf"                           ,
      "\033[0;31mModel (6): There were 1 divergent transitions.\033[0m",
      ""                                                               ,
      "Model-averaged group parameter estimates:"                      ,
      "          Mean Median   0.1   0.5   0.9"                        ,
      "mu[1]    1.691  1.689 0.264 1.689 3.148"                        ,
      "mu[2]    2.035  1.944 0.576 1.944 3.630"                        ,
      "sigma[1] 4.937  3.801 2.479 3.801 7.802"                        ,
      "sigma[2] 5.115  3.946 2.610 3.946 7.920"                        ,
      "\033[0;31mModel (6): There were 1 divergent transitions.\033[0m",
      ""                                                               ,
      "Conditional estimates:"                                         ,
      "       Mean Median    0.1   0.5   0.9"                          ,
      "delta 0.272  0.233 -0.261 0.233 0.862"                          ,
      "rho   0.549  0.558  0.276 0.558 0.805"                          ,
      "nu    3.054  2.754  2.136 2.754 4.362"                          ,
      "\033[0;31mModel (6): There were 1 divergent transitions.\033[0m"
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
      "     1       normal                S(0.5)                0.125       -24.39       0.175        1.490",
      "     2            t                S(0.5)     E(1)       0.125       -24.07       0.244        2.253",
      "     3       normal               B(1, 1)                0.125       -24.97       0.099        0.769",
      "     4            t               B(1, 1)     E(1)       0.125       -24.54       0.152        1.253",
      "     5       normal  C(0, 0.71)    S(0.5)                0.125       -24.82       0.115        0.911",
      "     6            t  C(0, 0.71)    S(0.5)     E(1)       0.125       -25.01       0.094        0.730",
      "     7       normal  C(0, 0.71)   B(1, 1)                0.125       -25.43       0.062        0.466",
      "     8            t  C(0, 0.71)   B(1, 1)     E(1)       0.125       -25.50       0.058        0.432",
      "\033[0;31mModel (6): There were 1 divergent transitions.\033[0m"  
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

  saved_files <- paste0("fit_", 1:4, ".RDS")
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
