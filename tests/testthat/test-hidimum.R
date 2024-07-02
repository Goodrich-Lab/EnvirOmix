# hidimum with Early, intermediate, and late integration ----
test_that("test hidimum", {
  simulated_data <- simulated_data
  # Define exposure and outcome name
  covars <- c("e3_sex_None", "hs_child_age_yrs_None")

  # Extract exposure and outcome data
  exposure <- simulated_data[["phenotype"]]$hs_hg_m_scaled
  outcome  <- simulated_data[["phenotype"]]$ck18_scaled

  # Get numeric matrix of covariates
  covs <- simulated_data[["phenotype"]][covars]
  covs$e3_sex_None <- ifelse(covs$e3_sex_None == "male", 1, 0)

  # create list of omics data
  omics_lst <- simulated_data[-which(names(simulated_data) == "phenotype")]

  # Check that the number of cores is not greater than the max allowed
  chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")

  if (nzchar(chk) && chk == "TRUE") {
    # use 2 cores in CRAN/Travis/AppVeyor
    num_workers <- 2L
  } else {
    # use all cores in devtools::test()
    num_workers <- parallel::detectCores()
  }

  ## Run Early Analysis ----
  result_hidimum_early <- hidimum(exposure = exposure,
                                  outcome = outcome,
                                  omics_lst = omics_lst,
                                  covs = covs,
                                  Y.family = "gaussian",
                                  M.family = "gaussian",
                                  integration = "early")

  ### Test Early -------
  testthat::expect_equal(object = ncol(result_hidimum_early), expected = 15)

  ## Run Intermediate Analysis ----
  set.seed(1234)
  result_hidimum_int <- hidimum(exposure = exposure,
                                outcome = outcome,
                                omics_lst = omics_lst,
                                covs = covs,
                                Y.family = "gaussian",
                                M.family = "gaussian",
                                integration = "intermediate",
                                n_boot = num_workers,
                                n_cores = num_workers)

  ### Test intermediate -------
  testthat::expect_equal(object = ncol(result_hidimum_int), expected = 16)
  testthat::expect_equal(object = nrow(result_hidimum_int), expected = 140)

  ## Run Late Analysis ----
  result_hidimum_int <- hidimum(exposure = exposure,
                                outcome = outcome,
                                omics_lst = omics_lst,
                                covs = covs,
                                Y.family = "gaussian",
                                M.family = "gaussian",
                                integration = "late")

  ### Test that late ran -------
  testthat::expect_equal(object = ncol(result_hidimum_int), expected = 14)
  testthat::expect_equal(object = nrow(result_hidimum_int), expected = 49)

  # Compare results to HIMA analysis run for a single layer
  hima_individual <- HIMA::hima(X = exposure,
                                Y = outcome,
                                M = omics_lst[[1]],
                                COV.XM = covs,
                                Y.family = "gaussian",
                                M.family = "gaussian",
                                max.iter = 100000,
                                scale = FALSE)

  # Subset late analysis to only include the methylome
  himum_l_m <- result_hidimum_int[result_hidimum_int$omic_layer == "methylome",]

  testthat::expect_equal(object = hima_individual$alpha,
                         expected = himum_l_m$alpha)

})
