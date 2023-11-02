# early integration hidimum ----
test_that("hidimum early", {
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

  # Run Analysis
  result_hima_early <- hidimum(exposure = exposure,
                               outcome = outcome,
                               omics_lst = omics_lst,
                               covs = covs,
                               Y.family = "gaussian",
                               M.family = "gaussian",
                               integration = "early")

  testthat::expect_equal(object = ncol(result_hima_early),
                         expected = 15)

})



# intermediate integration hidimum ----
test_that("hidimum intermediate", {
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

  # Run Analysis
  result_hima_int <- hidimum(exposure = exposure,
                             outcome = outcome,
                             omics_lst = omics_lst,
                             covs = covs,
                             Y.family = "gaussian",
                             M.family = "gaussian",
                             integration = "intermediate",
                             n_boot = 12)

  testthat::expect_equal(object = ncol(result_hima_int),
                         expected = 16)
  testthat::expect_gt(object = nrow(result_hima_int), expected = 1)
})


# late integration hidimum ----
test_that("hidimum late", {
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

  # Run Analysis
  result_hima_int <- hidimum(exposure = exposure,
                             outcome = outcome,
                             omics_lst = omics_lst,
                             covs = covs,
                             Y.family = "gaussian",
                             M.family = "gaussian",
                             integration = "late")

  testthat::expect_equal(object = ncol(result_hima_int),
                         expected = 14)
  testthat::expect_gt(object = nrow(result_hima_int), expected = 1)
})
