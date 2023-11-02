# early integration lafamum ----
test_that("lafamum early", {
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
  result_hidimum_early <- lafamum(exposure = exposure,
                               outcome  = outcome,
                               omics_lst = omics_lst,
                               covs = covs,
                               Y.family = "gaussian",
                               integration = "early")

  testthat::expect_equal(object = length(result_hidimum_early), expected = 3)
})



# intermediate integration lafamum ----
test_that("lafamum intermediate", {
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
  result_hidimum_int <- lafamum(exposure = exposure,
                             outcome = outcome,
                             omics_lst = omics_lst,
                             covs = covs,
                             Y.family = "gaussian",
                             integration = "intermediate")

  testthat::expect_equal(object = length(result_hidimum_int),
                         expected = 3)
})


# late integration lafamum ----
test_that("lafamum late", {
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
  result_hidimum_late <- lafamum(exposure = exposure,
                             outcome = outcome,
                             omics_lst = omics_lst,
                             covs = covs,
                             Y.family = "gaussian",
                             integration = "late")

  testthat::expect_equal(object = length(result_hidimum_late),
                         expected = 3)
})
