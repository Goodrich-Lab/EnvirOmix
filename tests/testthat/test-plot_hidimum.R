# early integration hidimum ----
test_that("plot_hidimum", {
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
  result_hidimum_early <- hidimum(exposure = exposure,
                               outcome = outcome,
                               omics_lst = omics_lst,
                               covs = covs,
                               Y.family = "gaussian",
                               M.family = "gaussian",
                               integration = "early")

  # Create Plot
  plotout <- plot_hidimum(result_hidimum_early)

  # Test that plot was created
  testthat::expect_equal(object = length(plotout), expected = 11)
})


