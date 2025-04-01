#' High Dimensional Multiomic Mediation (hidimum)
#'
#' Given exposure, outcome and multiple omics data, this
#' function runs high dimensional mediation with early, intermediate, or
#' late integration of the multiomics datasets.#'
#'
#' @param exposure A numeric vector for the exposure variable
#' @param outcome A numeric vector for the outcome variable
#' @param omics_lst A list of numeric matrices representing omics data
#' @param covs A numeric matrix representing the covariates
#' @param Y.family A character string indicating the family of the outcome
#' @param M.family A character string indicating the family of the mediator
#' @param integration A character string indicating the integration method (one
#' of early, intermediate, or late)
#' @param bh.fdr Bonferroni-Hochberg FDR correction threshold for selecting
#' significant features from HIMA for early and late integration. Defaults
#' to 0.05.
#' @param n_boot number indicating number of bootstrap estimates to perform
#' for calculating the se of the mediation effect (only applies for
#' intermediate integration, otherwise, ignored)
#' @param n_cores Optional, number of cores for parallelization for
#' bootstrapping for intermediate integration. If unspecified and
#' integration = intermediate, defaults to parallel::detectCores().
#'
#' @return A tidy dataframe summarizing the results of HIMA analysis
#'
#' @import dplyr
#' @import tibble
#' @import purrr
#' @import tidyr
#' @importFrom stringr str_detect
#' @importFrom stats gaussian binomial coef glm lm sd
#' @importFrom utils capture.output tail
#' @importFrom HIMA hima_classic
#' @importFrom xtune xtune estimateVariance
#' @importFrom boot boot
#' @importFrom parallel detectCores
#' @importFrom epiomics owas
#' @importFrom RMediation medci
#
#' @export
#'
#' @examples
#'  # Load Example Data
#'  data("simulated_data")
#'  covars <- c("e3_sex_None", "hs_child_age_yrs_None")
#'
#'  # Extract exposure and outcome data
#'  exposure <- simulated_data[["phenotype"]]$hs_hg_m_scaled
#'  outcome  <- simulated_data[["phenotype"]]$ck18_scaled
#'
#'  # Get numeric matrix of covariates
#'  covs <- simulated_data[["phenotype"]][covars]
#'  covs$e3_sex_None <- ifelse(covs$e3_sex_None == "male", 1, 0)
#'
#'  # create list of omics data
#'  omics_lst <- simulated_data[-which(names(simulated_data) == "phenotype")]
#'
#'  # High Dimensional Multiomic Mediation with Early integration
#'  result_hidimum_early <- hidimum(exposure = exposure,
#'                                  outcome = outcome,
#'                                  omics_lst = omics_lst,
#'                                  covs = covs,
#'                                  Y.family = "gaussian",
#'                                  M.family = "gaussian",
#'                                  integration = "early")
#'
#'  # High Dimensional Multiomic Mediation with Intermediate integration (Not Run)
#'  #     chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
#'  #     if (nzchar(chk) && chk == "TRUE") {
#'  #     # use 2 cores in CRAN/Travis/AppVeyor
#'  #     num_workers <- 2L
#'  #     } else {
#'  #     # use all cores in devtools::test()
#'  #     num_workers <- parallel::detectCores()
#'  #     }
#'  #
#'  #   result_hidimum_int <- hidimum(exposure = exposure,
#'  #                                 outcome = outcome,
#'  #                                 omics_lst = omics_lst,
#'  #                                 covs = covs,
#'  #                                 Y.family = "gaussian",
#'  #                                 M.family = "gaussian",
#'  #                                 integration = "intermediate",
#'  #                                 n_boot = num_workers,
#'  #                                 n_cores = num_workers)
#'  #
#'  # High Dimensional Multiomic Mediation with Late integration
#'  result_hidimum_late <- hidimum(exposure = exposure,
#'                                outcome = outcome,
#'                                omics_lst = omics_lst,
#'                                covs = covs,
#'                                Y.family = "gaussian",
#'                                M.family = "gaussian",
#'                                integration = "late")

hidimum <- function(exposure,
                    outcome,
                    omics_lst,
                    covs,
                    Y.family = "binomial",
                    M.family = "gaussian",
                    integration,
                    n_boot,
                    bh.fdr = 0.05,
                    n_cores = NULL) {

  # Set all of these variables to NULL to fix message that they are not found
  alpha_hat <- beta_hat <- IDE <- rimp <-`% Total Effect scaled` <- `% total effect` <- Alpha <- Beta <- BH.FDR <-
    `Correlation TME (%)` <- beta_bootstrap <- data <- estimate <- feature <-
    feature_name <- ftr_name <- in_ind_omic <- indirect <- lcl <- lf_named <-
    lf_num <- lf_numeric <- lf_ordered <- name <- omic_layer <- omic_num <-
    omic_pc <- pte <- res <- sig <- te_direction <- ucl <- value <- NULL


  # Give error if integration is not early, intermediate, or late
  if (!(integration %in% c("early", "intermediate", "late"))) {
    stop("Please provide a valid integration method: early, intermediate, or late")
  }

  # Give error if Y.family and M.family are not either binomial or gaussian
  if (!(Y.family %in% c("binomial", "gaussian") &
        M.family %in% c("binomial", "gaussian"))) {
    stop("Please provide a valid family for Y and M: binomial or gaussian")
  }

  # Give error if covs is NULL
  if (is.null(covs)) {
    stop("Currently, this method does not support analysis without covariates.
         Please provide covariates.")
  }

  # If method is intermediate, make sure that n_boot is not missing and is > 0
  if (integration == "intermediate") {
    if (missing(n_boot)) {
      stop("Please provide a number of bootstrap estimates to perform for calculating the se of the mediation effect")
    }
    if (n_boot <= 0) {
      stop("Please provide a number of bootstrap estimates to perform for calculating the se of the mediation effect")
    }
  }

  # Early integration ---------
  if(integration == "early"){

    # Combines omics data into one dataframe
    omics_lst_df <- purrr::map(omics_lst, ~tibble::as_tibble(.x, rownames = "name"))

    meta_df <- imap_dfr(omics_lst_df,
                        ~tibble(omic_layer = .y, ftr_name = names(.x))) |>
      filter(ftr_name != "name") |>
      mutate(omic_num = case_when(str_detect(omic_layer, "meth") ~ 1,
                                  str_detect(omic_layer, "transc") ~ 2,
                                  str_detect(omic_layer, "miR") ~ 3,
                                  str_detect(omic_layer, "pro") ~ 4,
                                  str_detect(omic_layer, "met") ~ 5))

    # Create data frame of omics data
    omics_df <- omics_lst_df  |>
      purrr::reduce(left_join, by = "name") |>
      column_to_rownames("name")

    # Run hima
    result_hidimum_early <- hima_classic(X = exposure,
                              Y = outcome,
                              M = omics_df,
                              COV.XM = covs,
                              COV.MY = covs,
                              Y.family = Y.family,
                              M.family = M.family,
                              verbose = FALSE,
                              max.iter = 100000,
                              scale = FALSE) |>
      as_tibble(rownames = "ftr_name")

    # Reorders the columns and adds the omics layer information
    result_hidimum_early <- result_hidimum_early |>
      dplyr::mutate(
        multiomic_mthd = "Early Integration",
        mediation_mthd = "HIMA") |>
      dplyr::select("multiomic_mthd", "mediation_mthd",
                    "ftr_name",
                    everything())
    # Filter to significant features only and scale % total effect to 100
    result_hidimum_early <- result_hidimum_early |>
      dplyr::mutate(pte = 100*rimp/sum(rimp),
             sig = if_else(pmax < bh.fdr, 1, 0)) |>
      rename(ie = 'IDE',
             `TME (%)` = pte,
             alpha = alpha_hat,
             beta = beta_hat)

    # Merge results with feature metadata
    result_hidimum_early <- result_hidimum_early |>
      left_join(meta_df, by = "ftr_name") |>
      mutate(integration = integration)


    # Return result
    return(result_hidimum_early)
  }

  # Intermediate integration ---------
  if(integration == "intermediate") {
    ## Change omics elements to dataframes
    omics_lst_df <- purrr::map(omics_lst, ~as_tibble(.x, rownames = "name"))

    meta_df <- imap_dfr(omics_lst_df, ~tibble(omic_layer = .y, ftr_name = names(.x)))|>
      filter(ftr_name != "name") |>
      mutate(omic_num = case_when(str_detect(omic_layer, "meth") ~ 1,
                                  str_detect(omic_layer, "transc") ~ 2,
                                  str_detect(omic_layer, "miR") ~ 3,
                                  str_detect(omic_layer, "pro") ~ 4,
                                  str_detect(omic_layer, "met") ~ 5))

    ## Create data frame of omics data
    omics_df <- omics_lst_df  |>
      purrr::reduce(left_join, by = "name") |>
      column_to_rownames("name")

    # Rename family for xtune function
    if(Y.family != "gaussian") {
      stop("Only continuous outcomes currently supported")
    }

    # Get dataframe of all data
    full_data <- tibble(outcome = outcome,
                        exposure = exposure) |>
      bind_cols(omics_df)

    # Add covs if not null
    if(!is.null(covs)) {full_data <- full_data |> bind_cols(covs)}

    # Get external information matrix
    # Convert each data frame to a long format and extract the unique column names
    # Convert each data frame to a long format and extract the unique column names
    df_list <- purrr::map(omics_lst, function(omic_data) {
      data.frame(column = colnames(omic_data), val = 1)
    })
    # Rename the 'val' column using names from the 'omics_lst'
    named_df_list <- purrr::map2(df_list, names(omics_lst),
                                 function(df, name) {
                                   names(df)[names(df) == "val"] <- name
                                   return(df)
                                 })

    # Reduce the list of data frames into a single data frame with full join
    combined_df <- purrr::reduce(named_df_list, function(df1, df2) {
      full_join(df1, df2, by = "column")
    })
    # Replace NA values with 0
    final_df <- replace(combined_df, is.na(combined_df), 0)
    # Set rownames using "column" column and drop the column
    external_info <- column_to_rownames(final_df, "column")


    ## 0) calculate gamma (x --> y) ----
    if(Y.family == "binary"){
      gamma_est <- coef(glm(outcome ~ exposure + as.matrix(covs),
                            family = binomial))[["exposure"]]
    } else if(Y.family == "gaussian"){
      gamma_est <- coef(lm(outcome ~ exposure ))[["exposure"]]
    }

    # 1) X --> M ----------
    # Model 1: x --> m
    x_m_reg <- epiomics::owas(df = full_data,
                              omics = rownames(external_info),
                              covars = colnames(covs),
                              var = "exposure",
                              var_exposure_or_outcome = "exposure") |>
      dplyr::select("feature_name", "estimate", "se") |>
      dplyr::rename("alpha" = "estimate",
                    "alpha_se" = "se")

    # 2) M--> Y: select features associated with the outcome using group lasso ----
    # x+M-->Y Glasso
    X = as.matrix(full_data[, colnames(omics_df)])
    Y = full_data$outcome
    Z = as.matrix(external_info)
    U = as.matrix(full_data[,"exposure"])
    if(is.null(covs)){
      U = as.matrix(full_data[,"exposure"])
    } else {
      U = as.matrix(full_data[,c(colnames(covs), "exposure")])
    }

    # Run xtune
    invisible(
      capture.output(
        xtune.fit_all_data <- xtune(X = X, Y = Y, Z = Z, U = U,
                                    c = 1,
                                    family = "linear")
      )
    )

    # Extract estimates
    xtune_betas_all_data <- as_tibble(as.matrix(xtune.fit_all_data$beta.est),
                                      rownames = "feature_name") |>
      left_join(meta_df, by = c("feature_name" = "ftr_name")) |>
      dplyr::filter(feature_name %in% colnames(omics_df))

    # 3) Calculate SE for Model 3: x+m to y reg -----
    # 3.1) boot function --------------------------------------------------------
    group_lasso_boot <- function(data, indices, external_info, covs = NULL) {
      X = as.matrix(data)[indices, rownames(external_info)]
      Y = data$outcome[indices]
      if(is.null(covs)){
        U = as.matrix(data[indices,"exposure"])
      } else {
        U = as.matrix(data[indices,c(colnames(covs), "exposure")])
      }
      # Run xtune with error catching function, sometimes xtune needs to run again
      success <- FALSE
      attempts <- 0
      while(!success & attempts < 10) {
        tryCatch({
          # Run xtune
          xtune.fit <- xtune(X = X, Y = Y, Z = as.matrix(external_info), U = U,
                             sigma.square = estimateVariance(X,Y),
                             c = 0,
                             family = "linear", message = FALSE)


          # If the xtune call is successful, proceed with the rest of the code
          # Select betas, drop intercept
          xtune_betas <- as_tibble(as.matrix(xtune.fit$beta.est),
                                   rownames = "feature_name") |>
            dplyr::filter(feature_name %in% colnames(omics_df)) |>
            dplyr::select("s1") |>
            as.matrix()
          # Fix issue where sometimes lasso returns a null matrix
          if(sum(dim(xtune_betas) == c(nrow(external_info), 1))==2){
            return(xtune_betas)
          } else {
            return(as.matrix(rep(0, nrow(external_info))))
          }

          success <- TRUE
        }, error = function(e) {
          attempts <- attempts + 1
          print(paste0("Encountered error: ", e$message))
          if (attempts < 10) {
            print("Retrying...")
          } else {
            print("Maximum number of attempts reached. Exiting...")
            return(NULL)
          }
        })
      }
    }

    # 3.2) Run Bootstrap analysis ----------------------
    # Set up for parallelization
    if (!is.null(n_cores)) {
      # Use the specified ncores value
      num_workers <- n_cores
    } else {
      # Use a default or handle as necessary
      num_workers <- parallel::detectCores()
    }

    # Run Bootstrap
    if(is.null(covs)){
      boot_out <- boot(data = full_data,
                       statistic = group_lasso_boot,
                       R = n_boot,
                       ncpus = num_workers,
                       parallel = "multicore",
                       external_info = external_info)
    } else {
      boot_out <- boot(data = full_data,
                       statistic = group_lasso_boot,
                       R = n_boot,
                       ncpus = num_workers,
                       parallel = "multicore",
                       external_info = external_info,
                       covs = covs)
    }

    # 3.3) Calculate SE for Model 3: x+m to y reg -----
    # Calculate percent of times feature was selected
    glasso_boot_results <- tibble(
      feature_name = rownames(external_info),
      beta_bootstrap = colMeans(replace_na(boot_out$t, 0)),
      beta_se = apply(replace_na(boot_out$t, 0), 2, sd)) |>
      left_join(meta_df, by = c("feature_name" = "ftr_name"))

    # 3.4) Join unpenalized results with glasso results ----
    int_med_coefs <- dplyr::inner_join(xtune_betas_all_data,
                                       glasso_boot_results,
                                       by = c("feature_name",
                                              "omic_layer", "omic_num")) |>
      dplyr::inner_join(x_m_reg, by = "feature_name")

    # Calculate confidence intervals -----
    # mu.x: a1 from reg m = a0 + a1*X
    # mu.y: b2 from reg y = b0 + b1*X + b2*M
    int_med_res <- int_med_coefs |>
      group_by(feature_name) |>
      nest() |>
      mutate(res = purrr::map(data,
                              ~RMediation::medci(mu.x = .x$alpha,
                                                 se.x = .x$alpha_se,
                                                 mu.y = .x$beta_bootstrap,
                                                 se.y = .x$beta_se,
                                                 type = "MC") |>
                                unlist() |> t() |> as_tibble())) |>
      unnest(c(res, data)) |>
      ungroup()

    # Modify results
    intermediate_int_res <- int_med_res |>
      rename_with(~c("lcl", "ucl")[seq_along(.)], tidyr::contains("CI.")) |>
      rename(indirect = "Estimate",
             ind_effect_se = "SE") |>
      mutate(gamma = gamma_est,
             pte = (indirect)/gamma,
             sig = if_else(lcl>0|ucl<0, 1, 0)) |>
      dplyr::select(-contains(" Error"))


    # Filter to significant features only and scale % total effect to 100
    intermediate_int_res <- intermediate_int_res |>
      mutate(pte = 100*pte/sum(pte))

    # Rename feature name
    intermediate_int_res <- intermediate_int_res |>
      dplyr::rename(ftr_name = feature_name,
                    ie = indirect,
                    beta = beta_bootstrap,
                    `TME (%)` = pte) |>
      mutate(integration = integration)

    # Return message if intermediate_int_res has zero rows
    if(nrow(intermediate_int_res) == 0){
      return(message("No significant features identified."))
    } else {
      return(intermediate_int_res)
    }
  }

  # Late Integration ---------------
  if(integration == "late"){
    # Get number of omics layers
    n_omics <- length(omics_lst)
    omics_name <- names(omics_lst)

    # Meta data
    omics_lst_df <- purrr::map(omics_lst, ~as_tibble(.x, rownames = "name"))

    meta_df <- imap_dfr(omics_lst_df, ~tibble(omic_layer = .y, ftr_name = names(.x)))|>
      filter(ftr_name != "name") |>
      mutate(omic_num = case_when(str_detect(omic_layer, "meth") ~ 1,
                                  str_detect(omic_layer, "transc") ~ 2,
                                  str_detect(omic_layer, "miR") ~ 3,
                                  str_detect(omic_layer, "pro") ~ 4,
                                  str_detect(omic_layer, "met") ~ 5))
    # Start the computation
    result_hidimum_late <- vector(mode = "list", length = n_omics)
    for(i in 1:n_omics) {
      # Run HIMA with input data
      result_hidimum_late[[i]] <- hima_classic(X = exposure,
                                    Y = outcome,
                                    M = omics_lst[[i]],
                                    COV.XM = covs,
                                    Y.family = Y.family,
                                    M.family = M.family,
                                    max.iter = 100000,
                                    scale = FALSE) |>
        as_tibble(rownames = "ftr_name")
    }

    # Assign omic names
    names(result_hidimum_late) <- names(omics_lst)

    # Concatenate the resulting data frames
    result_hidimum_late_df <- bind_rows(result_hidimum_late, .id = "omic_layer")

    # Add key details
    result_hidimum_late_df <- result_hidimum_late_df |>
      dplyr::mutate(
        multiomic_mthd = "Late Integration",
        mediation_mthd = "HIMA") |>
      dplyr::select("multiomic_mthd", "mediation_mthd",
                    "omic_layer", "ftr_name",
                    everything())

    # Filter to significant features only and scale % total effect to 100
    result_hidimum_late_df <- result_hidimum_late_df |>
      mutate(pte = 100*rimp/sum(rimp),
             sig = if_else(pmax < bh.fdr, 1, 0)) |>
      rename(ie = IDE,
             `TME (%)` = pte,
             alpha = alpha_hat,
             beta = beta_hat) |>
      mutate(integration = integration)

    # Return the final table
    return(result_hidimum_late_df)
  }
}

