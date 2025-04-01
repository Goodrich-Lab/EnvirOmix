#' Plot for the High Dimensional Multiomic Mediation Analysis
#'
#' Given a tidy dataframe summarizing the results of HIMA analysis
#'
#' @param result_hidimum  an output dataframe from the [hidimum] analysis
#'
#' @return a figure of the result of high dimensional multiomic mediation analysis in given integration
#'
#' @import dplyr
#' @import tidyr
#' @import forcats
#' @import ggplot2
#'
#' @export
#' @examples
#' # Load Example Data
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
#'  # Plot result
#'  fig_early <- plot_hidimum(result_hidimum_early)
#'
plot_hidimum <- function(result_hidimum) {
  name <- ftr_name <- value <- omic_layer <- NULL

  # Pivot longer for figure
  result_hidimum_long <- result_hidimum |>
    rename("Alpha" = "alpha",
           "Beta" = "beta") |>
    pivot_longer(cols = c("Alpha", "Beta","TME (%)"),
                 names_to = "name") |>
    mutate(name = factor(name, levels = c("Alpha", "Beta", "TME (%)")))

  # Plot features
  p <- ggplot(result_hidimum_long,
              aes(x = fct_inorder(ftr_name),
                  y = value,
                  fill = omic_layer)) +
    geom_bar(stat = "identity") +
    facet_grid(name ~ omic_layer,
               scales = "free",
               space = "free_x") +
    scale_fill_brewer(type = "qual", palette = 2) +
    geom_hline(yintercept = 0, linetype = 1, color = "grey50") +
    ylab(NULL) + xlab(NULL) +
    theme(
      strip.background = element_blank(),
      strip.text.x = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
      legend.title = element_blank(),
      legend.position = "bottom", # Place the legend at the bottom
      legend.justification = c(.5, 0))

  return(p)
}
