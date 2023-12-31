#' Plot for the High Dimensional Multiomic Mediation Analysis
#'
#' Given a tidy dataframe summarizing the results of HIMA analysis
#'
#' @param result_hidimum  an output dataframe from the hidimum analysis
#'
#' @return a ggplot figure.
#'
#' @import dplyr
#' @import tidyr
#' @import forcats
#' @import ggplot2
#'
#' @export
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
