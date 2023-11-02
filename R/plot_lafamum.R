## ---- plot_med_lf ----
#' Plot result of mediation analysis with latent factors
#'
#' Given a list including two tidy dataframes and one integration type vector
#' summarizing  mediation analysis result
#' function runs plotting and returns a plot
#'
#' @param lafamum_res A list including two tidy dataframes
#' summarizing the results of HIMA analysis
#' and one vector indicating integration type
#' @param n_features Number of top features from each latent factor to be
#' plotted. May not exactly be related to the total number of features plotted
#' due to overlapping features across latent factors.
#' @param plot_type One of "all", "mediation", or "features".
#' If "all" is selected (the default), the mediation bargraph and heatmap of
#' individual features is plotted together. If "mediation" or "features"
#' is selected, only the mediation or the features heatmap is plotted.
#'
#' @return a figure of the result of mediation with latent factors in given integration
#'
#' @import dplyr
#' @importFrom dplyr left_join
#' @importFrom ggplot2 ggplot
#' @importFrom cowplot plot_grid
#'
#' @export
#'
plot_lafamum <- function(lafamum_res,
                         n_features = 10,
                         plot_type = "all") {

  # Make sure that plot_type is one of all, mediation, or features
  if(!plot_type %in% c("all", "mediation", "features")) {
    stop("plot_type must be one of all, mediation, or features")
  }

  # Panel A Bargraph of mediation effects of latent factors --------------
  med_long <- lafamum_res[[1]] %>%
    pivot_longer(cols = c(Alpha, Beta, `TME (%)`))

  # Plot
  panel_a <- ggplot(med_long, aes(x = lf_ordered, y = value)) +
    geom_bar(stat = "identity", fill = "grey50") +
    geom_hline(yintercept = 0) +
    facet_grid(name ~ omic_num, scales = "free", space = "free_x", switch = "y") +
    ggh4x::facetted_pos_scales(
      y = list(name == "Alpha"  ~ scale_y_continuous(
        limits = c(-.45,.45), breaks = c(-0.4, 0, .4)),
        name == "Beta"   ~ scale_y_continuous(
          limits = c(-.45,.45), breaks = c(-0.4, 0, .4)),
        name == "TME (%)" ~ scale_y_continuous(
          limits = c(-1,55), n.breaks = 4))) +
    theme(axis.title = element_blank(),
          strip.placement = "outside",
          strip.text.x = element_blank(),
          strip.text.y = element_text(angle = 0, size = 9),
          axis.text.x = element_blank(),
          strip.background = element_blank(),
          axis.text.y = element_text(size = 8))

  # Panel B: heatmap of correlation of features vs PC's --------------

  # Select features for plots -------
  # Select Features
  lafamum_res_ftrs <- lafamum_res[[2]] %>%
    select(-omic_layer, -omic_num) %>%
    filter(!rowSums(is.na(.)) == ncol(.)) %>%
    column_to_rownames("feature")

  ## Select rows with values in the top 10 of their respective columns
  top <- apply(lafamum_res_ftrs, 2,
               function(x) x %in% tail(sort(abs(x)), n_features))

  ## Filter selected rows
  ftr_cor_sig_lf_top_ft <- lafamum_res[[2]][which(rowSums(top) > 0), ]

  # Pivot longer
  ftr_cor_sig_lf_top_ft_l <-
    ftr_cor_sig_lf_top_ft %>%
    pivot_longer(cols = all_of(lafamum_res[[1]]$lf_num),
                 names_to = "lf_num",
                 values_to = "Correlation")

  # Change features which are in wrong JIVE individual component to zero
  if(lafamum_res[[3]] == "Intermediate") {
    ftr_cor_sig_lf_top_ft_l <- ftr_cor_sig_lf_top_ft_l %>%
      mutate(
        in_ind_omic = str_sub(lf_num, 1, 5) == str_sub(omic_layer, 1, 5),
        Correlation = ifelse(in_ind_omic | str_detect(lf_num, "Joint"),
                             Correlation, NA))
  }

  # Join with long format of mediation effects of latent factors
  # to get the ordered latent factor
  panel_b_dat_top_ft <- left_join(ftr_cor_sig_lf_top_ft_l,
                                  med_long %>%
                                    filter(name == "Alpha") %>%
                                    dplyr::select(lf_num, lf_ordered),
                                  by = "lf_num") %>%
    mutate(omic_num2 = case_when(str_detect(lf_num, "meth") ~ 1,
                                 str_detect(lf_num, "transc") ~ 2,
                                 str_detect(lf_num, "miR") ~ 3,
                                 str_detect(lf_num, "pro") ~ 4,
                                 str_detect(lf_num, "met") ~ 5,
                                 TRUE ~ 0))

  panel_b <- ggplot(data = panel_b_dat_top_ft,
                    aes(y = feature,
                        x = lf_ordered,
                        fill = Correlation)) +
    geom_tile(color = "white") +
    facet_grid(omic_layer ~ omic_num2, scales = "free", space = "free") +
    scale_fill_gradient2(low  = "blue",
                         mid  = "white",
                         high = "red",
                         midpoint = 0,
                         limits = c(-1, 1),
                         breaks = c(-1, 0, 1),
                         na.value = "grey50") +
    theme(
      axis.text.x = element_text(size = 8,angle = 90, hjust = 1, vjust = .5),
      axis.text.y = element_text(size = 8),
      strip.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "none",
      text = element_text(size = 8))

  # Return plots
  if(plot_type == "mediation") return(panel_a)
  if(plot_type == "features") return(panel_b)
  if(plot_type == "all") {
  # Combine Figures
  p <- cowplot::plot_grid(
    NULL, panel_a,  NULL, panel_b,
    ncol = 1, align = "v", axis = "lr",
    rel_heights  = c(.05, .6, .1, 1.75),
    labels = c("a)","", "b) "))
  return(p)
  }
}
