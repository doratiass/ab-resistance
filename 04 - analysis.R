# -----------------------------------------------------------------------------
# Risk Difference Analysis and Visualization for Antibiotic Resistance Development
#
# This script performs a full analysis of risk differences (RD) and odds ratios (OR)
# for various antibiotic resistance outcomes. It includes:
#
#   1. The `rd_calc` function which:
#      - Iterates through a matched dataset.
#      - Fits logistic regression models (both multivariable and univariable) for each outcome.
#      - Computes average comparisons (OR and RD) using the `marginaleffects` package.
#      - Aggregates results (including confidence intervals, p-values, and number needed to harm)
#        into a summary table and saves it as a CSV file.
#
#   2. The `rd_fig` function which:
#      - Filters the summary results to a specific reference antibiotic group.
#      - Generates a ggplot figure illustrating the risk differences (with error bars)
#        for each antibiotic resistance outcome.
#      - Adds additional annotations (labels, segments, and custom axis breaks).
#      - Saves the final plot as an SVG file.
#
# Dependencies:
#   - tidyverse: Data manipulation and plotting
#   - foreach: Parallel processing for model fitting
#   - WeightIt: Weight estimation for censoring
#   - sandwich: Robust covariance estimation
#   - marginaleffects: Average marginal effects and comparisons
#   - RColorBrewer, ggbreak, patchwork: Enhanced plotting and layout features
#
# -----------------------------------------------------------------------------

# Load required packages
library(tidyverse)
library(foreach)
library(WeightIt)
library(sandwich)
library(marginaleffects)
library(RColorBrewer)
library(ggbreak)
library(patchwork)

# -----------------------------------------------------------------------------
# Risk Difference Calculation Function
# -----------------------------------------------------------------------------

rd_calc <- function(df,
                    outcomes = c("outcome_sample_ab_res_cef_1", "outcome_sample_ab_res_cef_2", 
                                 "outcome_sample_ab_res_cef_3", 
                                 "outcome_sample_ab_res_augmentin",
                                 "outcome_sample_ab_res_fluoroquinolones", 
                                 "outcome_sample_ab_res_nitrofurantoin")) {
  # Expand the outcomes into a long format for iteration
  matched_df <- df %>%
    mutate(out = list(outcomes)) %>%
    unnest_longer(out)
  
  # Set up a parallel processing cluster with 32 workers
  my.cluster <- parallel::makeCluster(
    32, 
    type = "PSOCK",
    outfile = ""
  )
  doParallel::registerDoParallel(cl = my.cluster)
  
  # Iterate over each matched cohort using foreach in parallel
  summary <- foreach(i = 1:length(matched_df$ref_ab), 
                     .combine = rbind,
                     .packages = c("tidyverse", "marginaleffects")) %dopar% {
                       # Identify subclass values that are missing the outcome of interest for exclusion
                       sub_r <- matched_df$df[[i]][is.na(matched_df$df[[i]][, matched_df$out[[i]]]), "subclass", drop = TRUE]
                       
                       # Subset the data:
                       # - Remove columns starting with "outcome_sample" and "start_sample"
                       # - Exclude id, distance, weights, and ab_sample_diff columns
                       # - Retain the outcome variable for analysis
                       df <- matched_df$df[[i]] %>%
                         select(-starts_with("outcome_sample"),
                                -starts_with("start_sample"),
                                -c(id, distance, weights, ab_sample_diff),
                                all_of(matched_df$out[[i]])) %>%
                         # Exclude rows with subclass values identified earlier
                         filter(!(subclass %in% sub_r)) %>%
                         # Remove columns with only a single unique value (no variability)
                         select(-where(~ n_distinct(.) == 1))
                       
                       # Define the multivariable regression formula for the outcome
                       eq <- paste0(matched_df$out[[i]], " ~ . - subclass - censor_weights")
                       # Define the univariable regression formula using only the exposure (start_ab)
                       eq_uni <- paste0(matched_df$out[[i]], " ~ start_ab")
                       
                       # Fit the multivariable logistic regression model
                       model <- glm(as.formula(eq), 
                                    data = df,
                                    family = binomial)
                       
                       # Fit the univariable logistic regression model (for comparison)
                       model_uni <- glm(as.formula(eq_uni), 
                                        data = df,
                                        family = binomial)
                       
                       # Calculate average comparisons (Odds Ratio) using the multivariable model.
                       # The vcov argument is set to cluster on 'subclass' and censor_weights are applied.
                       or <- avg_comparisons(model, 
                                             variables = list(start_ab = c(matched_df$ref_ab[[i]],
                                                                           matched_df$match_ab[[i]])),
                                             comparison = "lnor",
                                             transform = "exp",
                                             vcov = ~subclass,
                                             wts = "censor_weights")
                       
                       # Calculate average comparisons (Risk Difference) using the multivariable model.
                       rd <- avg_comparisons(model, 
                                             variables = list(start_ab = c(matched_df$ref_ab[[i]],
                                                                           matched_df$match_ab[[i]])),
                                             comparison = "difference",
                                             vcov = ~subclass,
                                             wts = "censor_weights")
                       
                       # Calculate average comparisons (Risk Difference) using the univariable model.
                       rd_uni <- avg_comparisons(model_uni, 
                                                 variables = list(start_ab = c(matched_df$ref_ab[[i]],
                                                                               matched_df$match_ab[[i]])),
                                                 comparison = "difference",
                                                 vcov = ~subclass,
                                                 wts = "censor_weights")
                       
                       # Create a summary table containing the effect estimates and additional statistics.
                       sum_table <- tibble(
                         time = matched_df$follow_up[[i]],
                         ref_ab = matched_df$ref_ab[[i]],
                         match_ab = matched_df$match_ab[[i]],
                         outcome = matched_df$out[[i]],
                         drop_ref = matched_df$ref_ab_drop[[i]],
                         drop_match = matched_df$match_ab_drop[[i]],
                         n = nrow(df) / 2,  # Number of observations per group
                         OR = or$estimate,
                         LL_r = or$conf.low,
                         UP_r = or$conf.high,
                         pVal_r = round(or$p.value, 2),
                         RD = rd$estimate,
                         RD_std = rd$std.error,
                         LL_d = rd$conf.low,
                         UP_d = rd$conf.high,
                         pVal_d = round(rd$p.value, 2), 
                         NNH = 1 / rd$estimate,
                         RD_uni = rd_uni$estimate,
                         RD_uni_std = rd_uni$std.error,
                         LL_d_uni = rd_uni$conf.low,
                         UP_d_uni = rd_uni$conf.high,
                         pVal_d_uni = round(rd_uni$p.value, 2), 
                         NNH_uni = 1 / rd_uni$estimate
                       )
                       
                       return(sum_table)
                     }
  
  # Write the summary table to a CSV file for export.
  write_csv(summary, 
            file = file.path("export", "rd_summary", 
                             paste0(str_remove(deparse(substitute(df)), "matched_df_"), ".csv")))
  
  # Assign the summary table to a variable in the global environment.
  assign(paste0("rd_summary_", str_remove(deparse(substitute(df)), "matched_df_")), 
         summary, envir = .GlobalEnv)
}

# -----------------------------------------------------------------------------
# Figure 1 Generation Function
# -----------------------------------------------------------------------------

rd_fig <- function(df,
                   ref_abs = "Ceph 1gen",
                   univar = FALSE,
                   x = "RD", ul = "UP_d", ll = "LL_d",
                   p_size = 2,
                   l_size = 1,
                   x_text_ref = -5,
                   x_text_match = 5,
                   y_text = 6.5,
                   arrow_size = 2,
                   color_pal = "Dark2") {
  # If univariate results are requested, update column names accordingly.
  if (univar) {
    x <- paste0(x, "_uni")
    ul <- paste0(ul, "_uni")
    ll <- paste0(ll, "_uni")
  }
  
  # (Optional) Define caption and title texts.
  cap <- ""
  plot_title <- ""
  
  # Prepare the data for plotting:
  # - Filter for the reference antibiotic group.
  # - Multiply the effect estimates and CI bounds by 100 (to convert to percentages).
  # - Create labels for exposures and outcomes.
  df %>%
    filter(ref_ab == ref_abs) %>%
    mutate_at(c(x, ul, ll), ~.*100) %>%
    mutate(
      ref_label = paste("Exposure to", 
                        case_when(
                          ref_ab == "Ceph 1gen" ~ "Cephalosporins 1st generation",
                          ref_ab == "Ceph 2gen" ~ "Cephalosporins 2nd generation",
                          TRUE ~ ref_ab)),
      match_label = paste("Exposure to", 
                          case_when(
                            match_ab == "Ceph 1gen" ~ "Cephalosporins 1st generation",
                            match_ab == "Ceph 2gen" ~ "Cephalosporins 2nd generation",
                            match_ab == "Ceph 3gen" ~ "Cephalosporins 3rd generation",
                            TRUE ~ match_ab)),
      match_label_fav = paste("Favors", 
                              case_when(
                                ref_ab == "Ceph 1gen" ~ "Cephalosporins 1st generation",
                                ref_ab == "Ceph 2gen" ~ "Cephalosporins 2nd generation",
                                TRUE ~ ref_ab)),
      ref_label_fav = paste("Favors", 
                            case_when(
                              match_ab == "Ceph 1gen" ~ "Cephalosporins 1st generation",
                              match_ab == "Ceph 2gen" ~ "Cephalosporins 2nd generation",
                              match_ab == "Ceph 3gen" ~ "Cephalosporins 3rd generation",
                              TRUE ~ match_ab)),
      x_text_ref = x_text_ref,
      x_text_match = x_text_match,
      y_text = y_text,
      outcome = str_remove_all(outcome, "outcome_sample_ab_res_"),
      outcome = factor(case_when(
        outcome == "cef_1" ~ "Cephalosporins 1st generation",
        outcome == "cef_2" ~ "Cephalosporins 2nd generation",
        outcome == "cef_3" ~ "Cephalosporins 3rd generation",
        outcome == "augmentin" ~ "Augmentin",
        outcome == "fluoroquinolones" ~ "Fluoroquinolones",
        outcome == "nitrofurantoin" ~ "Nitrofurantoin"
      ),
      levels = c("Cephalosporins 1st generation",
                 "Cephalosporins 2nd generation",
                 "Cephalosporins 3rd generation",
                 "Augmentin", "Fluoroquinolones", "Nitrofurantoin"))
    ) %>%
    ggplot(aes(!!sym(x), outcome, color = outcome)) +
    # Plot the point estimates with error bars
    geom_pointrange(aes(xmin = !!sym(ll), xmax = !!sym(ul)), linewidth = l_size) +
    # Add a vertical reference line at zero
    geom_vline(aes(xintercept = 0)) +
    # Expand the y-axis limits
    expand_limits(y = c(0, 8)) +
    # Customize x-axis scale and labels
    scale_x_continuous(limits = c(-10, 30),
                       breaks = c(-10, -5, 0, 5, 10, 15, 20, 25, 30),
                       labels = c("10", "5", "0", "5", "10", "15", "20", "25", "30")) +
    # Break the x-axis between 10 and 15 for better visualization
    scale_x_break(c(10, 15), scales = 0.1, expand = FALSE) +
    # Add arrows to indicate directionality
    geom_segment(
      x = 0, y = 7,
      xend = -9, yend = 7,
      lineend = "round",
      linejoin = "mitre",
      linewidth = arrow_size, 
      arrow = arrow(length = unit(arrow_size / 10, "cm"), type = 'closed'),
      colour = "snow3"
    ) + 
    geom_segment(
      x = 0, y = 7,
      xend = 9, yend = 7,
      lineend = "round",
      linejoin = "mitre",
      linewidth = arrow_size, 
      arrow = arrow(length = unit(arrow_size / 10, "cm"), type = 'closed'),
      colour = "snow4"
    ) + 
    # Add custom text labels for the reference and matching groups
    geom_text(aes(x_text_ref, y_text, label = ref_label),
              inherit.aes = FALSE, check_overlap = TRUE) +
    geom_text(aes(x_text_match, y_text, label = match_label),
              inherit.aes = FALSE, check_overlap = TRUE) +
    geom_text(aes(x_text_ref, 0.5, label = ref_label_fav),
              inherit.aes = FALSE, check_overlap = TRUE) +
    geom_text(aes(x_text_match, 0.5, label = match_label_fav),
              inherit.aes = FALSE, check_overlap = TRUE) +
    labs(x = "Risk Difference in Antibiotic Resistance Development (%)",
         y = "",
         color = "",
         title = "Future resistance to:",
         caption = str_wrap(cap, width = 150)) +
    facet_grid(match_ab ~ .) +
    theme_classic() +
    scale_color_brewer(palette = color_pal) + 
    theme(legend.position = "none",
          axis.text.y = element_text(size = 12),
          plot.title = element_text(hjust = 0.05, vjust = -10),
          axis.text.x = element_text(size = 10),
          axis.text.x.top = element_blank(),  
          axis.line.x.top = element_blank(),  
          axis.ticks.x.top = element_blank(),  
          strip.background = element_blank(),  
          strip.text = element_blank(),
          plot.caption = element_text(size = 12, hjust = 0))
  
  # Save and print the plot
  -> p
  
  print(p)
  
  # Generate a standardized file name for the plot based on input parameters
  plot_name <- str_to_lower(gsub("[^A-Za-z0-9_]", "_",
                                 paste("fig_1", ref_abs,
                                       str_remove(deparse(substitute(df)), "rd_summary_"),
                                       sep = "_")))
  file_name <- paste0(plot_name, ".svg")
  
  # Assign the plot object to the global environment
  assign(plot_name, p, envir = .GlobalEnv)
  
  # Save the plot as an SVG file in the export/figs directory
  ggsave(plot = p, width = 10, height = 8,
         filename = file.path("export", "figs", file_name))
}