# -----------------------------------------------------------------------------
# Matching Pipeline for Antibiotic Exposure Cohorts
#
# This script provides functions for performing propensity score matching 
# on antibiotic exposure cohorts and generating summary statistics and plots.
#
# The workflow includes:
#   1. Class functions to print, summarize, and plot matched cohort objects.
#   2. A matching function (match_ab) that performs propensity score matching 
#      based on specified exposure groups (e.g., different antibiotics) using 
#      the MatchIt package.
#   3. Helper functions nested within match_ab to generate file names, plot 
#      propensity score distributions and standardized mean differences (SMDs), 
#      and summarize inverse probability of treatment weights (IPTW).
#   4. A function (create_matched_df) to aggregate matched cohort results into 
#      a summary data frame.
#
# Dependencies: tidyverse, gtsummary, smd, MatchIt, survey
# Additional user functions are sourced from "Rscripts/05 - funcs.R"
#
# -----------------------------------------------------------------------------

# Load required packages
library(tidyverse)
library(gtsummary)
library(smd)
library(MatchIt)
library(survey)
source("Rscripts/05 - funcs.R")  # Load additional helper functions

# -----------------------------------------------------------------------------
# Class Functions for Matched Cohort Objects
# -----------------------------------------------------------------------------

# 'print.matched_cohort': Print method for a matched cohort object.
#   Displays follow-up days, control and match antibiotic groups, and sample count.
print.matched_cohort <- function(x) {
  cat(x$follow_up, "follow up days cohort\n")
  cat("Control AB:", x$ref_ab, "Match AB:", x$match_ab, "\n")
  cat("Contains", nrow(x$df), "samples\n")
}

# 'summary.matched_cohort': Summary method for a matched cohort object.
#   Prints a summary of the matching procedure including drop information.
summary.matched_cohort <- function(x) {
  cat("Match summary for", x$follow_up, "follow up days cohort\n")
  cat("Control AB:", x$ref_ab, "\n", x$match_summary$drop_txt[1], "\n")
  cat("Match AB:", x$match_ab, "\n", x$match_summary$drop_txt[2], "\n")
}

# 'plot.matched_cohort': Plot method for a matched cohort object.
#   Displays the standardized mean difference (SMD) plot.
plot.matched_cohort <- function(x) {
  print(x$smd_plot)
}

# -----------------------------------------------------------------------------
# Matching Function: match_ab
# -----------------------------------------------------------------------------

# 'match_ab': Performs propensity score matching for antibiotic exposure cohorts.
#   - cohort: The input cohort object with pre-matching data.
#   - short: If TRUE, prints a brief summary and SMD match plot; otherwise,
#            generates full plots and tables for pre- and post-match comparisons.
#   - caliper: Caliper for matching.
#   - std_caliper: Logical; if TRUE, uses standardized caliper.
#   - ref_ab: Control antibiotic group (must be one of the allowed values).
#   - match_ab: Matching antibiotic group (must be one of the allowed values).
#   - outcome_var: Outcome variable for exact matching; defaults to a ceftriaxone result.
#   - exact_vars: Additional variables to exactly match on.
#   - fig_w, fig_h: Figure width and height for saved plots.
match_ab <- function(cohort,
                     short = TRUE,
                     caliper = 0.2, std_caliper = FALSE,
                     ref_ab = c("Nitrofurantoin", "Augmentin", "Ceph 1gen",
                                "Ceph 2gen", "Fluoroquinolones", "Macrolides"),
                     match_ab = c("Nitrofurantoin", "Augmentin", "Ceph 1gen",
                                  "Ceph 2gen", "Fluoroquinolones", "Macrolides"),
                     outcome_var = NULL,
                     exact_vars = NULL,
                     fig_w = 30,
                     fig_h = 15
) {
  # Choose the control and match antibiotics from allowed values
  ref_ab <- match.arg(ref_ab)
  match_ab <- match.arg(match_ab)
  
  # Extract follow-up days from the cohort object
  follow_up <- cohort$follow_up
  
  # ---------------------------------------------------------------------------
  # Nested Helper Functions within match_ab
  # ---------------------------------------------------------------------------
  
  # 'create_file_name': Creates a standardized file name using follow-up, 
  # control AB, matching AB, a file identifier (x), and a suffix.
  create_file_name <- function(x, surfix,
                               fu = follow_up, ab_1 = ref_ab, 
                               ab_2 = match_ab) {
    ab_1 <- str_to_lower(gsub("[^A-Za-z0-9_]", "_", ab_1))
    ab_2 <- str_to_lower(gsub("[^A-Za-z0-9_]", "_", ab_2))
    return(paste0(fu, "_", ab_1, "_", ab_2, "_", x, surfix))
  }
  
  # 'ps_plot': Generates density plots of the propensity score distributions
  # before and after matching with a detailed caption.
  ps_plot <- function(df, df_matched, text) {
    cap <- paste0("Propensity Score Distributions Before (left panel) and After (right panel) Matching for Patients Receiving ",
                  text[1, "start_ab", drop = TRUE], " n = (", text[1, "all", drop = TRUE],
                  ") in red and ", text[2, "start_ab", drop = TRUE], " n = (" , text[1, "all", drop = TRUE],
                  ") in blue. After matching, ",
                  text[1, "dropped", drop = TRUE], " of ", text[1, "start_ab", drop = TRUE],
                  " are dropped compared to ", text[2, "dropped", drop = TRUE], " of ", text[2, "start_ab", drop = TRUE],
                  ". A total of ", text[1, "matched", drop = TRUE], " patients were matched in each group.")
    
    bind_rows(df %>% mutate(matched = "Pre match"),
              df_matched %>% mutate(matched = "Post match")) %>%
      ggplot(aes(distance, fill = start_ab)) +
      geom_density(alpha = 0.4) +
      geom_rug(aes(distance, color = start_ab),
               data = df_matched %>% filter(is.na(subclass)),
               alpha = 0.4) +
      facet_grid(~fct_relevel(matched, "Pre match", "Post match")) + 
      labs(
        x = "Propensity score",
        y = "Density",
        fill = "",
        caption = str_wrap(cap, width = fig_w * 5)
      ) +
      theme_classic() +
      theme(
        legend.position = "bottom",
        legend.title = element_blank()
      )
  }
  
  # 'smd_plot': Creates a plot of the absolute standardized mean differences 
  # (SMD) before and after matching.
  smd_plot <- function(df1, df2, line_size = 0.5) {
    bind_rows(
      df1$table_body %>% transmute(variable, estimate = abs(estimate), match = "pre-match"),
      df2$table_body %>% transmute(variable, estimate = abs(estimate), match = "post-match")
    ) %>%
      mutate(n = c(1:(nrow(.)/2), 1:(nrow(.)/2))) %>%
      drop_na() %>%
      filter(!str_detect(variable, "outcome|ab_sample_diff")) %>%
      ggplot(aes(estimate, reorder(variable, desc(n)), color = match)) +
      labs(
        y = "",
        x = "Absolute standardized\nmean difference",
        color = ""
      ) +
      geom_point() +
      geom_vline(xintercept = 0, linewidth = line_size) +
      geom_vline(xintercept = 0.05, linetype = "dashed", linewidth = line_size) +
      geom_vline(xintercept = 0.1, linewidth = line_size) +
      scale_y_discrete(labels = vars_label) +
      theme(legend.position = "bottom")
  }
  
  # 'iptw_table': Generates a table summarizing the IPTW-adjusted covariate balance 
  # using survey-weighted summaries.
  iptw_table <- function(df) {
    svydesign(~0, data = df %>% select(-c(id, subclass, distance, weights)), 
              weights = ~censor_weights) %>%
      tbl_svysummary(
        by = "start_ab",
        missing = "no",
        type = list(
          c(starts_with("start_sample"), "sex") ~ "dichotomous"
        ),
        value = list(
          sex ~ "Female",
          c(starts_with("start_sample_ab")) ~ "Resistant",
          c(starts_with("start_sample_gram")) ~ "Found",
          c(starts_with("start_sample_mic")) ~ "Present"
        ),
        include = -c(censor_weights)
      ) %>%
      add_n(statistic = "{n_miss} ({p_miss})") %>%
      add_overall() %>%
      add_difference(everything() ~ "smd") %>%
      modify_column_hide(ci) %>%
      modify_header(estimate = "**SMD**", n = "**Missing**") -> tbl
    
    tbl$table_body$label <- sapply(tbl$table_body$label, label_get, USE.NAMES = FALSE)
    
    return(tbl)
  }
  
  # ---------------------------------------------------------------------------
  # Create Antibiotic-Specific Data Frame for Matching
  # ---------------------------------------------------------------------------
  
  # Filter the pre-match data to include only records with the specified control
  # (ref_ab) and match (match_ab) antibiotic exposures.
  pre_match_ab <- cohort$pre_match_df %>%
    filter(start_ab %in% c(ref_ab, match_ab)) %>%
    mutate(start_ab = factor(start_ab, levels = c(ref_ab, match_ab))) %>%
    drop_na()
  
  # ---------------------------------------------------------------------------
  # Define Exact Matching Formula
  # ---------------------------------------------------------------------------
  
  # If an outcome variable is provided, include it in the exact matching formula.
  if (!is.null(outcome_var)) {
    exact_formula <- paste("~", outcome_var)
  } else {
    exact_formula <- "~ start_sample_ab_res_cef_3"
  }
  
  # If additional exact matching variables are specified, append them to the formula.
  if (!is.null(exact_vars)) {
    exact_formula <- as.formula(paste(exact_formula, "+", paste(exact_vars, collapse = "+")))
  } else {
    exact_formula <- as.formula(exact_formula)
  }
  
  # ---------------------------------------------------------------------------
  # Perform Propensity Score Matching using MatchIt
  # ---------------------------------------------------------------------------
  
  ab_match <- matchit(start_ab ~ . + poly(age, 3) - age - censor_weights,
                      data = pre_match_ab %>% select(-starts_with("outcome_sample"), -c(id, ab_sample_diff)),
                      estimand = "ATC",
                      method = "nearest",
                      distance = "glm",
                      exact = exact_formula,
                      caliper = caliper,
                      std.caliper = std_caliper)
  
  # Summarize the matching result
  sum_ab_match <- summary(ab_match)
  
  # Extract matched data with and without dropping unmatched records
  matched_df_all <- match.data(ab_match, data = pre_match_ab, drop.unmatched = FALSE)
  matched_df <- match.data(ab_match, data = pre_match_ab)
  
  # ---------------------------------------------------------------------------
  # Generate Matching Summary and Plots
  # ---------------------------------------------------------------------------
  
  # Create a summary table with match counts and drop percentages.
  plot_txt_ab <- tibble(
    start_ab = c(ref_ab, match_ab),
    all = c(sum_ab_match$nn[2, 1], sum_ab_match$nn[2, 2]),
    dropped = c(
      paste0(sum_ab_match$nn[5, 1], " (", round(100 * sum_ab_match$nn[5, 1] / sum_ab_match$nn[2, 1], 2), "%)"),
      paste0(sum_ab_match$nn[5, 2], " (", round(100 * sum_ab_match$nn[5, 2] / sum_ab_match$nn[2, 2], 2), "%)")
    ),
    matched = c(sum_ab_match$nn[4, 1], sum_ab_match$nn[4, 2]),
    drop_txt = case_when(
      start_ab == ref_ab ~ paste0(ref_ab, "\nAll ", sum_ab_match$nn[2, 1], "\nMatched ", sum_ab_match$nn[4, 1],
                                  "\nUnmatched ", sum_ab_match$nn[5, 1], " (", round(100 * sum_ab_match$nn[5, 1] / sum_ab_match$nn[2, 1], 2), "%)"),
      start_ab == match_ab ~ paste0(match_ab, "\nAll ", sum_ab_match$nn[2, 2], "\nMatched ", sum_ab_match$nn[4, 2],
                                    "\nUnmatched ", sum_ab_match$nn[5, 2], " (", round(100 * sum_ab_match$nn[5, 2] / sum_ab_match$nn[2, 2], 2), "%)")
    )
  )
  
  # If short matching summary is requested, simply print the summary and SMD plot.
  if (short) {
    print(plot_txt_ab)
    plot(sum_ab_match)
  } 
  # Otherwise, generate full match plots and tables.
  else {
    # Generate and save the propensity score distribution plot.
    ps_match_plot <- ps_plot(matched_df_all, matched_df, plot_txt_ab)
    
    ggsave(filename = file.path("export", deparse(substitute(cohort)), create_file_name("ps_plot_ab", ".svg")), 
           plot = ps_match_plot,
           width = fig_w, height = fig_h, dpi = 300, units = "cm", bg = "white")
    
    # Generate IPTW summary tables for pre-match and post-match data.
    iptw_pre_match_t_ab <- iptw_table(matched_df_all)
    iptw_post_match_t_ab <- iptw_table(matched_df)
    
    # Generate and save the standardized mean difference (SMD) plot.
    smd_plot <- smd_plot(iptw_pre_match_t_ab, iptw_post_match_t_ab)
    
    ggsave(filename = file.path("export", deparse(substitute(cohort)), create_file_name("smd_ab", ".svg")), 
           plot = ggplot2::last_plot(), 
           width = fig_w, height = fig_h, dpi = 300, units = "cm", bg = "white")
    
    # Merge pre-match and post-match IPTW tables and save as HTML.
    comp_iptw_ab <- tbl_merge(
      tbls = list(iptw_pre_match_t_ab, iptw_post_match_t_ab),
      tab_spanner = c("**Pre-Match**", "**Post-Match**")
    )
    
    gt::gtsave(as_gt(comp_iptw_ab), 
               file = file.path("export", deparse(substitute(cohort)), create_file_name("comp_iptw_ab", ".html")))
    
    # Create the matched cohort object with all relevant components.
    matched_cohort <- list(
      follow_up = follow_up, 
      ref_ab = ref_ab, 
      match_ab = match_ab, 
      ab_match = ab_match, 
      matched_df = matched_df, 
      match_summary = plot_txt_ab,
      ps_match_plot = ps_match_plot, 
      smd_plot = smd_plot, 
      match_table = comp_iptw_ab
    )
    
    class(matched_cohort) <- "matched_cohort"
    
    # Assign the matched cohort object to the global environment using a standardized name.
    assign(str_to_lower(gsub("[^A-Za-z0-9_]", "_", paste0("matched_", deparse(substitute(cohort)), "_", ref_ab, "_", match_ab))), 
           matched_cohort, envir = .GlobalEnv)
  }
}

# -----------------------------------------------------------------------------
# Function to Create a Summary Data Frame from Matched Cohort Objects
# -----------------------------------------------------------------------------

# 'create_matched_df': Aggregates information from a list of matched cohort objects
# into a single tibble containing follow-up days, counts of subjects, and drop percentages.
create_matched_df <- function(matched_list) {
  # For each matched cohort, compute the summary of matching (using MatchIt summary).
  sum_list <- lapply(matched_list, function(x) summary(x[[4]]))
  
  matched_df <- tibble(
    follow_up = sapply(matched_list, function(x) x[[1]]),
    df = lapply(matched_list, function(x) x[[5]]),
    ref_ab = sapply(matched_list, function(x) x[[2]]),
    ref_ab_total = sapply(sum_list, function(x) {x$nn[2, 1]}),
    ref_ab_matched = sapply(sum_list, function(x) {x$nn[4, 1]}),
    ref_ab_drop = 100 * (ref_ab_total - ref_ab_matched) / ref_ab_total,
    match_ab = sapply(matched_list, function(x) x[[3]]),
    match_ab_total = sapply(sum_list, function(x) {x$nn[2, 2]}),
    match_ab_matched = sapply(sum_list, function(x) {x$nn[4, 2]}),
    match_ab_drop = 100 * (match_ab_total - match_ab_matched) / match_ab_total
  )
  
  return(matched_df)
}