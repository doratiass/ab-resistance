# -----------------------------------------------------------------------------
# Cohort Creation and Analysis Pipeline for Antibiotic Exposure Studies
#
# This script defines a set of functions for managing and analyzing cohorts
# in studies investigating antibiotic exposure and microbiology outcomes.
#
# The pipeline includes:
#   1. Utility functions to print and summarize cohort objects.
#   2. A function to generate an example flowchart using DiagrammeR.
#   3. A comprehensive function (create_cohort) that processes raw data 
#      (antibiotic consumption, microbiology, hospitalization, and covariates) 
#      to assemble a final cohort. This includes applying exclusion criteria, 
#      matching microbiology and antibiotic dates, censoring based on hospitalization, 
#      and generating summary tables.
#
# Dependencies: tidyverse, gtsummary, WeightIt, foreach, DiagrammeR
#
# -----------------------------------------------------------------------------

# Load required packages
library(tidyverse)
library(gtsummary)
library(WeightIt)
library(foreach)
library(DiagrammeR)

# -----------------------------------------------------------------------------
# Utility Functions for Cohort Object
# -----------------------------------------------------------------------------

# 'print.cohort': Prints key information about a cohort (follow-up days and sample count)
print.cohort <- function(x) {
  cat(x$follow_up, "follow up days cohort\n")
  cat("Contains", nrow(x$df), "samples\n")
}

# 'summary.cohort': Prints a summary of the cohort details
summary.cohort <- function(x) {
  cat("Cohort summary for", x$follow_up, "follow up days\n")
  paste0(x$summary[1], x$summary[2],
         x$summary[3], x$summary[4],
         sep = "\n")
}

# -----------------------------------------------------------------------------
# Function to Generate a Visual Flowchart
# -----------------------------------------------------------------------------

# 'ex_graph': Creates an example flowchart using DiagrammeR to visualize cohort summary steps
ex_graph <- function(x) {
  grViz("digraph flowchart {

  node [fontname = Helvetica, shape = box, style = filled, fillcolor = lightblue]
  
  A [label = '@@1-1']
  int_1 [label = '']
  B [label = '@@1-2']
  C [label = '@@1-3']
  int_2 [label = '@@1-4']
  D [label = '@@2-1']
  E [label = '@@2-2']
  F [label = '@@2-3']
  G [label = '@@2-4']

  A -> int_1 -> C
  int_1 -> B
  C -> int_2 -> {D E F G}
  
  subgraph {rank=same; int_1; B}
  }
  
  [1]: c(x$summary[1], x$summary[2], x$summary[3], x$summary[4])
  [2]: c(x$summary[5], x$summary[6], x$summary[7], x$summary[8])")
}

# -----------------------------------------------------------------------------
# Data Creation and Cohort Assembly Function
# -----------------------------------------------------------------------------

# 'create_cohort': Generates a cohort object from raw data. The function:
#   - Applies exclusion criteria based on antibiotic washout and dual usage.
#   - Processes microbiology data (filtering by sample source, gram stain, etc.).
#   - Matches antibiotic and microbiology dates.
#   - Applies hospitalization-based censoring.
#   - Generates summary tables (Table 1 and Table 2) using gtsummary.
#   - Returns a cohort object with the final dataset and summary information.
create_cohort <- function(ab_washout = 90,
                          follow_up = 182,
                          ab_consp_all = ab_consumption_all,
                          micro = micro_df,
                          hosp_all = hospital_all,
                          outcome_abs = c("Cef 3"),
                          sample_sources = c("Urine"),
                          micro_names = NULL,
                          micro_grams = c("Negative"),
                          obj_name = NULL) {
  
  # Calculate the antibiotic washout date
  ab_washout_date <- as.Date("2017-01-01") + ab_washout
  
  # Create standardized outcome column name (e.g., outcome_sample_ab_res_cef_3)
  outcome_abs <- paste0("outcome_sample_ab_res_", str_to_lower(gsub("[^A-Za-z0-9_]", "_", outcome_abs)))
  
  # Generate a default object name if not provided
  obj_name <- ifelse(is.null(obj_name),
                     paste0("cohort_", follow_up, "_", str_to_lower(gsub("[^A-Za-z0-9_]", "_", paste(micro_names, collapse = "_")))),
                     obj_name)
  
  message(paste("Creating cohort for", follow_up, "follow up days"))
  
  # ---------------------------------------------------------------------------
  # Nested Functions for Data Processing within create_cohort
  # ---------------------------------------------------------------------------
  
  # 'sum_ab': Summarizes antibiotic consumption data ensuring sufficient washout between doses
  sum_ab <- function(ab_consumption) {
    ab_summary <- ab_consumption %>%
      group_by(id) %>%
      arrange(id, ab_date) %>%
      mutate(before = lag(ab_date),
             diff = difftime(ab_date, before, units = "days"),
             date_ok = (is.na(diff) | (diff > ab_washout))) %>%
      ungroup() %>%
      mutate(sumVar = rowSums(select(., contains("ab_first")))) %>%
      filter(sumVar == 1,
             date_ok,
             ab_date >= ab_washout_date) %>%
      select(-sumVar)
    
    return(ab_summary)
  }
  
  # 'return_ab_date': Returns valid antibiotic dates that are between 14 days and the follow-up period
  return_ab_date <- function(date, times) {
    df <- tibble(micro = date,
                 ab = times) %>%
      mutate(diff = difftime(micro, ab, units = "days"),
             date_ok = ((diff > 14) & (diff < follow_up))) %>%
      filter(date_ok)
    if(nrow(df) > 0) {
      return(df$ab)
    } else {
      return(NA)
    }
  }
  
  # 'return_hosp_date': Checks hospitalization timing relative to ab/micro dates
  return_hosp_date <- function(start_date, end_date, hosp_start, hosp_end) {
    df <- tibble(start_date = start_date,
                 end_date = end_date,
                 hosp_start = hosp_start,
                 hosp_end = hosp_end) %>%
      mutate(diff_ab_end = difftime(start_date, hosp_end, units = "days"),
             diff_mic_start = difftime(end_date, hosp_start, units = "days"),
             date_not_ok = !(diff_ab_end > 90 | diff_mic_start < 2)) 
    return(!any(df$date_not_ok, na.rm = TRUE))
  }
  
  # 'create_new_out_ab': Ensures outcome columns exist in the data; if missing, creates them as "Not-tested"
  create_new_out_ab <- function(df, abs = c("cephalexin", "cefazolin", "cephalothin",
                                            "cefuroxime", "cefoxitin",
                                            "ofloxacin", "ciprofloxacin", "levofloxacin",
                                            "nitrofurantoin", "gentamicin", "amoxicillin_clavul_a")) {
    abs_full <- paste0("outcome_sample_ab_res_", abs)
    for (i in 1:length(abs_full)) {
      if (!any(str_detect(colnames(df), abs_full[i]))) {
        df[, abs_full[i]] <- "Not-tested"
      }
    }
    
    return(df)
  }
  
  # 'summarise_resist': Summarizes a set of resistance test results into a single outcome
  summarise_resist <- function(sample) {
    case_when(
      any(sample == "Resistant", na.rm = TRUE) ~ "Resistant", 
      any(sample == "Susceptible", na.rm = TRUE) ~ "Susceptible", 
      any(sample == "Not-tested", na.rm = TRUE) ~ "Not-tested")
  }
  
  # ---------------------------------------------------------------------------
  # Begin Cohort Creation Process
  # ---------------------------------------------------------------------------
  
  # Process antibiotic consumption data to flag valid initiation (based on several antibiotic types)
  ab_consumption <- ab_consp_all %>%
    mutate(ab_ok = ab_first_fluoroquinolones |
             ab_first_cephalosporins_1st_gen |
             ab_first_cephalosporins_2nd_gen |
             ab_first_augmentin)
  
  # Identify unique patient IDs with valid antibiotic purchases
  start_ids <- unique(ab_consumption$id[ab_consumption$ab_ok])
  ex_start <- paste(length(start_ids),
                    "patients who purchased antibiotic between 2017-2019")
  message(ex_start)
  
  # ---------------------------------------------------------------------------
  # Exclusion Based on DMG Criteria (e.g., Organ Transplant, Nursing Home)
  # ---------------------------------------------------------------------------
  
  dmg_sw <- dmg_sw_all %>%
    filter(id %in% start_ids) %>%
    distinct(id, .keep_all = TRUE)
  
  ex_before_total <- paste("Among", nrow(dmg_sw), "patients")
  message(ex_before_total)
  
  ex_nurse <- paste(sum(dmg_sw$sw_LTCF), "reside in a nursing home")
  message(ex_nurse)
  
  ex_trans <- paste(sum(dmg_sw[!dmg_sw$sw_LTCF,]$sw_trans), "underwent organ transplant")
  message(ex_trans)
  
  # Exclude patients who either underwent an organ transplant or reside in a nursing home
  dmg_no_sw <- dmg_sw %>%
    filter(!sw_trans,
           !sw_LTCF)
  
  dmg_id <- dmg_no_sw$id
  
  # ---------------------------------------------------------------------------
  # Process Antibiotic Data: Apply Washout Criteria and Exclude Ineligible Patients
  # ---------------------------------------------------------------------------
  
  ab_id <- ab_consumption %>%
    filter(id %in% dmg_id) %>%
    select(id) %>%
    distinct()
  
  n <- nrow(ab_id) / 32
  nr <- nrow(ab_id)
  ab_con <- split(ab_id, rep(1:ceiling(nr/n), each = n, length.out = nr))
  
  # Set up parallel processing using 32 clusters
  my.cluster <- parallel::makeCluster(
    32, 
    type = "PSOCK",
    outfile = ""
  )
  
  doParallel::registerDoParallel(cl = my.cluster)
  
  # Summarize antibiotic consumption data in parallel
  ab_summary <- foreach(i = 1:length(ab_con), 
                        .combine = rbind,
                        .packages = "tidyverse") %dopar% {
                          ab_consumption %>%
                            filter(id %in% ab_con[[i]]$id) %>%
                            group_by(id) %>%
                            arrange(id, ab_date) %>%
                            mutate(before = lag(ab_date),
                                   diff = difftime(ab_date, before, units = "days"),
                                   date_ok = (is.na(diff) | (diff > ab_washout))) %>%
                            ungroup() %>%
                            mutate(sumVar = rowSums(select(., contains("ab_first")))) %>%
                            filter(sumVar == 1,
                                   date_ok,
                                   ab_ok,
                                   ab_date >= ab_washout_date) %>%
                            select(-sumVar)
                        }
  
  ab_wash_ids <- dmg_id[dmg_id %in% ab_summary$id]
  
  ex_ab_wash <- paste((length(dmg_id) - length(ab_wash_ids)),
                      "dropped due dual antibiotic or insufficient washout time")
  message(ex_ab_wash)
  gc()
  
  # ---------------------------------------------------------------------------
  # Process Microbiology Data: Filter Based on Available Data and Criteria
  # ---------------------------------------------------------------------------
  
  no_micro_id <- ab_wash_ids[ab_wash_ids %in% micro$id]
  ex_no_micro <- paste((length(ab_wash_ids) - length(no_micro_id)),
                       "dropped due to lack of micro")
  message(ex_no_micro)
  gc()
  
  # Filter microbiology records by outcome, sample source, gram stain, and organism group
  if (is.null(micro_names)) {
    micro_filtered <- micro %>%
      filter(id %in% no_micro_id) %>%
      filter_at(vars(all_of(outcome_abs)), all_vars(. != "Not-tested")) %>% 
      filter(sample_source %in% sample_sources,
             micro_gram %in% micro_grams,
             !(micro_group %in% c("Candida", "Cryptococcus"))) %>%
      droplevels() %>% 
      mutate_at(vars(matches("outcome_")), function(x) {
        factor(ifelse(is.na(x), "Not-tested", as.character(x)))
      })
  } else {
    micro_filtered <- micro %>%
      filter(id %in% no_micro_id) %>%
      filter_at(vars(all_of(outcome_abs)), all_vars(. != "Not-tested")) %>% 
      filter(sample_source %in% sample_sources,
             micro_gram %in% micro_grams,
             micro_name %in% micro_names,
             !(micro_group %in% c("Candida", "Cryptococcus"))) %>%
      droplevels() %>% 
      mutate_at(vars(matches("outcome_")), function(x) {
        factor(ifelse(is.na(x), "Not-tested", as.character(x)))
      })
  }
  
  mic_id <- micro_filtered %>%
    select(id) %>%
    distinct()
  
  n <- nrow(mic_id) / 32
  nr <- nrow(mic_id)
  mic_out <- split(mic_id, rep(1:ceiling(nr/n), each = n, length.out = nr))
  
  doParallel::registerDoParallel(cl = my.cluster)
  
  # Match microbiology sample dates to antibiotic dates in parallel
  ab_micro_match_all <- foreach(i = 1:length(mic_out), 
                                .combine = rbind,
                                .packages = "tidyverse") %dopar% {
                                  micro_filtered %>%
                                    filter(id %in% mic_out[[i]]$id) %>%
                                    distinct(id, sample_date) %>%
                                    group_by(id, sample_date) %>%
                                    reframe(ab_date = return_ab_date(sample_date, ab_summary[ab_summary$id == id,]$ab_date)) %>%
                                    ungroup() %>%
                                    filter(!is.na(ab_date))
                                }
  
  # Exclude patients with mismatched antibiotic and microbiology dates
  ab_micro_ids <- no_micro_id[no_micro_id %in% ab_micro_match_all$id]
  
  ex_ab_micro <- paste((length(no_micro_id) - length(ab_micro_ids)),
                       "dropped due to ab-micro miss-match")
  message(ex_ab_micro)
  gc()
  
  # ---------------------------------------------------------------------------
  # Process Hospitalization Data and Exclude Ineligible Records
  # ---------------------------------------------------------------------------
  
  hospital <- hosp_all %>%
    filter(id %in% ab_micro_ids)
  
  rem_hos_id <- ab_micro_match_all %>%
    select(id) %>%
    distinct()
  
  n <- nrow(rem_hos_id) / 32
  nr <- nrow(rem_hos_id)
  hosp_out <- split(rem_hos_id, rep(1:ceiling(nr/n), each = n, length.out = nr))
  
  doParallel::registerDoParallel(cl = my.cluster)
  
  # Check hospitalization criteria in parallel
  ab_micro_match_hosp <- foreach(i = 1:length(hosp_out), 
                                 .combine = rbind,
                                 .packages = "tidyverse") %dopar% {
                                   ab_micro_match_all %>%
                                     filter(id %in% hosp_out[[i]]$id) %>%
                                     distinct(id, sample_date, ab_date) %>%
                                     group_by(id, sample_date, ab_date) %>%
                                     reframe(hosp_ok = return_hosp_date(ab_date, sample_date, 
                                                                        hospital[hospital$id == id,]$hospital_start_date,
                                                                        hospital[hospital$id == id,]$hospital_end_date)) %>%
                                     ungroup()
                                 }
  
  hosp_ids <- ab_micro_ids[ab_micro_ids %in% ab_micro_match_hosp$id[ab_micro_match_hosp$hosp_ok]]
  ex_hosp <- paste((length(ab_micro_ids) - length(hosp_ids)),
                   "dropped due hospitalization during the period")
  message(ex_hosp)
  gc()
  
  # ---------------------------------------------------------------------------
  # Summarize Microbiology Culture Data
  # ---------------------------------------------------------------------------
  
  micro_summary <- ab_micro_match_hosp %>%
    filter(id %in% hosp_ids) %>%
    group_by(id) %>%
    arrange(id, ab_date) %>%
    mutate(before = lag(ab_date),
           diff = difftime(ab_date, before, units = "days")) %>%
    filter(sample_date == min(sample_date)) %>%
    filter(ab_date == max(ab_date)) %>%
    ungroup() %>%
    left_join(micro_filtered, by = c("id", "sample_date")) %>%
    droplevels() %>%
    distinct() %>%
    select(-c(before, diff, hosp_ok))
  
  micro_sense_outcome <- micro_summary %>%
    group_by(id, sample_date, ab_date) %>%
    summarise_at(vars(starts_with("outcome_sample_ab_res")), summarise_resist) %>%
    ungroup()
  
  micro_name_outcome <- micro_summary %>%
    mutate(eq = TRUE) %>%
    distinct(id, micro_name, eq) %>%
    pivot_wider(id_cols = id, 
                names_from = micro_name,
                values_from = eq,
                values_fill = FALSE,
                names_glue = "outcome_sample_mic_{micro_name}") %>%
    janitor::clean_names()
  
  micro_hosp_outcome <- micro_summary %>% group_by(id) %>%
    summarise(outcome_sample_hospital = any(sample_hospital))
  
  micro_outcome_combine <- micro_name_outcome %>%
    left_join(micro_sense_outcome, by = "id") %>%
    left_join(micro_hosp_outcome, by = "id")
  
  cohort_dates <- micro_summary %>% 
    distinct(id, sample_date, ab_date)
  
  # ---------------------------------------------------------------------------
  # Process Initial Microbiology Data for Cohort Start
  # ---------------------------------------------------------------------------
  
  micro_start <- micro %>%
    left_join(cohort_dates %>% select(id, ab_date), by = c("id")) %>%
    mutate(start_time = difftime(sample_date, ab_date, units = "days"),
           date_ok = start_time >= 0 & start_time <= 14) %>%
    filter(date_ok)
  
  micro_start_name <- micro_start %>%
    mutate(eq = TRUE) %>%
    distinct(id, micro_name, eq) %>%
    pivot_wider(id_cols = id, 
                names_from = micro_name,
                values_from = eq,
                values_fill = FALSE,
                names_glue = "start_sample_mic_{micro_name}") %>%
    janitor::clean_names()
  
  micro_start_gram <- micro_start %>%
    mutate(eq = TRUE) %>%
    distinct(id, micro_gram, eq) %>%
    drop_na(micro_gram) %>%
    pivot_wider(id_cols = id, 
                names_from = micro_gram,
                values_from = eq,
                values_fill = FALSE,
                names_glue = "start_sample_gram_{micro_gram}") %>%
    janitor::clean_names()
  
  micro_start_ab <- micro_start %>%
    group_by(id) %>%
    summarise_at(vars(starts_with("outcome_sample_ab_res")), summarise_resist) %>%
    ungroup() %>%
    rename_with(~ gsub("^outcome", "start", .), starts_with("outcome"))
  
  micro_start_combine <- micro_start_gram %>%
    left_join(micro_start_name, by = "id") %>%  
    left_join(micro_start_ab, by = "id")
  
  gc()
  
  # ---------------------------------------------------------------------------
  # Censoring Process: Define Censor Dates 
  # ---------------------------------------------------------------------------
  
  doParallel::registerDoParallel(cl = my.cluster)
  
  censor_start <- foreach(i = 1:length(hosp_out), 
                          .combine = rbind,
                          .packages = "tidyverse") %dopar% {
                            ab_micro_match_all %>%
                              filter(id %in% hosp_out[[i]]$id) %>%
                              distinct(id, sample_date, ab_date) %>%
                              left_join(hospital, by = "id") %>%
                              mutate(diff_ab_end = difftime(hospital_end_date, ab_date, units = "days"),
                                     diff_ab_start = difftime(hospital_start_date, ab_date, units = "days"),
                                     date_end_ok = diff_ab_end > follow_up | diff_ab_end < -90 | is.na(diff_ab_end),
                                     date_start_ok = diff_ab_start > follow_up | diff_ab_start < -90 | is.na(diff_ab_start),
                                     date_ok = date_start_ok & date_end_ok)
                          }
  
  censor_dates <- censor_start %>%
    filter(date_ok) %>%
    group_by(id) %>%
    summarise(min_date = min(ab_date)) %>%
    full_join(cohort_dates, by = "id") %>%
    mutate(start_date = case_when(
      is.na(ab_date) ~ min_date,
      TRUE ~ ab_date)) %>%
    select(-c(min_date, ab_date)) %>%
    distinct(id, .keep_all = TRUE)
  
  # Process microbiology data for censoring
  micro_censor <- micro %>%
    left_join(censor_dates %>% select(id, start_date), by = c("id")) %>%
    mutate(start_time = difftime(sample_date, start_date, units = "days"),
           date_ok = start_time >= 0 & start_time <= 14) %>%
    filter(date_ok)
  
  micro_censor_name <- micro_censor %>%
    mutate(eq = TRUE) %>%
    distinct(id, micro_name, eq) %>%
    pivot_wider(id_cols = id,
                names_from = micro_name,
                values_from = eq,
                values_fill = FALSE,
                names_glue = "start_sample_mic_{micro_name}") %>%
    janitor::clean_names()
  
  micro_censor_gram <- micro_censor %>%
    mutate(eq = TRUE) %>%
    distinct(id, micro_gram, eq) %>%
    drop_na(micro_gram) %>%
    pivot_wider(id_cols = id,
                names_from = micro_gram,
                values_from = eq,
                values_fill = FALSE,
                names_glue = "start_sample_gram_{micro_gram}") %>%
    janitor::clean_names()
  
  micro_censor_ab <- micro_censor %>%
    group_by(id) %>%
    summarise_at(vars(starts_with("outcome_sample_ab_res")), summarise_resist) %>%
    ungroup() %>%
    rename_with(~ gsub("^outcome", "start", .), starts_with("outcome"))
  
  micro_censor_source <- micro_censor %>%
    group_by(id) %>%
    summarise(sample_source_urine = any(sample_source == "Urine"),
              sample_source_blood = any(sample_source == "Blood"))
  
  # Combine censoring information with covariate data
  censor_df <- dmg_sw_all %>% 
    distinct(id, .keep_all = TRUE) %>%
    left_join(cov_df, by = "id") %>%  
    left_join(censor_dates, by = "id") %>%  
    left_join(micro_censor_name, by = "id") %>%  
    left_join(micro_censor_gram, by = "id") %>%  
    left_join(micro_censor_ab, by = "id") %>%  
    left_join(micro_censor_source, by = "id") %>%  
    ungroup() %>%
    mutate(age = as.numeric(difftime(start_date, birth_date, units = "days"))/365,
           start_date = factor(case_when(
             month(start_date) <= 3 ~ paste0("1_", year(start_date)),
             month(start_date) <= 6 ~ paste0("2_", year(start_date)),
             month(start_date) <= 9 ~ paste0("3_", year(start_date)),
             month(start_date) <= 13 ~ paste0("4_", year(start_date))
           )),
           censor = ifelse(is.na(sample_date), 1, 0)) %>%
    mutate_at(vars(matches("start_sample_mic")), function(x) {
      factor(case_when(
        x ~ "Found",
        !x ~ "not found",
        is.na(x) ~ "Not taken"
      ))
    }) %>%
    mutate_at(vars(matches("start_sample_gram")), function(x) {
      factor(case_when(
        x ~ "Found",
        !x ~ "not found",
        is.na(x) ~ "Not taken"
      ))
    }) %>%
    mutate_at(vars(matches("sample_source")), function(x) {
      factor(case_when(
        x ~ "Yes",
        !x ~ "No",
        is.na(x) ~ "Not taken"
      ))
    }) %>%
    mutate_at(vars(matches("start_sample_ab")), function(x) {
      factor(ifelse(is.na(x), "Not-tested", as.character(x)))
    })  %>%
    mutate_at(vars(matches("start_sample_mic")), ~relevel(.x, "Not taken"))  %>%
    mutate_at(vars(matches("start_sample_ab")), ~relevel(.x, "Not-tested"))  %>%
    mutate_at(vars(matches("sample_source")), ~relevel(.x, "Not taken")) %>%
    filter(!is.na(start_date)) %>%
    select(-c(birth_date, sample_date, sw_AB_RegishYr, 
              sw_trans, sw_LTCF, death_date)) %>%
    droplevels()
  
  # ---------------------------------------------------------------------------
  # Calculate Inverse Probability of Censoring Weights (IPCW)
  # ---------------------------------------------------------------------------
  
  censor_w <- weightit(censor ~ . + poly(age, 3) - id - age,
                       data = censor_df, 
                       stabilize = TRUE,
                       estimand = "ATE", 
                       method = "glm")
  
  censor_ip_weights <- tibble(
    id = censor_df$id,
    censor = censor_df$censor,
    ps = censor_w$ps,
    censor_weights = censor_w$weights
  )
  
  # ---------------------------------------------------------------------------
  # Assemble Final Cohort Data Frame
  # ---------------------------------------------------------------------------
  
  df <- dmg_no_sw %>% 
    filter(id %in% cohort_dates$id) %>%
    left_join(cov_df, by = "id") %>%  
    left_join(micro_start_combine, by = "id") %>%
    left_join(micro_outcome_combine, by = "id") %>%
    left_join(ab_summary %>% select(-c(before, diff, date_ok)), 
              by = c("id", "ab_date")) %>% 
    left_join(censor_ip_weights %>% select(id, censor_weights), by = "id") %>%
    ungroup() %>%
    mutate_if(is.character, as.factor) %>%
    droplevels()
  
  miss_cov_id <- df %>%
    filter(is.na(cov_dm) |
             is.na(cov_copd) |
             is.na(cov_onco) |
             is.na(cov_preg)) %>%
    pull(id)
  
  ex_cov <- paste(length(miss_cov_id),
                  "dropped due missing covs")
  message(ex_cov)
  
  pre_match_df <- df %>%
    filter(!(id %in% miss_cov_id)) %>%
    mutate(start_ab = factor(case_when(
      ab_first_augmentin ~ "Augmentin",
      ab_first_cephalosporins_1st_gen ~ "Ceph 1gen",
      ab_first_cephalosporins_2nd_gen ~ "Ceph 2gen",
      ab_first_fluoroquinolones ~ "Fluoroquinolones"
    ),
    levels = c("Ceph 1gen", "Ceph 2gen", "Augmentin",
               "Fluoroquinolones")),
    ab_sample_diff = as.numeric(difftime(sample_date, ab_date, units = "days")),
    start_q = factor(case_when(
      month(ab_date) <= 3 ~ paste0("1_", year(ab_date)),
      month(ab_date) <= 6 ~ paste0("2_", year(ab_date)),
      month(ab_date) <= 9 ~ paste0("3_", year(ab_date)),
      month(ab_date) <= 12 ~ paste0("4_", year(ab_date))
    ))) %>%
    select(id, age, sex, sector, cov_dm, cov_copd,
           cov_onco, cov_preg, start_q,
           start_ab, starts_with("start_sample_"),
           ab_sample_diff, censor_weights,
           starts_with("outcome_sample_")) %>%
    mutate_at(vars(matches("start_sample_mic")), function(x) {
      case_when(
        x ~ "Present",
        is.na(x) ~ "Not taken",
        !x ~ "Not present"
      )
    }) %>%
    mutate_at(vars(matches("start_sample_gram")), function(x) {
      factor(case_when(
        x ~ "Found",
        !x ~ "not found",
        is.na(x) ~ "Not taken"
      ))
    }) %>%
    mutate_at(vars(matches("start_sample_ab")), function(x) {
      ifelse(is.na(x), "Not-tested", as.character(x))
    }) %>%
    mutate_at(vars(matches("outcome_sample_ab_res")), function(x) {
      case_when(
        x == "Resistant" ~ TRUE,
        x == "Susceptible" ~ FALSE
      )
    }) %>%
    select(-c(start_sample_mic_other_gram_na,
              start_sample_ab_res_macrolides,
              outcome_sample_ab_res_macrolides)) %>%
    filter(!is.na(start_ab)) %>%
    filter(censor_weights < 10)
  
  # Output exclusion summary messages
  ex_total <- paste("A total of", 
                    (length(start_ids) - nrow(pre_match_df)),
                    "removed")
  
  ex_total_left <- paste(nrow(pre_match_df), "remained in final cohort")
  
  message(paste(ex_total, "\n", ex_total_left))
  
  # ---------------------------------------------------------------------------
  # Generate Summary Tables Using gtsummary
  # ---------------------------------------------------------------------------
  
  # Table 1: Demographic and sample characteristics stratified by antibiotic exposure
  pre_match_df %>%
    mutate_at(vars(starts_with("start_sample")),
              ~ifelse(. %in% c("Not-tested", "Not taken"), NA, .)) %>%
    select(-id) %>%
    mutate(start_ab = paste("Exposed to", start_ab)) %>%
    tbl_summary(
      by = "start_ab",
      type = list(
        c(starts_with("start_sample"), "sex") ~ "dichotomous"
      ),
      value = list(
        sex ~ "Female",
        c(starts_with("start_sample_ab")) ~ "Resistant",
        c(starts_with("start_sample_mic")) ~ "Present",
        c(starts_with("start_sample_gram")) ~ "1"
      ),
      missing = "no"
    ) %>%
    add_n(statistic = "{n_miss} ({p_miss})") %>%
    add_overall() %>%
    modify_header(n = "**Missing**") -> tbl_1
  
  tbl_1$table_body$label <- sapply(tbl_1$table_body$label, label_get, USE.NAMES = FALSE)
  
  gt::gtsave(as_gt(tbl_1), file = file.path("export", paste0(obj_name, "_tbl_1", ".html")))
  
  # Table 2: Summary of antibiotic resistance outcomes
  pre_match_df %>%
    mutate_at(vars(starts_with("start_sample")),
              ~ifelse(. %in% c("Not-tested", "Not taken"), NA, .)) %>%
    select(starts_with("start_sample_ab")) %>%
    tbl_summary(
      type = list(
        c(starts_with("start_sample")) ~ "dichotomous"
      ),
      value = list(
        c(starts_with("start_sample_ab")) ~ "Resistant"
      ),
      missing = "no"
    ) %>%
    add_n(statistic = "{n_miss} ({p_miss})") %>%
    modify_header(n = "**Missing**") -> tbl_2 
  
  tbl_2$table_body$label <- sapply(tbl_2$table_body$label, label_get, USE.NAMES = FALSE)
  
  gt::gtsave(as_gt(tbl_2), file = file.path("export", paste0(obj_name, "_tbl_2", ".html")))
  
  # ---------------------------------------------------------------------------
  # Create Final Cohort Summary and Return Cohort Object
  # ---------------------------------------------------------------------------
  
  final_ab_sum <- pre_match_df %>%
    group_by(start_ab) %>%
    summarise(n = n()) %>%
    mutate(p = round(100 * n / sum(n), 2))
  
  ex_text <- paste(ex_trans, ex_nurse, ex_no_micro, ex_ab_wash, ex_ab_micro, 
                   ex_hosp, ex_cov, ex_total, sep = "\n")
  
  ex_total_left_ab <- paste(nrow(pre_match_df), "exposed to relevant antibiotics")
  
  ab_a <- paste0(final_ab_sum[1, 2, drop = T], " (", final_ab_sum[1, 3, drop = T], "%) ", final_ab_sum[1, 1, drop = T])
  ab_b <- paste0(final_ab_sum[2, 2, drop = T], " (", final_ab_sum[2, 3, drop = T], "%) ", final_ab_sum[2, 1, drop = T])
  ab_c <- paste0(final_ab_sum[3, 2, drop = T], " (", final_ab_sum[3, 3, drop = T], "%) ", final_ab_sum[3, 1, drop = T])
  ab_d <- paste0(final_ab_sum[4, 2, drop = T], " (", final_ab_sum[4, 3, drop = T], "%) ", final_ab_sum[4, 1, drop = T])
  
  cohort_summary <- c(ex_start, ex_text, ex_total_left, ex_total_left_ab,
                      ab_a, ab_b, ab_c, ab_d)
  
  # Create the cohort object with class "cohort"
  cohort <- list(
    df = df,
    pre_match_df = pre_match_df,
    summary = cohort_summary, 
    follow_up = follow_up
  )
  
  class(cohort) <- "cohort"
  
  # Assign the cohort object to the global environment with the specified name
  assign(obj_name, cohort, envir = .GlobalEnv)
}