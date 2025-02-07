# -----------------------------------------------------------------------------
# Comprehensive Workflow for Cohort Creation, Matching, and Risk Difference Analysis 
# in Bacterial Antibiotic Exposure Studies
#
# This R script orchestrates a complete analysis pipeline to investigate the impact 
# of different antibiotic exposures on bacterial resistance. The workflow includes:
#
# 1. Setup: Loading required packages and sourcing custom functions.
# 2. Cohort Creation: Generating datasets for key microbial cohorts (e.g., 
#    Escherichia coli, and a combined cohort including Proteus mirabilis and Klebsiella pneumoniae)
#    with varying antibiotic washout periods and outcome definitions.
# 3. Matching: Applying propensity score matching (via the `match_ab` function) to 
#    compare antibiotic treatments (e.g., Cephalosporins of the 1st vs. 2nd generation, 
#    Augmentin, and Fluoroquinolones) while enforcing exact matching on relevant covariates 
#    (such as sex, sector, and start quarter) and optimizing caliper settings.
# 4. Data Aggregation: Compiling the matched datasets into summary data frames.
# 5. Analysis: Calculating risk differences using the `rd_calc` function and saving detailed summaries.
# 6. Visualization: Generating figures (with `rd_fig`) to graphically display the risk differences 
#    across various antibiotic comparisons and washout scenarios.
#
# This structured approach supports robust evaluation of antibiotic exposure effects, facilitating
# both statistical inference and clear visualization of the study findings.
# -----------------------------------------------------------------------------

# packages & functions --------------------------------------------------------
library(tidyverse)
load("bigfiles.RData")
source("Rscripts/02 - df_create.R")
source("Rscripts/03 - matching.R")
source("Rscripts/04 - analysis.R")
source("Rscripts/05 - funcs.R")

# Create DFs ------------------------------------------------------------------
create_cohort(micro_names = "Escherichia coli",
              obj_name = "ecoli_182")
create_cohort(micro_names = "Escherichia coli",
              ab_washout = 182,
              obj_name = "ecoli_182_wo")

create_cohort(micro_names = c("Escherichia coli","Proteus mirabilis","Klebsiella pneumoniae"),
              obj_name = "ecoli_prot_kleb_182")
create_cohort(micro_names = c("Escherichia coli","Proteus mirabilis","Klebsiella pneumoniae"),
              ab_washout = 182,
              obj_name = "ecoli_prot_kleb_wo_182")
create_cohort(micro_names = "Escherichia coli",
              outcome_abs = "Cef 2",
              obj_name = "ecoli_182_cef_2")

save(ecoli_182, 
     ecoli_182_wo,
     ecoli_prot_kleb_182, 
     ecoli_prot_kleb_wo_182,
     ecoli_182_cef_2,
     file = "cohorts.RData")

rm(ab_consumption_all, cov_df, dmg_sw_all, hospital_all, micro_df, micro_df_long)
gc()

# Matching --------------------------------------------------------------------
## E.coli ---------------------------------------------------------------------
### 3 months WO ---------------------------------------------------------------
match_ab(ecoli_182, ref_ab = "Ceph 1gen", match_ab = "Ceph 2gen", 
         exact_vars = "sex", caliper = 0.2, short = FALSE)
match_ab(ecoli_182, ref_ab = "Ceph 1gen", match_ab = "Augmentin", 
         exact_vars = c("sex"), caliper = 0.5, short = FALSE) #,"start_q"
match_ab(ecoli_182, ref_ab = "Ceph 1gen", match_ab = "Fluoroquinolones", 
         exact_vars =  c("sex"), caliper = 0.03, short = FALSE)
match_ab(ecoli_182, ref_ab = "Ceph 2gen", match_ab = "Augmentin", 
         caliper = 0.018, short = FALSE) 
match_ab(ecoli_182, ref_ab = "Ceph 2gen", match_ab = "Fluoroquinolones", 
         exact_vars = c("sex"), caliper = 0.005, short = FALSE)

### 6 mo wo -------------------------------------------------------------------
match_ab(ecoli_182_wo, ref_ab = "Ceph 1gen", match_ab = "Ceph 2gen", 
         exact_vars = "sex", caliper = 0.5, short = FALSE)
match_ab(ecoli_182_wo, ref_ab = "Ceph 1gen", match_ab = "Augmentin", 
         caliper = 1.5, short = FALSE) 
match_ab(ecoli_182_wo, ref_ab = "Ceph 1gen", match_ab = "Fluoroquinolones", 
         exact_vars = "sector", caliper = 0.1, short = FALSE)
match_ab(ecoli_182_wo, ref_ab = "Ceph 2gen", match_ab = "Augmentin", 
         exact_vars = "start_q", caliper = 0.02, short = FALSE) 
match_ab(ecoli_182_wo, ref_ab = "Ceph 2gen", match_ab = "Fluoroquinolones", 
         exact_vars = c("sex","start_q", "sector"), caliper = 0.01, short = FALSE)

## E.coli + Proteus + Klebsiella ----------------------------------------------
### 3 months WO ---------------------------------------------------------------
match_ab(ecoli_prot_kleb_182, 
         ref_ab = "Ceph 1gen", match_ab = "Ceph 2gen", 
         exact_vars = "sex", caliper = 0.2, short = FALSE)
match_ab(ecoli_prot_kleb_182, 
         ref_ab = "Ceph 1gen", match_ab = "Augmentin", 
         exact_vars = "sex", caliper = 0.5, short = FALSE)
match_ab(ecoli_prot_kleb_182, 
         ref_ab = "Ceph 1gen", match_ab = "Fluoroquinolones", 
         exact_vars = c("sex"), caliper = 0.04, short = FALSE)
match_ab(ecoli_prot_kleb_182, 
         ref_ab = "Ceph 2gen", match_ab = "Augmentin", 
         exact_vars = "sex", caliper = 0.01, short = FALSE) 
match_ab(ecoli_prot_kleb_182, 
         ref_ab = "Ceph 2gen", match_ab = "Fluoroquinolones", 
         exact_vars = "sex", caliper = 0.005, short = FALSE)

### 6 mo wo -------------------------------------------------------------------
match_ab(ecoli_prot_kleb_wo_182, 
         ref_ab = "Ceph 1gen", match_ab = "Ceph 2gen", 
         exact_vars = "sex", caliper = 0.2, short = FALSE)
match_ab(ecoli_prot_kleb_wo_182, 
         ref_ab = "Ceph 1gen", match_ab = "Augmentin", 
         exact_vars = "sex", caliper = 0.5, short = FALSE)
match_ab(ecoli_prot_kleb_wo_182, 
         ref_ab = "Ceph 1gen", match_ab = "Fluoroquinolones", 
         exact_vars = c("sex"), caliper = 0.03, short = FALSE)
match_ab(ecoli_prot_kleb_wo_182, 
         ref_ab = "Ceph 2gen", match_ab = "Augmentin", 
         exact_vars = "sex", caliper = 0.02, short = FALSE) 
match_ab(ecoli_prot_kleb_wo_182, 
         ref_ab = "Ceph 2gen", match_ab = "Fluoroquinolones", 
         exact_vars = "sex", caliper = 0.005, short = FALSE)

## Cef 2 gen ------------------------------------------------------------------
match_ab(ecoli_182_cef_2, ref_ab = "Ceph 1gen", match_ab = "Ceph 2gen", 
         outcome_var = "start_sample_ab_res_cef_2",
         exact_vars = "sex", caliper = 0.5, short = FALSE)
match_ab(ecoli_182_cef_2, ref_ab = "Ceph 1gen", match_ab = "Augmentin", 
         outcome_var = "start_sample_ab_res_cef_2",
         exact_vars = c("sex","start_q"), caliper = 1, short = FALSE) 
match_ab(ecoli_182_cef_2, ref_ab = "Ceph 1gen", match_ab = "Fluoroquinolones", 
         outcome_var = "start_sample_ab_res_cef_2",
         exact_vars = c("sex","sector"), caliper = 0.08, short = FALSE)
match_ab(ecoli_182_cef_2, ref_ab = "Ceph 2gen", match_ab = "Augmentin", 
         outcome_var = "start_sample_ab_res_cef_2",
         exact_vars = "start_q", caliper = 0.02, short = FALSE) 
match_ab(ecoli_182_cef_2, ref_ab = "Ceph 2gen", match_ab = "Fluoroquinolones", 
         outcome_var = "start_sample_ab_res_cef_2",
         exact_vars = c("sex","start_q", "sector"), caliper = 0.015, short = FALSE)

## save matched data ----------------------------------------------------------
matched_list_ecoli <- lapply(ls(pattern = "^matched_ecoli_182_ceph"), get)
matched_list_ecoli_wo <- lapply(ls(pattern = "^matched_ecoli_182_wo"), get)
matched_list_ecoli_prot_kleb <- lapply(ls(pattern = "^matched_ecoli_prot_kleb_182"), get)
matched_list_ecoli_prot_kleb_wo <- lapply(ls(pattern = "^matched_ecoli_prot_kleb_wo"), get)
matched_list_ecoli_cef_2 <- lapply(ls(pattern = "^matched_ecoli_182_cef_2"), get)

save(matched_list_ecoli, matched_list_ecoli_prot_kleb,
     matched_list_ecoli_wo, matched_list_ecoli_prot_kleb_wo,
     matched_list_ecoli_cef_2, 
     file = "matched_list.RData")

matched_df_ecoli <- create_matched_df(matched_list_ecoli)
matched_df_ecoli_wo <- create_matched_df(matched_list_ecoli_wo)
matched_df_ecoli_prot_kleb <- create_matched_df(matched_list_ecoli_prot_kleb)
matched_df_ecoli_prot_kleb_wo <- create_matched_df(matched_list_ecoli_prot_kleb_wo)
matched_df_ecoli_cef_2 <- create_matched_df(matched_list_ecoli_cef_2)

save(matched_df_ecoli, matched_df_ecoli_prot_kleb,
     matched_df_ecoli_wo,matched_df_ecoli_prot_kleb_wo,
     matched_df_ecoli_cef_2, 
     file = "matched_df.RData")

# Analyze matched data -------------------------------------------------------
rd_calc(matched_df_ecoli)
rd_calc(matched_df_ecoli_wo)
rd_calc(matched_df_ecoli_prot_kleb)
rd_calc(matched_df_ecoli_prot_kleb_wo)
rd_calc(matched_df_ecoli_cef_2)

save(rd_summary_ecoli, rd_summary_ecoli_prot_kleb,
     rd_summary_ecoli_wo,rd_summary_ecoli_prot_kleb_wo,
     rd_summary_ecoli_cef_2, 
     file = "rd_summary.RData")

## Fig 1 -----------------------------------------------------------------------
rd_fig(rd_summary_ecoli)
rd_fig(rd_summary_ecoli, ref_abs = "Ceph 2gen")
rd_fig(rd_summary_ecoli_wo)
rd_fig(rd_summary_ecoli_wo, ref_abs = "Ceph 2gen")
rd_fig(rd_summary_ecoli_prot_kleb)
rd_fig(rd_summary_ecoli_prot_kleb, ref_abs = "Ceph 2gen")
rd_fig(rd_summary_ecoli_prot_kleb_wo)
rd_fig(rd_summary_ecoli_prot_kleb_wo, ref_abs = "Ceph 2gen")
rd_fig(rd_summary_ecoli_cef_2)
rd_fig(rd_summary_ecoli_cef_2, ref_abs = "Ceph 2gen")
