# -----------------------------------------------------------------------------
# Data Preparation for Antibiotic Resistance Analysis in Hebrew Locale
#
# This script performs data preparation and transformation tasks for an
# analysis of antibiotic resistance. The following key steps are executed:
#
#   1. Set the locale to Hebrew for proper handling of Hebrew text and dates.
#   2. Load the tidyverse package for data manipulation.
#   3. Define a helper function 'summarise_resist' to summarize antibiotic
#      resistance test results into a single outcome.
#   4. Read and transform raw microbiology data from a CSV file into two data
#      frames:
#         a. micro_df_long: The raw long-format microbiology data.
#         b. micro_df: A cleaned and restructured version with additional derived
#            variables (e.g., bacterial groupings, sensitivity outcomes).
#   5. Read and process additional datasets:
#         - Antibiotic consumption data (ab_consumption_all)
#         - Comorbidity and demographic data (cov_df, dmg_sw_all)
#         - Hospitalization data (hospital_all)
#   6. Save all processed datasets into a single RData file for later use.
#
# Dependencies: tidyverse (for data manipulation), janitor (for cleaning names),
#               and base R functions for date conversion and locale setting.
#
# Author: [Your Name]
# Date: [Today's Date]
# -----------------------------------------------------------------------------

# Set locale to Hebrew to ensure proper handling of Hebrew characters and date names.
Sys.setlocale("LC_ALL", "Hebrew")

# Load tidyverse for data manipulation and piping functions.
library(tidyverse)

# -----------------------------------------------------------------------------
# Define Helper Functions
# -----------------------------------------------------------------------------

# summarise_resist: Summarizes a vector of antibiotic test results into a single
# outcome. If any test result in 'sample' is "Resistant", return "Resistant";
# if any are "Susceptible", return "Susceptible"; otherwise, if any are "Not-tested",
# return "Not-tested".
summarise_resist <- function(sample) {
  case_when(
    any(sample == "Resistant", na.rm = TRUE) ~ "Resistant", 
    any(sample == "Susceptible", na.rm = TRUE) ~ "Susceptible", 
    any(sample == "Not-tested", na.rm = TRUE) ~ "Not-tested"
  )
}

# -----------------------------------------------------------------------------
# Read and Transform Microbiology Data
# -----------------------------------------------------------------------------

# Read the raw microbiology CSV file into a long-format dataframe and transform it.
micro_df_long <- read_csv("raw_data/NNH_AB_Micro17_19ByTZ_Hash.csv",
                          show_col_types = FALSE) %>%
  transmute(
    id = PatientIDInternal,
    sample_source = factor(DgimaType),
    sample_hospital = ifelse(SentByBTH == 1, TRUE, FALSE),
    sample_date = as.Date(date_bdika_day_Hash, format = "%Y-%m-%d"),
    sample_num = mispar_dgima,
    # Create a new grouping variable based on HaidakGroup_Name; if missing or "Not_grouped",
    # extract a value from shem_haidak.
    new_group = ifelse((is.na(HaidakGroup_Name) | HaidakGroup_Name == "Not_grouped"),
                       str_extract(shem_haidak, regex('^\\S*')),
                       HaidakGroup_Name),
    # Derive a factor variable (micro_gram) to indicate Gram stain status based on group names.
    micro_gram = factor(case_when(
      new_group %in% c("Candida", "Enterococcus", "Staphylococcus",
                       "Staphylococcus aureus", "Streptococcus",
                       "Streptococcus pneumoniae", "Strep.", "Abiotrophia",
                       "Actinomyces", "Aerococcus", "Arcanobacter",
                       "Bifidobacterium", "Clostridium", "Corynebacterium",
                       "Kocuria", "Lactobacillus", "Micrococcus", 
                       "Pediococcus", "Peptostreptococcus", "Propionibacterium",
                       "Ruminococcus", "STOMATOCOCCUS", "St.", "Str.gallolityc.sp",
                       "Strep", "Veillonella", "Coryneb.striatum", "Listeria",
                       "Arcanobacterium", "Pasteurella", "Eubacterium",
                       "Leuconostoc", "Mycobacterium", "Gemella", "streptococcus",
                       "Peptococcus", "Rhodococcus", "LACTOCOCCUS", "Mycobacteium",
                       "Erysipelothrix", "Bacillus") ~ "Positive",
      new_group %in% c("Acinetobacter", "Bacteroides","Campylobacter",
                       "Citrobacter", "E.coli","Enterobacter",
                       "Haemophilus influenzae", "Klebsiella",
                       "Morganella", "Proteus", "Providencia",
                       "Pseudomonas", "Pseudomonas aeruginosa",
                       "Salmonela", "Serratia", "Stenotrophormonas",
                       "Alcaligenes", "Brucella", "Cardiobacterium",
                       "Chryseobacterium", "Ewingella", "Fusobacterium",
                       "Haemophilus", "Moraxella", "Mycoplasma", 
                       "Neisseria", "Neisseriae", "Pantoea", "Prevotella",
                       "Saccharomyces", "Salmonella", "Tatumella",
                       "Achromobacter", "Burkholderia", "Leclercia",
                       "Kluyvera", "Hafnia", "Flavobacterium",
                       "Cedecea", "Aeromonas", "Sphingomonas", "Comamonas",
                       "Agrobacterium", "Ochrobacterium", "Brevundimonas",
                       "Cryptococcus", "Oligella", "Sphingobacterium",
                       "Eikenella", "Vibrio", "Capnocytophaga", "Flavimonas",
                       "Escherichia", "Actinobacillus", "Yersinia") ~ "Negative",
      # Additional conditions based on shem_haidak description.
      shem_haidak %in% c("Anaerobic gram pos bacilli", "Anaerobic gram pos cocci",
                         "Bacilli gram positive", "Beta hemol. strept. not group A") ~ "Positive",
      shem_haidak %in% c("Anaerobic gram neg bacilli", "Anaerobic gram neg cocci",
                         "Cocci gram negative", "Diplococci gram minus ex cell") ~ "Negative"
    )),
    # Group the microbiology results into broader bacterial groups.
    micro_group = factor(case_when(
      new_group  == "Coryneb.striatum" ~ "Corynebacterium",
      new_group  == "Haemophilus influenzae" ~ "Haemophilus",
      new_group  == "LACTOCOCCUS" ~ "Lactococcus",
      new_group  == "Mycobacteium" ~ "Mycobacterium",
      new_group  == "Neisseria" ~ "Neisseriae",
      new_group  == "Pseudomonas aeruginosa" ~ "Pseudomonas",
      new_group  == "Salmonela" ~ "Salmonella",
      new_group  == "Staphylococcus aureus" ~ "Staphylococcus",
      new_group  == "STOMATOCOCCUS" ~ "Stomatococcus",
      new_group %in% c("St.", "Str.gallolityc.sp", "Strep",
                       "Strep.", "streptococcus","Streptococcus",
                       "Streptococcus pneumoniae",
                       "Beta") ~ "Streptococcus",
      new_group %in% c("Anaerobic", "Bacilli",
                       "CDC", "Cocci", "Diplococci",
                       "Nonferment.GN", "Nonfermenter") ~ NA,
      TRUE ~ new_group
    )),
    # Set the microbiology name as a factor.
    micro_name = factor(shem_haidak),
    # Create a temporary variable for sensitivity text (handle problematic strings)
    micro_sens_t = ifelse(
      str_detect(regishut, "×\u009c×\u0090 ×\u0099×\u0093×\u0095×¢"),
      NA,
      regishut
    ),
    # Convert sensitivity test result to a factor with custom mapping.
    micro_sens = factor(case_when(
      micro_sens_t == "Intermediate" ~ "Resistant",   
      micro_sens_t == "Normal" ~ "Susceptible",
      TRUE ~ as.character(micro_sens_t)
    )),
    # Map antibiotic sensitivity for a specific set of antibiotics.
    micro_sens_ab = factor(case_when(
      shem_antibiotic %in% c("Cefuroxime",
                             "Cefuroxime axetil",
                             "Cefuroxime sodium") ~ "Cefuroxime",
      # Uncomment and modify the next line if needed for additional antibiotics.
      # shem_antibiotic %in% c("MIC Ceftriaxone") ~ "Ceftriaxone",
      shem_antibiotic %in% c("MIC Ceftriaxone") ~ "Ceftriaxone",
      TRUE ~ shem_antibiotic
    ))
  ) %>%
  # Remove temporary variables no longer needed
  select(-c(new_group, micro_sens_t)) %>%
  janitor::clean_names()

# Define a vector of bacterial names in the Enterobacteriaceae family for later use.
enterobacteriaceae <- c("Citrobacter", "E.coli", "Proteus", "Klebsiella", 
                        "Enterobacter", "Leclercia", "Salmonella", "Providencia",
                        "Morganella", "Serratia")

# -----------------------------------------------------------------------------
# Create a Clean Microbiology Dataset (micro_df)
# -----------------------------------------------------------------------------
micro_df <- micro_df_long %>%
  # Lump rare levels in micro_name together; if lumped as "Other", append the gram stain info.
  mutate(
    micro_name = fct_lump_prop(micro_df_long$micro_name, p = 0.01),
    micro_name = ifelse(micro_name == "Other", paste0("Other_gram", micro_gram), as.character(micro_name)),
    micro_sens = as.character(micro_sens)
  ) %>%
  group_by(id, sample_date, sample_source, sample_hospital,
           micro_gram, micro_group, micro_name, micro_sens_ab) %>%
  # Summarise sensitivity test results using the summarise_resist function.
  summarise(micro_sens = summarise_resist(micro_sens)) %>%
  ungroup() %>%
  # Pivot the data wider to create separate columns for each antibiotic sensitivity outcome.
  pivot_wider(
    id_cols = c(id, sample_date, sample_source, sample_hospital, micro_gram, micro_group, micro_name), 
    names_from = micro_sens_ab, 
    names_glue = "outcome_sample_ab_res_{micro_sens_ab}",
    values_from = micro_sens,
    values_fill = "Not-tested"
  ) %>%
  janitor::clean_names() %>%
  # Apply recoding to specific antibiotic outcomes based on the micro_name.
  mutate(
    outcome_sample_ab_res_ceftriaxone = ifelse(micro_name == "Pseudomonas", "Resistant",
                                               as.character(outcome_sample_ab_res_ceftriaxone)),
    outcome_sample_ab_res_augmentin = ifelse(micro_name == "Pseudomonas", "Resistant",
                                             as.character(outcome_sample_ab_res_amoxicillin_clavul_a)),
    outcome_sample_ab_res_nitrofurantoin = ifelse(micro_name == "Pseudomonas", "Resistant",
                                                  as.character(outcome_sample_ab_res_nitrofurantoin)),
    outcome_sample_ab_res_cef_3 = case_when(
      outcome_sample_ab_res_ceftazidime == "Resistant" ~ "Resistant",
      outcome_sample_ab_res_ceftriaxone == "Resistant" ~ "Resistant",
      outcome_sample_ab_res_ceftazidime == "Susceptible" ~ "Susceptible",
      outcome_sample_ab_res_ceftriaxone == "Susceptible" ~ "Susceptible",
      TRUE ~ "Not-tested"
    )
  ) %>%
  rowwise() %>%
  # Compute composite outcomes using summarise_resist on groups of antibiotic tests.
  mutate(
    outcome_sample_ab_res_cef_1 = summarise_resist(c(outcome_sample_ab_res_cephalexin,
                                                     outcome_sample_ab_res_cefazolin,
                                                     outcome_sample_ab_res_cephalothin)),
    outcome_sample_ab_res_cef_2 = summarise_resist(c(outcome_sample_ab_res_cefuroxime,
                                                     outcome_sample_ab_res_cefoxitin)),
    outcome_sample_ab_res_cef_3_int = (outcome_sample_ab_res_ceftazidime == "Resistant") & (outcome_sample_ab_res_ceftriaxone == "Resistant"),
    outcome_sample_ab_res_fluoroquinolones = summarise_resist(c(outcome_sample_ab_res_ofloxacin,
                                                                outcome_sample_ab_res_ciprofloxacin,
                                                                outcome_sample_ab_res_levofloxacin)),
    outcome_sample_ab_res_macrolides = summarise_resist(c(outcome_sample_ab_res_azithromicin,
                                                          outcome_sample_ab_res_erythromycin,
                                                          outcome_sample_ab_res_clarithromycin))
  ) %>%
  ungroup() %>%
  mutate(
    outcome_sample_ab_res_cef_2 = case_when(
      (micro_name %in% enterobacteriaceae) & outcome_sample_ab_res_cef_3_int ~ "Resistant",
      micro_name == "Pseudomonas" ~ "Resistant",
      TRUE ~ as.character(outcome_sample_ab_res_cef_2)
    ),
    outcome_sample_ab_res_cef_1 = case_when(
      (micro_name %in% enterobacteriaceae) & outcome_sample_ab_res_cef_3_int ~ "Resistant",
      micro_name == "Pseudomonas" ~ "Resistant",
      TRUE ~ as.character(outcome_sample_ab_res_cef_1)
    )
  ) %>%
  # Select only the final set of variables needed for analysis.
  select(
    id, sample_date, sample_source, sample_hospital,
    micro_gram, micro_group, micro_name,
    outcome_sample_ab_res_cef_1, outcome_sample_ab_res_cef_2, outcome_sample_ab_res_cef_3,
    outcome_sample_ab_res_augmentin, outcome_sample_ab_res_fluoroquinolones,
    outcome_sample_ab_res_macrolides, outcome_sample_ab_res_nitrofurantoin
  )

# -----------------------------------------------------------------------------
# Process Antibiotic Consumption Data
# -----------------------------------------------------------------------------

ab_consumption_all <- read_csv("raw_data/NNH_AB_Consumption17_19ByTZcsv_Hash.csv",
                               show_col_types = FALSE) %>%
  # Exclude unwanted antibiotic types.
  filter(!(Ab_type %in% c("IMIDAZOLE", "NR", "MACROLIDES", "CLINDAMYCIN"))) %>%
  transmute(
    id = PatientIDInternal,
    age = gil,
    ab_date = as.Date(date_nipuk_Hash, format = "%Y-%m-%d"),
    # Recode antibiotic type; e.g., if type equals "ANTISTAPHYLOCOCCAL PENICILLINS", label it as "PENICILLIN"
    ab_type = factor(ifelse(Ab_type == "ANTISTAPHYLOCOCCAL PENICILLINS", "PENICILLIN", Ab_type)),
    ab_amount = !is.na(kamut)
  ) %>%
  # Reshape the data to have one column per antibiotic type with a logical indicator.
  pivot_wider(
    id_cols = c(id, age, ab_date), 
    names_from = ab_type, 
    values_from = ab_amount, 
    names_glue = "ab_first_{ab_type}",
    values_fill = FALSE
  ) %>%
  janitor::clean_names() %>%
  # Ensure that if both amoxicillin and augmentin are present, adjust the flag for amoxicillin.
  mutate(ab_first_amoxicillin = ifelse(ab_first_amoxicillin & ab_first_augmentin,
                                       !ab_first_amoxicillin,
                                       ab_first_amoxicillin))

# -----------------------------------------------------------------------------
# Process Comorbidity Data (Covariates)
# -----------------------------------------------------------------------------

cov_df <- read_csv("raw_data/NNH_NewVars4Exclusion_30Jun_2BHashed_hashed.csv",
                   show_col_types = FALSE) %>%
  transmute(
    id = ID,
    cov_dm = IsDIABBef,
    cov_copd = IsCOPDBef,
    cov_onco = IsOncoActive5YrBef,
    cov_preg = IsPregInDate
  ) %>%
  # Convert comorbidity flags to logical values (TRUE if 1, FALSE otherwise)
  mutate_at(c("cov_dm", "cov_copd", "cov_onco", "cov_preg"), ~. == 1)

# -----------------------------------------------------------------------------
# Process Hospitalization Data
# -----------------------------------------------------------------------------

hospital_all <- read_csv("raw_data/NNH_AB_Ishp17_19ByTZ__correected_Hash.csv",
                         show_col_types = FALSE,
                         locale = locale(date_names = "he", encoding = "utf8")) %>%
  transmute(
    id = PatientIDInternal,
    hospital_start_date = as.Date(date_start_Hash, format = "%Y-%m-%d"),
    hospital_end_date = as.Date(date_end_Hash, format = "%Y-%m-%d")
  )

# -----------------------------------------------------------------------------
# Process Demographic Data and Exclusions
# -----------------------------------------------------------------------------

dmg_sw_all <- read_csv("raw_data/NNH_AB_Demog17_19ByTZ_Hash.csv",
                       show_col_types = FALSE,
                       locale = locale(date_names = "he", encoding = "utf8")) %>%
  transmute(
    id = PatientIDInternal,
    birth_date = date_leida_Hash,
    death_date = date_ptira_Hash,
    sex = factor(case_when(
      str_detect(teur_min, "×§×\u0091×\u0094") ~ "Female",
      str_detect(teur_min, "×\u0096×\u009b×¨") ~ "Male"
    )),
    sector = factor(case_when(
      str_detect(teur_ifyun_mirpaa, "×\u009b×\u009c×\u009c×\u0099") ~ "General",
      str_detect(teur_ifyun_mirpaa, "×\u0094×\u0097×¨×\u0093×\u0099×\u009d") ~ "Haredic",
      str_detect(teur_ifyun_mirpaa, "×¢×¨×\u0091×\u0099") ~ "Arab"
    )),
    sw_AB_RegishYr = sw_AB_RegishYr,
    sw_trans = ifelse(is.na(sw_trans), FALSE, TRUE),
    sw_LTCF = ifelse(is.na(sw_LTCF), FALSE, TRUE)
  )

# -----------------------------------------------------------------------------
# Save Processed Data to File
# -----------------------------------------------------------------------------

# Save the processed datasets into a single RData file for later use.
save(micro_df, micro_df_long, ab_consumption_all, cov_df,
     hospital_all, dmg_sw_all, 
     file = "bigfiles.RData")