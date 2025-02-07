# -----------------------------------------------------------------------------
# Analysis of Antibiotic Consumption and Microbiology Data
# Supplemantary Tables
#
# This script performs several analyses on antibiotic consumption and
# microbiology datasets to:
#
#   1. Inspect the structure of the long-format microbiology data.
#   2. Identify and summarize cases where patients received two or more
#      antibiotics. The script computes the combination of antibiotics used,
#      counts the frequency and proportion of each combination, and prints the
#      top combinations.
#   3. For urine samples with a negative gram stain, determine the presence of
#      specific bacteria (Escherichia coli, Klebsiella pneumoniae, Proteus mirabilis)
#      for each patient and sample date, then summarize these frequencies.
#   4. For the same urine samples, assess whether antibiotic resistance testing
#      was performed (i.e. outcomes not labeled "Not-tested") for various antibiotics.
#      The results are summarized as counts and proportions.
#   5. Define a helper function (`overlap_table`) to calculate the overlap (or
#      complement) in testing outcomes between a reference antibiotic and other
#      antibiotic resistance outcomes.
#
# Dependencies: tidyverse (for data manipulation and piping)
#
# Author: [Your Name]
# Date: [Today's Date]
# -----------------------------------------------------------------------------

# Inspect the structure of the long-format microbiology data
glimpse(micro_df_long)

# -----------------------------------------------------------------------------
# Analyze Cases with Two or More Antibiotics
# -----------------------------------------------------------------------------
ab_consumption_all %>%
  # Create a new variable "sumVar" counting the number of antibiotic indicators
  mutate(sumVar = rowSums(select(., contains("ab_first")))) %>%
  # Keep only rows where more than one antibiotic was consumed
  filter(sumVar > 1) %>%
  # Remove columns not needed for the combination summary
  select(-c("id", "age", "ab_date", "sumVar")) %>%
  rowwise() %>%
  # Create a string listing all antibiotics consumed (removing the "ab_first_" prefix)
  mutate(
    ab_combination = paste(gsub("ab_first_", "", names(.)[c_across(starts_with("ab_first")) == TRUE]), collapse = ", ")
  ) %>%
  ungroup() -> duble_ab_sum

# Clean and summarize the antibiotic combinations
duble_ab_sum %>%
  # Remove any unwanted text (e.g., "ab_date, ") from the combination string
  mutate(ab_combination = str_remove_all(ab_combination, "ab_date, ")) %>%
  # Count the frequency of each antibiotic combination
  count(ab_combination) %>%
  # Calculate total records and proportion for each combination
  mutate(total = nrow(ab_consumption_all),
         proportion = 100 * n / nrow(ab_consumption_all)) %>%
  # Sort the results by decreasing proportion
  arrange(desc(proportion)) %>%
  # Print the top 50 combinations
  print(n = 50)

# -----------------------------------------------------------------------------
# Summarize Microbiology Data for Urine Samples (Bacterial Presence)
# -----------------------------------------------------------------------------
micro_df %>%
  # Filter to include only urine samples with a negative gram stain
  filter(micro_gram == "Negative",
         sample_source == "Urine") %>%
  # For each patient (id) and sample date, determine if specific bacteria are present
  group_by(id, sample_date) %>%
  summarise(
    ecoli = any(micro_name == "Escherichia coli"),
    kleb = any(micro_name == "Klebsiella pneumoniae"),
    prot = any(micro_name == "Proteus mirabilis")
  ) -> all_micro

# Summarize bacterial frequencies, including an 'other' category when none are present
all_micro %>%
  mutate(other = !ecoli & !kleb & !prot) %>%
  pivot_longer(-c(id, sample_date)) %>%  # Convert data to long format for each bacteria
  group_by(name) %>%
  summarise(
    n = sum(value, na.rm = TRUE),  # Count of positive cases per bacteria
    total = n(),                   # Total observations for that group
    p = round(100 * n / total, 2)  # Percentage of positive cases
  )

# -----------------------------------------------------------------------------
# Summarize Antibiotic Resistance Testing in Urine Samples
# -----------------------------------------------------------------------------
micro_df %>%
  # Filter to include only urine samples with a negative gram stain
  filter(micro_gram == "Negative",
         sample_source == "Urine") %>%
  group_by(id, sample_date) %>%
  # For each sample, check if resistance testing was performed (i.e. outcome not "Not-tested")
  summarise(
    cef_1 = any(outcome_sample_ab_res_cef_1 != "Not-tested"),
    cef_2 = any(outcome_sample_ab_res_cef_2 != "Not-tested"),
    cef_3 = any(outcome_sample_ab_res_cef_3 != "Not-tested"),
    augmentin = any(outcome_sample_ab_res_augmentin != "Not-tested"),
    fluoroquinolones = any(outcome_sample_ab_res_fluoroquinolones != "Not-tested"),
    nitrofurantoin = any(outcome_sample_ab_res_nitrofurantoin != "Not-tested")
  ) -> all_sens

# Summarize the frequency and percentage of tests performed for each antibiotic
all_sens %>%
  pivot_longer(-c(id, sample_date)) %>%
  group_by(name) %>%
  summarise(
    n = sum(value, na.rm = TRUE),
    total = n(),
    p = round(100 * n / total, 2)
  )

# -----------------------------------------------------------------------------
# Define a Function to Calculate Overlap Between Antibiotic Outcomes
# -----------------------------------------------------------------------------
overlap_table <- function(df, c_ab = "cef_3") {
  # Identify all antibiotic outcome columns, excluding id, sample_date, and the reference antibiotic (c_ab)
  antibiotic_cols <- colnames(df)[!colnames(df) %in% c("id", "sample_date", c_ab)]
  tibble(
    ab = antibiotic_cols,
    # For each antibiotic column, compute the frequency from a proportion table comparing
    # the reference outcome (c_ab) to the current antibiotic outcome.
    freq = sapply(antibiotic_cols, function(x) {
      prop_table <- table(df[[c_ab]], df[[x]]) %>%
        prop.table() %>%
        as.data.frame()
      # Extract a specific cell from the table (here, row 3, column 3) as the frequency.
      return(prop_table[3, 3])
    })
  ) %>%
    # Calculate the complement frequency (i.e., 100 - frequency in percentage)
    mutate(freq = 100 - round(100 * freq, 2))
}

# Apply the overlap_table function to the antibiotic resistance data
overlap_table(all_sens)