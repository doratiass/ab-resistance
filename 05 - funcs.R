# -----------------------------------------------------------------------------
# Variable Labeling Utility Functions
#
# This script loads a variable dictionary from a CSV file and defines several
# helper functions to map variable names to more descriptive labels and vice
# versa. These functions are useful for data cleaning, recoding, and for
# enhancing the readability of output tables and plots.
#
# Functions:
#   - label_get: Returns a human-readable label for a given variable name.
#   - var_get: Returns the original variable name given a human-readable label.
#   - level_get: Converts coded level values to a more interpretable format.
#   - label_all: Constructs a combined label from a variable name and an associated
#                level, handling special cases (e.g., "other_combined", polynomial terms).
#   - vars_label: Applies the label_get function to a vector of variable names.
#
# Dependencies: tidyverse (for string and data manipulation functions)
# -----------------------------------------------------------------------------

# Load the tidyverse package for data manipulation and string processing
library(tidyverse)

# Read the variable dictionary CSV file, which maps variable codes to descriptive names.
# The file "vars_dict.csv" should contain at least two columns: 'var' and 'name'.
vars_dict <- read_csv("vars_dict.csv", show_col_types = FALSE)

# -----------------------------------------------------------------------------
# label_get: Map a variable code to its descriptive label.
#
# If the input x is found in the 'var' column of vars_dict, return the corresponding
# 'name' (i.e., the human-readable label); otherwise, return x as-is.
label_get <- function(x) {
  ifelse(x %in% vars_dict$var,
         vars_dict[vars_dict$var == x, "name", drop = TRUE],
         x)
}

# -----------------------------------------------------------------------------
# var_get: Map a descriptive label back to its variable code.
#
# If the input x is found in the 'name' column of vars_dict, return the corresponding
# 'var' value; otherwise, return x as-is.
var_get <- function(x) {
  ifelse(x %in% vars_dict$name,
         vars_dict[vars_dict$name == x, "var", drop = TRUE],
         x)
}

# -----------------------------------------------------------------------------
# level_get: Convert coded level values to a more interpretable format.
#
# This function processes strings that represent level values:
#   - For values matching the pattern "X[0-9][0-9].[0-9][0-9]", it removes the leading "X"
#     and replaces dots with dashes.
#   - For values matching "X[0-9][0-9].", it replaces dots with pluses.
#   - For values containing an underscore, it replaces underscores with spaces.
#   - Otherwise, it returns the input unchanged.
level_get <- function(x) {
  if (str_detect(x, "X[0-9][0-9]\\.[0-9][0-9]")) {
    str_replace_all(str_remove(x, "X"), "\\.", "-")
  } else if (str_detect(x, "X[0-9][0-9]\\.")) {
    str_replace_all(str_remove(x, "X"), "\\.", "+")
  } else if (str_detect(x, "_")) {
    str_replace_all(x, "_", " ")
  } else {
    x
  }
}

# -----------------------------------------------------------------------------
# label_all: Create a combined label from a variable name and an associated level.
#
# This function handles different cases:
#   - If the input x contains the separator "_@_", it splits the string.
#     * If the second part equals "other_combined" or contains "TRUE", only the
#       label for the first part is returned.
#     * If the second part contains a dot, it is reformatted by replacing dots with spaces.
#     * Otherwise, the label for the first part is combined with the processed second part.
#   - If x contains "_poly_", it creates a label indicating a polynomial term.
#   - Otherwise, it simply applies label_get to x.
label_all <- function(x) {
  if (str_detect(x, "_@_")) {
    # Split the string on the delimiter "_@_"
    parts <- str_split(x, "_@_", simplify = TRUE)
    if (parts[,2] == "other_combined") {
      label_get(parts[,1])
    } else if (str_detect(parts[,2], "TRUE")) {
      label_get(parts[,1])
    } else if (str_detect(parts[,2], "\\.")) {
      paste0(label_get(parts[,1]), " - ",
             str_replace_all(parts[,2], "\\.", " "))
    } else {
      paste0(label_get(parts[,1]), " - ",
             level_get(parts[,2]))
    }
  } else if (str_detect(x, "_poly_")) {
    # Handle polynomial terms by splitting on "_poly_"
    parts <- str_split(x, "_poly_", simplify = TRUE)
    paste0(label_get(parts[,1]), " - poly ",
           label_get(parts[,2]))
  } else {
    label_get(x)
  }
}

# -----------------------------------------------------------------------------
# vars_label: Apply the label_get function to a vector of variable names.
#
# This function takes a vector x and returns a vector of labels by applying label_get
# to each element.
vars_label <- function(x) {
  sapply(x, label_get, USE.NAMES = FALSE)
}