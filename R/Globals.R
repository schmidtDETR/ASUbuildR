# globals.R
# This file defines global variables to avoid R CMD CHECK warnings
# These variables are typically column names used in data manipulation functions

# Suppress R CMD CHECK warnings for undefined global variables
utils::globalVariables(c(
  # ASU (Area of Substantial Unemployment) related variables
  "asunum",
  "original_asu",

  # Geographic identifiers
  "GEOID",
  "name",

  # Employment and labor force variables
  "tract_ASU_clf",      # Civilian Labor Force
  "tract_ASU_unemp",    # Unemployment
  "tract_ASU_emp",      # Employment
  "tract_ASU_urate",    # Unemployment Rate
  "tract_pop_cur",      # Current Population
  "Unemployment",       # Unemployment (capitalized version)
  "Employed",           # Employed (capitalized version)
  "unemp",             # Unemployment (abbreviated)
  "lf",                # Labor Force (abbreviated)
  "ur",                # Unemployment Rate (abbreviated)

  # Analysis variables
  "comp",              # Comparison variable
  "continuous",        # Continuous variable indicator
  "row_num"            # Row number variable
))
