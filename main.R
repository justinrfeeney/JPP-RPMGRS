# ================================================================= #
# MASTER ANALYSIS SCRIPT
# ================================================================= #
#
# This script runs the entire analysis pipeline for the project.
# It sources each step in order, from data cleaning to final analysis.
# To run the entire analysis, simply execute this file.
#
# ================================================================= #

# --- 1. Preamble -------------------------------------------------

# Globally set the CRAN mirror to avoid interactive prompts
options(repos = c(CRAN = "https://cran.rstudio.com/"))

# Set a global option to use the `here` package for path management# This makes the project portable and avoids `setwd()` issues.
if (!requireNamespace("here", quietly = TRUE)) {
  install.packages("here", dependencies = TRUE)
}
library(here)

# Record start time
start_time <- Sys.time()
message("Starting analysis pipeline at: ", start_time)

# Ensure output directories exist
dir.create(here::here("output"), recursive = TRUE, showWarnings = FALSE)
dir.create(here::here("output", "tables"), recursive = TRUE, showWarnings = FALSE)

# --- 2. Execute Analysis Steps ----------------------------------- 

# Each script prints messages about its progress and saves its output
# to the `/output` directory.

message("\n--- STEP 02: Loading and Cleaning Data ---")
source(here::here("R", "02_load_clean.R"))

message("\n--- STEP 03: ICC Analysis ---")
source(here::here("R", "03_icc_analysis.R"))

message("\n--- STEP 04: Confirmatory Factor Analysis (CFA) ---")
source(here::here("R", "04_cfa_analysis.R"))

message("\n--- STEP 05: Leniency Analyses ---")
source(here::here("R", "05_leniency.R"))

message("\n--- STEP 06: Correlation Analysis ---")
source(here::here("R", "06_correlations.R"))

# --- 3. Conclusion ----------------------------------------------- 

end_time <- Sys.time()
time_taken <- end_time - start_time

message("\n========================================================")
message("Analysis pipeline complete!")
message("Total execution time: ", format(time_taken))
message("All outputs have been saved to the /output directory.")
message("========================================================")
