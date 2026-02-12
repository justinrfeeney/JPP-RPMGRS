# ================================================================= #
# STEP 03: INTRACLASS CORRELATION (ICC) ANALYSIS
# ================================================================= #

# --- 1. Load Packages and Helpers --------------------------------
# Use here::here to find the project root and build paths reliably
source(here::here("R", "01_helpers.R"))
load_packages(c("here", "dplyr", "tidyr", "lme4", "knitr"))

# --- 2. Load Cleaned Data ----------------------------------------
my_data <- load_cleaned_data()

# --- 3. Define Analysis Parameters -------------------------------
# Variables for which to calculate ICCs
vars_to_analyze <- c("Self_Mean", "Peer_Mean", "Super_Mean")
# Clustering variables
clusters <- c("Group", "Experimenter")
# Conditions to loop through
conditions <- c("RPM", "GRS")

# --- 4. Perform ICC Analysis ------------------------------------- 
message("Running ICC analysis...")

# Create a grid of all combinations to analyze
icc_grid <- tidyr::expand_grid(
  Variable = vars_to_analyze,
  Cluster = clusters,
  Condition = conditions
)

# Calculate ICC for each row in the grid
icc_results <- icc_grid %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    ICC1 = get_icc1(Variable, Cluster, dplyr::filter(my_data, Condition == .data$Condition))
  ) %>%
  dplyr::ungroup()

# --- 5. Save and Display Results --------------------------------
output_file <- here::here("output", "tables", "03_icc_summary.csv")
message("Saving ICC results to: ", output_file)

# Save the results to a CSV file
write.csv(icc_results, output_file, row.names = FALSE)

# Display a rounded summary in the console and viewer
display_table(icc_results, "ICC(1) Summary", digits = 3)

message("Step 03: ICC analysis complete.")
