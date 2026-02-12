# ================================================================= #
# STEP 05: LENIENCY ANALYSES (ELEVATION & DISCREPANCY)
# ================================================================= #

# --- 1. Load Packages and Helpers --------------------------------
source(here::here("R", "01_helpers.R"))
load_packages(c("here", "dplyr", "tidyr", "lme4", "lmerTest", "glmmTMB", "nlme", "knitr", "pwr"))

# --- 2. Load Cleaned Data ----------------------------------------
my_data <- load_cleaned_data()

# --- 3. Analysis Part 1: Elevation Leniency ----------------------
# Does the mean rating (elevation) differ by Condition (RPM vs. GRS)?
message("Running Elevation Leniency analysis...")

# Define the outcomes to test
elevation_outcomes <- c("Self_Mean", "Peer_Mean", "Super_Mean")

elevation_results <- list()

for (outcome in elevation_outcomes) {
  # Define the model formula
  formula <- as.formula(paste0(outcome, " ~ Condition"))

  # Identify potential random effects using ICC-guided helper
  random_effects_candidates <- build_random_effect_candidates(
    vars = outcome,
    data = my_data
  )

  # Fit the best model using the helper function
  model_fit <- fit_best_model(formula, data = my_data, random_effects = random_effects_candidates)

  # Store the result for the 'ConditionGRS' coefficient
  elevation_results[[outcome]] <- model_fit %>%
    dplyr::filter(term == "ConditionGRS") %>%
    dplyr::mutate(outcome_variable = outcome, .before = 1)
}

elevation_summary <- dplyr::bind_rows(elevation_results)

# --- 4. Analysis Part 2: Discrepancy-Based Leniency --------------
# Is the discrepancy from the supervisor (Self/Peer - Super) different from zero?
message("Running Discrepancy-Based Leniency analysis...")

discrepancy_outcomes <- c("Self_Discrep", "Peer_Discrep")
discrepancy_results <- list()

for (outcome in discrepancy_outcomes) {
  for (cond in c("RPM", "GRS")) {
    # Filter data for the specific condition
    df_cond <- dplyr::filter(my_data, Condition == cond)

    # Model formula is an intercept-only model
    formula <- as.formula(paste0(outcome, " ~ 1"))

    # Check ICCs within the condition using helper
    random_effects_candidates <- build_random_effect_candidates(
      vars = outcome,
      data = df_cond
    )

    # Fit the best model
    model_fit <- fit_best_model(formula, data = df_cond, random_effects = random_effects_candidates)

    # Store the result for the intercept
    discrepancy_results[[paste(outcome, cond)]] <- model_fit %>%
      dplyr::filter(term == "(Intercept)") %>%
      dplyr::mutate(outcome_variable = outcome, condition = cond, .before = 1)
  }
}

discrepancy_summary <- dplyr::bind_rows(discrepancy_results)

# --- 5. Power Analysis: Elevation Leniency -------------------
message("Computing achieved power for Elevation Leniency (GRS vs. RPM)...")

elevation_power_summary <- elevation_outcomes %>%
  lapply(function(outcome) compute_condition_power(my_data, outcome)) %>%
  dplyr::bind_rows() %>%
  dplyr::mutate(power_estimate = round(power_estimate, 3))

# --- 6. Power Analysis: Discrepancy-Based Leniency -----------
message("Computing achieved power for Discrepancy Leniency (GRS vs. RPM)...")

discrepancy_power_summary <- tidyr::expand_grid(
  outcome = discrepancy_outcomes,
  condition = c("RPM", "GRS")
) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    result = list(compute_onesample_power(
      dplyr::filter(my_data, Condition == condition), outcome, condition
    ))
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(result) %>%
  tidyr::unnest(result) %>%
  dplyr::mutate(power_estimate = round(power_estimate, 3))

# --- 7. Save and Display Results --------------------------------

elevation_output_file <- here::here("output", "tables", "05_elevation_leniency_summary.csv")
discrepancy_output_file <- here::here("output", "tables", "05_discrepancy_leniency_summary.csv")
elevation_power_file <- here::here("output", "tables", "05_elevation_power_summary.csv")
discrepancy_power_file <- here::here("output", "tables", "05_discrepancy_power_summary.csv")

message("Saving elevation leniency results to: ", elevation_output_file)
write.csv(elevation_summary, elevation_output_file, row.names = FALSE)

message("Saving discrepancy leniency results to: ", discrepancy_output_file)
write.csv(discrepancy_summary, discrepancy_output_file, row.names = FALSE)

message("Saving elevation power results to: ", elevation_power_file)
write.csv(elevation_power_summary, elevation_power_file, row.names = FALSE)

message("Saving discrepancy power results to: ", discrepancy_power_file)
write.csv(discrepancy_power_summary, discrepancy_power_file, row.names = FALSE)

# Display tables to the console and viewer
display_table(elevation_summary, "Elevation Leniency Summary (Effect of GRS vs. RPM)", digits = 3)
display_table(discrepancy_summary, "Discrepancy Leniency Summary (Test vs. Zero)", digits = 3)
display_table(elevation_power_summary, "Power Analysis for Elevation Leniency", digits = 3)
display_table(discrepancy_power_summary, "Power Analysis for Discrepancy Leniency", digits = 3)

message("Step 05: Leniency analyses complete.")
