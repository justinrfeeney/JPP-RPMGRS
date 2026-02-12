# ================================================================= #
# STEP 06: CORRELATION ANALYSIS (CONVERGENT VALIDITY)
# ================================================================= #

# --- 1. Load Packages and Helpers --------------------------------
source(here::here("R", "01_helpers.R"))
load_packages(c("here", "dplyr", "tidyr", "correlation", "knitr", "pwr"))

# --- 2. Load Cleaned Data ----------------------------------------
my_data <- load_cleaned_data()

# --- 3. Define Analysis Parameters -------------------------------

# Define the pairs of variables to correlate
pairs_to_correlate <- list(
  "Self-Peer" = c("Self_Mean", "Peer_Mean"),
  "Self-Super" = c("Self_Mean", "Super_Mean"),
  "Peer-Super" = c("Peer_Mean", "Super_Mean")
)

# --- 4. Run Correlation Analysis --------------------------------
message("Running multilevel correlation analysis...")

correlation_results <- list()

# Loop over each pair of raters
for (pair_name in names(pairs_to_correlate)) {
  vars <- pairs_to_correlate[[pair_name]]
  v1 <- vars[1]
  v2 <- vars[2]

  # Loop over each condition
  for (cond in c("RPM", "GRS")) {
    df_cond <- my_data %>%
      dplyr::filter(Condition == cond) %>%
      dplyr::select(dplyr::all_of(c(v1, v2, "Group", "Experimenter"))) %>%
      tidyr::drop_na()

    # Determine random effects structure using helper to mirror primary analyses
    random_effect_candidates <- build_random_effect_candidates(
      vars = c(v1, v2),
      data = df_cond
    )
    random_effects <- if (length(random_effect_candidates) > 0) random_effect_candidates[1] else NULL

    # Calculate the correlation
    corr_test <- suppress_known_model_warnings(
      correlation::correlation(
        data = df_cond,
        vars = v1,
        vars2 = v2,
        multilevel = !is.null(random_effects),
        random_effects = random_effects
      )
    )

    # Store the results with achieved power for the correlation
    corr_tbl <- corr_test %>%
      tibble::as_tibble() %>%
      dplyr::select(r, p, n_Obs) %>%
      dplyr::mutate(
        Pair = pair_name,
        Condition = cond,
        power_condition = compute_correlation_power(r, n_Obs)
      ) %>%
      dplyr::relocate(Pair, Condition)

    correlation_results[[paste(pair_name, cond)]] <- corr_tbl
  }
}

correlation_summary <- dplyr::bind_rows(correlation_results)

# --- 5. Compare Correlations with Steiger's Z-test ------------
message("Comparing correlations between conditions...")

# Pivot the data to have RPM and GRS side-by-side
corr_wide <- correlation_summary %>%
  tidyr::pivot_wider(
    names_from = Condition,
    values_from = c(r, p, n_Obs, power_condition)
  )

# Apply the Steiger test and compute achieved power for each pair
comparison_results <- corr_wide %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    steiger_result = list(steiger_test(r_RPM, n_Obs_RPM, r_GRS, n_Obs_GRS)),
    steiger_z = steiger_result[["Z"]],
    steiger_p = steiger_result[["p"]],
    z_crit = qnorm(1 - 0.05 / 2),
    power_diff = 1 - pnorm(z_crit - abs(steiger_z)) + pnorm(-z_crit - abs(steiger_z))
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(-z_crit, -steiger_result)

# --- 6. Save and Display Results --------------------------------
output_file <- here::here("output", "tables", "06_correlation_summary.csv")
message("Saving correlation results to: ", output_file)
write.csv(comparison_results, output_file, row.names = FALSE)

# Display table to the console and viewer
display_table(comparison_results, "Correlation Summary and Comparison", digits = 3)

message("Step 06: Correlation analysis complete.")
