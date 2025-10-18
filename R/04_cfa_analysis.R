# ================================================================= #
# STEP 04: CONFIRMATORY FACTOR ANALYSIS (CFA)
# ================================================================= #

# --- 1. Load Packages and Helpers --------------------------------
source(here::here("R", "01_helpers.R"))
load_packages(c("here", "dplyr", "tidyr", "lavaan", "knitr"))

# --- 2. Load Cleaned Data ----------------------------------------
cleaned_data_path <- here::here("output", "cleaned_data.rds")
my_data <- readRDS(cleaned_data_path)

# --- 3. Define Model and Parameters ------------------------------

# Define the single-factor CFA model template.
# The indicators (dim1, dim2, etc.) will be replaced dynamically.
cfa_model_template <- '
  latent_factor =~ dim1 + dim2 + dim3 + dim4
'

# Define the sets of variables for each rater type
variable_sets <- list(
  Supervisor = c("Super_Organization", "Super_Physical", "Super_Visual", "Super_Vocal"),
  Peer = c("Peer_Organization", "Peer_Physical", "Peer_Visual", "Peer_Vocal")
)

# --- 4. Run Multilevel CFAs --------------------------------------- 
message("Running multilevel CFA...")

# Initialize lists to store results
all_fit_indices <- list()
all_loadings <- list()

# Loop over each condition (RPM, GRS)
for (cond in c("RPM", "GRS")) {
  df_cond <- dplyr::filter(my_data, Condition == cond)

  # Loop over each rater type (Supervisor, Peer)
  for (rater in names(variable_sets)) {
    vars <- variable_sets[[rater]]
    df_rater <- df_cond %>%
      dplyr::select(dplyr::all_of(c(vars, "Group", "Experimenter"))) %>%
      tidyr::drop_na()

    # Determine the appropriate cluster based on ICCs
    # We check if there is meaningful variance at the Group or Experimenter level
    icc_vals_group <- sapply(vars, get_icc1, cl = "Group", df = df_rater)
    icc_vals_experimenter <- sapply(vars, get_icc1, cl = "Experimenter", df = df_rater)

    icc_g <- if (all(is.na(icc_vals_group))) NA_real_ else max(icc_vals_group, na.rm = TRUE)
    icc_e <- if (all(is.na(icc_vals_experimenter))) NA_real_ else max(icc_vals_experimenter, na.rm = TRUE)

    cluster_var <- NULL
    if (!is.na(icc_g) && icc_g > 0.05) cluster_var <- "Group"
    # If Experimenter ICC is also high, you might reconsider, but for now, Group takes precedence
    else if (!is.na(icc_e) && icc_e > 0.05) cluster_var <- "Experimenter"

    # Define a unique model label
    model_label <- paste0(rater, "_", cond)

    # Build the specific model string for this iteration
    model_string <- cfa_model_template
    for (i in 1:length(vars)) {
      model_string <- gsub(paste0("dim", i), vars[i], model_string)
    }

    # Fit the CFA model using lavaan
    fit <- lavaan::cfa(
      model = model_string,
      data = df_rater,
      cluster = cluster_var, # `cluster` handles non-independence
      estimator = "MLR",       # Robust Maximum Likelihood for non-normality
      missing = "pairwise"   # Handle missing data
    )

    # --- Extract Results ---
    # 1. Fit Indices
    fit_measures <- lavaan::fitMeasures(fit, c("chisq.scaled", "df.scaled", "pvalue.scaled", "cfi.scaled", "tli.scaled", "rmsea.scaled", "srmr"))
    all_fit_indices[[model_label]] <- tibble::tibble(
        Model = model_label,
        Cluster = cluster_var %||% "None",
        chisq.scaled = fit_measures["chisq.scaled"],
        df.scaled = fit_measures["df.scaled"],
        pvalue.scaled = fit_measures["pvalue.scaled"],
        cfi.scaled = fit_measures["cfi.scaled"],
        tli.scaled = fit_measures["tli.scaled"],
        rmsea.scaled = fit_measures["rmsea.scaled"],
        srmr = fit_measures["srmr"]
    )

    # 2. Standardized Loadings
    solution <- lavaan::standardizedSolution(fit)
    all_loadings[[model_label]] <- solution %>%
      dplyr::filter(op == "=~") %>%
      dplyr::select(Indicator = rhs, Loading = est.std) %>%
      dplyr::mutate(Model = model_label, .before = 1)
  }
}

# Combine results from the loops into single dataframes
cfa_fit_summary <- dplyr::bind_rows(all_fit_indices)
cfa_loadings_summary <- dplyr::bind_rows(all_loadings)

# --- 5. Save and Display Results --------------------------------
fit_output_file <- here::here("output", "tables", "04_cfa_fit_summary.csv")
loadings_output_file <- here::here("output", "tables", "04_cfa_loadings_summary.csv")

message("Saving CFA fit indices to: ", fit_output_file)
write.csv(cfa_fit_summary, fit_output_file, row.names = FALSE)

message("Saving CFA loadings to: ", loadings_output_file)
write.csv(cfa_loadings_summary, loadings_output_file, row.names = FALSE)

# Display tables to the console and viewer
display_table(cfa_fit_summary, "CFA Fit Indices", digits = 3)
display_table(cfa_loadings_summary, "CFA Standardized Loadings", digits = 3)

message("Step 04: CFA analysis complete.")