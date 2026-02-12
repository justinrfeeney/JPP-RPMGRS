options(repos = c(CRAN = "https://cran.rstudio.com/"))

if (!requireNamespace("here", quietly = TRUE)) {
  install.packages("here", dependencies = TRUE)
}
library(here)

source(here::here("R", "01_helpers.R"))

load_packages(c(
  "here", "dplyr", "tidyr", "tibble"
))

start_time <- Sys.time()
message("Starting analysis pipeline at: ", start_time)

dir.create(here::here("output"), recursive = TRUE, showWarnings = FALSE)
dir.create(here::here("output", "tables"), recursive = TRUE, showWarnings = FALSE)

run_step <- function(label, path) {
  message("\n--- ", label, " ---")
  tryCatch(
    {
      source(here::here("R", path), echo = FALSE, max.deparse.length = Inf)
      message(label, " completed successfully.")
    },
    error = function(e) {
      message("ERROR in ", label, ": ", conditionMessage(e))
      stop(e)
    }
  )
}

run_step("STEP 02: Loading and Cleaning Data",          "02_load_clean.R")
run_step("STEP 03: ICC Analysis",                       "03_icc_analysis.R")
run_step("STEP 04: Confirmatory Factor Analysis (CFA)", "04_cfa_analysis.R")
run_step("STEP 05: Leniency Analyses",                  "05_leniency.R")
run_step("STEP 06: Correlation Analysis",               "06_correlations.R")

end_time   <- Sys.time()
time_taken <- end_time - start_time

message("\n========================================================")
message("Analysis pipeline complete!")
message("Total execution time: ", format(time_taken))
message("All outputs have been saved to the /output directory.")
message("========================================================")
