# ================================================================= #
# STEP 02: LOAD, CLEAN, AND PREPARE DATA
# ================================================================= #

# --- 1. Load Packages --------------------------------------------
# Source the helper script to get the package loader function
source(here::here("R", "01_helpers.R"))
load_packages(c("here", "haven", "dplyr", "tidyr"))

# --- 2. Define File Paths ----------------------------------------
# Using here() makes file paths robust to where the project is located
raw_data_path <- here::here("data", "data.sav")
output_path <- here::here("output", "cleaned_data.rds")

# --- 3. Load Raw Data --------------------------------------------
message("Loading raw data from: ", raw_data_path)
my_data_raw <- haven::read_sav(raw_data_path)

# --- 4. Clean and Prepare Data -----------------------------------
message("Cleaning and preparing data...")
my_data_clean <- my_data_raw %>%
  # Rename variables for clarity and consistency
  dplyr::rename(Group = GroupID) %>%
  # Mutate variables into the correct format
  dplyr::mutate(
    # Convert clustering variables to factors
    Group = factor(Group),
    Experimenter = factor(Experimenter),

    # Convert Condition to a factor with meaningful labels and a specific reference level
    Condition = factor(Condition, levels = c(1, 2), labels = c("RPM", "GRS")),
    Condition = relevel(Condition, ref = "RPM"),

    # --- Create Composite Scores (Means) ---
    # Calculate mean scores for each rater type across all dimensions
    Self_Mean = rowMeans(dplyr::across(dplyr::starts_with("Self_")), na.rm = TRUE),
    Peer_Mean = rowMeans(dplyr::across(dplyr::starts_with("Peer_")), na.rm = TRUE),
    Super_Mean = rowMeans(dplyr::across(dplyr::starts_with("Super_")), na.rm = TRUE),

    # --- Create Discrepancy Scores ---
    # Calculate the difference between self/peer ratings and supervisor ratings
    Peer_Discrep = Peer_Mean - Super_Mean,
    Self_Discrep = Self_Mean - Super_Mean
  )

# --- 5. Save Cleaned Data ----------------------------------------
message("Saving cleaned data to: ", output_path)
saveRDS(my_data_clean, file = output_path)

message("Step 02: Data loading and cleaning complete. ", nrow(my_data_clean), " rows processed.")