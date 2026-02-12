# ================================================================= #
# STEP 02: LOAD, CLEAN, AND PREPARE DATA
# ================================================================= #

source(here::here("R", "01_helpers.R"))
load_packages(c("here", "haven", "dplyr", "tidyr", "tibble"))

raw_data_path <- here::here("data", "data.sav")
output_path   <- here::here("output", "cleaned_data.rds")

message("Step 02: starting load and clean.")
message("Raw data path: ", raw_data_path)

if (!file.exists(raw_data_path)) {
  stop("Raw data file not found at: ", raw_data_path)
}

read_spss_safe <- function(path) {
  out <- tryCatch(
    {
      message("Reading SPSS file with haven::read_sav() ...")
      haven::read_sav(path)
    },
    error = function(e1) {
      message("haven::read_sav() failed: ", conditionMessage(e1))
      message("Falling back to foreign::read.spss() ...")
      load_packages("foreign")
      tryCatch(
        {
          foreign::read.spss(path, to.data.frame = TRUE, use.value.labels = TRUE)
        },
        error = function(e2) {
          stop(
            "Failed to read SPSS file with both haven::read_sav() and foreign::read.spss(). ",
            "First error: ", conditionMessage(e1), " | Second error: ", conditionMessage(e2)
          )
        }
      )
    }
  )
  tibble::as_tibble(out)
}

my_data_raw <- read_spss_safe(raw_data_path)

message("Raw data loaded. Rows: ", nrow(my_data_raw), " | Cols: ", ncol(my_data_raw))

my_data_clean <- my_data_raw %>%
  dplyr::rename(Group = GroupID) %>%
  dplyr::mutate(
    Group        = factor(Group),
    Experimenter = factor(Experimenter),
    Condition    = factor(Condition, levels = c(1, 2), labels = c("RPM", "GRS"))
  ) %>%
  dplyr::mutate(
    Condition    = stats::relevel(Condition, ref = "RPM"),
    Self_Mean    = rowMeans(dplyr::pick(Self_Organization, Self_Physical, Self_Visual, Self_Vocal),  na.rm = TRUE),
    Peer_Mean    = rowMeans(dplyr::pick(Peer_Organization, Peer_Physical, Peer_Visual, Peer_Vocal),  na.rm = TRUE),
    Super_Mean   = rowMeans(dplyr::pick(Super_Organization, Super_Physical, Super_Visual, Super_Vocal), na.rm = TRUE),
    Peer_Discrep = Peer_Mean - Super_Mean,
    Self_Discrep = Self_Mean - Super_Mean
  )

message("Cleaning complete. Rows: ", nrow(my_data_clean), " | Cols: ", ncol(my_data_clean))

dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
saveRDS(my_data_clean, file = output_path)

message("Step 02 complete. Cleaned data saved to: ", output_path)
