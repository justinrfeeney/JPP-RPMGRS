# ================================================================= #
# ANALYSIS HELPER FUNCTIONS
# ================================================================= #

# ---- Package Loading ----------------------------------------------
#' Load required packages quietly.
#'
#' Installs any missing packages from CRAN before attaching them
#' without emitting startup messages.
#'
#' @param pkgs Character vector of package names.
#' @return Invisibly returns TRUE when all packages are attached.
load_packages <- function(pkgs) {
  if (length(pkgs) == 0) {
    return(invisible(TRUE))
  }
  if (!is.character(pkgs)) {
    stop("`pkgs` must be a character vector.")
  }
  pkgs <- unique(pkgs)
  new_pkgs <- setdiff(pkgs, rownames(installed.packages()))
  if (length(new_pkgs) > 0) {
    message("Installing missing packages: ", paste(new_pkgs, collapse = ", "))
    install.packages(new_pkgs, dependencies = TRUE, quiet = TRUE, repos = "https://cran.rstudio.com/")
  }
  suppressPackageStartupMessages({
    invisible(lapply(pkgs, require, character.only = TRUE))
  })
  invisible(TRUE)
}

# ---- Warning Suppression -------------------------------------------

#' Suppress known, non-critical warnings from model fitting routines.
#'
#' This helper muffles boundary/singularity convergence warnings that are
#' expected in small samples but otherwise harmless for downstream use.
suppress_known_model_warnings <- function(expr) {
  withCallingHandlers(
    expr,
    warning = function(w) {
      msg <- conditionMessage(w)
      if (grepl("boundary (singular) fit", msg, fixed = TRUE) ||
        grepl("Model convergence problem; non-positive-definite Hessian matrix", msg, fixed = TRUE)) {
        invokeRestart("muffleWarning")
      }
    }
  )
}


# ---- ICC Calculation ------------------------------------------------

#' Calculate Intraclass Correlation Coefficient (ICC1)
#' @param var A string for the variable name.
#' @param cl A string for the clustering variable name.
#' @param df The dataframe.
#' @return The ICC1 value as a numeric scalar.
get_icc1 <- function(var, cl, df) {
  if (!is.character(var) || !is.character(cl) || !is.data.frame(df)) stop("Invalid input types.")
  if (!var %in% names(df) || !cl %in% names(df)) stop("Variable or cluster not in dataframe.")

  cl_values <- stats::na.omit(df[[cl]])
  var_values <- stats::na.omit(df[[var]])
  if (length(unique(cl_values)) < 2 || length(var_values) < 2L) {
    return(NA_real_)
  }

  formula <- stats::as.formula(paste0("`", var, "` ~ 1 + (1|`", cl, "`)"))
  fm <- tryCatch(
    suppress_known_model_warnings(
      lme4::lmer(formula, data = df, REML = TRUE)
    ),
    error = function(e) NULL
  )

  if (is.null(fm) || lme4::isSingular(fm)) return(NA_real_)

  vc <- as.data.frame(lme4::VarCorr(fm))
  clvar <- vc$vcov[vc$grp == cl]
  res <- vc$vcov[vc$grp == "Residual"]

  if (length(clvar) == 0 || length(res) == 0) return(NA_real_)
  denom <- clvar + res
  if (!is.finite(denom) || denom <= 0) return(NA_real_)
  as.numeric(clvar / denom)
}


# ---- Random Effects Utilities --------------------------------------

#' Identify clustering variables that warrant random intercepts.
#'
#' @param vars Character vector of variable names to evaluate.
#' @param data Data frame containing the variables.
#' @param clusters Candidate clustering variables. Defaults to c("Group", "Experimenter").
#' @param threshold Numeric ICC1 threshold beyond which a random intercept is retained.
#' @return Character vector of cluster names meeting the inclusion threshold.
identify_random_effects <- function(vars, data, clusters = c("Group", "Experimenter"), threshold = 0.05) {
  if (!is.data.frame(data)) stop("`data` must be a data frame.")
  vars <- unique(vars)
  clusters <- unique(clusters)
  if (!all(vars %in% names(data))) stop("All `vars` must be present in `data`.")

  clusters <- clusters[clusters %in% names(data)]
  if (length(clusters) == 0) {
    return(character())
  }

  keep <- character()
  for (cl in clusters) {
    icc_vals <- vapply(vars, function(v) get_icc1(v, cl, data), numeric(1), USE.NAMES = FALSE)
    if (any(icc_vals > threshold, na.rm = TRUE)) {
      keep <- c(keep, cl)
    }
  }
  unique(keep)
}

#' Compose candidate random-effects terms ordered by complexity.
#'
#' @param vars Character vector of variable names to evaluate.
#' @param data Data frame containing the variables.
#' @param clusters Candidate clustering variables. Defaults to c("Group", "Experimenter").
#' @param threshold Numeric ICC1 threshold beyond which a random intercept is retained.
#' @return Character vector of random-effect expressions ordered from most to least complex.
build_random_effect_candidates <- function(vars, data, clusters = c("Group", "Experimenter"), threshold = 0.05) {
  selected <- identify_random_effects(vars = vars, data = data, clusters = clusters, threshold = threshold)
  if (length(selected) == 0) {
    return(character())
  }
  if (length(selected) == 1) {
    return(sprintf("(1|%s)", selected))
  }
  combos <- paste(sprintf("(1|%s)", selected), collapse = " + ")
  c(combos, sprintf("(1|%s)", selected))
}


# ---- Statistical Tests ----------------------------------------------

#' Compare two independent correlations using Steiger's Z-test.
#' @param r1 Correlation coefficient of group 1.
#' @param n1 Sample size of group 1.
#' @param r2 Correlation coefficient of group 2.
#' @param n2 Sample size of group 2.
#' @return A list with the Z-statistic (`Z`) and p-value (`p`).
steiger_test <- function(r1, n1, r2, n2) {
  if (any(is.na(c(r1, r2, n1, n2)))) {
    return(list(Z = NA_real_, p = NA_real_))
  }
  if (any(c(n1, n2) <= 3)) {
    warning("Steiger's test requires sample sizes greater than 3 per group; returning NA.")
    return(list(Z = NA_real_, p = NA_real_))
  }

  clamp_r <- function(x) {
    pmin(0.999999, pmax(-0.999999, x))
  }
  r1 <- clamp_r(r1)
  r2 <- clamp_r(r2)

  z1 <- atanh(r1)
  z2 <- atanh(r2)
  se <- sqrt(1 / (n1 - 3) + 1 / (n2 - 3))
  if (!is.finite(se) || se <= 0) {
    return(list(Z = NA_real_, p = NA_real_))
  }
  z_stat <- (z1 - z2) / se
  p_val  <- 2 * (1 - pnorm(abs(z_stat)))
  list(Z = z_stat, p = p_val)
}


# ---- Model Fitting --------------------------------------------------

#' Fit the best-fitting mixed-effects model from a set of candidates.
#' @param formula A formula for the fixed effects.
#' @param data The dataframe.
#' @param random_effects A character vector of random effect terms.
#' @return A tibble summarizing the best-fitting model.
fit_best_model <- function(formula, data, random_effects = NULL) {
  if (!inherits(formula, "formula")) {
    if (is.character(formula) && length(formula) == 1) {
      formula <- stats::as.formula(formula)
    } else {
      stop("`formula` must be a formula or a length-one character string.")
    }
  }
  if (!is.data.frame(data)) stop("`data` must be a data frame.")
  if (!is.null(random_effects)) {
    random_effects <- unique(as.character(random_effects))
  }

  if (!is.null(random_effects) && length(random_effects) > 0) {
    for (re in random_effects) {
      full_formula <- stats::update(formula, paste0(". ~ . + ", re))

      fit_glmmTMB <- try(
        suppress_known_model_warnings(
          glmmTMB::glmmTMB(
            full_formula, data = data, family = stats::gaussian(), REML = FALSE,
            control = glmmTMB::glmmTMBControl(optimizer = stats::optim, optArgs = list(method = "BFGS"))
          )
        ),
        silent = TRUE
      )

      if (!inherits(fit_glmmTMB, "try-error") && isTRUE(fit_glmmTMB$sdr$pdHess)) {
        summary_fit <- summary(fit_glmmTMB)$coefficients$cond
        if (is.null(summary_fit)) {
          next
        }
        result <- as.data.frame(summary_fit) %>%
          tibble::rownames_to_column("term") %>%
          dplyr::select(term, estimate = "Estimate", std.error = "Std. Error", statistic = "z value", p.value = "Pr(>|z|)") %>%
          tibble::as_tibble()
        return(result %>% dplyr::mutate(model_engine = "glmmTMB", random_effect = re))
      }

      if (length(strsplit(re, " +", fixed = TRUE)[[1]]) == 1 && startsWith(re, "(1|") && endsWith(re, ")")) {
          grouping_var <- substring(re, 4, nchar(re) - 1)
          random_formula <- stats::as.formula(paste0("~ 1 | ", grouping_var))

          fit_lme <- try(
            nlme::lme(
              fixed = formula, random = random_formula, data = data,
              na.action = na.omit, method = "ML"
            ), silent = TRUE
          )

          if (!inherits(fit_lme, "try-error")) {
            summary_fit <- summary(fit_lme)$tTable
            result <- as.data.frame(summary_fit) %>%
              tibble::rownames_to_column("term") %>%
              dplyr::select(term, estimate = "Value", std.error = "Std.Error", statistic = "t-value", p.value = "p-value") %>%
              tibble::as_tibble()
            return(result %>% dplyr::mutate(model_engine = "nlme", random_effect = re))
          }
      }
    }
  }

  fit_lm <- stats::lm(formula, data = data, na.action = stats::na.omit)
  summary_fit <- summary(fit_lm)$coefficients
  result <- as.data.frame(summary_fit) %>%
    tibble::rownames_to_column("term") %>%
    dplyr::select(term, estimate = "Estimate", std.error = "Std. Error", statistic = "t value", p.value = "Pr(>|t|)") %>%
    tibble::as_tibble()
  return(result %>% dplyr::mutate(model_engine = "lm", random_effect = "None"))
}


# ---- Power Analysis Utilities --------------------------------------

#' Fit a Gaussian model suitable for simulation-based power analysis.
#'
#' Attempts to fit the supplied random-effect structures (from most to least
#' complex) with `lme4::lmer()`. Falls back to `stats::lm()` when none of the
#' structures converge without singularity.
#'
#' @param formula A model formula or character string describing the fixed effects.
#' @param data Data frame used for the fit.
#' @param random_effects Character vector of random-effect terms ordered by preference.
#' @param na_action NA handling function passed to the model fitting call.
#' @return A list with elements `model` (the fitted object) and `random_effect` (the term used).
fit_power_model <- function(formula, data, random_effects = character(), na_action = stats::na.omit) {
  if (!inherits(formula, "formula")) {
    formula <- stats::as.formula(formula)
  }
  if (!is.data.frame(data)) stop("`data` must be a data frame.")
  random_effects <- unique(as.character(random_effects))

  formula_text <- paste(deparse(formula), collapse = "")

  lmer_fn <- if (requireNamespace("lmerTest", quietly = TRUE)) lmerTest::lmer else lme4::lmer
  if (length(random_effects) > 0) {
    for (re in random_effects) {
      re <- trimws(re)
      full_formula <- stats::as.formula(paste(formula_text, "+", re))
      attempt <- try(
        suppress_known_model_warnings(
          lmer_fn(
            full_formula,
            data = data,
            REML = FALSE,
            na.action = na_action,
            control = lme4::lmerControl(
              optimizer = "bobyqa",
              optCtrl = list(maxfun = 200000),
              check.conv.singular = "ignore"
            )
          )
        ),
        silent = TRUE
      )
      if (!inherits(attempt, "try-error") && !lme4::isSingular(attempt, tol = 1e-6)) {
        return(list(model = attempt, random_effect = re))
      }
    }
  }

  fit_lm <- stats::lm(formula, data = data, na.action = na_action)
  list(model = fit_lm, random_effect = "None")
}

#' Convert a `simr::powerSim` result to a tidy tibble.
#'
#' @param power_object Output from `simr::powerSim()`.
#' @param outcome Label for the outcome variable.
#' @param effect Label for the effect being tested.
#' @param model_label Character description of the fitted model class.
#' @param random_effect Random-effect structure used for the fit.
#' @param test_method Hypothesis test used within `simr::fixed()`.
#' @param ... Additional named columns to append to the returned tibble.
#' @return A tibble summarising the simulation-based power estimate.
summarize_power_sim <- function(power_object, outcome, effect, model_label, random_effect, test_method, ...) {
  if (missing(power_object) || is.null(power_object)) {
    stop("`power_object` must be supplied.")
  }
  summ <- tryCatch(summary(power_object), error = function(e) stop("Unable to summarise `power_object`: ", e$message, call. = FALSE))
  power_col <- if ("power" %in% names(summ)) summ$power else summ$mean
  trials_col <- if ("nsim" %in% names(summ)) summ$nsim else summ$trials
  lower_col <- if ("lower" %in% names(summ)) summ$lower else NA_real_
  upper_col <- if ("upper" %in% names(summ)) summ$upper else NA_real_

  base_tbl <- tibble::tibble(
    outcome = outcome,
    effect = effect,
    model = model_label,
    random_effect = random_effect,
    test = test_method,
    nsim = as.integer(trials_col),
    power = as.numeric(power_col),
    lower = as.numeric(lower_col),
    upper = as.numeric(upper_col)
  )
  extra <- list(...)
  if (length(extra) > 0) {
    extra_tbl <- tibble::as_tibble(extra)
    return(dplyr::bind_cols(base_tbl, extra_tbl))
  }
  base_tbl
}

# ---- Reporting Utilities -----------------------------------------
display_table <- function(df, title, digits = 3) {
  if (!is.data.frame(df)) stop("`df` must be a data frame.")
  if (length(digits) != 1 || !is.numeric(digits)) stop("`digits` must be a single numeric value.")
  digits <- as.integer(digits)
  if (is.na(digits) || digits < 0) stop("`digits` must be a non-negative integer.")
  title <- as.character(title)[1]

  tbl <- df %>%
    tibble::as_tibble() %>%
    dplyr::mutate(dplyr::across(where(is.numeric), ~round(.x, digits)))

  cat("\n--- ", title, " ---\n", sep = "")
  print(knitr::kable(tbl, digits = digits))
  cat("\n")

  if (interactive()) {
    try(utils::View(tbl, title), silent = TRUE)
  }

  invisible(tbl)
}
