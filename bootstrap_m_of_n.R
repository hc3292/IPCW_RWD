################################################################################
# m-of-n Bootstrap for IPCW + IPTW Workflow
# 
# This script performs m-of-n bootstrap resampling to estimate standard errors
# for Hazard Ratio (HR) estimates from the combined IPTW + IPCW workflow.
#
# Workflow:
# 1. find_ps_legend_script.R - Creates PS and fits outcome model with IPTW
# 2. censoring_vars_legend_script.R - Fits Cox censoring model
# 3. find_kZ_script.R - Computes IPCW weights
# 4. weights_script.R - Combines weights and fits final models
#
# Also computes bootstrap-based RMST and survival curve confidence intervals
################################################################################

library(DatabaseConnector)
library(CohortMethod)
library(dplyr)
library(tidyr)
library(tibble)
library(data.table)
library(survival)
library(Cyclops)
library(progress)
library(digest)  # For hash computation

################################################################################
# Configuration
################################################################################

MASTER_SEED <- 20260123

# Bootstrap parameters
B <- 200                    # Number of bootstrap iterations
m <- 10000                    # Bootstrap sample size (NULL = use full sample size n)
alpha <- 0.05                # For confidence intervals (1-alpha)*100%

# RMST and survival curve parameters
tau <- 2500                  # Restriction time for RMST (days)
# time_grid <- seq(0, tau, by = 60)  # Time points for survival curves (every 60 days up to tau)
# for longer time for survival plotting
time_grid <- c(
  seq(0, tau, by = 60),
  seq(max(tau, tail(seq(0, tau, by = 60), 1)) + 200, 5000, by = 200)
)

# Data file paths
cohortMethodData_file <- "results/cohortMethodData_t1788868_c1788867_o1788866.zip"
cohortMethodData_allcovar_file <- "results/cohortMethodData_t1788868_c1788867_o1788866_allcovar.zip"

# Output directory for bootstrap results
output_dir <- "bootstrap_results"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Temporary file names (will be suffixed with bootstrap iteration)
temp_ps_file <- "ps_study_boot.rds"
temp_censoring_file <- "Cox_censoring_boot.rds"
temp_weights_file <- "survival_weights_boot.csv"


################################################################################
# Helper Functions
################################################################################

#' Compute RMST from survfit object
#' @param survfit_obj survfit object with two treatment groups
#' @param tau Restriction time
#' @return List with RMST for each group and difference
compute_rmst_from_survfit <- function(survfit_obj, tau) {
  tryCatch({
    # Extract summary with full time grid
    surv_summary <- summary(survfit_obj, extend = TRUE)
    
    # Identify strata labels
    strata_names <- unique(surv_summary$strata)
    if (length(strata_names) != 2) {
      return(list(rmst0 = NA, rmst1 = NA, rmst_diff = NA))
    }
    
    # Helper: extract time and survival for one stratum
    extract_group_data <- function(group_label) {
      idx <- which(surv_summary$strata == group_label)
      time <- surv_summary$time[idx]
      surv <- surv_summary$surv[idx]
      return(list(time = time, surv = surv))
    }
    
    # Compute RMST manually
    compute_rmst <- function(time, surv, tau) {
      idx <- which(time <= tau)
      if (length(idx) == 0) {
        return(0)  # No events before tau
      }
      time_trunc <- c(0, time[idx], tau)
      surv_trunc <- c(1, surv[idx], surv[max(idx)])
      sum(diff(time_trunc) * head(surv_trunc, -1))
    }
    
    # Extract and compute RMST for both groups
    group0 <- extract_group_data("treatment=0")
    group1 <- extract_group_data("treatment=1")
    
    rmst0 <- compute_rmst(group0$time, group0$surv, tau)
    rmst1 <- compute_rmst(group1$time, group1$surv, tau)
    
    return(list(rmst0 = rmst0, rmst1 = rmst1, rmst_diff = rmst1 - rmst0))
  }, error = function(e) {
    return(list(rmst0 = NA, rmst1 = NA, rmst_diff = NA))
  })
}

#' Extract survival curve at specific time points
#' @param survfit_obj survfit object
#' @param time_points Vector of time points to extract
#' @return Data frame with time, survival for treatment=0, survival for treatment=1
extract_survival_curve <- function(survfit_obj, time_points) {
  tryCatch({
    surv_summary <- summary(survfit_obj, times = time_points, extend = TRUE)
    
    # Identify strata
    strata_names <- unique(surv_summary$strata)
    if (length(strata_names) != 2) {
      return(data.frame(time = time_points, surv0 = NA, surv1 = NA))
    }
    
    # Extract survival for each treatment group
    extract_group_surv <- function(group_label) {
      idx <- which(surv_summary$strata == group_label)
      time <- surv_summary$time[idx]
      surv <- surv_summary$surv[idx]
      # Interpolate/extend to match time_points
      surv_interp <- approx(time, surv, xout = time_points, method = "constant", 
                           yleft = 1, yright = tail(surv, 1), rule = 2)$y
      return(surv_interp)
    }
    
    surv0 <- extract_group_surv("treatment=0")
    surv1 <- extract_group_surv("treatment=1")
    
    return(data.frame(time = time_points, surv0 = surv0, surv1 = surv1))
  }, error = function(e) {
    return(data.frame(time = time_points, surv0 = NA, surv1 = NA))
  })
}

#' Run single bootstrap iteration
#' @param b Bootstrap iteration number
#' @param studyPop Original study population to sample from
#' @param cohortMethodData Original CohortMethodData object
#' @param cohortMethodData_allcovar Original CohortMethodData with all covariates
#' @param m Bootstrap sample size
#' @return List containing bootstrap results (ps_boot, Cox_censoring_boot, etc.)
runBootstrapIteration <- function(b, studyPop, cohortMethodData, 
                                  cohortMethodData_allcovar, m) {
  
  set.seed(MASTER_SEED + b)
  
  # Sample m out of n WITHOUT replacement
  boot_idx <- sample(seq_len(nrow(studyPop)), size = m, replace = FALSE)
  boot_studyPop <- studyPop[boot_idx, ]
  boot_oldRowIds <- studyPop$rowId[boot_idx]
  
  cat("Iteration", b, "hash:", digest::digest(boot_oldRowIds), "\n")
  print(table(boot_studyPop$treatment))
  
  # Map each bootstrap draw to a NEW unique rowId
  # Since we're sampling without replacement, each oldRowId appears exactly once
  boot_map <- tibble(
    oldRowId = as.integer(boot_oldRowIds),
    newRowId = seq_len(length(boot_oldRowIds))  # 1..m unique "subjects"
  )
  
  # For subsetting Andromeda tables efficiently:
  # Since no replacement, boot_old_unique is the same as boot_map (no duplicates)
  boot_old_unique <- boot_map %>% distinct(oldRowId)
  
  
  ### BOOTSTRAP COHORTMETHODDATA
  cohorts_sub <- cohortMethodData$cohorts %>% collect() %>%
    semi_join(boot_old_unique, by = c("rowId" = "oldRowId"))
  
  covariates_sub <- cohortMethodData$covariates %>% collect() %>%
    semi_join(boot_old_unique, by = c("rowId" = "oldRowId"))
  
  outcomes_sub <- cohortMethodData$outcomes %>% collect() %>%
    semi_join(boot_old_unique, by = c("rowId" = "oldRowId"))
  
  # Since sampling without replacement, each oldRowId appears exactly once
  # So we can use regular left_join (no need for many-to-many)
  cohorts_boot <- boot_map %>%
    left_join(cohorts_sub, by = c("oldRowId" = "rowId")) %>%
    mutate(rowId = newRowId) %>%
    select(-oldRowId, -newRowId)
  
  covariates_boot <- boot_map %>%
    left_join(covariates_sub, by = c("oldRowId" = "rowId")) %>%
    transmute(
      rowId = newRowId,
      covariateId = covariateId,
      covariateValue = covariateValue
    )
  
  # outcomes are optional for PS creation, but keeping them consistent is good
  outcomes_boot <- boot_map %>%
    left_join(outcomes_sub, by = c("oldRowId" = "rowId")) %>%
    mutate(rowId = newRowId) %>%
    select(-oldRowId, -newRowId)
  
  population_boot <- boot_studyPop %>%
    as_tibble() %>%
    mutate(rowId = seq_len(n()))  # newRowId 1..m
  
  # sanity checks
  stopifnot(nrow(population_boot) == m)
  stopifnot(all(population_boot$rowId == seq_len(m)))
  
  
  cohortMethodData_boot <- cohortMethodData
  cohortMethodData_boot$cohorts <- cohorts_boot
  cohortMethodData_boot$covariates <- covariates_boot
  cohortMethodData_boot$outcomes <- outcomes_boot
  
  # Ensure CohortMethodData class is preserved
  class(cohortMethodData_boot) <- class(cohortMethodData)
  
  ps_boot <- createPs(
    cohortMethodData = cohortMethodData_boot,
    population = population_boot,
    prior = createPrior("laplace", exclude = c(0), useCrossValidation = TRUE),
    control = createControl(threads = 8, fold = 5)
  )
  
  #### STOP HERE ### next we fit bootstrap censoring models
  # Note: cohortMethodData_allcovar is passed as parameter, no need to reload
  
  #----------------------------
  # 1) Subset required tables from cohortMethodData_allcovar
  #    (only the sampled original subjects)
  #----------------------------
  cohorts_sub <- cohortMethodData_allcovar$cohorts %>%
    collect() %>%
    semi_join(boot_old_unique, by = c("rowId" = "oldRowId")) %>%
    collect() %>%
    as_tibble()
  
  # outcomes_sub: original outcome events (NOT censoring events)
  # Expecting columns: rowId, daysToEvent, outcomeId, etc.
  outcomes_sub <- cohortMethodData_allcovar$outcomes %>%
    collect() %>%
    semi_join(boot_old_unique, by = c("rowId" = "oldRowId")) %>%
    collect() %>%
    as_tibble()
  
  # covariates_sub: long format sparse covariates: rowId, covariateId, covariateValue
  covariates_sub <- cohortMethodData_allcovar$covariates %>%
    collect() %>%
    semi_join(boot_old_unique, by = c("rowId" = "oldRowId")) %>%
    collect() %>%
    as_tibble()
  
  #----------------------------
  # 2) Build censoring outcome (y/time) on the ORIGINAL rowIds
  #    y = 1 means "censored" (not in outcomes table)
  #    time = daysToObsEnd if censored else daysToEvent
  #----------------------------
  # One row per original subject with their event time (if any)
  event_time <- outcomes_sub %>%
    group_by(rowId) %>%
    summarise(daysToEvent = min(daysToEvent), .groups = "drop")
  
  outcomes_for_cyclops_old <- cohorts_sub %>%
    select(rowId, daysToObsEnd) %>%
    left_join(event_time, by = "rowId") %>%
    mutate(
      y = ifelse(is.na(daysToEvent), 1L, 0L),
      time = ifelse(y == 1L, daysToObsEnd, daysToEvent)
    ) %>%
    select(rowId, y, time)
  
  #----------------------------
  # 3) EXPAND to bootstrap replicates by joining through boot_map
  #    and replacing rowId with newRowId
  #----------------------------
  # Expand outcomes_for_cyclops (becomes length m, unique rowId)
  outcomes_for_cyclops_boot <- boot_map %>%
    left_join(outcomes_for_cyclops_old, by = c("oldRowId" = "rowId")) %>%
    transmute(
      rowId = newRowId,
      y = y,
      time = time
    )
  
  # Expand covariates: since no replacement, each oldRowId appears once
  covariates_boot <- boot_map %>%
    left_join(covariates_sub, by = c("oldRowId" = "rowId")) %>%
    transmute(
      rowId = newRowId,
      covariateId = covariateId,
      covariateValue = covariateValue
    )
  
  # (Optional) If you need cohorts_boot for anything else:
  cohorts_boot <- boot_map %>%
    left_join(cohorts_sub, by = c("oldRowId" = "rowId")) %>%
    mutate(rowId = newRowId) %>%
    select(-oldRowId, -newRowId)
  
  #----------------------------
  # 4) Sanity checks
  #----------------------------
  stopifnot(nrow(outcomes_for_cyclops_boot) == m)
  stopifnot(length(unique(outcomes_for_cyclops_boot$rowId)) == m)
  stopifnot(all(!is.na(outcomes_for_cyclops_boot$time)))
  stopifnot(all(outcomes_for_cyclops_boot$time >= 0))
  
  #----------------------------
  # 5) Convert to Cyclops + Fit regularized Cox censoring model
  #----------------------------
  censored_df <- convertToCyclopsData(
    outcomes = outcomes_for_cyclops_boot,
    covariates = covariates_boot,
    modelType = "cox",
    addIntercept = TRUE
  )
  
  lassoPrior <- Cyclops::createPrior(
    priorType = "laplace",
    useCrossValidation = TRUE
  )
  
  Cox_censoring_boot <- fitCyclopsModel(
    censored_df,
    prior = lassoPrior,
    control = createControl(threads = 8, fold = 5)
  )
  
  cat("Cox censoring model is fit!")
  
  #----------------------------
  # 6) Compute IPCW weights (from find_kZ_script.R)
  #----------------------------
  
  # Restrict CohortMethodData_allcovar_boot to PS analysis cohort
  cohort_ids <- ps_boot$rowId
  
  # Filter cohorts, covariates, and outcomes to PS analysis cohort
  cohorts_filtered <- cohorts_boot %>%
    filter(rowId %in% cohort_ids)
  covariates_filtered <- covariates_boot %>%
    filter(rowId %in% cohort_ids)
  outcomes_filtered <- outcomes_boot %>%
    filter(rowId %in% cohort_ids)
  
  # Build survival df from ps_boot
  outcomes_for_cyclops <- ps_boot %>%
    transmute(
      rowId = rowId,
      y = if_else(outcomeCount == 0, 1, 0),  # 1 = censored, 0 = event
      time = survivalTime,
      treatment = treatment, 
      iptw = iptw
    )
  
  # Extract non-zero covariates from censoring model
  coefs <- coef(Cox_censoring_boot)
  non_zero_coefs <- coefs[coefs != 0]
  
  # Filter to non-zero covariates only
  filtered_covariates <- covariates_filtered %>%
    filter(covariateId %in% names(non_zero_coefs))
  
  # Pivot longer -> wider (rowId = patient, columns = covariates)
  X_wide <- filtered_covariates %>%
    select(rowId, covariateId, covariateValue) %>%
    mutate(covariateId = as.character(covariateId)) %>%
    pivot_wider(
      names_from = covariateId,
      values_from = covariateValue,
      values_fill = 0
    )
  
  # Set rowId as rownames
  X_matrix <- X_wide %>%
    column_to_rownames(var = "rowId") %>%
    as.matrix()
  
  # Create coeff_vector matching X_matrix
  coeff_vector <- non_zero_coefs[colnames(X_matrix)]
  stopifnot(length(coeff_vector) == ncol(X_matrix))
  stopifnot(all(names(coeff_vector) == colnames(X_matrix)))
  
  # Align outcome rows to the covariate matrix row order
  outcomes_df <- data.table(outcomes_for_cyclops)
  outcomes_df <- outcomes_df %>%
    arrange(match(as.character(rowId), rownames(X_matrix)))
  
  # Baseline survival for censoring [numerator for stabilized IPCW]
  cox_baseline <- coxph(Surv(time, y) ~ 1, data = outcomes_df)
  baseline_surv <- survfit(cox_baseline)
  baseline_surv_fun <- stepfun(baseline_surv$time, c(1, baseline_surv$surv))
  
  # Breslow estimator for covariate-adjusted censoring model [denom]
  breslow_est <- function(time, status, X, B) {
    data <- data.frame(time, status, X)
    data <- data[order(data$time), ]
    t <- unique(data$time)
    k <- length(t)
    h <- rep(0, k)
    LP_indiv <- X %*% B
    
    for(i in 1:k) {
      lp <- (LP_indiv)[data$time >= t[i]]
      risk <- exp(lp)
      h[i] <- sum(data$status[data$time == t[i]]) / sum(risk)
    }
    res <- cumsum(h)
    return(res)
  }
  
  H0 <- breslow_est(time = outcomes_df$time, status = outcomes_df$y, 
                    X = X_matrix, B = coeff_vector)
  haz_step_fun <- stepfun(sort(unique(outcomes_df$time)), c(0, H0))
  
  # Transform to person-time (long) format with intervals
  transform.data <- function(data, cut.times) {
    data$Tstart <- 0
    data$ami <- 1 - data$y
    
    data.long <- survSplit(data = data,
                           cut = cut.times,
                           end = "time",
                           start = "Tstart",
                           event = "y")
    data.long <- data.table(data.long)
    data.long <- data.long[order(data.long$rowId, data.long$time),]
    
    data.long.cens <- survSplit(data,
                                cut = cut.times,
                                end = "time",
                                start = "Tstart",
                                event = "ami")
    data.long.cens <- data.long.cens[order(data.long.cens$rowId, data.long.cens$time),]
    data.long$ami <- data.long.cens$ami
    data.long$rowId <- as.numeric(data.long$rowId)
    return(data.long)
  }
  
  # Determine cut times
  dist <- summary(outcomes_df$time)
  cut.times = seq(from = 60, to = floor(max(outcomes_df$time) / 60) * 60, by = 60) 
  
  # Transform to long data
  outcomes_df.long <- transform.data(outcomes_df, cut.times)
  
  # Compute IPCW weights
  eta <- X_matrix %*% coeff_vector
  outcomes_df.long$H0_Tstart <- haz_step_fun(outcomes_df.long$Tstart)
  rowid_to_eta <- data.frame(rowId = as.numeric(rownames(X_matrix)), 
                             eta = as.numeric(eta))
  outcomes_df.long <- outcomes_df.long %>%
    left_join(rowid_to_eta, by = "rowId")
  outcomes_df.long <- outcomes_df.long %>%
    mutate(KZ = exp(-H0_Tstart * exp(eta)))
  outcomes_df.long$K0_ti <- baseline_surv_fun(outcomes_df.long$Tstart)
  outcomes_df.long$Unstab_ipcw <- 1/outcomes_df.long$KZ
  outcomes_df.long$Stab_ipcw <- outcomes_df.long$K0_ti/outcomes_df.long$KZ
  
  #----------------------------
  # 7) Combine weights and fit final models (from weights_script.R)
  #----------------------------
  
  # Truncate stabilized IPCW weights
  lo <- quantile(outcomes_df.long$Stab_ipcw, 0.01, na.rm = TRUE)
  hi <- quantile(outcomes_df.long$Stab_ipcw, 0.99, na.rm = TRUE)
  outcomes_df.long$Stab_ipcw_trunc <- pmin(pmax(outcomes_df.long$Stab_ipcw, lo), hi)
  
  # Truncate unstabilized IPCW weights
  lo <- quantile(outcomes_df.long$Unstab_ipcw, 0.01, na.rm = TRUE)
  hi <- quantile(outcomes_df.long$Unstab_ipcw, 0.99, na.rm = TRUE)
  outcomes_df.long$Unstab_ipcw_trunc <- pmin(pmax(outcomes_df.long$Unstab_ipcw, lo), hi)
  
  # Truncate IPTW weights
  outcomes_df.long <- outcomes_df.long %>%
    group_by(Tstart) %>%
    mutate(
      lower_iptw = quantile(iptw, 0.01, na.rm = TRUE),
      upper_iptw = quantile(iptw, 0.99, na.rm = TRUE),
      iptw_trunc = pmin(pmax(iptw, lower_iptw), upper_iptw)
    ) %>%
    ungroup()
  
  # Combined weights
  outcomes_df.long$comb <- outcomes_df.long$Stab_ipcw_trunc * outcomes_df.long$iptw_trunc
  
  #----------------------------
  # 8) Fit final models and extract coefficients (log(HR))
  #----------------------------
  
  results <- list()
  results$bootstrap_iteration <- b
  
  # Unadjusted
  tryCatch({
    fit_unadjusted <- coxph(Surv(Tstart, time, ami) ~ treatment, 
                            data = outcomes_df.long, id = rowId)
    results$unadjusted_coef <- coef(fit_unadjusted)
  }, error = function(e) {
    results$unadjusted_coef <- NA
  })
  
  # IPCW stabilized truncated
  tryCatch({
    fit_ipcw_Stab_trunc <- coxph(Surv(Tstart, time, ami) ~ treatment, 
                                 data = outcomes_df.long, 
                                 id = rowId, weights = Stab_ipcw_trunc)
    results$ipcw_stab_trunc_coef <- coef(fit_ipcw_Stab_trunc)
  }, error = function(e) {
    results$ipcw_stab_trunc_coef <- NA
  })
  
  # IPTW truncated
  tryCatch({
    fit_iptw_trunc <- coxph(Surv(Tstart, time, ami) ~ treatment, 
                            data = outcomes_df.long, 
                            id = rowId, weights = iptw_trunc)
    results$iptw_trunc_coef <- coef(fit_iptw_trunc)
  }, error = function(e) {
    results$iptw_trunc_coef <- NA
  })
  
  # Combined weights
  tryCatch({
    fit_comb <- coxph(Surv(Tstart, time, ami) ~ treatment, 
                      data = outcomes_df.long, 
                      id = rowId, weights = comb)
    results$combined_coef <- coef(fit_comb)
  }, error = function(e) {
    results$combined_coef <- NA
  })
  
  #----------------------------
  # 9) Compute survival curves and RMST for each model
  #----------------------------
  
  # Helper function to compute survfit and extract results
  compute_survival_rmst <- function(data, weight_var = NULL, model_name) {
    tryCatch({
      # Fit survfit
      if (is.null(weight_var)) {
        survfit_obj <- survfit(Surv(Tstart, time, ami) ~ treatment, 
                               data = data, id = rowId)
      } else {
        survfit_obj <- survfit(Surv(Tstart, time, ami) ~ treatment, 
                              data = data, id = rowId, weights = data[[weight_var]])
      }
      
      # Extract survival curve at time grid
      surv_curve <- extract_survival_curve(survfit_obj, time_grid)
      
      # Compute RMST
      rmst_results <- compute_rmst_from_survfit(survfit_obj, tau)
      
      return(list(
        surv_curve = surv_curve,
        rmst0 = rmst_results$rmst0,
        rmst1 = rmst_results$rmst1,
        rmst_diff = rmst_results$rmst_diff
      ))
    }, error = function(e) {
      return(list(
        surv_curve = data.frame(time = time_grid, surv0 = NA, surv1 = NA),
        rmst0 = NA,
        rmst1 = NA,
        rmst_diff = NA
      ))
    })
  }
  
  # Unadjusted
  unadj_results <- compute_survival_rmst(outcomes_df.long, weight_var = NULL, "unadjusted")
  results$unadjusted_rmst0 <- unadj_results$rmst0
  results$unadjusted_rmst1 <- unadj_results$rmst1
  results$unadjusted_rmst_diff <- unadj_results$rmst_diff
  results$unadjusted_surv_curve <- unadj_results$surv_curve
  
  # IPTW truncated
  iptw_results <- compute_survival_rmst(outcomes_df.long, weight_var = "iptw_trunc", "iptw")
  results$iptw_rmst0 <- iptw_results$rmst0
  results$iptw_rmst1 <- iptw_results$rmst1
  results$iptw_rmst_diff <- iptw_results$rmst_diff
  results$iptw_surv_curve <- iptw_results$surv_curve
  
  # IPCW stabilized truncated
  ipcw_results <- compute_survival_rmst(outcomes_df.long, weight_var = "Stab_ipcw_trunc", "ipcw")
  results$ipcw_rmst0 <- ipcw_results$rmst0
  results$ipcw_rmst1 <- ipcw_results$rmst1
  results$ipcw_rmst_diff <- ipcw_results$rmst_diff
  results$ipcw_surv_curve <- ipcw_results$surv_curve
  
  # Combined weights
  comb_results <- compute_survival_rmst(outcomes_df.long, weight_var = "comb", "combined")
  results$combined_rmst0 <- comb_results$rmst0
  results$combined_rmst1 <- comb_results$rmst1
  results$combined_rmst_diff <- comb_results$rmst_diff
  results$combined_surv_curve <- comb_results$surv_curve
  
  return(results)
}


################################################################################
# Main Bootstrap Procedure
################################################################################

cat("=== Starting m-of-n Bootstrap (WITHOUT replacement) ===\n")
cat(sprintf("Bootstrap iterations: %d\n", B))
cat(sprintf("Bootstrap sample size: m = %d\n", m))
cat("Note: Sampling WITHOUT replacement - each subject appears at most once per iteration\n")
cat(sprintf("RMST restriction time: tau = %d days\n", tau))
cat(sprintf("Survival curve time grid: every 60 days up to %d days\n", tau))

# Load original data ONCE (outside loop for efficiency)
cat("Loading original data...\n")
cohortMethodData <- loadCohortMethodData(cohortMethodData_file)
cohortMethodData_allcovar <- loadCohortMethodData(cohortMethodData_allcovar_file)

# Create study pop ONCE (this is what we will sample from)
cat("Creating study population...\n")
studyPop <- CohortMethod::createStudyPopulation(
  cohortMethodData = cohortMethodData, 
  outcomeId = 1788866, # AMI
  firstExposureOnly = TRUE,
  washoutPeriod = 365,
  removeDuplicateSubjects = "keep first",
  censorAtNewRiskWindow = FALSE,
  removeSubjectsWithPriorOutcome = TRUE,
  priorOutcomeLookback = 99999,
  riskWindowStart = 1,
  startAnchor = "cohort start",
  riskWindowEnd = 9999,
  endAnchor = "cohort end",
  minDaysAtRisk = 1
)

n <- nrow(studyPop)
cat(sprintf("Original sample size (from studyPop): n = %d\n", n))

# Check that m <= n for sampling without replacement
if (m > n) {
  stop(sprintf("Cannot sample m = %d subjects without replacement from n = %d subjects. Please set m <= n.", m, n))
}

# Check for existing results to resume from
results_file <- file.path(output_dir, "bootstrap_results_iterative.csv")
rmst_file <- file.path(output_dir, "bootstrap_rmst_iterative.csv")
survcurves_file <- file.path(output_dir, "bootstrap_survcurves_iterative.csv")
checkpoint_file <- file.path(output_dir, "bootstrap_results_checkpoint.rds")
completed_iters <- integer(0)

if (file.exists(results_file)) {
  cat("\nFound existing results file. Loading completed iterations...\n")
  existing <- read.csv(results_file)
  completed_iters <- existing$iteration
  cat(sprintf("Found %d completed iterations. Will resume from iteration %d.\n", 
              length(completed_iters), max(completed_iters) + 1))
} else if (file.exists(checkpoint_file)) {
  cat("\nFound checkpoint file. Loading...\n")
  bootstrap_results_checkpoint <- readRDS(checkpoint_file)
  completed_iters <- sapply(bootstrap_results_checkpoint, function(x) x$bootstrap_iteration)
  cat(sprintf("Found %d completed iterations in checkpoint. Will resume from iteration %d.\n", 
              length(completed_iters), max(completed_iters) + 1))
}

# Initialize results storage
bootstrap_results <- list()

# Progress bar
pb <- progress_bar$new(
  format = "  Progress [:bar] :percent in :elapsed, ETA: :eta",
  total = B,
  clear = FALSE,
  width = 60
)

# Run bootstrap iterations
cat("\nRunning bootstrap iterations...\n")
for (b in 1:B) {
  # Skip if already completed
  if (b %in% completed_iters) {
    pb$tick()
    next
  }
  
  pb$tick()
  
  cohortMethodData <- loadCohortMethodData(cohortMethodData_file)
  cohortMethodData_allcovar <- loadCohortMethodData(cohortMethodData_allcovar_file)
  
  # Create study pop ONCE (this is what we will sample from)
  cat("Creating study population...\n")
  studyPop <- CohortMethod::createStudyPopulation(
    cohortMethodData = cohortMethodData, 
    outcomeId = 1788866, # AMI
    firstExposureOnly = TRUE,
    washoutPeriod = 365,
    removeDuplicateSubjects = "keep first",
    censorAtNewRiskWindow = FALSE,
    removeSubjectsWithPriorOutcome = TRUE,
    priorOutcomeLookback = 99999,
    riskWindowStart = 1,
    startAnchor = "cohort start",
    riskWindowEnd = 9999,
    endAnchor = "cohort end",
    minDaysAtRisk = 1
  )
  
  
  tryCatch({
    results <- runBootstrapIteration(b, studyPop, cohortMethodData, 
                                     cohortMethodData_allcovar, m)
    bootstrap_results[[b]] <- results
    
    # Save this iteration's result immediately (coefficients only)
    res_row <- data.frame(
      iteration = results$bootstrap_iteration,
      unadjusted_coef = ifelse(is.null(results$unadjusted_coef) || is.na(results$unadjusted_coef), 
                               NA, results$unadjusted_coef),
      ipcw_stab_trunc_coef = ifelse(is.null(results$ipcw_stab_trunc_coef) || is.na(results$ipcw_stab_trunc_coef), 
                                    NA, results$ipcw_stab_trunc_coef),
      iptw_trunc_coef = ifelse(is.null(results$iptw_trunc_coef) || is.na(results$iptw_trunc_coef), 
                               NA, results$iptw_trunc_coef),
      combined_coef = ifelse(is.null(results$combined_coef) || is.na(results$combined_coef), 
                             NA, results$combined_coef)
    )
    
    # Append to CSV file
    if (!file.exists(results_file)) {
      write.csv(res_row, results_file, row.names = FALSE)
    } else {
      write.table(res_row, results_file, row.names = FALSE,
                  col.names = FALSE, sep = ",", append = TRUE)
    }
    
    # Save RMST results progressively to CSV
    rmst_row <- data.frame(
      iteration = b,
      unadjusted_rmst0 = ifelse(is.null(results$unadjusted_rmst0) || is.na(results$unadjusted_rmst0), 
                                NA, results$unadjusted_rmst0),
      unadjusted_rmst1 = ifelse(is.null(results$unadjusted_rmst1) || is.na(results$unadjusted_rmst1), 
                                NA, results$unadjusted_rmst1),
      unadjusted_rmst_diff = ifelse(is.null(results$unadjusted_rmst_diff) || is.na(results$unadjusted_rmst_diff), 
                                    NA, results$unadjusted_rmst_diff),
      iptw_rmst0 = ifelse(is.null(results$iptw_rmst0) || is.na(results$iptw_rmst0), 
                         NA, results$iptw_rmst0),
      iptw_rmst1 = ifelse(is.null(results$iptw_rmst1) || is.na(results$iptw_rmst1), 
                         NA, results$iptw_rmst1),
      iptw_rmst_diff = ifelse(is.null(results$iptw_rmst_diff) || is.na(results$iptw_rmst_diff), 
                             NA, results$iptw_rmst_diff),
      ipcw_rmst0 = ifelse(is.null(results$ipcw_rmst0) || is.na(results$ipcw_rmst0), 
                         NA, results$ipcw_rmst0),
      ipcw_rmst1 = ifelse(is.null(results$ipcw_rmst1) || is.na(results$ipcw_rmst1), 
                         NA, results$ipcw_rmst1),
      ipcw_rmst_diff = ifelse(is.null(results$ipcw_rmst_diff) || is.na(results$ipcw_rmst_diff), 
                             NA, results$ipcw_rmst_diff),
      combined_rmst0 = ifelse(is.null(results$combined_rmst0) || is.na(results$combined_rmst0), 
                             NA, results$combined_rmst0),
      combined_rmst1 = ifelse(is.null(results$combined_rmst1) || is.na(results$combined_rmst1), 
                             NA, results$combined_rmst1),
      combined_rmst_diff = ifelse(is.null(results$combined_rmst_diff) || is.na(results$combined_rmst_diff), 
                                 NA, results$combined_rmst_diff)
    )
    
    # Append RMST to CSV file
    if (!file.exists(rmst_file)) {
      write.csv(rmst_row, rmst_file, row.names = FALSE)
    } else {
      write.table(rmst_row, rmst_file, row.names = FALSE,
                  col.names = FALSE, sep = ",", append = TRUE)
    }
    
    # Save survival curves progressively to CSV
    # Each model's survival curve has multiple rows (one per time point)
    models_surv <- c("unadjusted", "iptw", "ipcw", "combined")
    
    for (model in models_surv) {
      # Map model names to result keys
      model_key_map <- list(
        "unadjusted" = "unadjusted_surv_curve",
        "iptw" = "iptw_surv_curve",
        "ipcw" = "ipcw_surv_curve",
        "combined" = "combined_surv_curve"
      )
      model_key <- model_key_map[[model]]
      
      surv_curve <- results[[model_key]]
      
      if (!is.null(surv_curve) && is.data.frame(surv_curve) && nrow(surv_curve) > 0) {
        surv_row <- surv_curve %>%
          mutate(iteration = b, model = model) %>%
          select(iteration, model, time, surv0, surv1)
        
        # Append survival curves to CSV file
        if (!file.exists(survcurves_file)) {
          write.csv(surv_row, survcurves_file, row.names = FALSE)
        } else {
          write.table(surv_row, survcurves_file, row.names = FALSE,
                      col.names = FALSE, sep = ",", append = TRUE)
        }
      }
    }
    
    # Save checkpoint RDS every 10 iterations
    if (b %% 10 == 0 || b == B) {
      # Load existing results if any
      if (file.exists(results_file)) {
        all_results <- read.csv(results_file)
        # Convert to list format for checkpoint
        checkpoint_list <- lapply(1:nrow(all_results), function(i) {
          list(
            bootstrap_iteration = all_results$iteration[i],
            unadjusted_coef = all_results$unadjusted_coef[i],
            ipcw_stab_trunc_coef = all_results$ipcw_stab_trunc_coef[i],
            iptw_trunc_coef = all_results$iptw_trunc_coef[i],
            combined_coef = all_results$combined_coef[i]
          )
        })
        saveRDS(checkpoint_list, checkpoint_file)
      }
    }
    
  }, error = function(e) {
    cat(sprintf("\nError in bootstrap iteration %d: %s\n", b, e$message))
    # Save error as NA values for coefficients
    res_row <- data.frame(
      iteration = b,
      unadjusted_coef = NA,
      ipcw_stab_trunc_coef = NA,
      iptw_trunc_coef = NA,
      combined_coef = NA
    )
    
    if (!file.exists(results_file)) {
      write.csv(res_row, results_file, row.names = FALSE)
    } else {
      write.table(res_row, results_file, row.names = FALSE,
                  col.names = FALSE, sep = ",", append = TRUE)
    }
    
    # Save error as NA values for RMST
    rmst_row <- data.frame(
      iteration = b,
      unadjusted_rmst0 = NA, unadjusted_rmst1 = NA, unadjusted_rmst_diff = NA,
      iptw_rmst0 = NA, iptw_rmst1 = NA, iptw_rmst_diff = NA,
      ipcw_rmst0 = NA, ipcw_rmst1 = NA, ipcw_rmst_diff = NA,
      combined_rmst0 = NA, combined_rmst1 = NA, combined_rmst_diff = NA
    )
    
    if (!file.exists(rmst_file)) {
      write.csv(rmst_row, rmst_file, row.names = FALSE)
    } else {
      write.table(rmst_row, rmst_file, row.names = FALSE,
                  col.names = FALSE, sep = ",", append = TRUE)
    }
    
    # Save error as NA values for survival curves (one row per model per time point)
    models_surv <- c("unadjusted", "iptw", "ipcw", "combined")
    for (model in models_surv) {
      surv_row <- data.frame(
        iteration = b,
        model = model,
        time = time_grid,
        surv0 = NA,
        surv1 = NA
      )
      
      if (!file.exists(survcurves_file)) {
        write.csv(surv_row, survcurves_file, row.names = FALSE)
      } else {
        write.table(surv_row, survcurves_file, row.names = FALSE,
                    col.names = FALSE, sep = ",", append = TRUE)
      }
    }
    
    bootstrap_results[[b]] <- list(
      bootstrap_iteration = b,
      error = e$message
    )
  })
}

cat("\n=== Bootstrap Complete ===\n")
cat(sprintf("Completed %d iterations\n", length(bootstrap_results)))

################################################################################
# Calculate Bootstrap Statistics for HR Coefficients
# 
# All bootstrap statistics are calculated on the log(HR) scale (coefficients).
# This is the correct scale for bootstrap inference:
# - Bootstrap SD of log(HR) is calculated directly
# - 95% CI is constructed using quantiles of log(HR)
# - HR values are transformed back (exp) only for display/interpretation
################################################################################

cat("\n=== Calculating Bootstrap Statistics for HR Coefficients ===\n")
cat("Working on log(HR) scale (coefficients) for all bootstrap statistics\n")

# Load results from saved file (progressive save) or use in-memory results
if (file.exists(results_file)) {
  cat("Loading results from progressive save file...\n")
  results_df <- read.csv(results_file)
} else {
  # Fallback to in-memory results if file doesn't exist
  cat("Using in-memory results...\n")
  results_df <- do.call(rbind, lapply(bootstrap_results, function(x) {
    if (is.null(x) || !is.null(x$error)) {
      return(data.frame(
        iteration = ifelse(is.null(x$bootstrap_iteration), NA, x$bootstrap_iteration),
        unadjusted_coef = NA,
        ipcw_stab_trunc_coef = NA,
        iptw_trunc_coef = NA,
        combined_coef = NA
      ))
    }
    data.frame(
      iteration = x$bootstrap_iteration,
      unadjusted_coef = ifelse(is.null(x$unadjusted_coef), NA, x$unadjusted_coef),
      ipcw_stab_trunc_coef = ifelse(is.null(x$ipcw_stab_trunc_coef), NA, x$ipcw_stab_trunc_coef),
      iptw_trunc_coef = ifelse(is.null(x$iptw_trunc_coef), NA, x$iptw_trunc_coef),
      combined_coef = ifelse(is.null(x$combined_coef), NA, x$combined_coef)
    )
  }))
}

# Save final bootstrap results (copy of progressive save for consistency)
write.csv(results_df, file.path(output_dir, "bootstrap_results.csv"), row.names = FALSE)

# Function to calculate bootstrap SE and CI on log(HR) scale (coefficients)
calcBootstrapStats <- function(coef_estimates, alpha = 0.05) {
  coef_estimates <- coef_estimates[!is.na(coef_estimates)]
  if (length(coef_estimates) == 0) {
    return(list(mean_logHR = NA, se_logHR = NA, median_logHR = NA, 
                ci_lower_logHR = NA, ci_upper_logHR = NA, n_valid = 0))
  }
  
  mean_logHR <- mean(coef_estimates)
  se_logHR <- sd(coef_estimates)  # Bootstrap SD of log(HR)
  median_logHR <- median(coef_estimates)
  
  # Percentile method CI on log(HR) scale
  ci_lower_logHR <- quantile(coef_estimates, alpha/2, na.rm = TRUE)
  ci_upper_logHR <- quantile(coef_estimates, 1 - alpha/2, na.rm = TRUE)
  
  return(list(
    mean_logHR = mean_logHR,
    se_logHR = se_logHR,
    median_logHR = median_logHR,
    ci_lower_logHR = as.numeric(ci_lower_logHR),
    ci_upper_logHR = as.numeric(ci_upper_logHR),
    n_valid = length(coef_estimates)
  ))
}

# Calculate statistics for each model
summary_stats <- list()

models <- c("unadjusted", "ipcw_stab_trunc", "iptw_trunc", "combined")
for (model in models) {
  coef_col <- paste0(model, "_coef")
  
  if (coef_col %in% names(results_df)) {
    # Calculate bootstrap statistics on log(HR) scale
    logHR_stats <- calcBootstrapStats(results_df[[coef_col]], alpha)
    
    # Transform back to HR scale for display (optional, but useful for interpretation)
    HR_from_mean_logHR <- exp(logHR_stats$mean_logHR)
    HR_from_median_logHR <- exp(logHR_stats$median_logHR)
    HR_CI_lower <- exp(logHR_stats$ci_lower_logHR)
    HR_CI_upper <- exp(logHR_stats$ci_upper_logHR)
    
    summary_stats[[model]] <- data.frame(
      Model = model,
      # Log(HR) scale statistics (primary)
      logHR_mean = logHR_stats$mean_logHR,
      logHR_se = logHR_stats$se_logHR,  # Bootstrap SD of log(HR)
      logHR_median = logHR_stats$median_logHR,
      logHR_CI_lower = logHR_stats$ci_lower_logHR,
      logHR_CI_upper = logHR_stats$ci_upper_logHR,
      # HR scale (transformed back for interpretation)
      HR_from_mean_logHR = HR_from_mean_logHR,
      HR_from_median_logHR = HR_from_median_logHR,
      HR_CI_lower = HR_CI_lower,
      HR_CI_upper = HR_CI_upper,
      n_valid = logHR_stats$n_valid
    )
  }
}

summary_df <- do.call(rbind, summary_stats)
write.csv(summary_df, file.path(output_dir, "bootstrap_summary.csv"), row.names = FALSE)

# Print summary
cat("\n=== Bootstrap Summary (HR Coefficients) ===\n")
print(summary_df)

################################################################################
# Calculate Bootstrap Statistics for RMST
################################################################################

cat("\n=== Calculating Bootstrap Statistics for RMST ===\n")

# Load RMST results from CSV
if (file.exists(rmst_file)) {
  rmst_df <- read.csv(rmst_file)
  
  # Calculate statistics for each model
  rmst_summary <- list()
  
  models_rmst <- c("unadjusted", "iptw", "ipcw", "combined")
  for (model in models_rmst) {
    rmst0_col <- paste0(model, "_rmst0")
    rmst1_col <- paste0(model, "_rmst1")
    rmst_diff_col <- paste0(model, "_rmst_diff")
    
    if (all(c(rmst0_col, rmst1_col, rmst_diff_col) %in% names(rmst_df))) {
      rmst0_vals <- rmst_df[[rmst0_col]][!is.na(rmst_df[[rmst0_col]])]
      rmst1_vals <- rmst_df[[rmst1_col]][!is.na(rmst_df[[rmst1_col]])]
      rmst_diff_vals <- rmst_df[[rmst_diff_col]][!is.na(rmst_df[[rmst_diff_col]])]
      
      if (length(rmst_diff_vals) > 0) {
        rmst_summary[[model]] <- data.frame(
          Model = model,
          RMST0_mean = mean(rmst0_vals),
          RMST0_se = sd(rmst0_vals),
          RMST0_CI_lower = quantile(rmst0_vals, alpha/2, na.rm = TRUE),
          RMST0_CI_upper = quantile(rmst0_vals, 1 - alpha/2, na.rm = TRUE),
          RMST1_mean = mean(rmst1_vals),
          RMST1_se = sd(rmst1_vals),
          RMST1_CI_lower = quantile(rmst1_vals, alpha/2, na.rm = TRUE),
          RMST1_CI_upper = quantile(rmst1_vals, 1 - alpha/2, na.rm = TRUE),
          RMST_diff_mean = mean(rmst_diff_vals),
          RMST_diff_se = sd(rmst_diff_vals),
          RMST_diff_CI_lower = quantile(rmst_diff_vals, alpha/2, na.rm = TRUE),
          RMST_diff_CI_upper = quantile(rmst_diff_vals, 1 - alpha/2, na.rm = TRUE),
          n_valid = length(rmst_diff_vals)
        )
      }
    }
  }
  
  if (length(rmst_summary) > 0) {
    rmst_summary_df <- do.call(rbind, rmst_summary)
    write.csv(rmst_summary_df, file.path(output_dir, "bootstrap_rmst_summary.csv"), row.names = FALSE)
    
    cat("\n=== Bootstrap RMST Summary ===\n")
    print(rmst_summary_df)
  }
} else {
  cat("No RMST results file found.\n")
}

################################################################################
# Calculate Bootstrap Statistics for Survival Curves
################################################################################

cat("\n=== Calculating Bootstrap Statistics for Survival Curves ===\n")

# Load survival curves from CSV
if (file.exists(survcurves_file)) {
  survcurves_df <- read.csv(survcurves_file)
  
  # Aggregate survival curves across bootstrap iterations
  # For each model and time point, compute mean and quantiles
  surv_summaries <- survcurves_df %>%
    group_by(model, time) %>%
    summarise(
      surv0_mean = mean(surv0, na.rm = TRUE),
      surv0_CI_lower = quantile(surv0, alpha/2, na.rm = TRUE),
      surv0_CI_upper = quantile(surv0, 1 - alpha/2, na.rm = TRUE),
      surv1_mean = mean(surv1, na.rm = TRUE),
      surv1_CI_lower = quantile(surv1, alpha/2, na.rm = TRUE),
      surv1_CI_upper = quantile(surv1, 1 - alpha/2, na.rm = TRUE),
      n_valid = sum(!is.na(surv0) & !is.na(surv1)),
      .groups = "drop"
    ) %>%
    arrange(model, time)
  
  # Save summary
  write.csv(surv_summaries, file.path(output_dir, "bootstrap_survival_curves_summary.csv"), row.names = FALSE)
  
  cat("\n=== Bootstrap Survival Curves Summary ===\n")
  cat(sprintf("Survival curves computed at %d time points (every 60 days up to %d days)\n", 
              length(time_grid), tau))
  cat(sprintf("Summary saved with point estimates and %d%% confidence intervals\n", (1-alpha)*100))
  cat(sprintf("Total rows in summary: %d (4 models × %d time points)\n", 
              nrow(surv_summaries), length(time_grid)))
} else {
  cat("No survival curves results file found.\n")
}

################################################################################
# Final Summary
################################################################################

cat(sprintf("\nResults saved to: %s/\n", output_dir))
cat("- bootstrap_results_iterative.csv: Progressive save (coefficients, updated after each iteration)\n")
cat("- bootstrap_results.csv: Final complete results (copy of iterative save)\n")
cat("- bootstrap_rmst_iterative.csv: RMST results (progressive save, updated after each iteration)\n")
cat("- bootstrap_rmst_summary.csv: RMST summary statistics (SEs and CIs)\n")
cat("- bootstrap_survcurves_iterative.csv: Survival curves (progressive save, updated after each iteration)\n")
cat("- bootstrap_survival_curves_summary.csv: Survival curves with bootstrap CIs\n")
cat("- bootstrap_results_checkpoint.rds: Checkpoint file (saved every 10 iterations)\n")
cat("- bootstrap_summary.csv: HR coefficient summary statistics (SEs and CIs)\n")
cat("\nNote: If the script crashes, it will automatically resume from existing results.\n")
cat("All results (coefficients, RMST, survival curves) are saved progressively after each iteration.\n")
