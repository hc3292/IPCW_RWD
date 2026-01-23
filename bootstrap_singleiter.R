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

################################################################################
# Configuration
################################################################################

# Bootstrap parameters
B <- 10                    # Number of bootstrap iterations
m <- 10000                    # Bootstrap sample size (NULL = use full sample size n)
alpha <- 0.05                # For confidence intervals (1-alpha)*100%

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
## writing 1 iteration of bootstrap
################################################################################

cohortMethodData <- loadCohortMethodData(cohortMethodData_file)

# create study pop -- this is what we will sample from
studyPop <- CohortMethod::createStudyPopulation(
  cohortMethodData = cohortMethodData, 
  outcomeId = 1788866, # AMI
  firstExposureOnly = TRUE,
  washoutPeriod = 0,
  removeDuplicateSubjects = "keep first",
  censorAtNewRiskWindow = TRUE,
  removeSubjectsWithPriorOutcome = TRUE,
  priorOutcomeLookback = 99999,
  minDaysAtRisk = 1,
  maxDaysAtRisk = 99999,
  riskWindowStart = 0,
  startAnchor = "cohort start",
  riskWindowEnd = 99999,
  endAnchor = "cohort end"
)

# Sample m out of n with replacement
boot_idx <- sample(seq_len(nrow(studyPop)), size = m, replace = TRUE)
boot_studyPop <- studyPop[boot_idx, ]
boot_oldRowIds <- studyPop$rowId[boot_idx]

# Map each bootstrap draw to a NEW unique rowId
boot_map <- tibble(
  oldRowId = as.integer(boot_oldRowIds),
  newRowId = seq_len(length(boot_oldRowIds))  # 1..m unique "subjects"
)

# For subsetting Andromeda tables efficiently:
boot_old_unique <- boot_map %>% distinct(oldRowId)


### BOOTSTRAP COHORTMETHODDATA
cohorts_sub <- cohortMethodData$cohorts %>% collect() %>%
  semi_join(boot_old_unique, by = c("rowId" = "oldRowId"))

covariates_sub <- cohortMethodData$covariates %>% collect() %>%
  semi_join(boot_old_unique, by = c("rowId" = "oldRowId"))

outcomes_sub <- cohortMethodData$outcomes %>% collect() %>%
  semi_join(boot_old_unique, by = c("rowId" = "oldRowId"))

cohorts_boot <- boot_map %>%
  left_join(cohorts_sub, by = c("oldRowId" = "rowId"), relationship = "many-to-many") %>%
  mutate(rowId = newRowId) %>%
  select(-oldRowId, -newRowId)

covariates_boot <- boot_map %>%
  left_join(covariates_sub, by = c("oldRowId" = "rowId"), relationship = "many-to-many") %>%
  transmute(
    rowId = newRowId,
    covariateId = covariateId,
    covariateValue = covariateValue
  )

# outcomes are optional for PS creation, but keeping them consistent is good
outcomes_boot <- boot_map %>%
  left_join(outcomes_sub, by = c("oldRowId" = "rowId"), relationship = "many-to-many") %>%
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

ps_boot <- createPs(
  cohortMethodData = cohortMethodData_boot,
  population = population_boot
)


#### STOP HERE ### next we fit bootstrap censoring models

cohortMethodData_allcovar <- loadCohortMethodData(cohortMethodData_allcovar_file)

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

# Expand covariates: replicate all covariate rows for each bootstrap draw
covariates_boot <- boot_map %>%
  left_join(covariates_sub, by = c("oldRowId" = "rowId"), relationship = "many-to-many") %>%
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
  prior = lassoPrior
)

# saveRDS(Cox_censoring, "Cox_censoring.rds")


################################################################################
# Helper Functions
################################################################################

#' Create bootstrap sample of CohortMethodData
#' @param cohortMethodData Original CohortMethodData object
#' @param boot_rowIds Vector of rowIds to include in bootstrap sample (with replacement)
#' @return Bootstrap CohortMethodData object
createBootstrapCohortMethodData <- function(cohortMethodData, boot_rowIds) {
  # Collect original data
  cohorts_orig <- collect(cohortMethodData$cohorts)
  covariates_orig <- collect(cohortMethodData$covariates)
  outcomes_orig <- collect(cohortMethodData$outcomes)
  
  # Create mapping: new_rowId (1, 2, ..., m) -> original_rowId (sampled with replacement)
  m <- length(boot_rowIds)
  rowId_map <- data.frame(
    new_rowId = 1:m,
    original_rowId = boot_rowIds
  )
  
  # Get cohorts for sampled rowIds and assign new rowIds
  cohorts_boot <- cohorts_orig %>%
    filter(rowId %in% boot_rowIds) %>%
    rename(original_rowId = rowId) %>%
    # Join with mapping, but need to handle duplicates
    inner_join(rowId_map, by = "original_rowId") %>%
    select(-original_rowId) %>%
    rename(rowId = new_rowId)
  
  # Get covariates for sampled rowIds
  covariates_boot <- covariates_orig %>%
    filter(rowId %in% boot_rowIds) %>%
    rename(original_rowId = rowId) %>%
    inner_join(rowId_map, by = "original_rowId") %>%
    select(-original_rowId) %>%
    rename(rowId = new_rowId)
  
  # Get outcomes for sampled rowIds
  outcomes_boot <- outcomes_orig %>%
    filter(rowId %in% boot_rowIds) %>%
    rename(original_rowId = rowId) %>%
    inner_join(rowId_map, by = "original_rowId") %>%
    select(-original_rowId) %>%
    rename(rowId = new_rowId)
  
  # Create new CohortMethodData object
  cohortMethodData_boot <- Andromeda::andromeda()
  cohortMethodData_boot$cohorts <- cohorts_boot
  cohortMethodData_boot$covariates <- covariates_boot
  cohortMethodData_boot$outcomes <- outcomes_boot
  class(cohortMethodData_boot) <- class(cohortMethodData)
  
  # Copy metaData if it exists
  if ("metaData" %in% names(cohortMethodData)) {
    cohortMethodData_boot$metaData <- cohortMethodData$metaData
  }
  
  return(cohortMethodData_boot)
}

#' Run single bootstrap iteration
#' @param b Bootstrap iteration number
#' @param original_rowIds Original rowIds from full dataset
#' @param m Bootstrap sample size
#' @return List of HR estimates
runBootstrapIteration <- function(b, original_rowIds, m) {
  cat(sprintf("\n=== Bootstrap Iteration %d/%d ===\n", b, B))
  
  # Sample m out of n with replacement
  boot_rowIds <- sample(original_rowIds, size = m, replace = TRUE)
  
  tryCatch({
    # Step 1: Load and bootstrap CohortMethodData for PS
    cohortMethodData <- CohortMethod::loadCohortMethodData(cohortMethodData_file)
    cohortMethodData_boot <- createBootstrapCohortMethodData(cohortMethodData, boot_rowIds)
    
    # Step 2: Create study population and fit PS model
    # Use same parameters as original studyPop
    studyPop_boot <- CohortMethod::createStudyPopulation(
      cohortMethodData = cohortMethodData_boot, 
      outcomeId = 1788866, # AMI
      firstExposureOnly = FALSE,
      restrictToCommonPeriod = FALSE, 
      washoutPeriod = 0, 
      removeDuplicateSubjects = "keep first", 
      removeSubjectsWithPriorOutcome = TRUE, 
      minDaysAtRisk = 1,
      riskWindowStart = 0,
      startAnchor = "cohort start",
      riskWindowEnd = 99999,
      endAnchor = "cohort end"
    )
    
    ps_boot <- createPs(cohortMethodData = cohortMethodData_boot, 
                        population = studyPop_boot)
    saveRDS(ps_boot, temp_ps_file)
    
    # Step 3: Load all-covariate data and fit censoring model
    cohortMethodData_allcovar <- CohortMethod::loadCohortMethodData(cohortMethodData_allcovar_file)
    cohortMethodData_allcovar_boot <- createBootstrapCohortMethodData(cohortMethodData_allcovar, boot_rowIds)
    
    # Build censored cohort and outcomes
    cohortMethodData_allcovar_boot$censored_cohort <- cohortMethodData_allcovar_boot$cohorts %>%
      anti_join(cohortMethodData_allcovar_boot$outcomes, by = "rowId")
    
    cohortMethodData_allcovar_boot$censored_outcomes <- cohortMethodData_allcovar_boot$censored_cohort %>%
      select(rowId, daysToEvent = daysToObsEnd) %>%
      mutate(outcomeId = 1) %>%
      select(rowId, outcomeId, daysToEvent)
    
    cohortMethodData_allcovar_boot$censored_outcomes <- cohortMethodData_allcovar_boot$censored_outcomes %>% 
      filter(outcomeId == 1)
    
    # Construct Cyclops dataset
    cohorts <- cohortMethodData_allcovar_boot$cohorts
    censored_outcomes <- cohortMethodData_allcovar_boot$censored_outcomes
    
    outcomes_for_cyclops <- cohorts %>%
      left_join(censored_outcomes, by = "rowId") %>%
      mutate(y = ifelse(!is.na(outcomeId), 1, 0)) %>%
      mutate(time = ifelse(!is.na(outcomeId), daysToEvent, daysToObsEnd)) %>%
      select(rowId, y, time)
    
    cohortMethodData_allcovar_boot$outcomes_for_cyclops <- outcomes_for_cyclops
    
    censored_df <- convertToCyclopsData(cohortMethodData_allcovar_boot$outcomes_for_cyclops, 
                                        cohortMethodData_allcovar_boot$covariates, 
                                        modelType = "cox",
                                        addIntercept = TRUE)
    
    # Fit censoring model
    lassoPrior <- Cyclops::createPrior(priorType = "laplace", useCrossValidation = TRUE)
    Cox_censoring_boot <- fitCyclopsModel(censored_df, prior = lassoPrior)
    saveRDS(Cox_censoring_boot, temp_censoring_file)
    
    # Step 4: Compute IPCW weights (simplified version of find_kZ_script.R)
    ps_boot <- readRDS(temp_ps_file)
    Cox_censoring_boot <- readRDS(temp_censoring_file)
    
    # Restrict CohortMethodData to PS analysis cohort
    # Note: ps_boot already has the correct rowIds from bootstrap sample
    cohort_ids <- ps_boot$rowId
    # Collect and filter (Andromeda objects need to be collected before filtering)
    cohorts_filtered <- collect(cohortMethodData_allcovar_boot$cohorts) %>%
      filter(rowId %in% cohort_ids)
    covariates_filtered <- collect(cohortMethodData_allcovar_boot$covariates) %>%
      filter(rowId %in% cohort_ids)
    outcomes_filtered <- collect(cohortMethodData_allcovar_boot$outcomes) %>%
      filter(rowId %in% cohort_ids)
    
    # Update Andromeda object (or use filtered data directly)
    cohortMethodData_allcovar_boot$cohorts <- cohorts_filtered
    cohortMethodData_allcovar_boot$covariates <- covariates_filtered
    cohortMethodData_allcovar_boot$outcomes <- outcomes_filtered
    
    # Build survival df
    outcomes_for_cyclops <- ps_boot %>%
      transmute(
        rowId = rowId,
        y = if_else(outcomeCount == 0, 1, 0),
        time = survivalTime,
        treatment = treatment, 
        iptw = iptw
      )
    
    # Extract non-zero covariates
    coefs <- coef(Cox_censoring_boot)
    non_zero_coefs <- coefs[coefs != 0]
    
    filtered_covariates <- collect(cohortMethodData_allcovar_boot$covariates) %>%
      filter(covariateId %in% names(non_zero_coefs))
    
    X_wide <- filtered_covariates %>%
      select(rowId, covariateId, covariateValue) %>%
      mutate(covariateId = as.character(covariateId)) %>%
      pivot_wider(names_from = covariateId, values_from = covariateValue, values_fill = 0)
    
    X_matrix <- X_wide %>%
      column_to_rownames(var = "rowId") %>%
      as.matrix()
    
    coeff_vector <- non_zero_coefs[colnames(X_matrix)]
    stopifnot(length(coeff_vector) == ncol(X_matrix))
    stopifnot(all(names(coeff_vector) == colnames(X_matrix)))
    
    # Align outcomes
    outcomes_df <- data.table(collect(outcomes_for_cyclops))
    outcomes_df <- outcomes_df %>%
      arrange(match(as.character(rowId), rownames(X_matrix)))
    
    # Baseline survival
    cox_baseline <- coxph(Surv(time, y) ~ 1, data = outcomes_df)
    baseline_surv <- survfit(cox_baseline)
    baseline_surv_fun <- stepfun(baseline_surv$time, c(1, baseline_surv$surv))
    
    # Breslow estimator (simplified - using existing function from find_kZ_script.R)
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
    
    # Transform to long format (simplified cut times)
    dist <- summary(outcomes_df$time)
    cut.times = seq(from = 60, to = floor(max(outcomes_df$time) / 60) * 60, by = 60) 
    
    transform.data <- function(data, cut.times) {
      data$Tstart <- 0
      data$ami <- 1 - data$y
      data.long <- survSplit(data = data, cut = cut.times, end = "time",
                             start = "Tstart", event = "y")
      data.long <- data.table(data.long)
      data.long <- data.long[order(data.long$rowId, data.long$time),]
      data.long.cens <- survSplit(data, cut = cut.times, end = "time",
                                  start = "Tstart", event = "ami")
      data.long.cens <- data.long.cens[order(data.long.cens$rowId, data.long.cens$time),]
      data.long$ami <- data.long.cens$ami
      data.long$rowId <- as.numeric(data.long$rowId)
      return(data.long)
    }
    
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
    
    # Step 5: Combine weights and fit models (from weights_script.R)
    # Truncate weights
    lo <- quantile(outcomes_df.long$Stab_ipcw, 0.01, na.rm = TRUE)
    hi <- quantile(outcomes_df.long$Stab_ipcw, 0.99, na.rm = TRUE)
    outcomes_df.long$Stab_ipcw_trunc <- pmin(pmax(outcomes_df.long$Stab_ipcw, lo), hi)
    
    lo <- quantile(outcomes_df.long$Unstab_ipcw, 0.01, na.rm = TRUE)
    hi <- quantile(outcomes_df.long$Unstab_ipcw, 0.99, na.rm = TRUE)
    outcomes_df.long$Unstab_ipcw_trunc <- pmin(pmax(outcomes_df.long$Unstab_ipcw, lo), hi)
    
    outcomes_df.long <- outcomes_df.long %>%
      group_by(Tstart) %>%
      mutate(
        lower_iptw = quantile(iptw, 0.01, na.rm = TRUE),
        upper_iptw = quantile(iptw, 0.99, na.rm = TRUE),
        iptw_trunc = pmin(pmax(iptw, lower_iptw), upper_iptw)
      ) %>%
      ungroup()
    
    outcomes_df.long$comb <- outcomes_df.long$Stab_ipcw_trunc * outcomes_df.long$iptw_trunc
    
    # Extract HRs
    results <- extractHRs(outcomes_df.long)
    results$bootstrap_iteration <- b
    
    # Clean up temporary files
    if (file.exists(temp_ps_file)) file.remove(temp_ps_file)
    if (file.exists(temp_censoring_file)) file.remove(temp_censoring_file)
    if (file.exists(temp_weights_file)) file.remove(temp_weights_file)
    
    return(results)
    
  }, error = function(e) {
    cat(sprintf("Error in bootstrap iteration %d: %s\n", b, e$message))
    return(list(
      unadjusted_HR = NA,
      unadjusted_coef = NA,
      ipcw_stab_trunc_HR = NA,
      ipcw_stab_trunc_coef = NA,
      iptw_trunc_HR = NA,
      iptw_trunc_coef = NA,
      combined_HR = NA,
      combined_coef = NA,
      bootstrap_iteration = b
    ))
  })
}


################################################################################
# Main Bootstrap Procedure
################################################################################

cat("=== Starting m-of-n Bootstrap ===\n")
cat(sprintf("Bootstrap iterations: %d\n", B))
cat(sprintf("Bootstrap sample size: m = %d\n", m))

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
  washoutPeriod = 0,
  removeDuplicateSubjects = "keep first",
  censorAtNewRiskWindow = TRUE,
  removeSubjectsWithPriorOutcome = TRUE,
  priorOutcomeLookback = 99999,
  minDaysAtRisk = 1,
  maxDaysAtRisk = 99999,
  riskWindowStart = 0,
  startAnchor = "cohort start",
  riskWindowEnd = 99999,
  endAnchor = "cohort end"
)

n <- nrow(studyPop)
cat(sprintf("Original sample size (from studyPop): n = %d\n", n))

# Check for existing results to resume from
results_file <- file.path(output_dir, "bootstrap_results_iterative.csv")
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
  
  tryCatch({
    results <- runBootstrapIteration(b, studyPop, cohortMethodData, 
                                     cohortMethodData_allcovar, m)
    bootstrap_results[[b]] <- results
    
    # Save this iteration's result immediately
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
    # Save error as NA values
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
    
    bootstrap_results[[b]] <- list(
      bootstrap_iteration = b,
      error = e$message
    )
  })
}

cat("\n=== Bootstrap Complete ===\n")
cat(sprintf("Completed %d iterations\n", length(bootstrap_results)))

################################################################################
# Calculate Bootstrap Statistics
# 
# All bootstrap statistics are calculated on the log(HR) scale (coefficients).
# This is the correct scale for bootstrap inference:
# - Bootstrap SD of log(HR) is calculated directly
# - 95% CI is constructed using quantiles of log(HR)
# - HR values are transformed back (exp) only for display/interpretation
################################################################################

cat("\n=== Calculating Bootstrap Statistics ===\n")
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
cat("\n=== Bootstrap Summary ===\n")
print(summary_df)

cat(sprintf("\nResults saved to: %s/\n", output_dir))
cat("- bootstrap_results_iterative.csv: Progressive save (updated after each iteration)\n")
cat("- bootstrap_results.csv: Final complete results (copy of iterative save)\n")
cat("- bootstrap_results_checkpoint.rds: Checkpoint file (saved every 10 iterations)\n")
cat("- bootstrap_summary.csv: Summary statistics (SEs and CIs)\n")
cat("\nNote: If the script crashes, it will automatically resume from existing results.\n")
