library(dplyr)
library(tidyr)
library(tibble)
library(data.table)
library(survival)
#library(CohortMethod)
#library(Cyclops)

# for progress bar
library(progress)

### -- 0) Read in data

ps = readRDS("results/ps_study.rds") # this reads in the propensity score results from find_ps_legend_script.R

# important: make sure you read in the cohort with all variables
cohortMethodData <- CohortMethod::loadCohortMethodData("results/cohortMethodData_t1788868_c1788867_o1788866_allcovar.zip")

# from cox censoring model that was fit:
Cox_censoring = readRDS("results/Cox_censoring.rds")

### -- 1) Restrict CohortMethodData to the PS analysis cohort --
# make sure cohort method data contains individuals included in PS analysis 

# use ps dataframe to get the relevant rowIds in cohortMethodData 
cohort_ids <- ps$rowId

# filter cohort 
cohortMethodData$cohorts <- collect(cohortMethodData$cohorts) %>%
  filter(rowId %in% cohort_ids)

# filter covariates
cohortMethodData$covariates <- collect(cohortMethodData$covariates) %>%
  filter(rowId %in% cohort_ids)

# filter outcomes 
cohortMethodData$outcomes <- collect(cohortMethodData$outcomes) %>%
  filter(rowId %in% cohort_ids)

### -- 2) build survival df from the ps table - # for sanity check, cumc had 42606 total
outcomes_for_cyclops <- ps %>%
  transmute(
    rowId = rowId,
    y = if_else(outcomeCount == 0, 1, 0),  # 1 = censored, 0 = event
    time = survivalTime, # survival time here = days till end of cohort
    treatment = treatment, 
    iptw = iptw
  )

### -- 3) extract the non-zero covariates used by the censoring model

coefs <- coef(Cox_censoring) # Extract coefficients
non_zero_coefs <- coefs[coefs != 0] # Keep only nonzero ones
non_zero_df <- data.frame(
  covariateId = names(non_zero_coefs),
  coefficient = as.numeric(non_zero_coefs)
) # Turn into a dataframe -- 1017 rows


# filter to non-zero covariates only
filtered_covariates <- collect(cohortMethodData$covariates) %>%
  filter(covariateId %in% names(non_zero_coefs))  # Use names() not non_zero_names if already extracted

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


# sanity checks:
# dim(X_matrix) # sanity check; should be 42606 ppl x 1017 covariates

# Create coeff_vector matching X_matrix
coeff_vector <- non_zero_coefs[colnames(X_matrix)]

# Check that lengths match
stopifnot(length(coeff_vector) == ncol(X_matrix))

# (more safety) Check names match
stopifnot(all(names(coeff_vector) == colnames(X_matrix)))

### -- 4) Align outcome rows to the covariate matrix row order
# outcomes_for_cyclops is the survival df, where y = censoring status
outcomes_df = data.table(collect(outcomes_for_cyclops))

# match it up with X_matrix
outcomes_df <- outcomes_df %>%
  arrange(match(as.character(rowId), rownames(X_matrix)))

# Check if the order is correct
all.equal(as.character(outcomes_df$rowId), rownames(X_matrix))

######################################
######################################
######################################
######################################
## now to estimate survival for ipcw # 
######################################
######################################
######################################
######################################


### -- 5) baseline survival for censoring [numerator for stabilized IPCW]

# non covar-adjusted survival (for stabilized weights)
cox_baseline <- coxph(Surv(time, y) ~ 1, data = outcomes_df)
baseline_surv <- survfit(cox_baseline) # median survival is 914
# Get baseline survival function
baseline_surv_fun <- stepfun(baseline_surv$time, c(1, baseline_surv$surv))

### -- 6) Compute a Breslow baseline cumulative hazard for covariate-adjusted censoring model [denom]

# fxn
breslow_est <- function(time, status, X, B){
  data <- data.frame(time, status, X)
  data <- data[order(data$time), ]
  t <- unique(data$time)
  k <- length(t)
  h <- rep(0,k)
  LP_indiv <- X %*% B # row wise X * B
  
  # Initialize progress bar
  pb <- progress::progress_bar$new(
    format = "  Progress [:bar] :percent in :elapsed, ETA: :eta",
    total = length(1:k),  # Set the total number of iterations
    clear = FALSE,  # Keeps progress bar visible after completion
    width = 60      # Customize the width of the bar
  )
  
  for(i in 1:k) {
    
    # Update the progress bar
    pb$tick()
    
    lp <- (LP_indiv)[data$time>=t[i]]
    risk <- exp(lp)
    h[i] <- sum(data$status[data$time==t[i]]) / sum(risk)
  }
  
  res <- cumsum(h)
  return(res)
}

# manual calculation
H0 <- breslow_est(time=outcomes_df$time, status=outcomes_df$y, 
                  X=X_matrix, B=coeff_vector)

##### interpolate for the times we want, which is datai$Tstart
haz_step_fun <- stepfun(sort(unique(outcomes_df$time)), c(0, H0))  # Include initial survival at t=0

## --7) Convert to person–time (long) format with intervals

# fxn
transform.data <- function(data, cut.times) 
{
  # Define Tstart and the indicator "censored" (y = 1 is censored):
  data$Tstart <- 0
  data$ami = 1-data$y
  
  # Times at which to split the intervals:
  # cut.times <- data$time
  
  # Split data with event = y: censoring = dependent var
  data.long <- survSplit(data = data,
                         cut = cut.times,
                         end = "time",
                         start = "Tstart", # all at 0 since anchored on entry date
                         event = "y") # dep var is censoring
  data.long %>% data.table()
  data.long <- data.long[order(data.long$rowId,
                               data.long$time),]
  
  # Split data with event = AMI: health outcome = dependent var
  # this ensures that we have the correct AMI indicator
  data.long.cens <- survSplit(data,
                              cut=cut.times,
                              end = "time",
                              start = "Tstart", # all at 0 since anchored on entry date
                              event = "ami")
  data.long.cens <- data.long.cens[order(data.long.cens$rowId,
                                         data.long.cens$time),]
  # Add "censored" indicator to long data format:
  data.long$ami <- data.long.cens$ami
  data.long$rowId <- as.numeric(data.long$rowId)
  # Return long data format:
  return(data.long)
}

# determining cut times and running the function 
# cut times should every time risk set changes, which is length(unique(outcomes_df$time))
# however that would mean lots of cut times for something like this

# NOTE: you might want to determine the vector of cut times so that it is more granular
# in the beginning (ex: every 10 days), then as the data is more sparse later/observations
# are more spaced out -- can space out the cut times

# hist(outcomes_df$time)
dist = summary(outcomes_df$time) # max is 9647 days = 26.47 years; mean is 8 years; median is 7.84 years

# one way to do cut times is space them out evenlly till the max
cut.times = seq(from = 60, to = floor(max(outcomes_df$time) / 60) * 60, by = 60) 
# another way is to do it more in the beginning then less till the end since dist of time is skewed
# cut.times = c(seq(from = 1, to = dist[2], by = 5),# go by 5 till the first quartile
#              seq(from = dist[2] + 5, to = dist[3], by = 10), # by 10 to the median
#              seq(from = dist[3] + 10, to = dist[5], by = 30), # by 30 to Q3
#              seq(from = dist[5] + 30, to = floor(max(outcomes_df$time) / 100) * 100, by = 100)) # by 100 rest of way

# cut.times = unique(outcomes_df$time) # I have tried this, but it makes the long data REALLY long 
# and doesn't really change the results

# transform to long data
outcomes_df.long <- transform.data(outcomes_df, cut.times) 


## -- 8)  Compute each person’s estimated probability of remaining uncensored up to each interval start
# MAIN IPCW CALCULATION STEP

# now to get the fitted curve for each person
# per person linear predictor
eta <- X_matrix %*% coeff_vector # this is a vector for every individual
# outcomes_df.long$KZ = NULL

# First: precompute baseline cumulative hazards at *all* Tstart (interval) times
outcomes_df.long$H0_Tstart <- haz_step_fun(outcomes_df.long$Tstart)

# Then: match eta values by rowId
rowid_to_eta <- data.frame(rowId = as.numeric(rownames(X_matrix)), eta = as.numeric(eta))

# Join eta onto outcomes_df.long
outcomes_df.long <- outcomes_df.long %>%
  left_join(rowid_to_eta, by = "rowId")

# Now: directly compute KZ (the conditional survival probability) row-by-row -- denominator for ipcw
outcomes_df.long <- outcomes_df.long %>%
  mutate(KZ = exp( - H0_Tstart * exp(eta) ))

# computer k0 for start time--need this for numerator for stabilized weights
### now for stabilized ipcw we need K0t
outcomes_df.long$K0_ti <- baseline_surv_fun(outcomes_df.long$Tstart) # compute baseline at each Tstart

## -- 9)  Compute stabilized and unstabilized weights
# some post process
outcomes_df.long$Unstab_ipcw = 1/outcomes_df.long$KZ
outcomes_df.long$Stab_ipcw = outcomes_df.long$K0_ti/outcomes_df.long$KZ

# save
write.csv(outcomes_df.long, "results/survival_weights_endObsDate.csv")
