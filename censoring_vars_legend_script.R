# library(DatabaseConnector)
# library(CohortMethod)

# install.packages("rJava", type="source")

# install.packages("DatabaseConnector") #6.3.2 was last version
#install.packages("CohortMethod") #5.3.0 was last version
#remotes::install_github("OHDSI/Cyclops", ref = 'v3.4.0') # 3.4.1 was last version

#########################

# This script constructs a censoring outcome from CohortMethod data and fits a
# penalized Cox proportional hazards model using Cyclops to estimate censoring
# risk as a function of baseline covariates, as a precursor to IPCW.

#########################

# library(knitr)
# library(tidyr)
# library(dplyr)
# library(Cyclops)


### 1) load cohort method data

# important: use the "all covar" version
cohortMethodData <- CohortMethod::loadCohortMethodData("results/cohortMethodData_t1788868_c1788867_o1788866_allcovar.zip")

### 2) build the censored_cohort table 

# identify who is censored (i.e. 1 - outcome)
# in the case of how cyclops data is written, we identify censored
# individuals via who is *not* in the outcomes table 
cohortMethodData$censored_cohort <- cohortMethodData$cohorts %>%
  anti_join(cohortMethodData$outcomes, by = "rowId")

### 3) now build the "censored_outcomes" table with daysToCohortEnd (or daysToObsEnd) column 

cohortMethodData$censored_outcomes <- cohortMethodData$censored_cohort %>%
  # 1) Select the needed column, renaming 'daysToObsEnd' to 'daysToEvent' -- "Event" here means "Censoring event"
  select(
    rowId,
    daysToEvent = daysToObsEnd # can use daysToCohortEnd or daysToObsEnd, depending on your original TAR defn
  ) %>%
  # 2) Add the constant column outcomeId
  mutate(
    outcomeId = 1 # because they are all censored; so censored outcome = 1
  ) %>%
  # 3) Reorder columns
  select(rowId, outcomeId, daysToEvent) 

# censored outcomes only -> add to cohortMethodData
cohortMethodData$censored_outcomes = cohortMethodData$censored_outcomes %>% filter(outcomeId == 1)

### 4) Construct Cyclops Dataset 

# Start by extracting the relevant tables from the Andromeda object
cohorts <- cohortMethodData$cohorts
tx = collect(cohortMethodData$cohorts)[,c("rowId", "treatment")]
censored_outcomes <- cohortMethodData$censored_outcomes

# Create the outcomes_for_cyclops table
outcomes_for_cyclops <- cohorts %>%
  # Perform a left join with censored_outcomes to check if the rowId exists in censored_outcomes
  left_join(censored_outcomes, by = "rowId") %>%
  # Create the "y" column: 1 if rowId is in censored_outcomes, otherwise 0
  mutate(y = ifelse(!is.na(outcomeId), 1, 0)) %>%
  # Create the "time" column based on the condition
  mutate(time = ifelse(!is.na(outcomeId), daysToEvent, daysToObsEnd)) %>%
  # Select only the required columns for the new table
  select(rowId, y, time)

# filter for only those in Studypop
outcomes_for_cyclops <- outcomes_for_cyclops %>%
  semi_join(studyPop %>% distinct(rowId), by = "rowId", copy = TRUE)

# Add the new table to the Andromeda object
cohortMethodData$outcomes_for_cyclops <- outcomes_for_cyclops

# now put it into Cyclops Data
censored_df = convertToCyclopsData(cohortMethodData$outcomes_for_cyclops, 
                                   cohortMethodData$covariates, 
                                   modelType = "cox",
                                   addIntercept = TRUE)

### 5) Fit Regularized Cox Model with Outcome = Censoring Status 
# L1 regularized Cox (cox is specified in previous function, convertToCyclopsData)
lassoPrior <- Cyclops::createPrior(
  priorType = "laplace", 
  useCrossValidation = TRUE
)

cat("\n running Cox censoring model")

Cox_censoring <- fitCyclopsModel(censored_df,
                                 prior = lassoPrior,
                                 control = createControl(threads = MAX_THREADS))

saveRDS(Cox_censoring, "results/Cox_censoring.rds")

