cohortMethodData_allcovar_file <- "results/cohortMethodData_t1788868_c1788867_o1788866_allcovar.zip"

### 1) load cohort method data
# important: use the "all covar" version
cohortMethodData_allcovar <- loadCohortMethodData(cohortMethodData_allcovar_file)

### 2) build the censored_cohort table 

# identify who is censored (i.e. 1 - outcome)
# in the case of how cyclops data is written, we identify censored
# individuals via who is *not* in the outcomes table 
cohortMethodData_allcovar$censored_cohort <- cohortMethodData_allcovar$cohorts %>%
  anti_join(cohortMethodData_allcovar$outcomes, by = "rowId")

### 3) now build the "censored_outcomes" table with daysToCohortEnd (or daysToObsEnd) column 

cohortMethodData_allcovar$censored_outcomes <- cohortMethodData_allcovar$censored_cohort %>%
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
cohortMethodData_allcovar$censored_outcomes = cohortMethodData_allcovar$censored_outcomes %>% filter(outcomeId == 1)

### 4) Construct Cyclops Dataset 

# Start by extracting the relevant tables from the Andromeda object
cohorts <- cohortMethodData_allcovar$cohorts
tx = collect(cohortMethodData_allcovar$cohorts)[,c("rowId", "treatment")]
censored_outcomes <- cohortMethodData_allcovar$censored_outcomes

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

# Add the new table to the Andromeda object
cohortMethodData_allcovar$outcomes_for_cyclops <- outcomes_for_cyclops

CohortMethod::saveCohortMethodData(cohortMethodData_allcovar, "cohortMethodData_bootstrap_allcovar_censoring.zip")

