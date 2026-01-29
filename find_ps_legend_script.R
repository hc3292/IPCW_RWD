## if necessary can run in terminal

library(DatabaseConnector)
library(CohortMethod)

# install.packages("rJava", type="source")

# install.packages("DatabaseConnector") #6.3.2 was last version
#install.packages("CohortMethod") #5.3.0 was last version
#remotes::install_github("OHDSI/Cyclops", ref = 'v3.4.0') # 3.4.1 was last version

# library(knitr)
# library(tidyr)
# library(dplyr)
# library(Cyclops)


cohortMethodData <- CohortMethod::loadCohortMethodData("results/cohortMethodData_t1788868_c1788867_o1788866.zip")

# For the LEGEND-HTN study, we pick target exposure thiazdes (atlas id 1788868), 
# comparator exposure ACEi (atlas id 1788867), and outcome  AMI (atlas id 1788866). 

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

getAttritionTable(studyPop)

# create propensity model
ps <- createPs(cohortMethodData = cohortMethodData, 
               population = studyPop,
               control = createControl(threads = MAX_THREADS))
saveRDS(ps, "results/ps_study.rds")

# outcome model
# ps = readRDS("ps_study.rds") # pick up where you left off if you need

# unadjusted outcome model
outcomeModel <- fitOutcomeModel(population = ps,
                                modelType = "cox",
                                inversePtWeighting = TRUE)

# outcome model but with matching
matchedPop <- matchOnPs(ps, caliper = 0.2)
# pop size after matching is 22154
outcomeModel_matching <- fitOutcomeModel(population = matchedPop,
                                         modelType = "cox")

bal = computeCovariateBalance(
  population = matchedPop,
  cohortMethodData,
  subgroupCovariateId = NULL,
  maxCohortSize = 250000,
  covariateFilter = NULL
)
