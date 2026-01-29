### 

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

getFollowUpDistribution(studyPop)

plotFollowUpDistribution(
  studyPop,
  targetLabel = "Thiazide Diuretics",
  comparatorLabel = "ACEi",
  yScale = "percent",
  logYScale = FALSE,
  dataCutoff = 0.95,
  title = NULL,
  fileName = NULL
)

drawAttritionDiagram(
  studyPop,
  targetLabel = "Target",
  comparatorLabel = "Comparator",
  fileName = NULL
)


covar_test = cohortMethodData$covariates %>% collect()
length(unique(covar_test$covariateId))


balance = computeCovariateBalance(
  studyPop,
  cohortMethodData,
  subgroupCovariateId = NULL,
  maxCohortSize = 250000,
  covariateFilter = NULL
)

# number of covariates
nrow(balance)

tab1 = createCmTable1(
  balance,
  specifications = getDefaultCmTable1Specifications(),
  beforeTargetPopSize = 22970,
  beforeComparatorPopSize = 24269,
  afterTargetPopSize = 22970,
  afterComparatorPopSize = 24269,
  beforeLabel = "Before matching",
  afterLabel = "After matching",
  targetLabel = "Thiazide Diuretics",
  comparatorLabel = "ACEi",
  percentDigits = 1,
  stdDiffDigits = 2
)

write.csv(tab1, "results/table_1.csv")

## get covariates 
covars = getPsModel(ps, cohortMethodData)
covars$coefficient = round(covars$coefficient, 3)

top200 <- covars[order(abs(covars$coefficient), decreasing = TRUE), ][1:51, ]


write.csv(covars, "results/covars.csv")
write.csv(top200, "results/top_covars.csv")



