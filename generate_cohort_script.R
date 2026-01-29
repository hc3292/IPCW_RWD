### THIS GENERATES THE COHORTS
# does covariate exclusion for thiazides and ace inhibitors

library(DatabaseConnector)
library(CohortMethod)

# library(knitr)
# library(tidyr)
# library(dplyr)

# library(survival)

# For the LEGEND-HTN study, we pick target exposure thiazdes (atlas id 1788868), 
# comparator exposure ACEi (atlas id 1788867), and outcome  AMI (atlas id 1788866). 
# The cohort generation codes are to make sure we are using the same data to run 
# the study. Feel free to skip the cohort generation part if you have already the 
# LEGEND-HTN data. 

# Again, by the end of cohort generation, we should have a CohortMethodData with 
# targetId 1788868, comparatorId 1788867, and outcomeId 1788866. 

################################################################################
######################## INput Database Information ############################
################################################################################

connectionDetails <- createConnectionDetails(dbms="sql server",
                                             server="")

conn <- DatabaseConnector::connect(connectionDetails)

cdmDatabaseSchema <- "" 
targetCohortTable <- "legend_results"
targetDatabaseSchema <- ""
tempEmulationSchema <- ""
vocabularyDatabaseSchema <- "" 
cdmSourceName <- "CUMC" 
cdmVersion <- "" 

################################################################################
###################### Generate LEGEND-HTN cohort ##############################
################################################################################

baseUrl <- "http://api.ohdsi.org:80/WebAPI" 

cohortIds <- c(1788868, # thiazides
               1788867, # ACEi
               1788866) # AMI

cohortDefinitionSet <- ROhdsiWebApi::exportCohortDefinitionSet(
  baseUrl = baseUrl,
  cohortIds = cohortIds,
  generateStats = TRUE
)

cohortTableNames <- CohortGenerator::getCohortTableNames(cohortTable = targetCohortTable)
cohortTableNames

CohortGenerator::createCohortTables(
  connectionDetails = connectionDetails,
  cohortTableNames = cohortTableNames,
  cohortDatabaseSchema = targetDatabaseSchema,
  incremental = T # if TRUE, then will not generate if these tables already exist
)

CohortGenerator::generateCohortSet(
  connectionDetails = connectionDetails,
  cdmDatabaseSchema = cdmDatabaseSchema,
  #tempEmulationSchema = tempEmulationSchema, 
  cohortDatabaseSchema = targetDatabaseSchema,
  cohortTableNames = cohortTableNames,
  cohortDefinitionSet = cohortDefinitionSet,
  incremental = FALSE, 
  stopOnError = FALSE
)

# the below is if we don't want to exclude anything
# covariateSettings and create cohortMethodData

covariateSettings <- FeatureExtraction::createDefaultCovariateSettings(
  excludedCovariateConceptIds = c(
    # thiazide
    907013,  # metolazone
    978555,  # indapamide
    974166,  # hydrochlorothiazide
    1395058, # chlorthalidone
    # ACEi
    1342439, # trandolapril			
    1334456, # ramipril		
    1331235, # quinapril			
    1373225, # perindopril			
    1310756, # moexipril				
    1308216, # lisinopril		
    1363749, # fosinopril			
    1341927, # enalapril			
    1340128  # captopril
  ), 
  addDescendantsToExclude = TRUE
)

# Create CohortMethodData. This step may take awhile
cohortMethodData <- CohortMethod::getDbCohortMethodData(
  connectionDetails = connectionDetails, 
  cdmDatabaseSchema = cdmDatabaseSchema,
  #tempEmulationSchema = tempEmulationSchema,
  targetId = 1788868,      # thiazides
  comparatorId = 1788867,  # ACEi
  outcomeIds = 1788866,    # AMI
  exposureDatabaseSchema = targetDatabaseSchema,
  exposureTable = targetCohortTable,
  outcomeDatabaseSchema = targetDatabaseSchema,
  outcomeTable = targetCohortTable,
  studyStartDate = '20000101',
  studyEndDate = '20241231',
  covariateSettings = covariateSettings
)
CohortMethod::saveCohortMethodData(cohortMethodData, "results/cohortMethodData_t1788868_c1788867_o1788866.zip")


# the below does covariate settings including ALL covariate
# for LSPS, we exclude the treatment drugs in the covariates
# but for IPCW, we want them all

covariateSettings <- FeatureExtraction::createDefaultCovariateSettings()

cohortMethodData_allcovar <- CohortMethod::getDbCohortMethodData(
  connectionDetails = connectionDetails, 
  cdmDatabaseSchema = cdmDatabaseSchema,
  #tempEmulationSchema = tempEmulationSchema,
  targetId = 1788868,      # thiazides
  comparatorId = 1788867,  # ACEi
  outcomeIds = 1788866,    # AMI
  exposureDatabaseSchema = targetDatabaseSchema,
  exposureTable = targetCohortTable,
  outcomeDatabaseSchema = targetDatabaseSchema,
  outcomeTable = targetCohortTable,
  studyStartDate = '20000101',
  studyEndDate = '20241231',
  covariateSettings = covariateSettings
)

CohortMethod::saveCohortMethodData(cohortMethodData_allcovar, "results/cohortMethodData_t1788868_c1788867_o1788866_allcovar.zip")
