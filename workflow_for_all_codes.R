### this is the code to run all codes
library(dplyr)
library(survival)
library(data.table)
# library(survminer)

MAX_THREADS = 4

# LSPS and IPTW workflow
source("generate_cohort_script.R") # this generates the cohorts for LSPS for treatment weights (excludes the treatment drugs)
source("find_ps_legend_script.R") # this runs LSPS + outcome model with IPTW weights

# IPCW workflow
source("censoring_vars_legend_script.R") # generates the Cox censoring model
source("find_kZ_script.R") # gets the IPCW weights and fits the outcome model with IPCW 

# combined weights workflow (IPTW + IPCW)
source("weights_script.R") # uses previous results to get the combined weights


### BOOTSTRAPPING PROCESS
# currently only m of n bootstrapping