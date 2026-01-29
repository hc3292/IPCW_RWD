# run after running find_kZ_script.R
library("dplyr")
library("data.table")
library("survival")

# outcomes_df.long = read.csv("survival_weights_endObsDate.csv") %>% data.table()

# distribution of weights
# summary(outcomes_df.long$Unstab_ipcw)
# summary(outcomes_df.long$Stab_ipcw)

# Truncate stabilized IPCW weights
lo <- quantile(outcomes_df.long$Stab_ipcw, 0.01, na.rm = TRUE)
hi <- quantile(outcomes_df.long$Stab_ipcw, 0.99, na.rm = TRUE)

outcomes_df.long$Stab_ipcw_trunc <- 
  pmin(pmax(outcomes_df.long$Stab_ipcw, lo), hi)


# truncate unstabilized IPCW weights
lo <- quantile(outcomes_df.long$Unstab_ipcw, 0.01, na.rm = TRUE)
hi <- quantile(outcomes_df.long$Unstab_ipcw, 0.99, na.rm = TRUE)

outcomes_df.long$Unstab_ipcw_trunc <- 
  pmin(pmax(outcomes_df.long$Unstab_ipcw, lo), hi)


# trucate iptw weights (alr stabilized)
outcomes_df.long <- outcomes_df.long %>%
  group_by(Tstart) %>%
  mutate(
    lower_iptw = quantile(iptw, 0.01, na.rm = TRUE),
    upper_iptw = quantile(iptw, 0.99, na.rm = TRUE),
    iptw_trunc = pmin(pmax(iptw, lower_iptw), upper_iptw)
  ) %>%
  ungroup()


# combined weights
outcomes_df.long$comb = outcomes_df.long$Stab_ipcw_trunc * outcomes_df.long$iptw_trunc


#####################
# FIT MODELS
#####################

# LEGEND THZ VALUE FOR AMI IS 0.84 (0.75-0.95)

# UNADJUSTED
fit_unadjusted <- coxph(Surv(Tstart, time, ami) ~ treatment, data=outcomes_df.long, 
                        id=rowId)
summary(fit_unadjusted) # coef = -0.73228  SE = 0.08192; HR = 2.07982
# -0.58


### IPCW STABILIZED 
# fit_ipcw_Stab <- coxph(Surv(Tstart, time, ami) ~ treatment , data=outcomes_df.long, 
#                       id=rowId, weights = Stab_ipcw) # unstab exp(coef) = 0.51 (1/0.51 = 1.95)
# summary(fit_ipcw_Stab) # coef = -0.77947, SE = 0.09738; HR = 0.45865

fit_ipcw_Stab_trunc <- coxph(Surv(Tstart, time, ami) ~ treatment , data=outcomes_df.long, 
                             id=rowId, weights = Stab_ipcw_trunc) # unstab exp(coef) = 0.51 (1/0.51 = 1.95)
summary(fit_ipcw_Stab_trunc) # coef = -0.78541, SE = 0.09659; HR = 0.45593

# -0.638, HR - 0.53

# IPTW weights ====== CONFOUNDING

#fit_iptw <- coxph(Surv(Tstart, time, ami) ~ treatment , data=outcomes_df.long, 
#                        id=rowId, weights = iptw) 
#summary(fit_iptw) # coef -0.3141; SE =  0.1408; HR = 0.7305 

fit_iptw_trunc <- coxph(Surv(Tstart, time, ami) ~ treatment , data=outcomes_df.long, 
                        id=rowId, weights = iptw_trunc) 
summary(fit_iptw_trunc) # coef = -0.2098; SE = 0.1022; HR = 0.8108

# -0.15, HR = 0.85

# combined
fit_comb <- coxph(Surv(Tstart, time, ami) ~ treatment , data=outcomes_df.long, 
                  id=rowId, weights = comb) 
summary(fit_comb) # coef = -0.08742 SE =  0.11319; exp: 0.91629

# -0.16, 0.85