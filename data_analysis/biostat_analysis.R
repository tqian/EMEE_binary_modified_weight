# Yihan Bao
# 2021.04.06

source("estimator_implementation/WCLS_modified.R")
source("estimator_implementation/WCLS_original.R")
source("estimator_implementation/GEE_estimators.R")

library(reshape)
library(kableExtra)
library(knitr)

# Please use both EMEE.mod and EMEE to analyze the data.
# Two estimands: one setting with S_t = 1 (only intercept), the other setting with S_t = (1, days since download)
# Control variables for both estimands and both estimators: 
#   the continuous variables of age, AUDIT score (continuous), days since download, 
#   the categorical variables of gender and employment type, the time-varying variables 
#   “did the use user the app before 8 pm that day?” and 
#   “did the users use the app after 9 pm the day before?”

dta <- readRDS("FINAL Dataset_A.rds")
dta <- as.data.frame(dta)
dta$days_since_download <- dta$days_since_download - 1
dta$after_9pm_day_before_before_8pm <- dta$after_9pm_day_before* dta$before_8pm
  
head(dta)
control_vars = c('age','AUDIT_score','days_since_download',
                 'gender','employment_type','before_8pm','after_9pm_day_before',
                 'after_9pm_day_before_before_8pm', 'treatment_day_before')

estimators <- function(Delta, moderator_var){ # , control_vars = control_vars
  ### creating features needed for WCLS function ###
  total_T = 30 # number of decision points for each individual
  sample_size = 349 # number of unique subject
  for (t in 1:total_T) {
    row_index <- seq(from = t, by = total_T, length = sample_size)
    for (j in 1:length(row_index)) {
      # notice here that when t+ Delta - 1 > total_T, we take the product up to total_T.
      end_index = min((row_index[j]+Delta-1),j*total_T)
      dta$proximal_outcome[row_index[j]] <- 1 - prod(1-dta$outcome_24hour[row_index[j]:end_index])
      # get the index of R's when Y = 1 for the first time, 
      # if Y = 0, assign k = Delta - 1.
      dta$k[row_index[j]] <- ifelse(dta$proximal_outcome[row_index[j]] == 0, (Delta-1),
                                    min(which(dta$outcome_24hour[row_index[j]:end_index] == 1))-1)
      
      dta$treatment_day_before[row_index[j]] <- ifelse(dta$days_since_download[row_index[j]] == 0, 
                                                       0, dta$treatment[row_index[j]-1])
    }
  }
  
  t_quantile <- qt(0.975, 349 - 1 - 2) 
  
  ### original EMEE estimator ###
  fit_wcls <- weighted_centered_least_square_withDelta(
    dta = dta,
    id_varname = "ID",
    decision_time_varname = "days_since_download",
    treatment_varname = "treatment",
    outcome_varname = "proximal_outcome",
    control_varname = control_vars,
    moderator_varname = moderator_var,
    rand_prob_varname = "prob_A",
    rand_prob_tilde_varname = NULL,
    rand_prob_tilde = 0.2,
    estimator_initial_value = NULL,
    Delta = Delta
  )
   
  ### modified-EMEE estimator ###
  fit_wcls_new <- weighted_centered_least_square_withDelta_new(
    dta = dta,
    id_varname = "ID",
    decision_time_varname = "days_since_download",
    treatment_varname = "treatment",
    outcome_varname = "proximal_outcome",
    control_varname = control_vars,
    moderator_varname = moderator_var,
    rand_prob_varname = "prob_A",
    rand_prob_tilde_varname = NULL,
    rand_prob_tilde = 0.2,
    estimator_initial_value = NULL,
    Delta = Delta
  )
  
  lci <- fit_wcls$beta_hat - t_quantile * fit_wcls$beta_se_adjusted
  rci <- fit_wcls$beta_hat + t_quantile * fit_wcls$beta_se_adjusted
    
  lci_new <- fit_wcls_new$beta_hat - t_quantile * fit_wcls_new$beta_se_adjusted
  rci_new <- fit_wcls_new$beta_hat + t_quantile * fit_wcls_new$beta_se_adjusted
  re <- (fit_wcls$beta_se_adjusted/fit_wcls_new$beta_se_adjusted)^2
  
  if (is.null(moderator_var)){
    results <- c(Delta, fit_wcls$beta_hat, fit_wcls$beta_se_adjusted,
                 rci, lci, "IPW", "beta0", 1)
    results <- rbind(results,
                     c(Delta, fit_wcls_new$beta_hat, fit_wcls_new$beta_se_adjusted,
                       rci_new, lci_new, "pd_IPW", "beta0", re))
  }else{
    results <- c(Delta, fit_wcls$beta_hat[1], fit_wcls$beta_se_adjusted[1],
                 rci[1], lci[1], "IPW", "beta0", 1)
    results <- rbind(results,
                     c(Delta, fit_wcls_new$beta_hat[1], fit_wcls_new$beta_se_adjusted[1],
                       rci_new[1], lci_new[1], "pd_IPW", "beta0", re[1]))
    results <- rbind(results,
                     c(Delta, fit_wcls$beta_hat[2], fit_wcls$beta_se_adjusted[2],
                       rci[2], lci[2], "IPW", "beta1", 1))
    results <- rbind(results,
                     c(Delta, fit_wcls_new$beta_hat[2], fit_wcls_new$beta_se_adjusted[2],
                       rci_new[2], lci_new[2], "pd_IPW", "beta1", re[2]))
  }
  results
  }


Deltas = 1:7
moderator_vars <- c("days_since_download", "treatment_day_before", 
                    "after_9pm_day_before_before_8pm", NULL)

for (moderator_var in moderator_vars) {
  results = c()
  for (Delta in Deltas) {
    result = estimators(Delta = Delta, moderator_var = moderator_var)
    results = rbind(results, result)
  }
  saveRDS(results, paste0("moderator=", moderator_var, ".RDS"))
}






