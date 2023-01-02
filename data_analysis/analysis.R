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
head(dta)

###### data engineering (creating other features needed for WCLS function)######

total_T = 30 # number of decision points for each individual
sample_size = 349 # number of unique subject
Delta = 3

  # "outcome_24hour" is "R_t" in our setting, we need to calculate the "proximal_outcome"
  # which is Y_t = max(R_{t+1}, R_{t+2}, R_{t+3})
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
    }
}  

dta$days_since_download <- dta$days_since_download - 1




###################
# 1. marginal excursion effect

control_vars = c('age','AUDIT_score','days_since_download',
                 'gender','employment_type','before_8pm','after_9pm_day_before')
moderator_vars = NULL

# 1.1 modified-EMEE estimator
fit_wcls_new <- weighted_centered_least_square_withDelta_new(
  dta = dta,
  id_varname = "ID",
  decision_time_varname = "days_since_download",
  treatment_varname = "treatment",
  outcome_varname = "proximal_outcome",
  control_varname = control_vars,
  moderator_varname = moderator_vars,
  rand_prob_varname = "prob_A",
  rand_prob_tilde_varname = NULL,
  rand_prob_tilde = 0.2,
  estimator_initial_value = NULL,
  Delta = 3
)
fit_wcls_new$beta_hat # 0.127397
fit_wcls_new$beta_se_adjusted # 0.0268428 
# 4939 data points used

  # estimator, SE, 95% CI, p-value
t_quantile <- qt(0.975, 349 - 1 - 2) 
fit_wcls_new$beta_hat # 0.127397 
fit_wcls_new$beta_se_adjusted # 0.0268428 
rbind(fit_wcls_new$beta_hat - t_quantile * fit_wcls_new$beta_se_adjusted,
      fit_wcls_new$beta_hat + t_quantile * fit_wcls_new$beta_se_adjusted)
# (0.07460137,  0.18019258)
2 * pt(abs(fit_wcls_new$beta_hat) / fit_wcls_new$beta_se_adjusted, 349 - 1 - 2, lower.tail = FALSE)
# 3.040731e-06

# fit_wcls_new_result <- c(fit_wcls_new$beta_hat, fit_wcls_new$beta_se_adjusted, 
#                          rbind(fit_wcls_new$beta_hat - t_quantile * fit_wcls_new$beta_se_adjusted,
#                                fit_wcls_new$beta_hat + t_quantile * fit_wcls_new$beta_se_adjusted),
#                          2 * pt(abs(fit_wcls_new$beta_hat) / fit_wcls_new$beta_se_adjusted, 349 - 1 - 2, lower.tail = FALSE))

# 1.2 original EMEE estimator
fit_wcls <- weighted_centered_least_square_withDelta(
  dta = dta,
  id_varname = "ID",
  decision_time_varname = "days_since_download",
  treatment_varname = "treatment",
  outcome_varname = "proximal_outcome",
  control_varname = control_vars,
  moderator_varname = moderator_vars,
  rand_prob_varname = "prob_A",
  rand_prob_tilde_varname = NULL,
  rand_prob_tilde = 0.2,
  estimator_initial_value = NULL,
  Delta = 3
)

  # estimator, SE, 95% CI, p-value
t_quantile <- qt(0.975, 349 - 1 - 2) 

fit_wcls$beta_hat # 0.07867746 
fit_wcls$beta_se_adjusted # 0.0557727
rbind(fit_wcls$beta_hat - t_quantile * fit_wcls$beta_se_adjusted,
      fit_wcls$beta_hat + t_quantile * fit_wcls$beta_se_adjusted)
# (-0.03101873,  0.18837365)
2 * pt(abs(fit_wcls$beta_hat) / fit_wcls$beta_se_adjusted, 349 - 1 - 2, lower.tail = FALSE)
# 0.1592373 


##################
# 2. with moderator 


moderator_vars = "days_since_download"
control_vars = c('age','AUDIT_score','days_since_download',
                 'gender','employment_type','before_8pm','after_9pm_day_before')

# modified-EMEE estimator
fit_wcls_newmod <- weighted_centered_least_square_withDelta_new(
  dta = dta,
  id_varname = "ID",
  decision_time_varname = "days_since_download",
  treatment_varname = "treatment",
  outcome_varname = "proximal_outcome",
  control_varname = control_vars,
  moderator_varname = moderator_vars,
  rand_prob_varname = "prob_A",
  rand_prob_tilde_varname = NULL,
  rand_prob_tilde = 0.2,
  estimator_initial_value = NULL,
  Delta = 3
)
fit_wcls_newmod$beta_hat # 0.195212334        -0.005523186  
fit_wcls_newmod$beta_se # 0.045330790         0.003225971 
# 4939 data points used

# estimator, SE, 95% CI, p-value
t_quantile <- qt(0.975, 349 - 1 - 2) 

fit_wcls_newmod$beta_hat # 0.195212334        -0.005523186  
fit_wcls_newmod$beta_se_adjusted # 0.045972948         0.003269156 
rbind(fit_wcls_newmod$beta_hat - t_quantile * fit_wcls_newmod$beta_se_adjusted,
      fit_wcls_newmod$beta_hat + t_quantile * fit_wcls_newmod$beta_se_adjusted)
# (0.1047907, 0.2856339)        
# (-0.0119531063, 0.0009067337)
2 * pt(abs(fit_wcls_newmod$beta_hat) / fit_wcls_newmod$beta_se_adjusted, 349 - 1 - 2, lower.tail = FALSE)
# 2.796715e-05        9.202757e-02 



# original EMEE estimator
fit_wclsmod <- weighted_centered_least_square_withDelta(
  dta = dta,
  id_varname = "ID",
  decision_time_varname = "days_since_download",
  treatment_varname = "treatment",
  outcome_varname = "proximal_outcome",
  control_varname = control_vars,
  moderator_varname = moderator_vars,
  rand_prob_varname = "prob_A",
  rand_prob_tilde_varname = NULL,
  rand_prob_tilde = 0.2,
  estimator_initial_value = NULL,
  Delta = 3
)

fit_wclsmod$beta_hat # 0.28392612         -0.01579804 
fit_wclsmod$beta_se # 0.080722926         0.005009648 

# estimator, SE, 95% CI, p-value
t_quantile <- qt(0.975, 349 - 1 - 2) 

fit_wclsmod$beta_hat # 0.284         -0.016 
fit_wclsmod$beta_se_adjusted # 0.081         0.005
rbind(fit_wclsmod$beta_hat - t_quantile * fit_wclsmod$beta_se_adjusted,
      fit_wclsmod$beta_hat + t_quantile * fit_wclsmod$beta_se_adjusted)
# (0.1226349, 0.4452174)        
# (-0.025801483, -0.005794598)
2 * pt(abs(fit_wclsmod$beta_hat) / fit_wclsmod$beta_se_adjusted, 349 - 1 - 2, lower.tail = FALSE)
# 0.0006027207        0.0020519879


####### for debugging use #######
dta = dta
id_varname = "ID"
decision_time_varname = "days_since_download"
treatment_varname = "treatment"
outcome_varname = "proximal_outcome"
control_varname = control_vars
moderator_varname = moderator_vars
rand_prob_varname = "prob_A"
rand_prob_tilde_varname = NULL
rand_prob_tilde = 0.2
estimator_initial_value = NULL
avail_varname = NULL
dta[1:100,c("ID","treatment", "primary_outcome", "k")]


######### Supplementary Materials #########
# 1. marginal excursion effect

# estimator, SE, 95% CI, p-value
t_quantile <- qt(0.975, 349 - 1 - 2) 
fit_wcls_new$alpha_hat 
fit_wcls_new$alpha_se_adjusted 
rbind(fit_wcls_new$alpha_hat - t_quantile * fit_wcls_new$alpha_se_adjusted,
      fit_wcls_new$alpha_hat + t_quantile * fit_wcls_new$alpha_se_adjusted)
2 * pt(abs(fit_wcls_new$alpha_hat) / fit_wcls_new$alpha_se_adjusted, 349 - 1 - 2, lower.tail = FALSE)

ctrl_names <- rep(names(fit_wcls_new$alpha_hat), each = 2)
ctrl_est <- round(c(rbind(fit_wcls_new$alpha_hat, fit_wcls$alpha_hat)), 3)
ctrl_se <- round(c(rbind(fit_wcls_new$alpha_se_adjusted, fit_wcls$alpha_se_adjusted)), 3)
# ctr_ci <- round(rbind(cbind(fit_wcls_new$alpha_hat - t_quantile * fit_wcls_new$alpha_se_adjusted,
#                       fit_wcls_new$alpha_hat + t_quantile * fit_wcls_new$alpha_se_adjusted), 
#                 cbind(fit_wcls$alpha_hat - t_quantile * fit_wcls$alpha_se_adjusted,
#                       fit_wcls$alpha_hat + t_quantile * fit_wcls$alpha_se_adjusted)), 3)

ctr_cilower <- round(c(rbind(fit_wcls_new$alpha_hat - t_quantile * fit_wcls_new$alpha_se_adjusted, 
              fit_wcls$alpha_hat - t_quantile * fit_wcls$alpha_se_adjusted)),3)
ctr_ciupper <- round(c(rbind(fit_wcls_new$alpha_hat + t_quantile * fit_wcls_new$alpha_se_adjusted, 
                             fit_wcls$alpha_hat + t_quantile * fit_wcls$alpha_se_adjusted)),3)
ctr_ci <- cbind(ctr_cilower, ctr_ciupper)
ctr_ci <- apply(ctr_ci,1, function(x) paste("(", paste(x[1],x[2], sep = ","), ")", sep = ""))
ctr_pval <- round(c(rbind(2 * pt(abs(fit_wcls_new$alpha_hat) / fit_wcls_new$alpha_se_adjusted, 
                                 349 - 1 - 2, lower.tail = FALSE),
                          2 * pt(abs(fit_wcls$alpha_hat) / fit_wcls$alpha_se_adjusted, 
                                 349 - 1 - 2, lower.tail = FALSE))), 3)
ctr_esmd <- rep(c("EMEE.mod","EMEE"), 8)
ctr_nomod <- cbind(ctrl_names, ctr_esmd, ctrl_est, ctrl_se, ctr_ci, ctr_pval)
colnames(ctr_nomod) <- c("Control Variables", "Estimator", "Estimate", "SE", "CI", "p-value")
mycaption <- "caption for simulation 1"
latex_code <- kable(ctr_nomod, format = "latex", booktabs = T, align = "c", caption = mycaption) %>%
  collapse_rows(columns = 1, latex_hline = "major")
print(latex_code)

# 2. excursion effect of moderation

# estimator, SE, 95% CI, p-value
t_quantile <- qt(0.975, 349 - 1 - 2) 
fit_wcls_newmod$alpha_hat 
fit_wcls_newmod$alpha_se_adjusted 
rbind(fit_wcls_newmod$alpha_hat - t_quantile * fit_wcls_newmod$alpha_se_adjusted,
      fit_wcls_newmod$alpha_hat + t_quantile * fit_wcls_newmod$alpha_se_adjusted)
2 * pt(abs(fit_wcls_newmod$alpha_hat) / fit_wcls_newmod$alpha_se_adjusted, 349 - 1 - 2, lower.tail = FALSE)

ctrl_names <- rep(names(fit_wcls_newmod$alpha_hat), each = 2)
ctrl_est <- round(c(rbind(fit_wcls_newmod$alpha_hat, fit_wclsmod$alpha_hat)), 3)
ctrl_se <- round(c(rbind(fit_wcls_newmod$alpha_se_adjusted, fit_wclsmod$alpha_se_adjusted)), 3)
# ctr_ci <- round(rbind(cbind(fit_wcls_new$alpha_hat - t_quantile * fit_wcls_new$alpha_se_adjusted,
#                       fit_wcls_new$alpha_hat + t_quantile * fit_wcls_new$alpha_se_adjusted), 
#                 cbind(fit_wcls$alpha_hat - t_quantile * fit_wcls$alpha_se_adjusted,
#                       fit_wcls$alpha_hat + t_quantile * fit_wcls$alpha_se_adjusted)), 3)

ctr_cilower <- round(c(rbind(fit_wcls_newmod$alpha_hat - t_quantile * fit_wcls_newmod$alpha_se_adjusted, 
                             fit_wclsmod$alpha_hat - t_quantile * fit_wclsmod$alpha_se_adjusted)),3)
ctr_ciupper <- round(c(rbind(fit_wcls_newmod$alpha_hat + t_quantile * fit_wcls_newmod$alpha_se_adjusted, 
                             fit_wclsmod$alpha_hat + t_quantile * fit_wclsmod$alpha_se_adjusted)),3)
ctr_ci <- cbind(ctr_cilower, ctr_ciupper)
ctr_ci <- apply(ctr_ci,1, function(x) paste("(", paste(x[1],x[2], sep = ","), ")", sep = ""))
ctr_pval <- round(c(rbind(2 * pt(abs(fit_wcls_newmod$alpha_hat) / fit_wcls_newmod$alpha_se_adjusted, 
                                 349 - 1 - 2, lower.tail = FALSE),
                          2 * pt(abs(fit_wclsmod$alpha_hat) / fit_wclsmod$alpha_se_adjusted, 
                                 349 - 1 - 2, lower.tail = FALSE))), 3)

ctr_esmd <- rep(c("EMEE.mod","EMEE"), 8)
ctr_mod <- cbind(ctrl_names, ctr_esmd, ctrl_est, ctrl_se, ctr_ci, ctr_pval)
colnames(ctr_mod) <- c("Control Variables", "Estimator", "Estimate", "SE", "CI", "p-value")
mycaption <- "caption for simulation 2"
latex_code <- kable(ctr_mod, format = "latex", booktabs = T, align = "c", caption = mycaption) %>%
  collapse_rows(columns = 1, latex_hline = "major")
print(latex_code)




