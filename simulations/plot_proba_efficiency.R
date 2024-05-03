rm(list = ls())

compute_result_beta <- function(beta_true, beta, beta_se, beta_se_adjusted, moderator_vars, control_vars, significance_level,
                                na.rm = FALSE) {
  
  beta_true_array <- array(NA, dim = dim(beta), dimnames = dimnames(beta))
  for (ind1 in 1:dim(beta_true_array)[1]) {
    for (ind3 in 1:dim(beta_true_array)[3]) {
      beta_true_array[ind1, , ind3] <- beta_true
    }
  }
  
  p <- length(moderator_vars) + 1
  q <- length(control_vars) + 1
  
  bias <- apply(beta - beta_true_array, c(1,2), mean, na.rm = na.rm)
  sd <- apply(beta, c(1,2), sd, na.rm = na.rm)
  rmse <- apply(beta - beta_true_array, c(1,2), function(v) sqrt(mean(v^2, na.rm = na.rm)))
  
  critical_factor <- qnorm(1 - significance_level/2)
  ci_left <- beta - critical_factor * beta_se
  ci_right <- beta + critical_factor * beta_se
  coverage_prob <- apply((ci_left < beta_true_array) & (ci_right > beta_true_array),
                         c(1,2), mean, na.rm = na.rm)
  
  critical_factor_adj <- qt(1 - significance_level/2, df = sample_size - 1 - q)
  ci_left_adj <- beta - critical_factor_adj * beta_se_adjusted
  ci_right_adj <- beta + critical_factor_adj * beta_se_adjusted
  coverage_prob_adj <- apply((ci_left_adj < beta_true_array) & (ci_right_adj > beta_true_array),
                             c(1,2), mean, na.rm = na.rm)
  
  return(list(bias = bias, sd = sd, rmse = rmse, coverage_prob = coverage_prob, coverage_prob_adjusted = coverage_prob_adj))
}

source("wrapper_implementation/dgm_simulation.R")
source("wrapper_implementation/helper_functions.R")
source("wrapper_implementation/EMEE_ImprovedEffByProjection_measurableWeightNotTakenOut.R")
source("estimator_implementation/WCLS_modified.R")
source("estimator_implementation/WCLS_original.R")

data_generating_process <- dgm_binary_categorical_covariate_new

library(rootSolve)
library(tidyverse)
library(foreach)
library(doMC)
library(doRNG)

max_cores <- 16
registerDoMC(min(detectCores() - 1, max_cores))
sample_sizes <- c(30, 50, 100)
nsim <- 1000

total_T <- 100 
Delta <- 3 

control_vars <- "S"
moderator_vars <- c()

set.seed(1)

pdEMEE_efficiency_proba = data.frame()
pdEMEE2_efficiency_proba = data.frame()
result_table = data.frame()

beta_true_marginal <- beta_true_marginal_generalDelta(Delta)
for (prob_a in (1:8)/10) {
  result_df_collected <- data.frame()
  for (i_ss in 1:length(sample_sizes)) {
    
    sample_size <- sample_sizes[i_ss]
    result <- foreach(isim = 1:nsim, .combine = "c") %dorng% {
      if (isim %% 10 == 0) {
        cat(paste("Starting iteration",isim,"\n"))
      }
 
      dta <- data_generating_process(sample_size, total_T, Delta = Delta, prob_a = prob_a)
      
      try({
        # improved EMEE estimator
        fit_wcls_improved <- EMEE_ImprovedEffByProjection_measurableWeightNotTakenOut_wrapper(
          dta = dta,
          id_varname = "userid",
          decision_time_varname = "day",
          treatment_varname = "A",
          outcome_varname = "Y",
          control_varname = control_vars,
          moderator_varname = moderator_vars,
          rand_prob_varname = "prob_A",
          rand_prob_tilde_varname = NULL,
          rand_prob_tilde = 0.2,
          estimator_initial_value = NULL,
          Delta = Delta,
          link_function_for_nuisance_prediction = "poisson"
        )
        
        # new EMEE estimator
        fit_wcls_new <- weighted_centered_least_square_withDelta_new(
          dta = dta,
          id_varname = "userid",
          decision_time_varname = "day",
          treatment_varname = "A",
          outcome_varname = "Y",
          control_varname = control_vars,
          moderator_varname = moderator_vars,
          rand_prob_varname = "prob_A",
          rand_prob_tilde_varname = NULL,
          rand_prob_tilde = 0.2,
          estimator_initial_value = NULL,
          Delta = Delta
        )
        
        # old EMEE estimator
        fit_wcls <- weighted_centered_least_square_withDelta(
          dta = dta,
          id_varname = "userid",
          decision_time_varname = "day",
          treatment_varname = "A",
          outcome_varname = "Y",
          control_varname = control_vars,
          moderator_varname = moderator_vars,
          rand_prob_varname = "prob_A",
          rand_prob_tilde_varname = NULL,
          rand_prob_tilde = 0.2,
          estimator_initial_value = NULL,
          Delta = Delta
        )
        
        output <- list(list(fit_wcls_improved = fit_wcls_improved, fit_wcls_new = fit_wcls_new, fit_wcls = fit_wcls))
      }, silent = TRUE)
    }
    #sink()
    
    ee_names <- c("improved-EMEE", "modified-EMEE", "EMEE")
    beta_names <- c("Intercept", moderator_vars)
    num_estimator <- length(ee_names)
    
    beta <- simplify2array(lapply(result, function(l) matrix(c(l$fit_wcls_improved$beta_hat, l$fit_wcls_new$beta_hat, l$fit_wcls$beta_hat),
                                                             nrow = length(ee_names), byrow = TRUE, dimnames = list(ee_names, beta_names))))
    beta_se <- simplify2array(lapply(result, function(l) matrix(c(l$fit_wcls_improved$beta_se, l$fit_wcls_new$beta_se, l$fit_wcls$beta_se),
                                                                nrow = length(ee_names), byrow = TRUE, dimnames = list(ee_names, beta_names))))
    beta_se_adjusted <- simplify2array(lapply(result, function(l) matrix(c(l$fit_wcls_improved$beta_se_adjusted, l$fit_wcls_new$beta_se_adjusted, l$fit_wcls$beta_se_adjusted),
                                                                         nrow = length(ee_names), byrow = TRUE, dimnames = list(ee_names, beta_names))))
    
    result <- compute_result_beta(beta_true_marginal, beta, beta_se, beta_se_adjusted, moderator_vars, control_vars, significance_level = 0.05)
    result_df <- data.frame(ss = rep(sample_size, num_estimator),
                            est = ee_names,
                            bias = result$bias,
                            sd = result$sd,
                            rmse = result$rmse,
                            cp.unadj = result$coverage_prob,
                            cp.adj = result$coverage_prob_adjusted)
    names(result_df) <- c("ss", "est", "bias", "sd", "rmse", "cp.unadj", "cp.adj")
    rownames(result_df) <- NULL
    
    result_df_collected <- rbind(result_df_collected, result_df)
    relative_efficiency_EMEE = (result_df_collected$sd[c(3,6,9)]/result_df_collected$sd[c(2,5,8)])^2
    relative_efficiency_EMEE2 = (result_df_collected$sd[c(3,6,9)]/result_df_collected$sd[c(1,4,7)])^2
  }
  pdEMEE_efficiency_proba = rbind(pdEMEE_efficiency_proba, c(prob_a, relative_efficiency_EMEE))
  pdEMEE2_efficiency_proba = rbind(pdEMEE2_efficiency_proba, c(prob_a, relative_efficiency_EMEE2))
  result_df_collected$prob_a = prob_a
  result_table = rbind(result_table, result_df_collected)
}

saveRDS(pdEMEE_efficiency_proba, file = "pdEMEE_efficiency_proba.RDS")
saveRDS(pdEMEE2_efficiency_proba, file = "pdEMEE2_efficiency_proba.RDS")
saveRDS(result_table, file = "result_table_proba.RDS")

