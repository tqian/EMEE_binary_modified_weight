# code to generate plot for efficiency vs various Delta in the paper

# Tianchen Qian
# 2018.08.12

# simulation part is copied from simulation.R
# table creation part is copied from U-stat paper

# update on 2019.02.06: to include eif_modified_weight

# update on 2019.03.29: to include GEE with exchangeable correlation structure

# Yihan Bao
# update on 2021.07.31: to include modified-EMEE

# update on 2021.09.15: implement visualization

##### simulation part #####

rm(list = ls())

source("GEE_estimators.R")

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

source("simulations/dgm_simulation.R")
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

control_vars <- "S"
moderator_vars <- c()

set.seed(1)

efficiency_table = data.frame()
result_table = data.frame()

for (Delta in 1:10) {
    beta_true_marginal <- beta_true_marginal_generalDelta(Delta)
    result_df_collected <- data.frame()
    for (i_ss in 1:length(sample_sizes)) {
        
        sample_size <- sample_sizes[i_ss]
        result <- foreach(isim = 1:nsim, .combine = "c") %dorng% {
            if (isim %% 10 == 0) {
                cat(paste("Starting iteration",isim,"\n"))
            }
            # set probability of accepting treatment 0.2 
            dta <- data_generating_process(sample_size, total_T, Delta = Delta, prob_a = 0.2)
            
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
            
            output <- list(list(fit_wcls_new = fit_wcls_new, fit_wcls = fit_wcls))
        }
        #sink()
        
        ee_names <- c("modified-EMEE", "EMEE")
        alpha_names <- c("Intercept", control_vars)
        beta_names <- c("Intercept", moderator_vars)
        num_estimator <- length(ee_names)
        
        alpha <- simplify2array(lapply(result, function(l) matrix(c(l$fit_wcls_new$alpha_hat, l$fit_wcls$alpha_hat),
                                                                  nrow = length(ee_names), byrow = TRUE, dimnames = list(ee_names, alpha_names))))
        alpha_se <- simplify2array(lapply(result, function(l) matrix(c(l$fit_wcls_new$alpha_se, l$fit_wcls$alpha_se),
                                                                     nrow = length(ee_names), byrow = TRUE, dimnames = list(ee_names, alpha_names))))
        alpha_se_adjusted <- simplify2array(lapply(result, function(l) matrix(c(l$fit_wcls_new$alpha_se_adjusted, l$fit_wcls$alpha_se_adjusted),
                                                                              nrow = length(ee_names), byrow = TRUE, dimnames = list(ee_names, alpha_names))))
        beta <- simplify2array(lapply(result, function(l) matrix(c(l$fit_wcls_new$beta_hat, l$fit_wcls$beta_hat),
                                                                 nrow = length(ee_names), byrow = TRUE, dimnames = list(ee_names, beta_names))))
        beta_se <- simplify2array(lapply(result, function(l) matrix(c(l$fit_wcls_new$beta_se, l$fit_wcls$beta_se),
                                                                    nrow = length(ee_names), byrow = TRUE, dimnames = list(ee_names, beta_names))))
        beta_se_adjusted <- simplify2array(lapply(result, function(l) matrix(c(l$fit_wcls_new$beta_se_adjusted, l$fit_wcls$beta_se_adjusted),
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
        relative_efficiency = (result_df_collected$sd[c(2,4,6)]/result_df_collected$sd[c(1,3,5)])^2
    }
    efficiency_table = rbind(efficiency_table, c(Delta,relative_efficiency))
    result_df_collected$delta = Delta
    result_table = rbind(result_table, result_df_collected)
}

#saveRDS(efficiency_table, file = "efficiency_table_delta.RDS")
saveRDS(efficiency_table, file = "efficiency_table_delta(with sd).RDS")
saveRDS(result_table, file = "result_table_delta(with sd).RDS")

##### create tables for paper #####

efficiency_table <- readRDS("efficiency_table_delta(with sd).RDS")

efficiency_means <- rowMeans(efficiency_table[,2:4])
efficiency_table$efficiency_means <- efficiency_means
colnames(efficiency_table) <- c("Delta", "SS = 30", "SS = 50", "SS = 100", "efficiency_means")
ggplot(efficiency_table) +
    geom_line(aes(x = Delta, y = efficiency_means)) +
    geom_point(aes(x = Delta, y = efficiency_means)) +
    theme_bw() +
    ylab("Relative Efficiency") +
    xlab("Delta") +
    theme(legend.position="bottom") 
