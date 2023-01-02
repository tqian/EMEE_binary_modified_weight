rm(list = ls())

source("estimator_implementation/GEE_estimators.R")

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

beta_true_marginalK <- function(Delta, K){
  beta_0 <- 0.1
  beta_1 <- 0.2
  pa <- 0.2
  
  C <- (0.5^(0.5/Delta) + 1 + 0.5^(-0.5/Delta))
  prob_S_weight <- c(0.5^(-0.5/Delta)/C, 1/C, 0.5^(0.5/Delta)/C)
  E_S = 3*0.5^(1/Delta)/C
  prob_R_0 <-  c(0.5^((1.5-0.5*0)/Delta),0.5^((1.5-0.5*1)/Delta),0.5^((1.5-0.5*2)/Delta))
  exp_h <- c(exp(beta_0+beta_1*0),exp(beta_0+beta_1*1),exp(beta_0+beta_1*2))
  prob_R_1 <- (1 - (1 - prob_R_0 * E_S^(Delta - 1)) * exp_h)/E_S^(Delta - 1)
    
  Q <- (C - sum(exp_h * (C * prob_S_weight)) + 0.5^(1/Delta) * E_S^(Delta - 1) * sum(exp_h))/
    (E_S^(Delta - 1) * C)
  Q_weighted <- (3/ C * 0.5^(1/ Delta)) * (1 - pa) + Q * pa
  
  numerator <- sum( ((1- prob_R_1*E_S^K * Q_weighted^(Delta-K-1))) * prob_S_weight)
  denominator <-  sum(((1- prob_R_0*E_S^K * Q_weighted^(Delta-K-1))) * prob_S_weight)
  
  beta_true_marginal <- log(numerator / denominator)
  return(beta_true_marginal)
}

source("simulations/dgm_simulation.R")
source("estimator_implementation/WCLS_modified.R")
source("simulations/simulation_forK.R")

data_generating_process <- dgm_binary_categorical_covariate_new

library(rootSolve)
library(tidyverse)
library(foreach)
library(doMC)
library(doRNG)

max_cores <- 16
registerDoMC(min(detectCores() - 1, max_cores))
sample_sizes <- c(30, 50, 100)
nsim <- 20 #1000

total_T <- 100 

control_vars <- "S"
moderator_vars <- c()

set.seed(1)

efficiency_table = data.frame()
result_table = data.frame()
Delta = 10

for (Kref in 0:(Delta-1)) {
  #beta_true_marginal <- beta_true_marginal_generalDelta(Delta)
  beta_true_marginal <- beta_true_marginalK(Delta, Kref)
  result_df_collected <- data.frame()
  
  for (i_ss in 1:length(sample_sizes)) {
    # progbar = txtProgressBar(min=1, max=nsim, style=3)
    sample_size <- sample_sizes[i_ss]
    result <- foreach(isim = 1:nsim, .combine = "c") %dorng% {
      if (isim %% 10 == 0) {
        cat(paste("Starting iteration", Kref, i_ss, isim, "\n"), file="log.txt", append=TRUE)
      }
      
      # set probability of accepting treatment 0.2 
      dta <- data_generating_process(sample_size, total_T, Delta = Delta, prob_a = 0.2)
      
      # new EMEE estimator
      fit_wcls_new <- weighted_centered_least_square_forK_new(
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
        Kref = Kref
      )
      
      # old EMEE estimator
      fit_wcls <- weighted_centered_least_square_forK(
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
        Kref = Kref
      )
      
      output <- list(list(fit_wcls_new = fit_wcls_new, fit_wcls = fit_wcls))
      # setTxtProgressBar(progbar, isim)
    }
    # close(progbar)
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
  efficiency_table = rbind(efficiency_table, c(Kref,relative_efficiency))
  result_df_collected$Kref = Kref
  result_table = rbind(result_table, result_df_collected)

}

saveRDS(efficiency_table, file = "efficiency_table_Kref(true beta).RDS")
saveRDS(result_table, file = "result_table_delta_Kref(true beta).RDS")

efficiency_table <- readRDS("efficiency_table_Kref.RDS.RDS")

efficiency_means <- rowMeans(efficiency_table[,2:4])
efficiency_table$efficiency_means <- efficiency_means
colnames(efficiency_table) <- c("Delta", "SS = 30", "SS = 50", "SS = 100", "efficiency_means")
ggplot(efficiency_table) +
  geom_line(aes(x = Delta, y = efficiency_means)) +
  geom_point(aes(x = Delta, y = efficiency_means)) +
  theme_bw() +
  ylab("Relative Efficiency") +
  xlab("K") +
  theme(legend.position="bottom") 
