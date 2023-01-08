library(reshape)
library(kableExtra)
library(knitr)
library(ggplot2)
source("estimator_implementation/GEE_estimators.R")

#######
dgm_baserate <- function(sample_size, total_T, Delta, gamma, prob_a = 0.2) {
  
  beta_0 <- 0.1
  beta_1 <- 0.2
  
  df_names <- c("userid", "day", "A", "S", "S2","prob_A","prob_R_0", "prob_R","R","Y","k") 
  # prob_R_0 is the probability of R = 0 given A = 0.
  # prob_R is the probability of R = 0 given A = 1.
  
  dta <- data.frame(matrix(NA, nrow = sample_size * total_T, ncol = length(df_names)))
  names(dta) <- df_names
  
  dta$userid <- rep(1:sample_size, each = total_T)
  dta$day <- rep(1:total_T, times = sample_size)
  
  C = (gamma^(-0.5/Delta) + 1 + gamma^(0.5/Delta))
  prob_S_weight = c(gamma^(-0.5/Delta)/C, 1/C, gamma^(0.5/Delta)/C)
  E_S = 3*gamma^(1/Delta)/C
  
  for (t in 1:total_T) {
    # row index for the rows corresponding to day t for every subject
    row_index <- seq(from = t, by = total_T, length = sample_size)
    dta$S[row_index] <- sample(c(0,1,2), sample_size, prob = prob_S_weight, replace = TRUE)
    dta$S2[row_index] <- ifelse(dta$S[row_index] == 2, 1, 0) 
    dta$prob_A[row_index] <- rep(prob_a, sample_size)
    dta$A[row_index] <- rbinom(sample_size, 1, dta$prob_A[row_index])
    
    dta$prob_R_0[row_index] <- gamma^((1.5-0.5*dta$S[row_index])/Delta)
    dta$prob_R[row_index] <- ( 1- (1-dta$prob_R_0[row_index]*E_S^(Delta-1)) * 
                                 exp(beta_0+beta_1*dta$S[row_index]) ) / E_S^(Delta-1)
    dta$R[row_index] <- rbinom(sample_size, 1, 
                               ifelse(dta$A[row_index] == 0, 1 - dta$prob_R_0[row_index], 1 - dta$prob_R[row_index]))
  }
  
  # after getting R_t's, calculate Y_t's respectively
  for (t in 1:total_T) {
    row_index <- seq(from = t, by = total_T, length = sample_size)
    for (j in 1:length(row_index)) {
      # notice here that when t+ Delta - 1 > total_T, we take the product up to total_T.
      dta$Y[row_index[j]] <- 1 - prod(1-dta$R[row_index[j]:min((row_index[j]+Delta-1),j*total_T)])
      
      # get the index of R's when Y = 1 for the first time, 
      # if Y = 0, assign k = Delta-1, since all A's must be 0 
      dta$k[row_index[j]] <- ifelse(dta$Y[row_index[j]] == 0, (Delta-1),
                                    min(which(dta$R[row_index[j]:min((row_index[j]+Delta-1),j*total_T)] == 1))-1)
      
    }
  }
  return(dta)
}

beta_true_marginal_baserate <- function(Delta, gamma){
  beta_0 <- 0.1
  beta_1 <- 0.2
  
  C = (gamma^(-0.5/Delta) + 1 + gamma^(0.5/Delta))
  prob_S_weight = c(gamma^(-0.5/Delta)/C, 1/C, gamma^(0.5/Delta)/C)
  E_S = 3*gamma^(1/Delta)/C
  prob_R_0 <-  c(gamma^((1.5-0.5*0)/Delta),gamma^((1.5-0.5*1)/Delta),gamma^((1.5-0.5*2)/Delta))
  exp_h <- c(exp(beta_0+beta_1*0),exp(beta_0+beta_1*1),exp(beta_0+beta_1*2))
  
  numerator <- sum( ((1- prob_R_0*E_S^(Delta-1)) * exp_h) * prob_S_weight)
  denominator <-  sum( ((1- prob_R_0*E_S^(Delta-1))) * prob_S_weight)
  beta_true_marginal <- log(numerator / denominator)
  return(beta_true_marginal)
}

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

data_generating_process <- dgm_baserate

library(tidyverse)
library(foreach)
library(doMC)
library(doRNG)

max_cores <- 16
registerDoMC(min(detectCores() - 1, max_cores))
sample_sizes <- c(30, 50, 100)
nsim <- 10
Delta <- 3
total_T <- 100 

control_vars <- "S"
moderator_vars <- c()

set.seed(1)

efficiency_table = data.frame()
result_table = data.frame()

for (gamma in (4:9)/10) {
  beta_true_marginal <- beta_true_marginal_baserate(Delta, gamma)
  result_df_collected <- data.frame()
  for (i_ss in 1:length(sample_sizes)) {
    
    sample_size <- sample_sizes[i_ss]
    result <- foreach(isim = 1:nsim, .combine = "c") %dorng% {
      if (isim %% 10 == 0) {
        cat(paste("Starting iteration",isim,"\n"))
      }
      # set probability of accepting treatment 0.2 
      dta <- data_generating_process(sample_size, total_T, Delta = Delta, gamma = gamma, prob_a = 0.2)
      
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
  efficiency_table = rbind(efficiency_table, c(gamma,relative_efficiency))
  result_df_collected$gamma = gamma
  result_table = rbind(result_table, result_df_collected)
}

saveRDS(efficiency_table, file = "efficiency_table_baserate.RDS")
saveRDS(result_table, file = "result_table_baserate.RDS")

# base rate
base_rate_perDeltagamma <- function(Delta, gamma) {
  C = (gamma^(-0.5/Delta) + 1 + gamma^(0.5/Delta))
  prob_S_weight = c(gamma^(-0.5/Delta)/C, 1/C, gamma^(0.5/Delta)/C)
  
  base_rate_vals = c(1 - gamma^(1 +0.5/Delta)* (3/C)^(Delta - 1),
                     1 - gamma^(1 +0  /Delta)* (3/C)^(Delta - 1),
                     1 - gamma^(1 -0.5/Delta)* (3/C)^(Delta - 1))
  e_base_rate = sum(prob_S_weight * base_rate_vals)
  return(e_base_rate)
}

base_rate_perDeltagamma(3, 0.5)

df_delta <- as.data.frame(rbind(cbind("Relative Efficiency", base_rates, efficiency_df_delta$efficiency_means),
                                cbind("Variance of pd-EMEE", base_rates, efficiency_df_delta$mod_var*200),
                                cbind("Variance of EMEE", base_rates, efficiency_df_delta$ori_var*200)))
colnames(df_delta) <- c("grps", "Base Rate","Relative Efficiency")
df_delta$`Base Rate` <- as.numeric(df_delta$`Base Rate`)
df_delta$`Relative Efficiency` <- as.numeric(df_delta$`Relative Efficiency`)

p1 <- ggplot(data=df_delta, aes(x=`Base Rate`, y=`Relative Efficiency`, group=grps)) +
  geom_line(aes(linetype=grps))+
  geom_point(aes(size=grps))+
  scale_size_manual(values=c(1, 0, 0))+
  #labs(linetype="",x="Delta",y="Relative Efficiency") + 
  scale_y_continuous(sec.axis = sec_axis(~./200, name = "Variance of Estimators"))+
  theme(plot.title = element_text(size = 8, hjust = 0.5),
        legend.title = element_blank(),
        legend.position="bottom",
        legend.text = element_text(size = 6), 
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 8))+ 
  ggtitle("Relative Efficiency of pd-EMEE to EMEE over Different Values of Base Rates")

# result_df_delta <- readRDS("result_table_delta(with sd).RDS")
# efficiency_df_delta <- readRDS("efficiency_table_delta(with sd).RDS")
# 
