# Yihan Bao, 2021.05.19
# Estimating Equation function for modified-EMEE
# Based on the estimating equation function for EMEE (in Qian et al. 2021)

# dta = dta
# id_varname = "userid"
# decision_time_varname = "day"
# treatment_varname = "A"
# outcome_varname = "Y"
# control_varname = control_vars
# moderator_varname = moderator_vars
# rand_prob_varname = "prob_A"
# rand_prob_tilde_varname = NULL
# rand_prob_tilde = 0.2
# estimator_initial_value = NULL
# Delta = 3
# avail_varname = NULL
# rand_prob_tilde_varname = NULL
# rand_prob_tilde = NULL
# estimator_initial_value = NULL

# Assume we already have k
# 1. estimate $E(Y_{it, \Delta}W_{u+1:t+\Delta-1} |H_{iu}, A_{iu}, R_{i,t+1}, ..., R_u = 0 )$
# for any u = 1, ... delta-1, if A_{iu} = 1, run a regression
# for any u in 1: (Delta-1)

weighted_centered_least_square_withDelta_improved <- function(
  dta,
  id_varname,
  decision_time_varname,
  treatment_varname,
  outcome_varname,
  control_varname,
  moderator_varname,
  rand_prob_varname,
  avail_varname = NULL,
  rand_prob_tilde_varname = NULL, # \tilde{p}_t(1|H_t) in WCLS (variable name in the data set)
  rand_prob_tilde = NULL,         # \tilde{p}_t(1|H_t) in WCLS (numeric number or vector)
  estimator_initial_value = NULL,
  Delta
)
{
  
  sample_size <- length(unique(dta[, id_varname]))
  total_person_decisionpoint <- nrow(dta)
  id_names <- unique(dta[,id_varname])
  
  if (is.null(avail_varname)) {
    avail <- rep(1, total_person_decisionpoint)
  } else {
    avail <- dta[, avail_varname]
  }
  
  A <- dta[, treatment_varname]
  # checking for NA in treatment indicator    
  if (any(is.na(A[avail == 1]))) {
    stop("Treatment indicator is NA where availability = 1.")
  }
  A[avail == 0] <- 0
  
  p_t <- dta[, rand_prob_varname]
  cA <- A - p_t # centered A
  Y <- dta[, outcome_varname]
  
  # X (moderator) design matrix, intercept added
  Xdm <- as.matrix( cbind( rep(1, nrow(dta)), dta[, moderator_varname] ) )
  # Z (control) design matrix, intercept added
  Zdm <- as.matrix( cbind( rep(1, nrow(dta)), dta[, control_varname] ) ) 
  
  if (is.null(rand_prob_tilde_varname) & is.null(rand_prob_tilde)) {
    p_t_tilde <- rep(0.5, nrow(dta))
  } else if (is.null(rand_prob_tilde_varname)) {
    if (length(rand_prob_tilde) == 1) {
      p_t_tilde <- rep(rand_prob_tilde, total_person_decisionpoint)
    } else if (length(rand_prob_tilde) == total_person_decisionpoint) {
      p_t_tilde <- rand_prob_tilde
    } else {
      stop("rand_prob_tilde is of incorrect length.")
    }
  } else {
    p_t_tilde <- dta[, rand_prob_tilde_varname]
  }
  cA_tilde <- A - p_t_tilde
  
  WCLS_weight_inverseProb <- c()
  if (Delta == 1){
    WCLS_weight_inverseProb = rep(1, total_person_decisionpoint)
  }else{
    for (i in 1:sample_size){
      dta_perid <- dta[dta[,id_varname] == id_names[i],]
      dta_timepoint <- length(dta_perid[,decision_time_varname])
      
      #append with A_{T+1} ... A_{T+\delta-1} = 0
      inverseProb_numerator <- 1 - c(dta_perid[,treatment_varname],rep(0,Delta-1)) 
      #append with p_{T+1} ... p_{T+\delta-1} = 0
      inverseProb_denominator <- 1- c(dta_perid[,rand_prob_varname],rep(0,Delta-1)) 
      
      for (j in 1:dta_timepoint) {
        # the first time R_{t+k} = 1
        first_R_1 <- dta_perid$k[j] 
        
        if (first_R_1 == 0){
          # if R_t = 1 right after A_t, set the weight to 1
          WCLS_weight_inverseProb <- c(WCLS_weight_inverseProb,1) 
          
        }else{ 
          # in other words, if A = 0 appears before k, we need to set it to 0, o.w keep this term
          WCLS_weight_inverseProb <- c(WCLS_weight_inverseProb,
                                       prod(inverseProb_numerator[(j+1):(j+first_R_1)]/
                                              inverseProb_denominator[(j+1):(j+first_R_1)]))
        }
      }
    }
  }
  
  p <- length(moderator_varname) + 1 # dimension of beta

  
  #######
  residuals_per_u <- function(u){
    WCLS_weight_till_u <- c()
    WCLS_weight_start_u <- c()
    A_u <- c()
    S_u <- c()
    if (Delta == 1){
      WCLS_weight_start_u = rep(1, total_person_decisionpoint)
      WCLS_weight_till_u = rep(1, total_person_decisionpoint)
      A_u <- A
      S_u <- dta$S
    }else{
      for (i in 1:sample_size){
        dta_perid <- dta[dta[,id_varname] == id_names[i],]
        dta_timepoint <- length(dta_perid[,decision_time_varname])
        
        #append with A_{T+1} ... A_{T+\delta-1} = 0
        inverseProb_numerator <- 1 - c(dta_perid[,treatment_varname],rep(0,Delta-1)) 
        #append with p_{T+1} ... p_{T+\delta-1} = 0
        inverseProb_denominator <- 1- c(dta_perid[,rand_prob_varname],rep(0,Delta-1))  
        
        for (j in 1:dta_timepoint) {
          A_u <- c(A_u, ifelse(j+u <= dta_timepoint, dta_perid$A[j+u], 0))
          S_u <- c(S_u, ifelse(j+u <= dta_timepoint, dta_perid$S[j+u], sample(c(0,1,2))))
          
        
          WCLS_weight_till_u <- c(WCLS_weight_till_u,
                                  prod(inverseProb_numerator[(j+1):(j+min(u, dta_perid$k[j]))]/
                                         inverseProb_denominator[(j+1):(j+min(u, dta_perid$k[j]))]))
          if (dta_perid$k[j] < u + 1){
            WCLS_weight_start_u <- c(WCLS_weight_start_u, 1)
          }else{
            WCLS_weight_start_u <- c(WCLS_weight_start_u,
                                     prod(inverseProb_numerator[(j+u+1):(j+ dta_perid$k[j])]/
                                            inverseProb_denominator[(j+u+1):(j+ dta_perid$k[j])]))}
        }
      }
    }
    dta$till_u <- WCLS_weight_till_u
    dta$start_u <- WCLS_weight_start_u
    dta$A_u <- A_u
    dta$S_u <- S_u
   
    Y_reg <- dta$Y * WCLS_weight_start_u
    reg <- lm(Y_reg ~ A_u + S_u, data = dta)
    Y_u_pred <- coef(reg)[1] + coef(reg)[2]*A_u + coef(reg)[3]*S_u
    reg0 <- lm(Y_reg ~ S_u, data = dta)
    Y_u_pred_mean <- coef(reg0)[1] + coef(reg0)[2] * S_u
    Y_u_residual <- WCLS_weight_till_u*(Y_u_pred - Y_u_pred_mean)
    
    return(Y_u_residual)
  }
  
  M <- ifelse(A, p_t_tilde / p_t, (1 - p_t_tilde) / (1 - p_t)) 
  dta$M <- M
  
  # 2. estimate E(Y_{it, \Delta} W_it|H_t, A_t)
  # regress on H_t
  dta$WCLS <- WCLS_weight_inverseProb
  Y_reg <- dta$Y * dta$WCLS
  
  # regress on moderator
  # if(is.null(moderator_varname)){
  #   reg <- lm(Y_reg ~ A, data = dta)
  #   Yt_pred_1 <- coef(reg)[1] + coef(reg)[2]*1 
  #   Yt_pred_0 <- coef(reg)[1] + coef(reg)[2]*0 
  #   Yt_pred <- coef(reg)[1] + coef(reg)[2]*dta$A 
  # }else{
  reg <- lm(Y_reg ~ dta$A + dta$S, data = dta)
  Yt_pred_1 <- coef(reg)[1] + coef(reg)[2]*1 + coef(reg)[3]*dta$S
  Yt_pred_0 <- coef(reg)[1] + coef(reg)[2]*0 + coef(reg)[3]*dta$S
  #Y_t_pred <- Yt_pred_1 * mean(A) +  Yt_pred_0 * (1 - mean(A)) 
  Yt_pred <- coef(reg)[1] + coef(reg)[2]*dta$A + coef(reg)[3]*dta$S
  #}
  
  # for (a in 0:1) {
  #   dta_a = dta[dta$A == a, ]
  #   Y_reg <- dta_a$Y * dta_a$WCLS
  #   # regress on moderator
  #   reg <- lm(Y_reg ~ dta_a$S, data = dta_a)
  #   Y_pred <- coef(reg)[1] + coef(reg)[2]*dta_a$S
  # }
  
  # Now calculating the estimating equation
  # Notice here we ignore the control variables
  ### 2. estimation ###
  
  estimating_equation <- function(beta) {
    beta <- as.matrix(beta)
    exp_AXdm_beta <- exp(A * (Xdm %*% beta))
    weight <- exp_AXdm_beta^(-1)
    residual_u <- c()
    for (u in 1:(Delta-1)) {
      residual_u <- cbind(residual_u, residuals_per_u(u))
    }
    residual_u_sum <- apply(residual_u, 1, sum)  
    
    ef <- rep(NA, length(beta)) # value of estimating function
    for (i in 1:p) {
      # Xdm is S_t, the moderator variables with intercept included
      error_term1 <- avail * weight * M * cA_tilde * Xdm[,i] * 
        (Y * WCLS_weight_inverseProb - residual_u_sum - Yt_pred)
      error_term2 <- avail * Xdm[,i] * p_t_tilde * (1-p_t_tilde) * 
        (exp(-1 *Xdm %*% beta) * Yt_pred_1 - Yt_pred_0)
      ef[i] <- sum(error_term1 + error_term2)
    }
    
    ef <- ef / sample_size
    return(ef)
  }
  
  if (is.null(estimator_initial_value)) {
    estimator_initial_value <- rep(0, length = p)
  }
  
  solution <- tryCatch(
    {
      multiroot(estimating_equation, estimator_initial_value, maxiter = 200)
    },
    error = function(cond) {
      message("\nCatched error in multiroot inside weighted_centered_least_square():")
      message(cond)
      return(list(root = rep(NaN, p), msg = cond,
                  f.root = rep(NaN, p)))
    })
  
  beta_hat <- as.vector(solution$root)
  
  
  ### 3. asymptotic variance ###
  ### 3.1 Compute M_n matrix (M_n is the empirical expectation of the derivative of the estimating function) ###
  Mn_summand <- array(NA, dim = c(total_person_decisionpoint, p, p))
  
  beta_hat <- as.matrix(beta_hat)
  exp_AXdm_beta_hat <- exp(A * (Xdm %*% beta_hat))
  weight_hat <- exp_AXdm_beta_hat^(-1)
  M2_multiplier <- exp(-1 *Xdm %*% beta_hat) 
  
  
  residual_u <- c()
  for (u in 1:(Delta-1)) {
    residual_u <- cbind(residual_u, residuals_per_u(u))
  }
  residual_u_sum <- apply(residual_u, 1, sum)  
  M1_error <- Y * WCLS_weight_inverseProb - residual_u_sum - Yt_pred
  
  M1 <- matrix(NA, nrow = p, ncol = total_person_decisionpoint)
  M2 <- matrix(NA, nrow = p, ncol = total_person_decisionpoint)
  M3 <- matrix(NA, nrow = p, ncol = total_person_decisionpoint)
  
  for (it in 1:total_person_decisionpoint) {
    # this is to make R code consistent whether X_it, Z_it contains more entries or is just the intercept.        
    if (p == 1) {
      Xbeta <- Xdm[it, ] * beta_hat
    } else {
      Xbeta <- as.numeric(Xdm[it, ] %*% beta_hat)
    }
    
    # derivatives
    M1_deriv <- avail[it] * weight_hat[it] * M[it] * (A[it] - p_t_tilde[it]) *
      M1_error[it] * as.matrix(Xdm[it, ]) %*% (-A[it] * t(Xdm[it, ]))
    M2_deriv <- avail[it] * p_t_tilde[it] * (1 - p_t_tilde[it]) * M2_multiplier[it] *
      Yt_pred_1[it] * as.matrix(Xdm[it, ])  %*% (-1 * t(Xdm[it, ]))
    M3_deriv <- 0
    Mn_summand[it, , ] <- M1_deriv + M2_deriv + M3_deriv
    
    M1[,it] <- avail[it] * weight_hat[it] *M[it] * (A[it] - p_t_tilde[it]) *
      M1_error[it] * as.matrix(Xdm[it, ])
    M2[,it] <- avail[it] * p_t_tilde[it] * (1 - p_t_tilde[it]) * M2_multiplier[it] *
      Yt_pred_1[it] * as.matrix(Xdm[it, ]) 
    M3[,it] <- -1 * avail[it] * p_t_tilde[it] * (1 - p_t_tilde[it]) * Yt_pred_0[it] * 
      as.matrix(Xdm[it, ]) 
    
  }
  Mn <- apply(Mn_summand, c(2,3), sum) / sample_size
  Mn_inv <- solve(Mn)
  
  ### 3.2 Compute \Sigma_n matrix (\Sigma_n is the empirical variance of the estimating function) ###
  
  Sigman_summand <- array(NA, dim = c(sample_size, p, p))
  # Sigman_summand is  \sum_{t=1}^T ( D^{(t),T} r^(t) )^{\otimes 2}
  # See note 2018.08.06 about small sample correction
  
  person_first_index <- c(find_change_location(dta[, id_varname]), total_person_decisionpoint + 1)
  M_matrix <- M1 + M2 + M3
  
  for (i in 1:sample_size) {
    if(p == 1){
      M_matrix_i <- t(M_matrix[, person_first_index[i] : (person_first_index[i+1] - 1)])
    }else{
      M_matrix_i <- M_matrix[, person_first_index[i] : (person_first_index[i+1] - 1)]
    }
    Sigman_summand[i, , ] <- M_matrix_i %*% t(M_matrix_i)
  }
  Sigman <- apply(Sigman_summand, c(2,3), sum) / sample_size
  
  varcov <- Mn_inv %*% Sigman %*% t(Mn_inv) / sample_size
  beta_se <- sqrt(diag(varcov)[1:p])

  return(list(beta_hat = beta_hat, beta_se = beta_se))
}


compute_result_beta_woadjusted <- function(beta_true, beta, beta_se, moderator_vars, control_vars, significance_level,
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
  
  # critical_factor_adj <- qt(1 - significance_level/2, df = sample_size - 1 - q)
  # ci_left_adj <- beta - critical_factor_adj * beta_se_adjusted
  # ci_right_adj <- beta + critical_factor_adj * beta_se_adjusted
  # coverage_prob_adj <- apply((ci_left_adj < beta_true_array) & (ci_right_adj > beta_true_array),
  #                            c(1,2), mean, na.rm = na.rm)
  
  return(list(bias = bias, sd = sd, rmse = rmse, coverage_prob = coverage_prob))
}





weighted_centered_least_square_withDelta_improved_old <- function(
  dta,
  id_varname,
  decision_time_varname,
  treatment_varname,
  outcome_varname,
  control_varname,
  moderator_varname,
  rand_prob_varname,
  avail_varname = NULL,
  rand_prob_tilde_varname = NULL, # \tilde{p}_t(1|H_t) in WCLS (variable name in the data set)
  rand_prob_tilde = NULL,         # \tilde{p}_t(1|H_t) in WCLS (numeric number or vector)
  estimator_initial_value = NULL,
  Delta
)
{
  
  sample_size <- length(unique(dta[, id_varname]))
  total_person_decisionpoint <- nrow(dta)
  id_names <- unique(dta[,id_varname])
  
  if (is.null(avail_varname)) {
    avail <- rep(1, total_person_decisionpoint)
  } else {
    avail <- dta[, avail_varname]
  }
  
  A <- dta[, treatment_varname]
  # checking for NA in treatment indicator    
  if (any(is.na(A[avail == 1]))) {
    stop("Treatment indicator is NA where availability = 1.")
  }
  A[avail == 0] <- 0
  
  p_t <- dta[, rand_prob_varname]
  cA <- A - p_t # centered A
  Y <- dta[, outcome_varname]
  
  # X (moderator) design matrix, intercept added
  Xdm <- as.matrix( cbind( rep(1, nrow(dta)), dta[, moderator_varname] ) )
  # Z (control) design matrix, intercept added
  Zdm <- as.matrix( cbind( rep(1, nrow(dta)), dta[, control_varname] ) ) 
  
  if (is.null(rand_prob_tilde_varname) & is.null(rand_prob_tilde)) {
    p_t_tilde <- rep(0.5, nrow(dta))
  } else if (is.null(rand_prob_tilde_varname)) {
    if (length(rand_prob_tilde) == 1) {
      p_t_tilde <- rep(rand_prob_tilde, total_person_decisionpoint)
    } else if (length(rand_prob_tilde) == total_person_decisionpoint) {
      p_t_tilde <- rand_prob_tilde
    } else {
      stop("rand_prob_tilde is of incorrect length.")
    }
  } else {
    p_t_tilde <- dta[, rand_prob_tilde_varname]
  }
  cA_tilde <- A - p_t_tilde
  
  WCLS_weight_inverseProb <- c()
  if (Delta == 1){
    WCLS_weight_inverseProb = rep(1, total_person_decisionpoint)
  }else{
    for (i in 1:sample_size){
      dta_perid <- dta[dta[,id_varname] == id_names[i],]
      dta_timepoint <- length(dta_perid[,decision_time_varname])
      
      #append with A_{T+1} ... A_{T+\delta-1} = 0
      inverseProb_numerator <- 1 - c(dta_perid[,treatment_varname],rep(0,Delta-1)) 
      #append with p_{T+1} ... p_{T+\delta-1} = 0
      inverseProb_denominator <- 1- c(dta_perid[,rand_prob_varname],rep(0,Delta-1)) 
      
      for (j in 1:dta_timepoint) {
        # the first time R_{t+k} = 1
        first_R_1 <- dta_perid$k[j] 
        
        if (first_R_1 == 0){
          # if R_t = 1 right after A_t, set the weight to 1
          WCLS_weight_inverseProb <- c(WCLS_weight_inverseProb,1) 
          
        }else{ 
          # in other words, if A = 0 appears before k, we need to set it to 0, o.w keep this term
          WCLS_weight_inverseProb <- c(WCLS_weight_inverseProb,
                                       prod(inverseProb_numerator[(j+1):(j+first_R_1)]/
                                              inverseProb_denominator[(j+1):(j+first_R_1)]))
        }
      }
    }
  }
  
  p <- length(moderator_varname) + 1 # dimension of beta
  
  ######
  # TODO: fix this
  # Get unique values and their frequencies
  # unique_values <- unique(vector)
  # frequencies <- table(vector)
  # 
  # # Calculate probabilities
  # probs <- frequencies / sum(frequencies)
  # 
  # # Sample from the vector with probabilities proportional to the frequencies
  # sampled_value <- sample(unique_values, size = 1, prob = probs)
  # 
  
  #######
  residuals_per_u <- function(u){
    WCLS_weight_till_u <- c()
    WCLS_weight_start_u <- c()
    A_u <- c()
    S_u <- c()
    if (Delta == 1){
      WCLS_weight_start_u = rep(1, total_person_decisionpoint)
      WCLS_weight_till_u = rep(1, total_person_decisionpoint)
      A_u <- A
      S_u <- dta$S
    }else{
      for (i in 1:sample_size){
        dta_perid <- dta[dta[,id_varname] == id_names[i],]
        dta_timepoint <- length(dta_perid[,decision_time_varname])
        
        #append with A_{T+1} ... A_{T+\delta-1} = 0
        inverseProb_numerator <- 1 - c(dta_perid[,treatment_varname],rep(0,Delta-1)) 
        #append with p_{T+1} ... p_{T+\delta-1} = 0
        inverseProb_denominator <- 1- c(dta_perid[,rand_prob_varname],rep(0,Delta-1))  
        
        for (j in 1:dta_timepoint) {
          A_u <- c(A_u, ifelse(j+u <= dta_timepoint, dta_perid$A[j+u], 0))
          #S_u <- c(S_u, ifelse(j+u <= dta_timepoint, dta_perid$S[j+u], 0))
          S_u <- c(S_u, ifelse(j+u <= dta_timepoint, dta_perid$S[j+u], dta_perid$S[j]))
          
          WCLS_weight_till_u <- c(WCLS_weight_till_u,
                                  prod(inverseProb_numerator[(j+1):(j+u)]/
                                         inverseProb_denominator[(j+1):(j+u)]))
          if (u == Delta -1){
            WCLS_weight_start_u <- c(WCLS_weight_start_u, 1)
          }else{
            WCLS_weight_start_u <- c(WCLS_weight_start_u,
                                     prod(inverseProb_numerator[(j+u+1):(j+Delta-1)]/
                                            inverseProb_denominator[(j+u+1):(j+Delta-1)]))}
        }
      }
    }
    dta$till_u <- WCLS_weight_till_u
    dta$start_u <- WCLS_weight_start_u
    dta$A_u <- A_u
    dta$S_u <- S_u
    # Y_reg <- dta$Y * dta$start_u
    # reg <- lm(Y_reg ~ dta$A_u + dta$S, data = dta)
    # Y_u_pred1 <- coef(reg)[1] + coef(reg)[2]*1 + coef(reg)[3]*dta$S
    # Y_u_pred0 <- coef(reg)[1] + coef(reg)[2]*0 + coef(reg)[3]*dta$S
    # Y_u_pred <- coef(reg)[1] + coef(reg)[2]*A_u + coef(reg)[3]*dta$S
    # 
    ##Y_u_pred_mean <- Y_u_pred1 * mean(A_u) + Y_u_pred0 * (1-mean(A_u))
    # reg0 <- lm(Y_reg ~ dta$S, data = dta)
    # Y_u_pred_mean <- coef(reg0)[1] + coef(reg0)[2] * dta$S
    # Y_u_residual <- WCLS_weight_till_u*(Y_u_pred - Y_u_pred_mean)
    
    Y_reg <- dta$Y * WCLS_weight_inverseProb
    # if(is.null(moderator_varname)){
    #   reg <- lm(Y_reg ~ A_u, data = dta)
    #   Y_u_pred <- coef(reg)[1] + coef(reg)[2]*A_u
    #   reg0 <- lm(Y_reg ~ 1, data = dta)
    #   Y_u_pred_mean <- coef(reg0)[1] 
    # }else{
    reg <- lm(Y_reg ~ A_u + S_u, data = dta)
    Y_u_pred <- coef(reg)[1] + coef(reg)[2]*A_u + coef(reg)[3]*S_u
    reg0 <- lm(Y_reg ~ S_u, data = dta)
    Y_u_pred_mean <- coef(reg0)[1] + coef(reg0)[2] * S_u
    #}
    Y_u_residual <- Y_u_pred - Y_u_pred_mean
    
    return(Y_u_residual)
  }
  
  M <- ifelse(A, p_t_tilde / p_t, (1 - p_t_tilde) / (1 - p_t)) 
  dta$M <- M
  
  # 2. estimate E(Y_{it, \Delta} W_it|H_t, A_t)
  # regress on H_t
  dta$WCLS <- WCLS_weight_inverseProb
  Y_reg <- dta$Y * dta$WCLS
  
  # regress on moderator
  # if(is.null(moderator_varname)){
  #   reg <- lm(Y_reg ~ A, data = dta)
  #   Yt_pred_1 <- coef(reg)[1] + coef(reg)[2]*1 
  #   Yt_pred_0 <- coef(reg)[1] + coef(reg)[2]*0 
  #   Yt_pred <- coef(reg)[1] + coef(reg)[2]*dta$A 
  # }else{
  reg <- lm(Y_reg ~ dta$A + dta$S, data = dta)
  Yt_pred_1 <- coef(reg)[1] + coef(reg)[2]*1 + coef(reg)[3]*dta$S
  Yt_pred_0 <- coef(reg)[1] + coef(reg)[2]*0 + coef(reg)[3]*dta$S
  #Y_t_pred <- Yt_pred_1 * mean(A) +  Yt_pred_0 * (1 - mean(A)) 
  Yt_pred <- coef(reg)[1] + coef(reg)[2]*dta$A + coef(reg)[3]*dta$S
  #}
  
  # for (a in 0:1) {
  #   dta_a = dta[dta$A == a, ]
  #   Y_reg <- dta_a$Y * dta_a$WCLS
  #   # regress on moderator
  #   reg <- lm(Y_reg ~ dta_a$S, data = dta_a)
  #   Y_pred <- coef(reg)[1] + coef(reg)[2]*dta_a$S
  # }
  
  # Now calculating the estimating equation
  # Notice here we ignore the control variables
  ### 2. estimation ###
  
  estimating_equation <- function(beta) {
    beta <- as.matrix(beta)
    exp_AXdm_beta <- exp(A * (Xdm %*% beta))
    weight <- exp_AXdm_beta^(-1)
    residual_u <- c()
    for (u in 1:(Delta-1)) {
      residual_u <- cbind(residual_u, residuals_per_u(u))
    }
    residual_u_sum <- apply(residual_u, 1, sum)  
    
    ef <- rep(NA, length(beta)) # value of estimating function
    for (i in 1:p) {
      # Xdm is S_t, the moderator variables with intercept included
      error_term1 <- avail * weight * M * cA_tilde * Xdm[,i] * 
        (Y * WCLS_weight_inverseProb - residual_u_sum - Yt_pred)
      error_term2 <- avail * Xdm[,i] * p_t_tilde * (1-p_t_tilde) * 
        (exp(-1 *Xdm %*% beta) * Yt_pred_1 - Yt_pred_0)
      ef[i] <- sum(error_term1 + error_term2)
    }
    
    ef <- ef / sample_size
    return(ef)
  }
  
  if (is.null(estimator_initial_value)) {
    estimator_initial_value <- rep(0, length = p)
  }
  
  solution <- tryCatch(
    {
      multiroot(estimating_equation, estimator_initial_value)
    },
    error = function(cond) {
      message("\nCatched error in multiroot inside weighted_centered_least_square():")
      message(cond)
      return(list(root = rep(NaN, p), msg = cond,
                  f.root = rep(NaN, p)))
    })
  
  beta_hat <- as.vector(solution$root)
  
  
  ### 3. asymptotic variance ###
  ### 3.1 Compute M_n matrix (M_n is the empirical expectation of the derivative of the estimating function) ###
  Mn_summand <- array(NA, dim = c(total_person_decisionpoint, p, p))
  
  beta_hat <- as.matrix(beta_hat)
  exp_AXdm_beta_hat <- exp(A * (Xdm %*% beta_hat))
  weight_hat <- exp_AXdm_beta_hat^(-1)
  M2_multiplier <- exp(-1 *Xdm %*% beta_hat) 
  
  
  residual_u <- c()
  for (u in 1:(Delta-1)) {
    residual_u <- cbind(residual_u, residuals_per_u(u))
  }
  residual_u_sum <- apply(residual_u, 1, sum)  
  M1_error <- Y * WCLS_weight_inverseProb - residual_u_sum - Yt_pred

  M1 <- matrix(NA, nrow = p, ncol = total_person_decisionpoint)
  M2 <- matrix(NA, nrow = p, ncol = total_person_decisionpoint)
  M3 <- matrix(NA, nrow = p, ncol = total_person_decisionpoint)
  
  for (it in 1:total_person_decisionpoint) {
    # this is to make R code consistent whether X_it, Z_it contains more entries or is just the intercept.        
    if (p == 1) {
      Xbeta <- Xdm[it, ] * beta_hat
    } else {
      Xbeta <- as.numeric(Xdm[it, ] %*% beta_hat)
    }

    # derivatives
    M1_deriv <- avail[it] * weight_hat[it] * M[it] * (A[it] - p_t_tilde[it]) *
      M1_error[it] * as.matrix(Xdm[it, ]) %*% (-A[it] * t(Xdm[it, ]))
    M2_deriv <- avail[it] * p_t_tilde[it] * (1 - p_t_tilde[it]) * M2_multiplier[it] *
      Yt_pred_1[it] * as.matrix(Xdm[it, ])  %*% (-1 * t(Xdm[it, ]))
    M3_deriv <- 0
    Mn_summand[it, , ] <- M1_deriv + M2_deriv + M3_deriv
    
    M1[,it] <- avail[it] * weight_hat[it] *M[it] * (A[it] - p_t_tilde[it]) *
      M1_error[it] * as.matrix(Xdm[it, ])
    M2[,it] <- avail[it] * p_t_tilde[it] * (1 - p_t_tilde[it]) * M2_multiplier[it] *
      Yt_pred_1[it] * as.matrix(Xdm[it, ]) 
    M3[,it] <- -1 * avail[it] * p_t_tilde[it] * (1 - p_t_tilde[it]) * Yt_pred_0[it] * 
      as.matrix(Xdm[it, ]) 
      
  }
  Mn <- apply(Mn_summand, c(2,3), sum) / sample_size
  Mn_inv <- solve(Mn)
  
  ### 3.2 Compute \Sigma_n matrix (\Sigma_n is the empirical variance of the estimating function) ###
  
  Sigman_summand <- array(NA, dim = c(sample_size, p, p))
  # Sigman_summand is  \sum_{t=1}^T ( D^{(t),T} r^(t) )^{\otimes 2}
  # See note 2018.08.06 about small sample correction
  
  person_first_index <- c(find_change_location(dta[, id_varname]), total_person_decisionpoint + 1)
  M_matrix <- M1 + M2 + M3
  
  for (i in 1:sample_size) {
    if(p == 1){
      M_matrix_i <- t(M_matrix[, person_first_index[i] : (person_first_index[i+1] - 1)])
    }else{
    M_matrix_i <- M_matrix[, person_first_index[i] : (person_first_index[i+1] - 1)]
    }
    Sigman_summand[i, , ] <- M_matrix_i %*% t(M_matrix_i)
  }
  Sigman <- apply(Sigman_summand, c(2,3), sum) / sample_size
  
  varcov <- Mn_inv %*% Sigman %*% t(Mn_inv) / sample_size
  beta_se <- sqrt(diag(varcov)[1:p])
  
  
  ### 4. small sample correction ###
  r_term_collected <- rep(NA, total_person_decisionpoint)
  D_term_collected <- matrix(NA, nrow = p, ncol = total_person_decisionpoint)
  partialr_partialtheta_collected <- matrix(NA, nrow = total_person_decisionpoint, ncol = p)
  #Mn_summand_small <- array(NA, dim = c(total_person_decisionpoint, p, p))
  
  for (it in 1:total_person_decisionpoint) {
    # this is to make R code consistent whether X_it, Z_it contains more entries or is just the intercept.        
    if (p == 1) {
      Xbeta <- Xdm[it, ] * beta_hat
    } else {
      Xbeta <- as.numeric(Xdm[it, ] %*% beta_hat)
    }
    
    r_term1 <- exp(-1 * A[it] * Xbeta) * M[it] * cA_tilde[it] * M1_error[it]
    r_term2 <- p_t_tilde[it] * (1 - p_t_tilde[it]) * (M2_multiplier[it] * Yt_pred_1[it])
    r_term3 <- p_t_tilde[it] * (1 - p_t_tilde[it]) * (- Yt_pred_0[it])
    # r_term = r^(t) (scalar)
    r_term <- r_term1 + r_term2 + r_term3
    r_term_collected[it] <- r_term
    
    # D_term = D^{(t),T} (dim = p * 1)
    D_term <- as.matrix(avail[it] * Xdm[it, ])
    D_term_collected[, it] <- D_term
    
    # partialr_partialtheta = \frac{\partial r^(t)}{\partial \theta^T}
    partialr_partialtheta <- r_term1  %*% (-A[it] * t(Xdm[it, ])) + r_term2  %*% (-1 * t(Xdm[it, ]))
    partialr_partialtheta_collected[it, ] <- partialr_partialtheta
    
    # Mn_summand_small[it, , ] <- D_term %*% partialr_partialtheta
  }
  
  Sigman_tilde <- 0
  for (i in 1:sample_size) {
    D_term_i <- D_term_collected[, person_first_index[i] : (person_first_index[i+1] - 1)]
    if (p == 1){D_term_i <- t(D_term_i)}
    r_term_i <- matrix(r_term_collected[person_first_index[i] : (person_first_index[i+1] - 1)], ncol = 1)
    partialr_partialtheta_i <- partialr_partialtheta_collected[person_first_index[i] : (person_first_index[i+1] - 1), ]
    H_ii <- partialr_partialtheta_i %*% Mn_inv %*% D_term_i / sample_size
    Ii_minus_Hii_inv <- solve(diag(nrow(H_ii)) - H_ii)
    
    Sigman_tilde <- Sigman_tilde + D_term_i %*% Ii_minus_Hii_inv %*% r_term_i %*% t(r_term_i) %*% t(Ii_minus_Hii_inv) %*% t(D_term_i)
  }
  Sigman_tilde <- Sigman_tilde / sample_size
  
  varcov_adjusted <- Mn_inv %*% Sigman_tilde %*% t(Mn_inv) / sample_size
  beta_se_adjusted <- sqrt(diag(varcov_adjusted)[1:p])
  
  ### 5. return the result with variable names ###
  return(list(beta_hat = beta_hat, 
              beta_se = beta_se, beta_se_adjusted = beta_se_adjusted,
              varcov = varcov, varcov_adjusted = varcov_adjusted,
              dims = list(p = p),
              f.root = solution$f.root))
}


