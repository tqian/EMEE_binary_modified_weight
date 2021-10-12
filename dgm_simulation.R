# Yihan Bao, 2021.05.19
# In this code, I constructed simulation dataset.

dgm_binary_categorical_covariate_new <- function(sample_size, total_T, Delta, prob_a = 0.2) {

  beta_0 <- 0.1
  beta_1 <- 0.2
  
  df_names <- c("userid", "day", "A", "S", "S2","prob_A","prob_R_0", "prob_R","R","Y","k") 
  # prob_R_0 is the probability of R = 0 given A = 0.
  # prob_R is the probability of R = 0 given A = 1.
  
  dta <- data.frame(matrix(NA, nrow = sample_size * total_T, ncol = length(df_names)))
  names(dta) <- df_names
  
  dta$userid <- rep(1:sample_size, each = total_T)
  dta$day <- rep(1:total_T, times = sample_size)
  
  C = (0.5^(0.5/Delta) + 1 + 0.5^(-0.5/Delta))
  prob_S_weight = c(0.5^(-0.5/Delta)/C, 1/C, 0.5^(0.5/Delta)/C)
  E_S = 3*0.5^(1/Delta)/C
  
  for (t in 1:total_T) {
    # row index for the rows corresponding to day t for every subject
    row_index <- seq(from = t, by = total_T, length = sample_size)
    dta$S[row_index] <- sample(c(0,1,2), sample_size, prob = prob_S_weight, replace = TRUE)
    dta$S2[row_index] <- ifelse(dta$S[row_index] == 2, 1, 0) 
    dta$prob_A[row_index] <- rep(prob_a, sample_size)
    dta$A[row_index] <- rbinom(sample_size, 1, dta$prob_A[row_index])
    
    dta$prob_R_0[row_index] <- 0.5^((1.5-0.5*dta$S[row_index])/Delta)
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

## true beta for (Intercept, S)
beta_true <- c(0.1, 0.2)

## true beta marginal
beta_true_marginal_generalDelta <- function(Delta){
  beta_0 <- 0.1
  beta_1 <- 0.2
  
  C <- (0.5^(0.5/Delta) + 1 + 0.5^(-0.5/Delta))
  prob_S_weight <- c(0.5^(-0.5/Delta)/C, 1/C, 0.5^(0.5/Delta)/C)
  E_S = 3*0.5^(1/Delta)/C
  prob_R_0 <-  c(0.5^((1.5-0.5*0)/Delta),0.5^((1.5-0.5*1)/Delta),0.5^((1.5-0.5*2)/Delta))
  exp_h <- c(exp(beta_0+beta_1*0),exp(beta_0+beta_1*1),exp(beta_0+beta_1*2))
  
  numerator <- sum( ((1- prob_R_0*E_S^(Delta-1)) * exp_h) * prob_S_weight)
  denominator <-  sum( ((1- prob_R_0*E_S^(Delta-1))) * prob_S_weight)
  beta_true_marginal <- log(numerator / denominator)
  return(beta_true_marginal)
}

# compute marginal beta_true
if (0){
  total_T = 30
  Delta = 3
  beta_0 <- 0.1
  beta_1 <- 0.2
  prob_a <- 0.2
  ### analytically, we have ###
  C = (0.5^(0.5/Delta) + 1 + 0.5^(-0.5/Delta))
  prob_S_weight = c(0.5^(-0.5/Delta)/C, 1/C, 0.5^(0.5/Delta)/C)
  E_S = 3*0.5^(1/Delta)/C
  
  prob_R_0 <-  c(0.5^((1.5-0.5*0)/Delta),0.5^((1.5-0.5*1)/Delta),0.5^((1.5-0.5*2)/Delta))
  exp_h <- c(exp(beta_0+beta_1*0),exp(beta_0+beta_1*1),exp(beta_0+beta_1*2))
  
  numerator <- sum( ((1- prob_R_0*E_S^(Delta-1)) * exp_h) * prob_S_weight)
  denominator <-  sum( ((1- prob_R_0*E_S^(Delta-1))) * prob_S_weight)
  
  beta_true_marginal <- log(numerator / denominator)
  print(beta_true_marginal) #0.2827493 when Delta = 3
}


