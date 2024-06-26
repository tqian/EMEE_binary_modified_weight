# Tianchen Qian
# 2023.12.01

# EMEE estimator with user-input IPW term.
# Note that this IPW term is the prod_{s=1}^{\Delta-1} term, not the stabilized weight term.

EMEE_userinput_ipw <- function(
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
        excursion_ipw_varname # the variabler name 
)
{
    ############## description ###############
    ##
    ## This function estimates the moderated treatment effect for binary outcome,
    ## and provides variance estimate.
    ##
    ## It incorporates two methods for small sample correction:
    ## 1) the usage of "Hat" matrix in the variance estimate (as in Mancl & DeRouen 2001)
    ## 2) the usage of t-distribution or F-distribution critical value with corrected degrees of freedom
    ##    (as in Liao et al. 2015)
    
    ############## arguments ###############
    ##
    ## dta...................the data set in long format
    ## id_varname............variable name for subject id (to distinguish between subjects in dta)
    ## decision_time_varname.....variable name for decision points in study
    ## outcome_varname.......variable name for outcome variable
    ## control_varname.......vector of variable names used to reduce noise (Z in the model),
    ##                       could be NULL (no control covariates)
    ## moderator_varname.....vector of variable names as effect modifiers (X in the model),
    ##                       could be NULL (no effect modifier)
    ## rand_prob_varname.....variable name for randomization probability (a column in dta)
    ## avail_varname.........variable name for availability (a column in dta)
    ##                       default to NULL (available at all decision points)
    ## rand_prob_tilde_varname.....variable name for \tilde{p}_t(1|H_t) (a column in dta)
    ##                             this is the arbitrary weight used in WCLS
    ##                             default to NULL (in which case \tilde{p}_t(1|H_t) is set to 0.5)
    ## rand_prob_tilde.............a numeric vector of the same length as dta
    ##                             this is another way to specify \tilde{p}_t(1|H_t)
    ##                             default to NULL (in which case \tilde{p}_t(1|H_t) is set to 0.5)
    ## estimator_initial_value.....initial value for the estimator in the root finding algorithm
    ##                             length is len(control_varname) + len(moderator_varname) + 2
    ##                             default to NULL (in which case the intial value = all 0's)
    ## Delta ................ the length of time window
    ############## return value ###############
    ##
    ## This function returns a list of the following components:
    ##
    ## beta_hat..............estimated beta
    ## alpha_hat.............estimated alpha
    ## beta_se...............standard error for beta_hat
    ## alpha_se..............standard error for alpha_hat
    ## beta_se_adjusted......standard error for beta_hat, with small sample correction (hat matrix)
    ## alpha_se_adjusted.....standard error for alpha_hat, with small sample correction (hat matrix)
    ## varcov................estimated variance-covariance matrix for (alpha_hat, beta_hat)
    ## varcov_adjusted.......estimated variance-covariance matrix for (alpha_hat, beta_hat), with small sample correction (hat matrix)
    ## f.root................value of the estimating function at (alpha_hat, beta_hat)
    
    ##########################################
    
    ### 1. preparation ###
    
    sample_size <- length(unique(dta[[id_varname]]))
    total_person_decisionpoint <- nrow(dta)
    id_names <- unique(dta[,id_varname])
    
    if (is.null(avail_varname)) {
        avail <- rep(1, total_person_decisionpoint)
    } else {
        avail <- dta[[avail_varname]]
    }
    
    A <- dta[[treatment_varname]]
    # checking for NA in treatment indicator    
    if (any(is.na(A[avail == 1]))) {
        stop("Treatment indicator is NA where availability = 1.")
    }
    A[avail == 0] <- 0
    
    p_t <- dta[[rand_prob_varname]]
    cA <- A - p_t # centered A
    Y <- dta[[outcome_varname]]
    Xdm <- as.matrix( cbind( rep(1, nrow(dta)), dta[, moderator_varname] ) ) # X (moderator) design matrix, intercept added
    Zdm <- as.matrix( cbind( rep(1, nrow(dta)), dta[, control_varname] ) ) # Z (control) design matrix, intercept added
    
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
        p_t_tilde <- dta[[rand_prob_tilde_varname]]
    }
    cA_tilde <- A - p_t_tilde
    
    excursion_ipw <- dta[[excursion_ipw_varname]]
    
    WCLS_weight_withDelta <- ifelse(A, p_t_tilde / p_t, (1 - p_t_tilde) / (1 - p_t)) * excursion_ipw
    
    p <- length(moderator_varname) + 1 # dimension of beta
    q <- length(control_varname) + 1 # dimension of alpha
    
    Xnames <- c("Intercept", moderator_varname)
    Znames <- c("Intercept", control_varname)
    
    ### 2. estimation ###
    
    estimating_equation <- function(theta) {
        alpha <- as.matrix(theta[1:q])
        beta <- as.matrix(theta[(q+1):(q+p)])
        
        exp_Zdm_alpha <- exp(Zdm %*% alpha)
        exp_AXdm_beta <- exp(A * (Xdm %*% beta))
        
        residual <- Y - exp_Zdm_alpha * exp_AXdm_beta
        weight <- exp_AXdm_beta^(-1)
        
        ef <- rep(NA, length(theta)) # value of estimating function
        for (i in 1:q) {
            ef[i] <- sum( weight * residual * avail * WCLS_weight_withDelta * Zdm[, i]) #Zdm is g(H_t), the control variables with interept included
        }
        for (i in 1:p) {
            ef[q + i] <- sum( weight * residual * avail * WCLS_weight_withDelta * cA_tilde * Xdm[, i]) #Xdm is S_t, the moderator variables with intercept included
        }
        
        ef <- ef / sample_size
        return(ef)
    }
    
    if (is.null(estimator_initial_value)) {
        estimator_initial_value <- rep(0, length = p + q)
    }
    
    solution <- tryCatch(
        {
            multiroot(estimating_equation, estimator_initial_value)
        },
        error = function(cond) {
            message("\nCatched error in multiroot inside weighted_centered_least_square():")
            message(cond)
            return(list(root = rep(NaN, p + q), msg = cond,
                        f.root = rep(NaN, p + q)))
        })
    
    estimator <- get_alpha_beta_from_multiroot_result(solution, p, q)
    alpha_hat <- as.vector(estimator$alpha)
    beta_hat <- as.vector(estimator$beta)
    
    ### 3. asymptotic variance ###
    
    ### 3.1 Compute M_n matrix (M_n is the empirical expectation of the derivative of the estimating function) ###
    
    Mn_summand <- array(NA, dim = c(total_person_decisionpoint, p+q, p+q))
    # Mn_summand is \frac{\partial D^{(t),T}}{\partial \theta^T} r^(t) + D^{(t),T} \frac{\partial r^(t)}{\partial \theta^T}
    # See note 2018.08.06 about small sample correction
    
    r_term_collected <- rep(NA, total_person_decisionpoint)
    D_term_collected <- matrix(NA, nrow = p+q, ncol = total_person_decisionpoint)
    partialr_partialtheta_collected <- matrix(NA, nrow = total_person_decisionpoint, ncol = p+q)
    
    for (it in 1:total_person_decisionpoint) {
        # this is to make R code consistent whether X_it, Z_it contains more entries or is just the intercept.        
        if (p == 1) {
            Xbeta <- Xdm[it, ] * beta_hat
        } else {
            Xbeta <- as.numeric(Xdm[it, ] %*% beta_hat)
        }
        if (q == 1) {
            Zalpha <- Zdm[it, ] * alpha_hat
        } else {
            Zalpha <- as.numeric(Zdm[it, ] %*% alpha_hat)
        }
        
        pre_multiplier <- exp(- A[it] * Xbeta) * WCLS_weight_withDelta[it]
        
        # partialD_partialtheta = \frac{\partial D^{(t),T}}{\partial \theta^T}, matrix of dim (p+q)*(p+q)
        partialD_partialtheta <- matrix(NA, nrow = p + q, ncol = p + q)
        partialD_partialtheta[1:q, 1:q] <- 0
        partialD_partialtheta[1:q, (q+1):(q+p)] <- - pre_multiplier * A[it] * (Zdm[it, ] %o% Xdm[it, ])
        partialD_partialtheta[(q+1):(q+p), 1:q] <- 0
        partialD_partialtheta[(q+1):(q+p), (q+1):(q+p)] <- - pre_multiplier * A[it] * cA_tilde[it] * (Xdm[it, ] %o% Xdm[it, ])
        
        # r_term = r^(t) (scalar)
        r_term <- (Y[it] - exp(Zalpha + A[it] * Xbeta)) * avail[it]
        r_term_collected[it] <- r_term
        
        # D_term = D^{(t),T} (dim = (p+q) * 1)
        D_term <- pre_multiplier * c(Zdm[it, ], cA_tilde[it] * Xdm[it, ])
        D_term_collected[, it] <- D_term
        
        # partialr_partialtheta = \frac{\partial r^(t)}{\partial \theta^T}
        partialr_partialtheta <- - exp(Zalpha + A[it] * Xbeta) * c(Zdm[it, ], A[it] * Xdm[it, ]) * avail[it]
        partialr_partialtheta_collected[it, ] <- partialr_partialtheta
        
        Mn_summand[it, , ] <- partialD_partialtheta * r_term + D_term %o% partialr_partialtheta
    }
    Mn <- apply(Mn_summand, c(2,3), sum) / sample_size
    Mn_inv <- solve(Mn)
    
    ### 3.2 Compute \Sigma_n matrix (\Sigma_n is the empirical variance of the estimating function) ###
    
    Sigman_summand <- array(NA, dim = c(sample_size, p+q, p+q))
    # Sigman_summand is  \sum_{t=1}^T ( D^{(t),T} r^(t) )^{\otimes 2}
    # See note 2018.08.06 about small sample correction
    
    person_first_index <- c(find_change_location(dta[[id_varname]]), total_person_decisionpoint + 1)
    
    for (i in 1:sample_size) {
        D_term_i <- D_term_collected[, person_first_index[i] : (person_first_index[i+1] - 1)]
        r_term_i <- matrix(r_term_collected[person_first_index[i] : (person_first_index[i+1] - 1)], ncol = 1)
        
        Sigman_summand[i, , ] <- D_term_i %*% r_term_i %*% t(r_term_i) %*% t(D_term_i)
    }
    Sigman <- apply(Sigman_summand, c(2,3), sum) / sample_size
    
    varcov <- Mn_inv %*% Sigman %*% t(Mn_inv) / sample_size
    alpha_se <- sqrt(diag(varcov)[1:q])
    beta_se <- sqrt(diag(varcov)[(q+1):(q+p)])
    
    
    ### 4. small sample correction ###
    
    Sigman_tilde <- 0
    for (i in 1:sample_size) {
        D_term_i <- D_term_collected[, person_first_index[i] : (person_first_index[i+1] - 1)]
        r_term_i <- matrix(r_term_collected[person_first_index[i] : (person_first_index[i+1] - 1)], ncol = 1)
        partialr_partialtheta_i <- partialr_partialtheta_collected[person_first_index[i] : (person_first_index[i+1] - 1), ]
        H_ii <- partialr_partialtheta_i %*% Mn_inv %*% D_term_i / sample_size
        Ii_minus_Hii_inv <- solve(diag(nrow(H_ii)) - H_ii)
        
        Sigman_tilde <- Sigman_tilde + D_term_i %*% Ii_minus_Hii_inv %*% r_term_i %*% t(r_term_i) %*% t(Ii_minus_Hii_inv) %*% t(D_term_i)
    }
    Sigman_tilde <- Sigman_tilde / sample_size
    
    varcov_adjusted <- Mn_inv %*% Sigman_tilde %*% t(Mn_inv) / sample_size
    alpha_se_adjusted <- sqrt(diag(varcov_adjusted)[1:q])
    beta_se_adjusted <- sqrt(diag(varcov_adjusted)[(q+1):(q+p)])
    
    
    ### 6. return the result with variable names ###
    
    names(alpha_hat) <- names(alpha_se) <- names(alpha_se_adjusted) <- Znames
    names(beta_hat) <- names(beta_se) <- names(beta_se_adjusted) <- Xnames
    
    return(list(beta_hat = beta_hat, alpha_hat = alpha_hat,
                beta_se = beta_se, alpha_se = alpha_se,
                beta_se_adjusted = beta_se_adjusted, alpha_se_adjusted = alpha_se_adjusted,
                varcov = varcov,
                varcov_adjusted = varcov_adjusted,
                dims = list(p = p, q = q),
                f.root = solution$f.root))
}

get_alpha_beta_from_multiroot_result <- function(root, p, q)
{
    if (p == 1) {
        beta_root <- root$root[q+1]
    } else {
        beta_root <- as.matrix(root$root[(q+1) : (q+p)])
    }
    if (q == 1) {
        alpha_root <- root$root[1]
    } else {
        alpha_root <- as.matrix(root$root[1:q])
    }
    return(list(alpha = alpha_root, beta = beta_root))
}

find_change_location <- function(v){
    n <- length(v)
    if (n <= 1) {
        stop("The vector need to have length > 1.")
    }
    return(c(1, 1 + which(v[1:(n-1)] != v[2:n])))
}
# examples
# v <- c("a", "a", "b", "c", "c"); find_change_location(v)
# [1] 1 3 4
