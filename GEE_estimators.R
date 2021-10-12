# Yihan Bao, 2021.05.19
# Estimators to compare with modified-EMEE for paper
# - log linear GEE
# - brm (by Richardson et al. 2017 JASA)
# Codes are adopted from the original EMEE paper. (in Qian et al. 2021)

library(rootSolve) # for solver function multiroot()
library(geepack) # for fitting GEE using package

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

log_linear_GEE <- function(
    dta,
    id_varname,
    decision_time_varname,
    treatment_varname,
    outcome_varname,
    control_varname,
    moderator_varname,
    avail_varname = NULL,
    estimator_initial_value = NULL
)
{
    ### 1. preparation ###
    
    sample_size <- length(unique(dta[, id_varname]))
    total_person_decisionpoint <- nrow(dta)
    
    A <- dta[, treatment_varname]
    Y <- dta[, outcome_varname]
    Xdm <- as.matrix( cbind( rep(1, nrow(dta)), dta[, moderator_varname] ) ) # X (moderator) design matrix, intercept added
    Zdm <- as.matrix( cbind( rep(1, nrow(dta)), dta[, control_varname] ) ) # Z (control) design matrix, intercept added
    
    if (is.null(avail_varname)) {
        avail <- rep(1, total_person_decisionpoint)
    } else {
        avail <- dta[, avail_varname]
    }
    
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
        residual <- Y - exp_AXdm_beta * exp_Zdm_alpha
        weight <- (1 - exp_AXdm_beta * exp_Zdm_alpha)^(-1)
        
        ef <- rep(NA, length(theta)) # value of estimating function
        for (i in 1:q) {
            ef[i] <- sum( weight * residual * avail * Zdm[, i])
        }
        for (i in 1:p) {
            ef[q + i] <- sum( weight * residual * avail * A * Xdm[, i])
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
            message("\nCatched error in multiroot inside log_linear_GEE():")
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
        
        pre_multiplier <- (1 - exp(Zalpha + A[it] * Xbeta))^(-2) * exp(Zalpha + A[it] * Xbeta)
        
        # partialD_partialtheta = \frac{\partial D^{(t),T}}{\partial \theta^T}, matrix of dim (p+q)*(p+q)
        partialD_partialtheta <- matrix(NA, nrow = p + q, ncol = p + q)
        partialD_partialtheta[1:q, 1:q] <- pre_multiplier * (Zdm[it, ] %o% Zdm[it, ])
        partialD_partialtheta[1:q, (q+1):(q+p)] <- pre_multiplier * A[it] * (Zdm[it, ] %o% Xdm[it, ])
        partialD_partialtheta[(q+1):(q+p), 1:q] <- pre_multiplier * A[it] * (Xdm[it, ] %o% Zdm[it, ])
        partialD_partialtheta[(q+1):(q+p), (q+1):(q+p)] <- pre_multiplier * A[it] * (Xdm[it, ] %o% Xdm[it, ])
        
        # r_term = r^(t) (scalar)
        r_term <- (Y[it] - exp(Zalpha + A[it] * Xbeta)) * avail[it]
        r_term_collected[it] <- r_term
        
        # D_term = D^{(t),T} (dim = (p+q) * 1)
        D_term <- (1 - exp(Zalpha + A[it] * Xbeta))^(-1) * c(Zdm[it, ], A[it] * Xdm[it, ])
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
    
    person_first_index <- c(find_change_location(dta[, id_varname]), total_person_decisionpoint + 1)
    
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
                # test_result_t = test_result_t,
                # test_result_f = test_result_f,
                varcov = varcov,
                varcov_adjusted = varcov_adjusted,
                dims = list(p = p, q = q),
                f.root = solution$f.root))
}



log_linear_GEE_geepack <- function(
    dta,
    id_varname,
    decision_time_varname,
    treatment_varname,
    outcome_varname,
    control_varname,
    moderator_varname,
    estimator_initial_value = NULL,
    corstr = "independence" # could also use, e.g., "exchangeable"
){
    p <- length(moderator_varname) + 1 # dimension of beta
    q <- length(control_varname) + 1 # dimension of alpha
    
    Xnames <- c("Intercept", moderator_varname)
    Znames <- c("Intercept", control_varname)
    
    control_summed <- paste0(control_varname, collapse = " + ")
    if (p > 1) {
        moderator_summed <- paste0("* (", paste0(moderator_varname, collapse = " + "), ")")
    } else {
        moderator_summed <- ""
    }
    
    gee_formula <- as.formula(paste(outcome_varname, "~", control_summed, "+", treatment_varname, moderator_summed))
    fit_geepack <- relRisk(gee_formula, data = dta, corstr = corstr, id = dta[, id_varname])
    
    alpha_hat <- fit_geepack$beta[1:q]
    beta_hat <- fit_geepack$beta[(q+1):(q+p)]
    varcov <- fit_geepack$vbeta
    alpha_se <- sqrt(diag(varcov)[1:q])
    beta_se <- sqrt(diag(varcov)[(q+1):(q+p)])
    alpha_se_adjusted <- alpha_se
    beta_se_adjusted <- beta_se
    
    names(alpha_hat) <- names(alpha_se) <- names(alpha_se_adjusted) <- Znames
    names(beta_hat) <- names(beta_se) <- names(beta_se_adjusted) <- Xnames
    
    return(list(beta_hat = beta_hat, alpha_hat = alpha_hat,
                beta_se = beta_se, alpha_se = alpha_se,
                beta_se_adjusted = beta_se_adjusted, alpha_se_adjusted = alpha_se_adjusted,
                # test_result_t = test_result_t,
                # test_result_f = test_result_f,
                varcov = varcov,
                dims = list(p = p, q = q),
                f.root = rep(-1, p+q)))
}




brm_MLE <- function(
    dta,
    id_varname,
    decision_time_varname,
    treatment_varname,
    outcome_varname,
    control_varname,
    moderator_varname,
    rand_prob_varname,
    estimator_initial_value = NULL
) {
    
    sample_size <- length(unique(dta[, id_varname]))
    total_person_decisionpoint <- nrow(dta)
    
    A <- dta[, treatment_varname]
    Y <- dta[, outcome_varname]
    Xdm <- as.matrix( cbind( rep(1, nrow(dta)), dta[, moderator_varname] ) ) # X (moderator) design matrix, intercept added
    Zdm <- as.matrix( cbind( rep(1, nrow(dta)), dta[, control_varname] ) ) # Z (control) design matrix, intercept added
    
    p <- length(moderator_varname) + 1 # dimension of beta
    q <- length(control_varname) + 1 # dimension of alpha
    
    Xnames <- c("Intercept", moderator_varname)
    Znames <- c("Intercept", control_varname)
    
    est <- brm(Y, A, va = Xdm, vb = Zdm, param = "RR")
    
    return(list(beta_hat = est$point.est[1:p], alpha_hat = est$point.est[(p+1):(p+q)],
                beta_se = est$se.est[1:p], alpha_se = est$se.est[(p+1):(p+q)],
                beta_se_adjusted = est$se.est[1:p], alpha_se_adjusted = est$se.est[(p+1):(p+q)],
                # test_result_t = test_result_t,
                # test_result_f = test_result_f,
                varcov = NA,
                # varcov_ssa = asymp_varcov_ssa / sample_size,
                dims = list(p = p, q = q),
                f.root = rep(-1, p+q)))
}

brm_DR <- function(
    dta,
    id_varname,
    decision_time_varname,
    treatment_varname,
    outcome_varname,
    control_varname,
    moderator_varname,
    rand_prob_varname,
    estimator_initial_value = NULL
) {
    
    sample_size <- length(unique(dta[, id_varname]))
    total_person_decisionpoint <- nrow(dta)
    
    A <- dta[, treatment_varname]
    Y <- dta[, outcome_varname]
    Xdm <- as.matrix( cbind( rep(1, nrow(dta)), dta[, moderator_varname] ) ) # X (moderator) design matrix, intercept added
    Zdm <- as.matrix( cbind( rep(1, nrow(dta)), dta[, control_varname] ) ) # Z (control) design matrix, intercept added
    
    p <- length(moderator_varname) + 1 # dimension of beta
    q <- length(control_varname) + 1 # dimension of alpha
    
    Xnames <- c("Intercept", moderator_varname)
    Znames <- c("Intercept", control_varname)
    
    est <- brm(Y, A, va = Xdm, vb = Zdm, param = "RR", est.method = "DR")
    
    return(list(beta_hat = est$point.est[1:p], alpha_hat = est$point.est[(p+1):(p+q)],
                beta_se = est$se.est[1:p], alpha_se = est$se.est[(p+1):(p+q)],
                beta_se_adjusted = est$se.est[1:p], alpha_se_adjusted = est$se.est[(p+1):(p+q)],
                # test_result_t = test_result_t,
                # test_result_f = test_result_f,
                varcov = NA,
                # varcov_ssa = asymp_varcov_ssa / sample_size,
                dims = list(p = p, q = q),
                f.root = rep(-1, p+q)))
}
