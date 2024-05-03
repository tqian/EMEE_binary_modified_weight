if (0) {
    rm(list = ls())
    setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
    source("dgm_simulation.R")
    source("helper_functions.R")
    library(rootSolve)


    Delta <- 3
    control_vars <- c("S", "S2")
    moderator_vars <- c()

    set.seed(123)
    dta <- dgm_binary_categorical_covariate_new(sample_size = 30, total_T = 30, Delta = Delta)

    fit <- EMEE_ImprovedEffByProjection_measurableWeightNotTakenOut_wrapper(
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
    
    c(fit$beta_hat, fit$beta_se_adjusted)
}


EMEE_ImprovedEffByProjection_measurableWeightNotTakenOut_wrapper <- function(
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
        rand_prob_tilde = NULL, # \tilde{p}_t(1|H_t) in WCLS (numeric number or vector)
        estimator_initial_value = NULL,
        Delta,
        link_function_for_nuisance_prediction = c("binomial", "poisson")
) {
    
    # compute excursion_ipw_matrix
    
    dta$ipw_pd <- compute_excursion_ipw(
        dta,
        id_varname,
        decision_time_varname,
        treatment_varname,
        outcome_varname,
        rand_prob_varname,
        avail_varname = NULL,
        Delta
    )
    
    fit <- EMEE_ImprovedEffByProjection_userinput_ipw_measurableWeightNotTakenOut(
        dta,
        id_varname,
        decision_time_varname,
        treatment_varname,
        outcome_varname,
        control_varname,
        moderator_varname,
        rand_prob_varname,
        avail_varname,
        rand_prob_tilde_varname,
        rand_prob_tilde,
        estimator_initial_value,
        Delta,
        excursion_ipw_varname = "ipw_pd",
        link_function_for_nuisance_prediction = link_function_for_nuisance_prediction
    )
    
    return(fit)
}



expit <- function(x) {
    1 / (1 + exp(-x))
}


EMEE_ImprovedEffByProjection_userinput_ipw_measurableWeightNotTakenOut <- function(
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
        rand_prob_tilde = NULL, # \tilde{p}_t(1|H_t) in WCLS (numeric number or vector)
        estimator_initial_value = NULL,
        Delta,
        excursion_ipw_varname, # the variable name
        link_function_for_nuisance_prediction = c("binomial", "poisson")
) {
    sample_size <- length(unique(dta[[id_varname]]))
    total_person_decisionpoint <- nrow(dta)
    id_names <- unique(dta[[id_varname]])
    
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
    
    S <- dta[, moderator_varname]
    
    p_t <- dta[[rand_prob_varname]]
    cA <- A - p_t # centered A
    Y <- dta[[outcome_varname]]
    
    # X (moderator) design matrix, intercept added
    Xdm <- as.matrix(cbind(rep(1, nrow(dta)), dta[, moderator_varname]))
    
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
    
    M <- ifelse(A, p_t_tilde / p_t, (1 - p_t_tilde) / (1 - p_t)) # marginalization weight
    dta$M <- M
    
    W <- dta[[excursion_ipw_varname]] # excursion weight
    YW <- Y * W
    
    p <- length(moderator_varname) + 1 # dimension of beta
    
    nuisance_predictions <- build_nuisance_predictions(
        dta = dta,
        id_varname = id_varname,
        decision_time_varname = decision_time_varname,
        treatment_varname = treatment_varname,
        outcome_varname = outcome_varname,
        control_varname = control_varname,
        rand_prob_varname = rand_prob_varname,
        avail_varname = avail_varname,
        Delta = Delta,
        excursion_ipw_varname = excursion_ipw_varname,
        link_function = link_function_for_nuisance_prediction
    )
    
    E_YWit_given_HAit <- nuisance_predictions$E_YWit_given_HAit
    E_YWit_given_HAit_Aeq1 <- nuisance_predictions$E_YWit_given_HAit_Aeq1
    E_YWit_given_HAit_Aeq0 <- nuisance_predictions$E_YWit_given_HAit_Aeq0
    E_YWit_given_HAiu <- nuisance_predictions$E_YWit_given_HAiu
    E_YWit_given_Hiu <- nuisance_predictions$E_YWit_given_Hiu
    
    ### 2. estimation ###
    
    # residual_term1 is the residual term in first two rows of the updated 
    # estimating equation that uses projection to improve efficiency,
    # which is eq (4.8) in the paper dated 2024.04.15, right before Algorithm 1
    
    # residual_term1 = Y_{it,\Delta}W_it
    #.      - \sum_{u=t+1}^{t+\Delta-1} [ E(Y_{it,\Delta}W_it | H_{iu}, A_{iu}) - E(Y_{it,\Delta}W_it | H_{iu}) ] 
    #.      - E(Y_{it,\Delta}W_it | H_{it}, A_{it})
    
    if (Delta >= 3) {
        residual_term1 <- as.numeric(YW) - rowSums(E_YWit_given_HAiu) + 
            rowSums(E_YWit_given_Hiu) - as.numeric(E_YWit_given_HAit)
    } else if (Delta == 2) {
        residual_term1 <- as.numeric(YW) - as.numeric(E_YWit_given_HAiu) + 
            as.numeric(E_YWit_given_Hiu) - as.numeric(E_YWit_given_HAit)
    } else if (Delta == 1) {
        residual_term1 <- as.numeric(YW) - as.numeric(E_YWit_given_HAit)
    } else {
        stop("Delta must be an integer >= 1.")
    }
    
    estimating_equation <- function(beta) {
        beta <- as.matrix(beta, ncol = 1)
        exp_AXdm_beta <- exp(A * (Xdm %*% beta))
        blip_down_term <- exp_AXdm_beta^(-1)
        
        # residual_term2 is the residual term in last (third) rows of the updated 
        # estimating equation that uses projection to improve efficiency,
        # which is eq (4.8) in the paper dated 2024.04.15, right before Algorithm 1
        
        # residual_term2 = e^{-S_{it}^T\beta} E(Y_{it,\Delta}W_it | H_{it}, A_{it} = 1) - E(Y_{it,\Delta}W_it | H_{it}, A_{it} = 0)
        
        residual_term2 <- as.numeric(exp(-1 * Xdm %*% beta)) * as.numeric(E_YWit_given_HAit_Aeq1) - as.numeric(E_YWit_given_HAit_Aeq0)
        
        ef <- rep(NA, length(beta)) # value of estimating function
        for (i in 1:p) {
            # Xdm is S_t, the moderator variables with intercept included
            ef_term1 <- avail * blip_down_term * M * cA_tilde * Xdm[, i] * residual_term1
            ef_term2 <- avail * Xdm[, i] * p_t_tilde * (1 - p_t_tilde) * residual_term2
            ef[i] <- sum(ef_term1 + ef_term2)
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
            return(list(
                root = rep(NaN, p), msg = cond,
                f.root = rep(NaN, p)
            ))
        }
    )
    
    beta_hat <- as.vector(solution$root)
    
    
    ### 3. asymptotic variance ###
    ### 3.1 Compute M_n matrix (M_n is the empirical expectation of the derivative of the estimating function) ###
    
    Mn_summand <- array(NA, dim = c(total_person_decisionpoint, p, p))
    
    beta_hat <- as.matrix(beta_hat)
    exp_AXdm_beta_hat <- exp(A * (Xdm %*% beta_hat))
    weight_hat <- exp_AXdm_beta_hat^(-1)
    M2_multiplier <- exp(-1 * Xdm %*% beta_hat)
    
    M1_error <- residual_term1
    
    # M1 is the first two lines of the estimating function in Eq. (4.8)
    # M2 is the first term in the last line of the estimating function in Eq. (4.8)
    # M3 is the second term in the last line of the estimating function in Eq. (4.8)
    
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
        
        M1[, it] <- avail[it] * weight_hat[it] * M[it] * (A[it] - p_t_tilde[it]) *
            M1_error[it] * as.matrix(Xdm[it, ])
        
        M2[, it] <- avail[it] * p_t_tilde[it] * (1 - p_t_tilde[it]) * M2_multiplier[it] *
            E_YWit_given_HAit_Aeq1[it] * as.matrix(Xdm[it, ])
        
        M3[, it] <- - 1 * avail[it] * p_t_tilde[it] * (1 - p_t_tilde[it]) * E_YWit_given_HAit_Aeq0[it] *
            as.matrix(Xdm[it, ])
        
        # derivatives
        M1_deriv <- avail[it] * weight_hat[it] * M[it] * (A[it] - p_t_tilde[it]) *
            M1_error[it] * as.matrix(Xdm[it, ]) %*% (-A[it] * t(Xdm[it, ]))
        
        M2_deriv <- avail[it] * p_t_tilde[it] * (1 - p_t_tilde[it]) * M2_multiplier[it] *
            E_YWit_given_HAit_Aeq1[it] * as.matrix(Xdm[it, ]) %*% (-1 * t(Xdm[it, ]))
        
        M3_deriv <- 0
        
        Mn_summand[it, , ] <- M1_deriv + M2_deriv + M3_deriv
    }
    Mn <- apply(Mn_summand, c(2, 3), sum) / sample_size
    Mn_inv <- solve(Mn)
    
    ### 3.2 Compute \Sigma_n matrix (\Sigma_n is the empirical variance of the estimating function) ###
    
    Sigman_summand <- array(NA, dim = c(sample_size, p, p))
    # Sigman_summand is  \sum_{t=1}^T ( D^{(t),T} r^(t) )^{\otimes 2}
    # See note 2018.08.06 about small sample correction
    
    person_first_index <- c(find_change_location(dta[[id_varname]]), total_person_decisionpoint + 1)
    M_matrix <- M1 + M2 + M3
    
    for (i in 1:sample_size) {
        # if(p == 1){
        #   M_matrix_i <- t(M_matrix[, person_first_index[i] : (person_first_index[i+1] - 1)])
        # }else{
        # M_matrix_i <- M_matrix[, person_first_index[i] : (person_first_index[i+1] - 1)]
        # }
        if (p == 1) {
            M_matrix_i <- sum(M_matrix[, person_first_index[i]:(person_first_index[i + 1] - 1)])
        } else {
            M_matrix_i <- rowSums(M_matrix[, person_first_index[i]:(person_first_index[i + 1] - 1)])
        }
        Sigman_summand[i, , ] <- M_matrix_i %*% t(M_matrix_i)
    }
    Sigman <- apply(Sigman_summand, c(2, 3), sum) / sample_size
    
    varcov <- Mn_inv %*% Sigman %*% t(Mn_inv) / sample_size
    beta_se <- sqrt(diag(varcov)[1:p])
    
    
    ### 4. small sample correction ###
    r_term_collected <- rep(NA, total_person_decisionpoint)
    D_term_collected <- matrix(NA, nrow = p, ncol = total_person_decisionpoint)
    partialr_partialtheta_collected <- matrix(NA, nrow = total_person_decisionpoint, ncol = p)
    # Mn_summand_small <- array(NA, dim = c(total_person_decisionpoint, p, p))
    
    for (it in 1:total_person_decisionpoint) {
        # this is to make R code consistent whether X_it, Z_it contains more entries or is just the intercept.
        if (p == 1) {
            Xbeta <- Xdm[it, ] * beta_hat
        } else {
            Xbeta <- as.numeric(Xdm[it, ] %*% beta_hat)
        }
        
        r_term1 <- exp(-1 * A[it] * Xbeta) * M[it] * cA_tilde[it] * M1_error[it]
        r_term2 <- p_t_tilde[it] * (1 - p_t_tilde[it]) * (M2_multiplier[it] * E_YWit_given_HAit_Aeq1[it])
        r_term3 <- p_t_tilde[it] * (1 - p_t_tilde[it]) * (-E_YWit_given_HAit_Aeq0[it])
        # r_term = r^(t) (scalar)
        r_term <- r_term1 + r_term2 + r_term3
        r_term_collected[it] <- r_term
        
        # D_term = D^{(t),T} (dim = p * 1)
        D_term <- as.matrix(avail[it] * Xdm[it, ])
        D_term_collected[, it] <- D_term
        
        # partialr_partialtheta = \frac{\partial r^(t)}{\partial \theta^T}
        partialr_partialtheta <- r_term1 %*% (-A[it] * t(Xdm[it, ])) + r_term2 %*% (-1 * t(Xdm[it, ]))
        partialr_partialtheta_collected[it, ] <- partialr_partialtheta
        
        # Mn_summand_small[it, , ] <- D_term %*% partialr_partialtheta
    }
    
    Sigman_tilde <- 0
    for (i in 1:sample_size) {
        D_term_i <- D_term_collected[, person_first_index[i]:(person_first_index[i + 1] - 1)]
        if (p == 1) {
            D_term_i <- t(D_term_i)
        }
        r_term_i <- matrix(r_term_collected[person_first_index[i]:(person_first_index[i + 1] - 1)], ncol = 1)
        partialr_partialtheta_i <- partialr_partialtheta_collected[person_first_index[i]:(person_first_index[i + 1] - 1), ]
        H_ii <- partialr_partialtheta_i %*% Mn_inv %*% D_term_i / sample_size
        Ii_minus_Hii_inv <- solve(diag(nrow(H_ii)) - H_ii)
        
        Sigman_tilde <- Sigman_tilde + D_term_i %*% Ii_minus_Hii_inv %*% r_term_i %*% t(r_term_i) %*% t(Ii_minus_Hii_inv) %*% t(D_term_i)
    }
    Sigman_tilde <- Sigman_tilde / sample_size
    
    varcov_adjusted <- Mn_inv %*% Sigman_tilde %*% t(Mn_inv) / sample_size
    beta_se_adjusted <- sqrt(diag(varcov_adjusted)[1:p])
    
    ### 5. return the result with variable names ###
    return(list(
        beta_hat = beta_hat,
        beta_se = beta_se, beta_se_adjusted = beta_se_adjusted,
        varcov = varcov, varcov_adjusted = varcov_adjusted,
        dims = list(p = p),
        f.root = solution$f.root
    ))
}



compute_excursion_ipw <- function(
        dta,
        id_varname,
        decision_time_varname,
        treatment_varname,
        outcome_varname,
        rand_prob_varname,
        avail_varname = NULL,
        Delta
) {
    # This is copy-pasted from Yihan's function weighted_centered_least_square_withDelta_improved_nonmem().
    
    sample_size <- length(unique(dta[[id_varname]]))
    total_person_decisionpoint <- nrow(dta)
    id_names <- unique(dta[[id_varname]])
    
    if (is.null(avail_varname)) {
        avail <- rep(1, total_person_decisionpoint)
    } else {
        avail <- dta[[avail_varname]]
    }
    
    A <- dta[, treatment_varname]
    # checking for NA in treatment indicator
    if (any(is.na(A[avail == 1]))) {
        stop("Treatment indicator is NA where availability = 1.")
    }
    A[avail == 0] <- 0
    
    p_t <- dta[[rand_prob_varname]]
    Y <- dta[[outcome_varname]]
    
    excursion_ipw <- c()
    if (Delta == 1) {
        excursion_ipw <- rep(1, total_person_decisionpoint)
    } else {
        for (i in 1:sample_size) {
            dta_perid <- dta[dta[[id_varname]] == id_names[i], ]
            dta_timepoint <- length(dta_perid[, decision_time_varname])
            
            # append with A_{T+1} ... A_{T+\delta-1} = 0
            inverseProb_numerator <- 1 - c(dta_perid[, treatment_varname], rep(0, Delta - 1))
            # append with p_{T+1} ... p_{T+\delta-1} = 0
            inverseProb_denominator <- 1 - c(dta_perid[, rand_prob_varname], rep(0, Delta - 1))
            
            for (j in 1:dta_timepoint) {
                # the first time R_{t+k} = 1
                first_R_1 <- dta_perid$k[j]
                
                if (first_R_1 == 0) {
                    # if R_t = 1 right after A_t, set the weight to 1
                    excursion_ipw <- c(excursion_ipw, 1)
                } else {
                    # in other words, if A = 0 appears before k, we need to set it to 0, o.w keep this term
                    excursion_ipw <- c(
                        excursion_ipw,
                        prod(inverseProb_numerator[(j + 1):(j + first_R_1)] /
                                 inverseProb_denominator[(j + 1):(j + first_R_1)])
                    )
                }
            }
        }
    }
    return(excursion_ipw)
}



build_nuisance_predictions <- function(
        dta,
        id_varname,
        decision_time_varname,
        treatment_varname,
        outcome_varname,
        control_varname,
        rand_prob_varname,
        avail_varname = NULL,
        Delta,
        excursion_ipw_varname, # the variable name
        link_function = c("binomial", "poisson")
) {
    # This does not take out the measurable part of the weight term in fitting
    # the conditional expectation E(Y_{it,\Delta} W_{it} | H_{iu}, A_{iu}) and E(Y_{it,\Delta} W_{it} | H_{iu})
    
    # Predict the following quantities for i = 1, ..., n and t = 1, ..., T
    # E(Y_{it,\Delta} W_{it} | H_{it}, A_{it})
    # E(Y_{it,\Delta} W_{it} | H_{it}, A_{it} = 1)
    # E(Y_{it,\Delta} W_{it} | H_{it}, A_{it} = 0)
    # E(Y_{it,\Delta} W_{it} | H_{iu}, A_{iu}) for u = t+1, ..., t+\Delta-1
    # E(Y_{it,\Delta} W_{it} | H_{iu}) for u = t+1, ..., t+\Delta-1
    
    link_function <- match.arg(link_function)
    
    sample_size <- length(unique(dta[[id_varname]]))
    total_person_decisionpoint <- nrow(dta)
    
    id_names <- unique(dta[[id_varname]])
    
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
    Zdm <- as.matrix(cbind(rep(1, nrow(dta)), dta[, control_varname])) # Z (control) design matrix, intercept added
    
    q <- ncol(Zdm)
    
    W <- dta[[excursion_ipw_varname]]
    
    ### 1. Fitting E(Y_{it,\Delta} W_{it} | H_{it}, A_{it}),
    ###.   E(Y_{it,\Delta} W_{it} | H_{it}, A_{it} = 1),
    ###.   E(Y_{it,\Delta} W_{it} | H_{it}, A_{it} = 0)
    
    # rescale Y_{it,\Delta} W_{it} so that it stays within [0,1], then use logistic regression
    
    YW <- Y * W
    YW_max <- max(YW)
    YW_rescaled <- YW / YW_max
    
    if (link_function == "binomial") {
        fit <- glm.fit(x = cbind(Zdm, A), y = YW_rescaled, family = binomial())
        fit_coef <- as.matrix(fit$coefficients, ncol = 1)
        
        E_YWit_given_HAit <- expit(cbind(Zdm, A) %*% fit_coef) * YW_max
        E_YWit_given_HAit_Aeq1 <- expit(cbind(Zdm, 1) %*% fit_coef) * YW_max
        E_YWit_given_HAit_Aeq0 <- expit(cbind(Zdm, 0) %*% fit_coef) * YW_max
    } else if (link_function == "poisson") {
        fit <- glm.fit(x = cbind(Zdm, A), y = YW, family = poisson())
        fit_coef <- as.matrix(fit$coefficients, ncol = 1)
        
        E_YWit_given_HAit <- exp(cbind(Zdm, A) %*% fit_coef)
        E_YWit_given_HAit_Aeq1 <- exp(cbind(Zdm, 1) %*% fit_coef)
        E_YWit_given_HAit_Aeq0 <- exp(cbind(Zdm, 0) %*% fit_coef)
    }
    
    if (Delta >= 2) {
        
        ### 2. Fitting  E(Y_{it,\Delta} W_{it} | H_{iu}, A_{iu}) for u = t+1, ..., t+\Delta-1
        
        # The prediction will be a list of \Delta-1 vectors:
        # one vector for each of lead = 1, 2, ..., \Delta-1
        
        person_first_index <- c(find_change_location(dta[[id_varname]])[1:sample_size], total_person_decisionpoint + 1)
        feature_mat <- cbind(Zdm, A)
        
        for (this_lead in 1:(Delta-1)) {
            for (i in 1:sample_size) {
                feature_mat_i <- feature_mat[person_first_index[i]:(person_first_index[i + 1] - 1), ]
                
                # fill in the last this_lead rows by LOCF for Z and 0 for A
                last_row_filler <- matrix(c(feature_mat_i[nrow(feature_mat_i), 1:(ncol(feature_mat_i)-1)], 0), nrow = 1)
                
                # the first rows of feature_mat_i_led are led versions of feature_mat_i
                # the last this_lead rows are 
                feature_mat_i_led <- feature_mat_i[(1+this_lead):nrow(feature_mat_i), ]
                for (j in 1:this_lead) {
                    feature_mat_i_led <- rbind(feature_mat_i_led, last_row_filler)
                }
                if (i == 1) {
                    feature_mat_led <- feature_mat_i_led
                } else {
                    feature_mat_led <- rbind(feature_mat_led, feature_mat_i_led)
                }
            }
            
            if (link_function == "binomial") {
                fit <- glm.fit(x = feature_mat_led, y = YW_rescaled, family = binomial())
                fit_coef <- as.matrix(fit$coefficients, ncol = 1)
                
                E_YWit_given_HAiu_this_lead <- expit(feature_mat_led %*% fit_coef) * YW_max
            } else if (link_function == "poisson") {
                fit <- glm.fit(x = feature_mat_led, y = YW, family = poisson())
                fit_coef <- as.matrix(fit$coefficients, ncol = 1)
                
                E_YWit_given_HAiu_this_lead <- exp(feature_mat_led %*% fit_coef)
            }
            
            
            if (this_lead == 1) {
                E_YWit_given_HAiu <- E_YWit_given_HAiu_this_lead
            } else {
                E_YWit_given_HAiu <- cbind(E_YWit_given_HAiu, E_YWit_given_HAiu_this_lead)
            }
        }
        
        colnames(E_YWit_given_HAiu) <- paste0("lead=", 1:(Delta-1))
        
        ### 3. Fitting  E(Y_{it,\Delta} W_{it} | H_{iu}) for u = t+1, ..., t+\Delta-1
        
        # The prediction will be a matrix of total_person_decisionpoint rows and \Delta-1 columns:
        # one vector for each of lead = 1, 2, ..., \Delta-1
        
        person_first_index <- c(find_change_location(dta[[id_varname]])[1:sample_size], total_person_decisionpoint + 1)
        feature_mat <- Zdm
        
        
        for (this_lead in 1:(Delta-1)) {
            for (i in 1:sample_size) {
                feature_mat_i <- matrix(feature_mat[person_first_index[i]:(person_first_index[i + 1] - 1), ],
                                        ncol = q)
                # fill in the last this_lead rows by LOCF for Z
                last_row_filler <- matrix(feature_mat_i[nrow(feature_mat_i), ], nrow = 1)
                
                # the first rows of feature_mat_i_led are led versions of feature_mat_i
                # the last this_lead rows are 
                feature_mat_i_led <- matrix(feature_mat_i[(1+this_lead):nrow(feature_mat_i), ], ncol = q)
                for (j in 1:this_lead) {
                    feature_mat_i_led <- rbind(feature_mat_i_led, last_row_filler)
                }
                if (i == 1) {
                    feature_mat_led <- feature_mat_i_led
                } else {
                    feature_mat_led <- rbind(feature_mat_led, feature_mat_i_led)
                }
            }
            
            if (link_function == "binomial") {
                fit <- glm.fit(x = feature_mat_led, y = YW_rescaled, family = binomial())
                fit_coef <- as.matrix(fit$coefficients, ncol = 1)
                
                E_YWit_given_Hiu_this_lead <- expit(feature_mat_led %*% fit_coef) * YW_max
            } else if (link_function == "poisson") {
                fit <- glm.fit(x = feature_mat_led, y = YW, family = poisson())
                fit_coef <- as.matrix(fit$coefficients, ncol = 1)
                
                E_YWit_given_Hiu_this_lead <- exp(feature_mat_led %*% fit_coef)
            }
            
            
            if (this_lead == 1) {
                E_YWit_given_Hiu <- E_YWit_given_Hiu_this_lead
            } else {
                E_YWit_given_Hiu <- cbind(E_YWit_given_Hiu, E_YWit_given_Hiu_this_lead)
            }
        }
        
        colnames(E_YWit_given_Hiu) <- paste0("lead=", 1:(Delta-1))
        
    } else {
        # Delta == 1
        E_YWit_given_HAiu <- NA
        E_YWit_given_Hiu <- NA
    }
    
    return(list(E_YWit_given_HAit = E_YWit_given_HAit,
                E_YWit_given_HAit_Aeq1 = E_YWit_given_HAit_Aeq1,
                E_YWit_given_HAit_Aeq0 = E_YWit_given_HAit_Aeq0,
                E_YWit_given_HAiu = E_YWit_given_HAiu,
                E_YWit_given_Hiu = E_YWit_given_Hiu))
}


