# Tianchen Qian
# 2024.03.02

rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(tidyverse)

# publicly available HeartStepsV1 data downloaded from https://github.com/klasnja/HeartStepsV1

jbslot <- read.csv("HeartStepsV1-main/data_files/jbsteps.csv")
gfslot <- read.csv("HeartStepsV1-main/data_files/gfsteps.csv")
users <- read.csv("HeartStepsV1-main/data_files/users.csv")
suggest <- read.csv("HeartStepsV1-main/data_files/suggestions.csv")

CREATE_PLOTS <- FALSE # if TRUE, will make sanity check plots
PRINT_SANITY_CHECK <- TRUE # if TRUE, will print sanity check results

jbslot <- as_tibble(jbslot)

suggest <- as_tibble(suggest)

# remove the baseline covariates (decision.index = NA) -------------------------------------------
jbslot <- jbslot %>% filter(!is.na(decision.index))
# the publicly available dataset has baseline information recorded (i.e. at decision.index = NA)
# we will remove these to keep the form of the datasets to be consistent with our analysis.


# Remove travel dates in suggest
var_from_suggest <- c("user.index", "decision.index.nogap",  "decision.index", "avail", "send",
                      "send.active", "send.sedentary", "sugg.decision.utime", 
                      "dec.dec.location.category", "jbsteps30pre", "jbsteps30pre.zero")
suggest_selectedvar <- suggest[, var_from_suggest]

suggest_selectedvar <- suggest_selectedvar %>% 
    mutate(sugg.decision.udate = date(ymd_hms(suggest_selectedvar$sugg.decision.utime)))
tmp_df <- data.frame()

for (i in 1:nrow(users)) {
    
    user_i_dta <- filter(suggest_selectedvar, user.index == users$user.index[i])
    if (users$travel.start[i] != "") {
        nrow_before <- nrow(user_i_dta)
        user_i_dta <- filter(user_i_dta, (sugg.decision.udate < users$travel.start[i]) | (sugg.decision.udate > users$travel.end[i]))
        nrow_after <- nrow(user_i_dta)
        cat(paste0("user ", users$user.index[i], ": travel.start ", users$travel.start[i], ", travel.end ", users$travel.end[i]), "\n")
        cat(paste0("        nrow_before ", nrow_before, ", nrow_after ", nrow_after, ", nrow_removed ",
                   nrow_before - nrow_after, "\n"))
    }
    tmp_df <- rbind(tmp_df, user_i_dta)
}
suggest_tdr <- tmp_df

# user 1: travel.start 2015-08-12, travel.end 2015-08-31 
# nrow_before 278, nrow_after 178, nrow_removed 100
# user 3: travel.start 2015-08-13, travel.end 2015-08-20 
# nrow_before 255, nrow_after 215, nrow_removed 40
# user 6: travel.start 2015-08-10, travel.end 2015-08-15 
# nrow_before 212, nrow_after 182, nrow_removed 30
# user 13: travel.start 2015-09-01, travel.end 2015-09-05 
# nrow_before 215, nrow_after 192, nrow_removed 23
# user 14: travel.start 2015-10-12, travel.end 2015-10-22 
# nrow_before 265, nrow_after 210, nrow_removed 55
# user 16: travel.start 2015-09-20, travel.end 2015-09-22 
# nrow_before 214, nrow_after 199, nrow_removed 15
# user 31: travel.start 2015-12-15, travel.end 2016-01-08 
# nrow_before 309, nrow_after 184, nrow_removed 125

# Construct a new decision.index.nogap variable in suggest dataset, then paste it into jbslot -------------------------------------------
## Issue: Why do we need to do this? b/c somehow there are observations that have decision.index.nogap skipped
## e.g. user 13's orginal decision.index.nogap has labeled wrong at j = 23, 24 (they got skipped)
## Solution: construct a new decision.index.nogap variable

suggest_tdr <- suggest_tdr %>% 
    group_by(user.index) %>% 
    mutate(max.decision.index.nogap = n(),
           decision.index.nogap.new = 0:(unique(max.decision.index.nogap) - 1))
a <- which(suggest_tdr$decision.index.nogap != suggest_tdr$decision.index.nogap.new)
b <- unique(suggest_tdr[a, ]) 

# Paste the newly constructed decision.index.nogap.new into jbslot
jbslot <- suggest_tdr %>% 
    dplyr::select("user.index", "decision.index", "decision.index.nogap.new", 
                  "sugg.decision.utime", "sugg.decision.udate") %>% 
    right_join(jbslot, by = c("user.index", "decision.index"))


# remove travel dates in jbslot-------------------------------------------

## Issue: At i  = 1, j = 203, the suggest dataset has notification sent when they are in travel dates. 
## But the steps are recorded. And the steps recorded time is not in travel dates (the morning they got back). 
## So there are NAs when joining these two datasets. We can remove these observations since their notifications are sent at travel dates.
## Solution: we use the sugg.decision.udate from suggest dataset (to filter the observations in travel dates) to avoid discrepancy between suggest and jbslot.

# travel day is stored in "users" data.frame, with variable names "travel.start" and "travel.end"
tmp_df <- data.frame()
for (i in 1:nrow(users)) {
    
    user_i_dta <- filter(jbslot, user.index == users$user.index[i])
    if (users$travel.start[i] != "") {
        nrow_before <- nrow(user_i_dta)
        user_i_dta <- filter(user_i_dta, (sugg.decision.udate < users$travel.start[i]) | (sugg.decision.udate > users$travel.end[i]))
        nrow_after <- nrow(user_i_dta)
        cat(paste0("user ", users$user.index[i], ": travel.start ", users$travel.start[i], ", travel.end ", users$travel.end[i]), "\n")
        cat(paste0("        nrow_before ", nrow_before, ", nrow_after ", nrow_after, ", nrow_removed ",
                   nrow_before - nrow_after, "\n"))
    }
    tmp_df <- rbind(tmp_df, user_i_dta)
}
jbslot_tdr <- tmp_df

# user 1: travel.start 2015-08-12, travel.end 2015-08-31 
# nrow_before 7965, nrow_after 7542, nrow_removed 423
# user 3: travel.start 2015-08-13, travel.end 2015-08-20 
# nrow_before 4611, nrow_after 4441, nrow_removed 170
# user 6: travel.start 2015-08-10, travel.end 2015-08-15 
# nrow_before 7364, nrow_after 7094, nrow_removed 270
# user 13: travel.start 2015-09-01, travel.end 2015-09-05 
# nrow_before 3389, nrow_after 2866, nrow_removed 523
# user 14: travel.start 2015-10-12, travel.end 2015-10-22 
# nrow_before 9199, nrow_after 7828, nrow_removed 1371
# user 16: travel.start 2015-09-20, travel.end 2015-09-22 
# nrow_before 7041, nrow_after 6496, nrow_removed 545
# user 31: travel.start 2015-12-15, travel.end 2016-01-08 
# nrow_before 5357, nrow_after 5108, nrow_removed 249


#  construct min.after.decision variable in the dataset -------------------------------------------
jbslot_tdr <- jbslot_tdr %>% 
    mutate(steps.utime = ymd_hms(steps.utime),
           sugg.decision.utime = ymd_hms(sugg.decision.utime),
           min.after.decision = floor(as.numeric(difftime(steps.utime, sugg.decision.utime, units = 'mins'))) + 1
    )

jbslot_tdr <- jbslot_tdr[order(jbslot_tdr$user.index, 
                               jbslot_tdr$decision.index.nogap.new, 
                               jbslot_tdr$min.after.decision), ]


if (PRINT_SANITY_CHECK) {
    # make sure primary keys are unique
    jbslot_tdr %>% 
        count(user.index, decision.index.nogap.new, min.after.decision) %>% 
        filter(n > 1) # no duplicate!
}

if (PRINT_SANITY_CHECK) {
    # make sure all user has decision.index.nogap from 0 to max consecutively, no gap in between
    print(jbslot_tdr %>% group_by(user.index) %>%
              summarise(n_dp = length(decision.index.nogap.new),
                        dp_min = min(decision.index.nogap.new),
                        dp_max = max(decision.index.nogap.new),
                        no_gap_in_dp = ((dp_max - dp_min + 1) == n_dp)), n = Inf)
}

#Note: User.index = 4 has no observation recorded at decision.index 0, because his/her connection failed.
# # A tibble: 37 × 5
#   user.index  n_dp dp_min dp_max no_gap_in_dp
# <int> <int>  <int>  <int> <lgl>       
# 1          1  7542      0    169 FALSE       
# 2          2  5020      0    207 FALSE       
# 3          3  4441      0    213 FALSE       
# 4          4 11715      1    217 FALSE       
# 5          5  4202      0    215 FALSE       
# 6          6  7094      0    181 FALSE       
# 7          7  6841      0    214 FALSE       
# 8          8  7898      0    219 FALSE       
# 9          9  5296      0    205 FALSE       
# 10         10  6937      0    213 FALSE       
# 11         11  6488      0    216 FALSE       
# 12         12  7744      0    242 FALSE       
# 13         13  2866      0    184 FALSE       
# 14         14  7828      0    208 FALSE       
# 15         15  7426      0    254 FALSE       
# 16         16  6496      0    197 FALSE       
# 17         17  7618      0    208 FALSE       
# 18         18  5226      0    210 FALSE       
# 19         19  6168      0    208 FALSE       
# 20         20  4779      0    233 FALSE       
# 21         21  8558      0    233 FALSE       
# 22         22  4733      0    213 FALSE       
# 23         23  5198      0    223 FALSE       
# 24         24  6215      0    173 FALSE       
# 25         25 10346      0    222 FALSE       
# 26         26  6289      0    211 FALSE       
# 27         27  3897      0    210 FALSE       
# 28         28  8327      0    208 FALSE       
# 29         29  2575      0     83 FALSE       
# 30         30  7494      0    208 FALSE       
# 31         31  5108      0    182 FALSE       
# 32         32  3922      0    223 FALSE       
# 33         33  7813      0    226 FALSE       
# 34         34  3039      0    211 FALSE       
# 35         35  5797      0    208 FALSE       
# 36         36  9155      0    201 FALSE       
# 37         37  5451      0    228 FALSE       

# Construct new template data set with "user", "decision.index.nogap", "min.after.decision" -----------------------------------
## make a combined keyword, then separate out (because expand.grid only accepts vector arguments)
user.di_unique <- unite(suggest_tdr, "user.di", "user.index", "decision.index.nogap.new")$user.di

dta_template <- expand.grid(user.di = user.di_unique, min.after.decision = 0:(total_T - 1)) %>%
    separate(user.di, into = c("user", "decision.index.nogap")) %>%
    arrange(user, decision.index.nogap, min.after.decision)
dta_template$user <- as.numeric(dta_template$user)
dta_template$decision.index.nogap <- as.numeric(dta_template$decision.index.nogap)
dta_template <- as_tibble(dta_template)

## sanity checks on the number of decision points and the number of minutes
if (PRINT_SANITY_CHECK) {
    summary(dta_template)
    print(dta_template %>% group_by(user) %>% summarise(n.dp = length(decision.index.nogap) / total_T), n = Inf)
    summary(dta_template %>% group_by(user, decision.index.nogap) %>%
                summarise(n.minute = length(min.after.decision),
                          minute.max = max(min.after.decision),
                          minute.min = min(min.after.decision),
                          n.unique.minute = length(unique(min.after.decision))))
}

# # A tibble: 37 × 2
# user  n.dp
# <dbl> <dbl>
#   1     1   178
# 2     2   209
# 3     3   215
# 4     4   219
# 5     5   217
# 6     6   182
# 7     7   216
# 8     8   221
# 9     9   207
# 10    10   215
# 11    11   220
# 12    12   244
# 13    13   192
# 14    14   210
# 15    15   256
# 16    16   199
# 17    17   210
# 18    18   211
# 19    19   210
# 20    20   235
# 21    21   235
# 22    22   214
# 23    23   225
# 24    24   175
# 25    25   223
# 26    26   212
# 27    27   211
# 28    28   210
# 29    29   211
# 30    30   210
# 31    31   184
# 32    32   225
# 33    33   228
# 34    34   212
# 35    35   213
# 36    36   202
# 37    37   230

which(suggest_tdr$sugg.decision.utime == "") # 603
print(suggest_tdr[600:630, "sugg.decision.utime"], n = Inf)
# decision points are around 11:00, 15:30, 18:30, 20:30, 23:00
# Therefore, we impute the missing utime for user4 by "2015-07-27 15:31:00"
suggest_tdr[603, "sugg.decision.utime"] <- "2015-07-27 15:31:00"

suggest_tdr$sugg.decision.utime <- as.POSIXct(suggest_tdr$sugg.decision.utime, tz = "UTC")

suggest_tdr$send <- as.numeric(suggest_tdr$send == "True")
suggest_tdr$send.active <- as.numeric(suggest_tdr$send.active == "True")
suggest_tdr$send.sedentary <- as.numeric(suggest_tdr$send.sedentary == "True")
suggest_tdr$avail <- as.numeric(suggest_tdr$avail == "True")


# Following is new code by Tianchen on 2023.12.01

# Input: Delta (number of minutes, the same number as above), a step count threshold

# construct a decision-point-level data set for analysis
    # For each decision point, determine whether the step count threshold is reached within the Delta window.
        # If no, the outcome is 0, the time window to compute the pd-IPW weight is the Delta window.
        # If yes, the outcome is 1, the time window to compute the pd-IPW weight is the first minute that the total step count reached the threshold.
    # Based on the pd-IPW weight window, compute the weight using treatment and availability in that window.

unique_users <- unique(suggest_tdr$user.index)
suggest_tdr_by_user <- c()
jbslot_tdr_by_user <- c()
for (this_user in unique_users) {
    suggest_tdr_by_user <- c(suggest_tdr_by_user,
                             list(filter(suggest_tdr, user.index == this_user)))
    jbslot_tdr_by_user <- c(jbslot_tdr_by_user,
                            list(filter(jbslot_tdr, user.index == this_user)))
}

delta_window_minutes <- 60 * 24 - 5 # 1440-5 minutes (24 hours - 5min)
# subtracting 5 minutes to avoid potentially including the decision point in exactly 24 hours
# and thus messing with the ipw_pd computation.


thresholds <- seq(from = 5000, to = 12000, by = 1000)

for (step_count_threshold in thresholds) {
    # Construct the following columns to suggest_tdr
    suggest_tdr$step_count_in_window <- NA
    suggest_tdr$step_reaching_threshold <- NA
    suggest_tdr$ipw_pd <- NA
    suggest_tdr$ipw <- NA
    
    for (irow in 1:nrow(suggest_tdr)) {
        if (irow %% 1000 == 0) {
            print(irow)
        }
        
        user_index <- which(unique_users == suggest_tdr$user.index[irow])
        this_dp_nogap <- suggest_tdr$decision.index.nogap[irow]
        this_dp_utime <- suggest_tdr$sugg.decision.utime[irow]
        
        this_dp_delta_window_end <- this_dp_utime + 60 * delta_window_minutes
        
        ### 1. compute step_count_in_window, step_reaching_threshold, ipw_pd_window_endtime
        
        jbslot_this_user_and_dp_window <- jbslot_tdr_by_user[[user_index]] %>%
            filter(steps.utime >= this_dp_utime,
                   steps.utime < this_dp_delta_window_end)
        
        if (nrow(jbslot_this_user_and_dp_window) == 0) {
            step_count_in_window <- 0
            step_reaching_threshold <- 0
            ipw_pd_window_endtime <- this_dp_delta_window_end
        } else {
            step_count_in_window <- sum(jbslot_this_user_and_dp_window$steps)
            step_reaching_threshold <- as.numeric(step_count_in_window >= step_count_threshold)
            step_count_cumsum <- cumsum(jbslot_this_user_and_dp_window$steps)
            ipw_pd_window_endtime <- jbslot_this_user_and_dp_window$steps.utime[which(step_count_cumsum >= step_count_threshold)[1]]
        }
        
        ### 2. compute ipw_pd
        
        suggest_this_user_and_pdipw_window <- suggest_tdr_by_user[[user_index]] %>%
            filter(sugg.decision.utime > this_dp_utime + 60,
                   sugg.decision.utime <= ipw_pd_window_endtime)
        A <- as.numeric(suggest_this_user_and_pdipw_window$send)
        avail <- as.numeric(suggest_this_user_and_pdipw_window$avail)
        rand_prob <- ifelse(avail, 0.6, 0)
        ipw_pd <- prod((1-A) / (1-rand_prob))
        
        ### 3. compute ipw
        
        suggest_this_user_and_dp_window <- suggest_tdr_by_user[[user_index]] %>%
            filter(sugg.decision.utime > this_dp_utime + 60,
                   sugg.decision.utime <= this_dp_delta_window_end)
        A <- as.numeric(suggest_this_user_and_dp_window$send)
        avail <- as.numeric(suggest_this_user_and_dp_window$avail)
        rand_prob <- ifelse(avail, 0.6, 0)
        ipw <- prod((1-A) / (1-rand_prob))
        
        ### 4. assign variable values to suggest_tdr
        suggest_tdr$step_count_in_window[irow] <- step_count_in_window
        suggest_tdr$step_reaching_threshold[irow] <- step_reaching_threshold
        suggest_tdr$ipw_pd[irow] <- ipw_pd
        suggest_tdr$ipw[irow] <- ipw
    }
    
    dir.create("data_processed", showWarnings = FALSE)
    saveRDS(suggest_tdr, file = paste0("data_processed/threshold=", step_count_threshold, ".RDS"))
}




# 2. Analysis ----------------------------------------------------------------


library(rootSolve)
source("HeartSteps_public/EMEE_userinput_ipw.R")
source("estimator_implementation/improve_efficiency.R")

moderators <- c("NULL", "decision.index.nogap.new", "at_home_or_work", "is_weekday")

controls <- c("decision.index.nogap.new", "jbsteps30pre.log",
                     "at_home_or_work", "is_weekday")
plot_est_collected <- c()
plot_re_collected <- c()

for (moderator in moderators) {
    
    if (moderator == "NULL") {
        moderator <- NULL
    }
    
    thresholds <- seq(from = 5000, to = 12000, by = 1000)
    
    beta0_ipw <- data.frame(threshold = thresholds, est = NA, se = NA)
    beta0_pdipw <- beta1_ipw <- beta1_pdipw <- beta0_ipw
    
    ### conduct analysis
    
    for (step_count_threshold in thresholds) {
        
        suggest_tdr <- readRDS(paste0("HeartSteps_public/data_processed/threshold=", step_count_threshold, ".RDS"))
        suggest_tdr$rand_prob <- 0.6
        suggest_tdr$at_home_or_work <- as.numeric(suggest_tdr$dec.location.category %in% c("home", "work"))
        
        suggest_tdr$is_weekday <- as.numeric(weekdays(suggest_tdr$sugg.decision.udate) %in%
                                                 c("Monday", "Tuesday", "Wednesday", "Thursday", "Friday"))
        suggest_tdr$jbsteps30pre.log <- log(suggest_tdr$jbsteps30pre + 0.5)
        suggest_tdr$jbsteps30pre.log[is.na(suggest_tdr$jbsteps30pre.log)] <- log(0.5)
        
        fit_ipw <- EMEE_userinput_ipw(
            dta = suggest_tdr,
            id_varname = "user.index",
            decision_time_varname = "decision.index.nogap.new",
            treatment_varname = "send",
            outcome_varname = "step_reaching_threshold",
            control_varname = controls,
            moderator_varname = moderator,
            rand_prob_varname = "rand_prob",
            avail_varname = "avail",
            rand_prob_tilde_varname = NULL, # \tilde{p}_t(1|H_t) in WCLS (variable name in the data set)
            rand_prob_tilde = 0.6,         # \tilde{p}_t(1|H_t) in WCLS (numeric number or vector)
            estimator_initial_value = NULL,
            excursion_ipw_varname = "ipw" # the variable name 
        )
        
        fit_pdipw <- EMEE_userinput_ipw(
            dta = suggest_tdr,
            id_varname = "user.index",
            decision_time_varname = "decision.index.nogap.new",
            treatment_varname = "send",
            outcome_varname = "step_reaching_threshold",
            control_varname = controls,
            moderator_varname = moderator,
            rand_prob_varname = "rand_prob",
            avail_varname = "avail",
            rand_prob_tilde_varname = NULL, # \tilde{p}_t(1|H_t) in WCLS (variable name in the data set)
            rand_prob_tilde = 0.6,         # \tilde{p}_t(1|H_t) in WCLS (numeric number or vector)
            estimator_initial_value = NULL,
            excursion_ipw_varname = "ipw_pd" # the variable name 
        )
        
        fit_pd2ipw <- weighted_centered_least_square_withDelta_improved(
          dta = suggest_tdr,
          id_varname = "user.index",
          decision_time_varname = "decision.index.nogap.new",
          treatment_varname = "send",
          outcome_varname = "step_reaching_threshold",
          control_varname = controls,
          moderator_varname = moderator,
          rand_prob_varname = "rand_prob",
          avail_varname = "avail",
          rand_prob_tilde_varname = NULL, # \tilde{p}_t(1|H_t) in WCLS (variable name in the data set)
          rand_prob_tilde = 0.6,         # \tilde{p}_t(1|H_t) in WCLS (numeric number or vector)
          estimator_initial_value = NULL
        )
        
        if (is.null(moderator)) {
            beta0_ipw[beta0_ipw$threshold == step_count_threshold, "est"] <- fit_ipw$beta_hat
            beta0_ipw[beta0_ipw$threshold == step_count_threshold, "se"] <- fit_ipw$beta_se_adjusted
            beta0_pdipw[beta0_pdipw$threshold == step_count_threshold, "est"] <- fit_pdipw$beta_hat
            beta0_pdipw[beta0_pdipw$threshold == step_count_threshold, "se"] <- fit_pdipw$beta_se_adjusted
        } else {
            beta0_ipw[beta0_ipw$threshold == step_count_threshold, "est"] <- fit_ipw$beta_hat[1]
            beta0_ipw[beta0_ipw$threshold == step_count_threshold, "se"] <- fit_ipw$beta_se_adjusted[1]
            beta0_pdipw[beta0_pdipw$threshold == step_count_threshold, "est"] <- fit_pdipw$beta_hat[1]
            beta0_pdipw[beta0_pdipw$threshold == step_count_threshold, "se"] <- fit_pdipw$beta_se_adjusted[1]
            beta1_ipw[beta1_ipw$threshold == step_count_threshold, "est"] <- fit_ipw$beta_hat[2]
            beta1_ipw[beta1_ipw$threshold == step_count_threshold, "se"] <- fit_ipw$beta_se_adjusted[2]
            beta1_pdipw[beta1_pdipw$threshold == step_count_threshold, "est"] <- fit_pdipw$beta_hat[2]
            beta1_pdipw[beta1_pdipw$threshold == step_count_threshold, "se"] <- fit_pdipw$beta_se_adjusted[2]
        }
    }

    ### make plots
    
    beta0_ipw <- beta0_ipw %>% 
        mutate(rci = est + 1.96 * se,
               lci = est - 1.96 * se,
               method = "IPW", estimand = "beta0")
    beta0_ipw$re <- 1
    beta0_pdipw <- beta0_pdipw %>% 
        mutate(rci = est + 1.96 * se,
               lci = est - 1.96 * se,
               method = "pd-IPW", estimand = "beta0")
    beta0_pdipw$re <- (beta0_ipw$se / beta0_pdipw$se)^2
    
    if (!is.null(moderator)) {
        beta1_ipw <- beta1_ipw %>% 
            mutate(rci = est + 1.96 * se,
                   lci = est - 1.96 * se,
                   method = "IPW", estimand = "beta1")
        beta1_ipw$re <- 1
        beta1_pdipw <- beta1_pdipw %>% 
            mutate(rci = est + 1.96 * se,
                   lci = est - 1.96 * se,
                   method = "pd-IPW", estimand = "beta1")
        beta1_pdipw$re <- (beta1_ipw$se / beta1_pdipw$se)^2
        
        beta_collected <- rbind(beta0_ipw, beta0_pdipw, beta1_ipw, beta1_pdipw)
    } else {
        beta_collected <- rbind(beta0_ipw, beta0_pdipw)
    }
    
    if (is.null(moderator)) {
        moderator_text <- "NULL"
    } else {
        moderator_text <- moderator
    }
    
    dir.create("analysis_result/", showWarnings = FALSE)
    saveRDS(beta_collected, file = paste0("analysis_result/moderator=", moderator_text, ".RDS"))
}

