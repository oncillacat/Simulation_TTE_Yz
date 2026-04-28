# =============================================================================
# pipeline.R  —  consolidated & bootstrap-ready version of your 5 Rmds
# =============================================================================

# ── STAGE 0: Setup ────────────────────────────────────────────────────────────

library(data.table)
library(dplyr)
library(FNN)
library(fastglm)
library(survival)
library(broom)
library(tidyr)

logit_ps <- function(p) log(p / (1 - p))   # used in PSM caliper

MAX_TIME <- 60    # end of study (month 59 = last time point)
N_BOOT   <- 100   # number of bootstrap iterations


# ── STAGE 1: Data generation ──────────────────────────────────────────────────

generate_data <- function(
    nsim     = 100000,
    K        = 60,
    A_params = c(-1.5, 0, -4.8, 0),       # beta_0_A, beta_L_A, beta_t_A, beta_t2_A
    Y_params = c(-6, 0.01, 0, 0, 0, 0),   # beta_0_Y .. beta_t2A_Y
    seed     = 42
) {
  beta_0_A  <- A_params[1]; beta_L_A  <- A_params[2]
  beta_t_A  <- A_params[3]; beta_t2_A <- A_params[4]
  beta_0_Y  <- Y_params[1]; beta_t_Y  <- Y_params[2]
  beta_t2_Y <- Y_params[3]; beta_A_Y  <- Y_params[4]
  beta_tA_Y <- Y_params[5]; beta_t2A_Y<- Y_params[6]

  dt <- CJ(id = 1:nsim, time = 0:(K - 1))

  indiv <- data.table(id = 1:nsim, L = runif(nsim, 0, 1))
  indiv[, A0 := rbinom(nsim, 1, plogis(beta_0_A + beta_L_A * L))]
  dt <- merge(dt, indiv, by = "id"); rm(indiv)

  dt[, `:=`(t2 = time^2, tstart = time, tstop = time + 1)]
  dt[, Astar := fifelse(A0 == 0 & time > 0,
                        rbinom(.N, 1, plogis(beta_t_A)), 0)]
  dt[, trt_time := ifelse(A0 == 1, 0, which(Astar == 1)[1]), by = id]
  dt[, A := fifelse(A0 == 1, 1,
             fifelse(!is.na(trt_time) & time >= trt_time, 1, 0))]
  dt[, c("Astar", "trt_time") := NULL]
  dt <- dt %>% mutate(across(A, ~replace_na(.x, 0)))

  dt[, Ystar := rbinom(nsim * K, 1, plogis(
    beta_0_Y + beta_t_Y * time + beta_t2_Y * t2 +
      beta_A_Y * A + beta_tA_Y * time * A + beta_t2A_Y * t2 * A))]
  dt[, event_time := which(Ystar == 1)[1], by = id]
  dt <- dt[time < event_time | is.na(event_time)][, event_time := NULL]

  setnames(dt,
           old = c("time", "t2",  "A",   "Ystar"),
           new = c("cal_time", "cal_timesqr", "trt", "outcome"))

  return(dt)
}


# ── STAGE 2: Propensity score matching ────────────────────────────────────────
run_psm <- function(dt) {

  # flag new initiators
  df <- dt %>%
    dplyr:: group_by(id) %>% arrange(cal_time) %>%
    mutate(new_initiator = as.integer(trt == 1 & lag(trt, 1, default = 0) == 0)) %>%
    ungroup()

  # fit PS model on untreated + new initiators
  ps_data <- df %>% filter(trt == 0 | new_initiator == 1) #filter all new iniator 
  ps_model <- glm(new_initiator ~ L + cal_time + cal_timesqr,
                  data = ps_data, family = binomial(link = "logit"))
  ps_data$ps <- predict(ps_model, newdata = ps_data, type = "response")
  df <- df %>% left_join(ps_data %>% select(id, cal_time, ps),
                         by = c("id", "cal_time"))

  # sequential matching loop — unchanged logic, now inside function scope
  matched_list <- list()
  pair_counter <- 0
  #skip to next month condtion
    #* 1. no treated in the month 
    #* 2. no fitted control in the month 
    #* 3. no eligiable vairance in the month
    #* 4. no enough valid treated in the month 
  for (t in sort(unique(df$cal_time))) {
    dat     <- df %>% filter(cal_time == t)
    treated <- dat %>% filter(new_initiator == 1)
    if (nrow(treated) == 0) next
    control <- dat %>% filter(trt == 0, is.na(outcome) | outcome == 0)
    if (nrow(control) == 0) next
    all_ps   <- c(treated$ps, control$ps)
    sd_logit <- sd(logit_ps(all_ps))
    if (is.na(sd_logit) || sd_logit == 0) next
    caliper <- 0.2 * sd_logit    # Austin (2011)

    treated_ps <- as.matrix(treated$ps)
    control_ps <- as.matrix(control$ps)
    nn_result  <- get.knnx(data = control_ps, query = treated_ps, k = 1)
    # 1 to 1 matching 

    matched_ctrl_ps <- as.vector(control_ps[nn_result$nn.index, , drop = TRUE])
    distances       <- abs(as.vector(treated_ps) - matched_ctrl_ps)
    valid <- distances <= caliper
    if (sum(valid) == 0) next

    matched_treated <- treated[valid, ]
    matched_control <- control[nn_result$nn.index[valid], ]
    n_pairs      <- nrow(matched_treated)
    new_pair_ids <- paste0("pair_", (pair_counter + 1):(pair_counter + n_pairs))
    pair_counter <- pair_counter + n_pairs

    matched_treated$pair_id    <- new_pair_ids
    matched_control$pair_id    <- new_pair_ids
    matched_treated$index_date <- t
    matched_control$index_date <- t

    matched_list[[as.character(t)]] <- rbind(matched_treated, matched_control)
  }

  ps_matched <- bind_rows(matched_list)
  return(list(ps_matched = ps_matched,df = df))   # return df and ps_matched 
}


# ── STAGE 3: Cox model ────────────────────────────────────────────────────────
run_cox <- function(ps_matched, df) {

  # build lookups inside the function (no global df needed)
  event_lookup <- df %>%
    dplyr::group_by(id) %>%
    summarise(
      event_time = min(cal_time[outcome == 1], na.rm = TRUE),
      had_event  = any(outcome == 1)
    ) %>%
    mutate(event_time = ifelse(is.infinite(event_time), NA_real_, event_time))

  trt_lookup <- df %>%
    dplyr::group_by(id) %>%
    summarise(trt_start_time = ifelse(
      any(trt == 1), min(cal_time[trt == 1]), NA_real_))

  ps_matched <- ps_matched %>%
    left_join(event_lookup, by = "id") %>%
    left_join(trt_lookup,   by = "id")

  max_time <- MAX_TIME

  ps_matched <- ps_matched %>%
    mutate(
      pp_censor_time = ifelse(trt == 0, trt_start_time, NA_real_),
      end_time = pmin(
        ifelse(had_event, event_time, max_time),
        ifelse(!is.na(pp_censor_time), pp_censor_time, max_time),
        na.rm = TRUE),
      end_ITT    = ifelse(had_event, event_time, max_time),
      follow_up  = end_time - index_date + 1,
      follow_ITT = end_ITT  - index_date + 1,
      outcome_flag = as.integer(
        had_event &
          (trt == 1 | is.na(trt_start_time) | event_time <= trt_start_time)),
      outcome_ITT  = as.integer(had_event)
    )

  # per-protocol Cox
  cox_pp <- coxph(
    Surv(follow_up, outcome_flag) ~ trt + strata(pair_id) + cluster(id),
    data = ps_matched, ties = "efron")

  # ITT Cox
  cox_itt <- coxph(
    Surv(follow_ITT, outcome_ITT) ~ trt + strata(pair_id) + cluster(id),
    data = ps_matched, ties = "efron")

  # KM curves at each time point (for bootstrap CI bands)
  km <- survfit(Surv(follow_up, outcome_flag) ~ trt, data = ps_matched)
  km_summary <- summary(km, times = 0:(max_time - 1))

  return(list(
    hr_pp    = exp(coef(cox_pp)),           # single HR for bootstrap collection
    hr_itt   = exp(coef(cox_itt)),
    cox_pp   = cox_pp,                      # full model if needed
    cox_itt  = cox_itt,
    km_times = km_summary$time,
    km_surv  = km_summary$surv,
    km_strata= km_summary$strata
  ))
}


# ── STAGE 4: IPW ──────────────────────────────────────────────────────────────
run_ipw_boot <- function(boot_dt) {

  dat <- boot_dt   # IMPORTANT: data.table modifies by reference without this
  
  dt_cleaned <- dat[, trt_time := which(trt == 1)[1], by = id]
  dt_cleaned <- dt_cleaned[cal_time <= trt_time -1 | is.na(trt_time) ]

  A         <- dt_cleaned$trt
  mat_num   <- model.matrix(~cal_time + cal_timesqr, dt_cleaned)
  mat_denom <- model.matrix(~L + cal_time + cal_timesqr, dt_cleaned)

  num_mod   <- fastglm(mat_num,   A, family = binomial("logit"))
  denom_mod <- fastglm(mat_denom, A, family = binomial("logit"))

  dt_cleaned[, ipw_num   := predict(num_mod,   mat_num,   type = "response")]
  dt_cleaned[, ipw_denom := predict(denom_mod, mat_denom, type = "response")]
  dt_cleaned[, wt := ifelse(trt == 1,
                      ipw_num / ipw_denom,
                      (1 - ipw_num) / (1 - ipw_denom))]
  df <- merge(dat,dt_cleaned[, .(boot_id,cal_time,wt)],by = c("cal_time","boot_id"), all.x = T, allow.cartesian = T)
  trt_wt <- dt_cleaned[cal_time == trt_time-1, .(boot_id, wt_trt = wt)]
  df <- merge(df, trt_wt, by = "boot_id", all.x = T, allow.cartesian = T)
  df[cal_time > trt_time-1 & is.na(wt), wt := wt_trt]
  df[ , wt_trt := NULL]

  return(df)
}

run_ipw <- function(dt) {
  
  dat <- copy(dt)  
  
  dt_cleaned <- dat[, trt_time := which(trt == 1)[1], by = id]
  dt_cleaned <- dt_cleaned[cal_time <= trt_time -1 | is.na(trt_time) ]
  
  A         <- dt_cleaned$trt
  mat_num   <- model.matrix(~cal_time + cal_timesqr, dt_cleaned)
  mat_denom <- model.matrix(~L + cal_time + cal_timesqr, dt_cleaned)
  
  num_mod   <- fastglm(mat_num,   A, family = binomial("logit"))
  denom_mod <- fastglm(mat_denom, A, family = binomial("logit"))
  
  dt_cleaned[, ipw_num   := predict(num_mod,   mat_num,   type = "response")]
  dt_cleaned[, ipw_denom := predict(denom_mod, mat_denom, type = "response")]
  dt_cleaned[, wt := ifelse(trt == 1,
                            ipw_num / ipw_denom,
                            (1 - ipw_num) / (1 - ipw_denom))]
  df <- merge(dat,dt_cleaned[, .(id,cal_time,wt)],by = c("cal_time","id"), all.x = T, allow.cartesian = T)
  
  #bc the model is built with new initators, we need to merge it back with all trt =1
  trt_wt <- dt_cleaned[cal_time == trt_time-1, .(id, wt_trt = wt)]
  df <- merge(df, trt_wt, by = "id", all.x = T, allow.cartesian = T)
  df[cal_time > trt_time-1 & is.na(wt), wt := wt_trt]
  df[ , wt_trt := NULL]
  
  return(df)
}

# ── STAGE 5: Pooled logistic regression ───────────────────────────────────────
run_pooled_lr <- function(dt_ipw) {

  frm <- ~trt + cal_time + cal_timesqr + trt:cal_time + trt:cal_timesqr
  mat <- model.matrix(frm, dt_ipw)

  model <- fastglm(x       = mat,
                   y       = dt_ipw$outcome,
                   weights = dt_ipw$wt,
                   family  = binomial(link = "logit"))
  
  tp    <- 0:(MAX_TIME - 1)
  notrt <- data.table(trt = 0L, cal_time = tp, cal_timesqr = tp^2)
  trtd  <- data.table(trt = 1L, cal_time = tp, cal_timesqr = tp^2)
  
  notrt[, risk := predict(model, model.matrix(frm, notrt), type = "response")]
  trtd[,  risk := predict(model, model.matrix(frm, trtd),  type = "response")]
  
  notrt_r <- notrt[, .(cal_time, cuminc_notrt = cumsum(risk) * 100)]
  trt_r   <- trtd[,  .(cal_time, cuminc_trt   = cumsum(risk) * 100)]
  
  res <- merge(trt_r, notrt_r, by = "cal_time")
  res[, RD := cuminc_trt - cuminc_notrt]
  res[, RR := cuminc_trt / cuminc_notrt]
  
  return(res)
}


# =============================================================================
# QUICK SANITY CHECK — run this to confirm the pipeline works before bootstrap
# =============================================================================
# dt         <- generate_data(seed = 42)
# ps_matched <- run_psm(dt)
# cox_res    <- run_cox(ps_matched, dt)
# dt_ipw     <- run_ipw(dt)
# lr_res     <- run_pooled_lr(dt_ipw)
#
# cat("HR (PP):", round(cox_res$hr_pp,  3), "\n")
# cat("RD at t59:", round(lr_res[cal_time==59, RD], 3), "\n")
# =============================================================================
