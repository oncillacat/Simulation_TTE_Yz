# Helper function to calculate IPW weights
calculate_ipw_weights <- function(data_boot_base, mechanism) {
  mat_boot_numA <- model.matrix(~ 1, data_boot_base)
  A <- data_boot_base$trt
  #W <- data_boot_base$nb_select
  
  numprobA <- fastglm(x = mat_boot_numA,
                      y = A,
                      #weights = W,
                      family = binomial(link = "logit"))
  

    mat_boot_denomA <- model.matrix(~ L, data_boot_base)
    
  denomprobA <- fastglm(x = mat_boot_denomA,
                        y = A,
                        weights = W,
                        family = binomial(link = "logit"))
  
  data_boot_base$ipw.num <- predict(numprobA, mat_boot_numA, type = "response")
  data_boot_base$ipw.denom <- predict(denomprobA, mat_boot_denomA, type = "response")
  
  data_boot_base$wt <- ifelse(data_boot_base$trt == 1,
                              data_boot_base$ipw.num / data_boot_base$ipw.denom,
                              (1 - data_boot_base$ipw.num) / (1 - data_boot_base$ipw.denom))
  
  
  return(data_boot_base[, c("id", "wt")])
}

# Function to expand follow-up for total effect estimation
#expand_followup_total <- function(data_boot) {
#  dead <- setDT(subset(data_boot, death == 1))
  
#  if (nrow(dead) == 0) return(data_boot)
  
#  dead$r <- pmax(0, 48 - dead$tstop)
#  dead <- dead[r > 0]
  
#  if (nrow(dead) == 0) return(data_boot)
  
#  dead <- dead[rep(seq_len(nrow(dead)), times = dead$r)]
#  dead$r <- NULL
  
#  dead[, `:=`(
#    cal_time = cal_time + seq_len(.N),
#    tstart = tstart + seq_len(.N),
#    tstop = tstop + seq_len(.N),
#    cal_timesqr = (cal_time + seq_len(.N))^2
#  ), by = id]
  
#  dat_exp <- rbind(data_boot, dead)
#  dat_exp$outcome <- ifelse(dat_exp$death == 1 & dat_exp$outcme != 1, 0, dat_exp$outcome)
  
#  return(dat_exp)
#}

# Function to expand follow-up for death outcome estimation
#expand_followup_death <- function(data_boot) {
#  event <- setDT(subset(data_boot, outcome == 1))
  
#  if (nrow(event) == 0) return(data_boot)
#  event$r <- pmax(0, 48 - event$tstop)
#  event <- event[r > 0]
  
#  if (nrow(event) == 0) return(data_boot)
#  event <- event[rep(seq_len(nrow(event)), times = event$r)]
#  event$r <- NULL
  
#  event[, `:=`(
#    cal_time = cal_time + seq_len(.N),
#    tstart = tstart + seq_len(.N),
#    tstop = tstop + seq_len(.N),
#    cal_timesqr = (cal_time + seq_len(.N))^2
#  ), by = id]
  
#  dat_exp <- rbind(data_boot, event)
#  dat_exp$death <- ifelse(dat_exp$outcome == 1 & dat_exp$death != 1, 0, dat_exp$death)
  
# return(dat_exp)
#}

# Function to prepare data based on analysis type
#prepare_analysis_data <- function(data_boot, type) {
  
#  if (type == "total") {
#    dat_exp <- expand_followup_total(data_boot)
#  } else if (type == "death") {
#    dat_exp <- expand_followup_death(data_boot)
#  } 
#  return(dat_exp)
#}

# Function to fit outcome models
fit_outcome_model <- function(dat_exp, type, weights = NULL) {
  mat_boot <- model.matrix(~ trt + cal_time + cal_timesqr + trt:cal_time + trt:cal_timesqr,
                           dat_exp)
  Y <- switch(type,
              "total" = dat_exp$outcome)
  
  W <- dat_exp$wt
  model <- fastglm(x = mat_boot,
                   y = Y,
                   weights = W,
                   family = binomial(link = "logit"))
  
  return(model)
}

# Function to calculate risk estimates
calculate_risk_estimates <- function(model, time_points = 0:47) {
  notrt_data <- data.table(
    trt = 0L,
    cal_time = time_points,
    cal_timesqr = time_points^2
  )
  
  trt_data <- data.table(
    trt = 1L,
    cal_time = time_points,
    cal_timesqr = time_points^2
  )
  mat_notrt <- model.matrix(~ trt + cal_time + cal_timesqr + trt:cal_time + trt:cal_timesqr,
                            notrt_data)
  mat_trt <- model.matrix(~ trt + cal_time + cal_timesqr + trt:cal_time + trt:cal_timesqr,
                          trt_data)
  
  notrt_data$risk <- predict(model, mat_notrt, type = "response")
  trt_data$risk <- predict(model, mat_trt, type = "response")
  
  notrt_risk <- notrt_data[, .(risk_notrt = mean(risk)), by = cal_time
  ][, cuminc_notrt := cumsum(risk_notrt) * 100]
  
  trt_risk <- trt_data[, .(risk_trt = mean(risk)), by = cal_time
  ][, cuminc_trt := cumsum(risk_trt) * 100]
  
  risk_results <- merge(trt_risk, notrt_risk, by = "cal_time")
  risk_results[, `:=`(
    risk_ratio = cuminc_trt / cuminc_notrt,
    risk_diff = cuminc_trt - cuminc_notrt
  )
  ][, c("risk_trt", "risk_notrt") := NULL]
  
  return(risk_results)
}

#process_bootstrap_iteration <- function(data_m, i, type, mechanism) {
#  if (i == 0) {
#    data_boot <- data_m
#    data_boot$nb_select <- 1L
#  } else {
#    data_boot <- boot(data = data_m, bsample = i)
#  }
  
#    data_boot_base <- subset(data_boot, cal_time == 0)
#    weight_data <- calculate_ipw_weights(data_boot_base, mechanism)
#    data_boot <- merge(data_boot, weight_data, by = "id")
#    weights <- data_boot$wt
  
#  dat_exp <- prepare_analysis_data(data_boot, type)
  
#  outcome_model <- fit_outcome_model(dat_exp, type)
  
#  risk_results <- calculate_risk_estimates(outcome_model)
#  risk_results[, `:=`(bsample = i)]
  
#  return(list(risk = risk_results))
#}

#' Parallel wrapper for bootstrap iterations
#' @param data_m Dataset for single simulation
#' @param nboot Number of bootstrap samples
#' @param type Analysis type
#' @param mechanism Confounding mechanism
#' @param parallel Whether to use parallel processing
#' @return Combined bootstrap results
#process_bootstrap <- function(data_m, nboot, type, mechanism, parallel = TRUE) {
  
#  if (parallel && nboot > 0) {
    # Parallel execution
#    bootstrap_results <- foreach::foreach(
#      i = 0:nboot,
#      .packages = c("data.table", "fastglm"),
#      .combine = function(...) list(
#        risk = rbindlist(lapply(list(...), function(x) x$risk))
#      ),
#      .options.future = list(seed = TRUE)
#    ) %dorng% {
#      process_bootstrap_iteration(data_m, i, type, mechanism)
#    }
#  } else {
#    # Sequential execution
#    bootstrap_results <- vector("list", nboot + 1)
#    for (i in 0:nboot) {
#      bootstrap_results[[i + 1]] <- process_bootstrap_iteration(
#        data_m, i, type, mechanism
#      )
#    }
    
    # Combine results
#    risk_list <- lapply(bootstrap_results, function(x) x$risk)
#    bootstrap_results <- list(
#      risk = rbindlist(risk_list)
#    )
#  }
  
#  return(bootstrap_results)
#}

#' Estimate Controlled Direct Effects, Total Effects, and Competing Risks
#'
#' @param data A data.table containing longitudinal data generated by the datagen function with columns including id, dataset, trt, outcome, death, cal_time, U1, U2, etc.
#' @param type Character string specifying the type of effect to estimate:
#'   \itemize{
#'     \item "total": Total effect of the treatment on outcome risk (expands follow-up for death)
#'     \item "death": Effect on competing risk of death (expands follow-up for outcome)
#'   }
#' @param boot Logical indicating whether to perform bootstrap resampling
#' @param nboot Integer number of bootstrap samples (ignored if boot = FALSE)
#' @param mechanism Integer specifying confounding mechanism for IPW:
#' @param parallel Logical indicating whether to use parallel processing
#' @importFrom data.table rbindlist
#' @export
#total_death <- function(data, type, boot, nboot, mechanism, parallel = TRUE) {
  
#  valid_types <- c("total", "death")
#  if (!type %in% valid_types) {
#    stop("Invalid type. Must be one of: ", paste(valid_types, collapse = ", "))
#  }
  
#  if (!boot) nboot <- 0
  
#  risk_results_list <- vector("list", length(unique(data$dataset)))
  
#  for (m in seq_along(unique(data$dataset))) {
#    data_m <- subset(data, dataset == m)
    
#    bootstrap_combined <- process_bootstrap(data_m, nboot, type, mechanism, parallel)
    
#    # Add dataset identifier
#    bootstrap_combined$risk[, dataset := m]
#    risk_results_list[[m]] <- bootstrap_combined$risk
#  }
  
#    return(list(final_risk_results))
#}