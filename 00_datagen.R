#' Improved Data generation function
#'
#' @param K max number of time points
#' @param nsim number of simulated datasets
#' @param A vector of treatment parameters c(beta_0_A, beta_U1_A) OR individual beta_0_A parameter
#' @param Y vector of outcome parameters c(beta_0_Y, beta_U2_Y, beta_t_Y, beta_t2_Y, beta_A_Y, beta_tA_Y, beta_t2A_Y) OR individual parameters
#' @param beta_0_A,beta_U1_A Treatment A parameters (used when A is not provided)
#' @param beta_0_Y,beta_U2_Y,beta_t_Y,beta_t2_Y,beta_A_Y,beta_tA_Y,beta_t2A_Y Outcome Y parameters (used when Y is not provided)
#' @param seed Random seed
#' @import data.table
#' @export

datagen_stage <- function(
    K,
    nsim,
    
    # Vector inputs (preferred by scenario functions)
    A = NULL,
    Y = NULL,
    
  # Individual parameter inputs (backward compatibility)
    beta_0_A = NULL, beta_L_A = NULL,beta_t_A = NULL,
    beta_0_Y = NULL, beta_t_Y = NULL, beta_t2_Y = NULL,
    beta_A_Y = NULL, beta_tA_Y = NULL, beta_t2A_Y = NULL,
    seed = NULL
) {

  
  # Handle vector vs individual parameter inputs
  if (!is.null(A)) {
    stopifnot(length(A) == 4)
    beta_0_A <- A[1]
    beta_L_A <- A[2]
    beta_t_A <- A[3]
    beta_t2_A <- A[4]
  }
  
  if (!is.null(Y)) {
    stopifnot(length(Y) == 6)
    beta_0_Y <- Y[1]
    beta_t_Y <- Y[2]
    beta_t2_Y <- Y[3]
    beta_A_Y <- Y[4]
    beta_tA_Y <- Y[5]
    beta_t2A_Y <- Y[6]
  }
  
  
  # Validate all required parameters are provided
  required_params <- c("beta_0_A", "beta_L_A","beta_t_A", "beta_0_Y", "beta_t_Y",
                       "beta_t2_Y", "beta_A_Y", "beta_tA_Y", "beta_t2A_Y"
                       )
  
  for (param in required_params) {
    if (is.null(get(param))) {
      stop(paste("Parameter", param, "must be provided either individually or through vector input"))
    }
  }
  
  #set.seed(seed)
  dt <- CJ(id = 1:nsim, time = 0:(K - 1))
  
  indiv <- data.table(
    id = 1:nsim,
    L = runif(nsim, min = 0, max = 1)
  )
  
  indiv[, A0 := rbinom(nsim, 1, plogis(beta_0_A + beta_L_A * L))]
  
  # Step 2: Merge in individual-level covariates 
  dt <- merge(dt, indiv, by = "id")
  rm(indiv)
  
  # Step 1: Expand to full panel (id x time)
  dt[, `:=`(
    t2 = time^2,
    tstart = time,
    tstop = time + 1)  # Step 3: Add squared time and other derived variables
  ][, Astar := fifelse(A0 == 0 & time>0, rbinom(.N, 1, plogis(beta_t_A)), 0)][, trt_time := ifelse(A0 == 1, 0, which(Astar == 1)[1]), by = id 
   ] [,A := fifelse(A0 == 1, 1, fifelse(!is.na(trt_time) & time >= trt_time, 1, 0))][, Astar := NULL][, trt_time := NULL]
  #here if there is patients that never started the treatment would return NA in all trt_time and mass up with rest of analysis
  #so we replace all NA A to  0
  dt = dt %>% 
    mutate(across(A, ~ replace_na(.x, 0)))
  #here we save percentage of A change with time to data frame 
  percent_A <- dt %>% select(A, time)%>% group_by(time)%>% summarise(percent_A = sum(A)/nsim) 
  # Step 2: Merge in individual-level covariates (U1, U2, A)
  dt[, Ystar := rbinom(nsim*K, 1, plogis(
    beta_0_Y  +
      beta_t_Y * time + beta_t2_Y * t2 +
      beta_A_Y * A + beta_tA_Y * time * A +
      beta_t2A_Y * t2 * A))
  ][, event_time := which(Ystar == 1)[1], by = id # Step 6: Censor at first Ystar per individual
  ]
  #save culmulative percent Y as data frame 
  culmulative_percent_Y <- dt %>%  
    filter(!is.na(event_time))%>% 
    select(id, event_time,time)%>%
    group_by(event_time)%>% summarise(n=n()/60)%>% mutate(per_Y = cumsum(n)/nsim)
  dt <- dt[time < event_time | is.na(event_time) == T][, event_time := NULL]
  setnames(dt, old = c("time", "t2", "A", "Ystar"), new = c("cal_time", "cal_timesqr","trt", "outcome"))
 
  return(list(dt = dt ,
              culmulative_percent_Y = culmulative_percent_Y,
              percent_A))
}
  
#to use those dataset later extract them out of list 
