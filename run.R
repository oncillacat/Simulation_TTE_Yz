rm(list =ls())
library(data.table)
library(dplyr)
library(FNN)
library(fastglm)
library(survival)
library(broom)
library(tidyr)
library(ggplot2)
library(parallel)
library(pbmcapply)  # install.packages("pbmcapply")

source("pipeline.R") 
logit_ps <- function(p) log(p / (1 - p))   # used in PSM caliper
MAX_TIME <- 60    # end of study (month 59 = last time point)
N_BOOT   <- 100   # number of bootstrap iterations
nsim <- 100000
dt <- read.csv('T100528A.csv')
dt<-dt[,-1]
setDT(dt)
#method traditional
#----method traditional----
dt <- generate_data( nsim     = 100000,
                     K        = 60,
                     beta_0_A = -1.9, beta_L_A = 0.5,beta_t_A = -4.8, beta_t2_A =0,
                     beta_0_Y = -6, beta_t_Y = 0.01, beta_L_Y = 0.5, 
                     beta_t2_Y = 0,beta_A_Y = 0, beta_tA_Y = 0, beta_t2A_Y = 0)
# %Y by time 
Per_Y <- dt[,
   .(new_events = sum(outcome, na.rm = TRUE)),
   by = cal_time
][order(cal_time)
][, cum_events := cumsum(new_events)][
  , per_Y := cum_events/nsim
]
p1 <- Per_Y%>%
  ggplot(aes(x = cal_time, y = per_Y))+ geom_point(shape = 1, size = 2) + theme_bw()+
  scale_y_continuous(
    limits = c(0, 1),
    labels = scales::number_format(accuracy = 0.1)
  ) +
  labs(x = "Time after eligible (Month)",
       y =" Percentage of cumulated outcome (%)",
       title = "Percentage of cumulative outcome by time")

#%A by time
dt[, start_treatment := as.integer(trt == 1 & shift(trt, fill = 0) == 0), by = id]

Per_A <- dt[,
   .(new_trts = sum(start_treatment, na.rm = TRUE)),
   by = cal_time
][order(cal_time)
][, cum_trts := cumsum(new_trts)][
  , per_A := cum_trts/nsim
]

p2 <- Per_A %>% 
  ggplot(aes(x = cal_time, y = per_A))+ geom_point(shape = 1, size = 2) + theme_bw()+
  scale_y_continuous(limits = c(0, 1),labels = scales::number_format(accuracy = 0.1)) +
  labs(title = "Cumulative Proportion of \n Treatment Initiation Overtime ",
       y =" Percentage of Initiation (%)",
       x = "Time after eligible (Month)", ylim = c(0,1))
write.csv(dt,"T100528A.csv")

#load("dt.RData") 
#data saved in 03232026 used in all presentations with RR = 1
#pipeline create all 5 run and mark max_60 = 60 + N_boot =100 
#Note to self: delete max_60 = 60 + N_boot =100 if we don't use them in the end 
#* @param dt this is the only input

ps_res<- run_psm(dt)
ps_matched <- ps_res$ps_matched
df_ps <- ps_res$df
cox_res <- run_cox(ps_matched,df_ps) # a list of hr, km curve, include both ITT & AT treatment 

res_lr <- data.table(bsample = NULL,
                     RD_t59_itt  = NULL,
                     RR_t59_itt  = NULL,
                     RD_t59_pp  = NULL,
                     RD_t59_pp  = NULL)
one_boot <- function(i) {
  unique_ids  <- unique(dt$id)
  sampled_ids <- sample(unique_ids, size = length(unique_ids), replace = TRUE)
  id_map      <- data.table(id = sampled_ids, boot_id = seq_along(sampled_ids))
  boot_dt     <- dt[id_map, on = "id", allow.cartesian = TRUE]
  
  res_IPW_boot <- run_ipw_boot(boot_dt)
  lr_res_itt   <- run_pooled_lr(res_IPW_boot$df_ipw_itt)
  lr_res_pp    <- run_pooled_lr(res_IPW_boot$df_ipw_pp)
  
  data.table(
    bsample    = i,
    RD_t59_itt = lr_res_itt[cal_time == 59, RD],
    RR_t59_itt = lr_res_itt[cal_time == 59, RR],
    RR_t59_pp  = lr_res_pp[cal_time == 59, RR],
    RD_t59_pp  = lr_res_pp[cal_time == 59, RD]
  )
}

results_list <- pbmclapply(
  1:100,
  one_boot,                                    
  mc.cores = 2
  
)

res_lr <- rbindlist(results_list)

#point esitimator

res_ipw <- run_ipw(dt)

df_ipw_pp <- res_ipw$df_ipw_pp
lr_res_pp <- run_pooled_lr(df_ipw_pp)
df_ipw_itt<- res_ipw$df_ipw_itt
lr_res_itt <- run_pooled_lr(df_ipw_itt)

#---- result extract -----
mod_pp <- cox_res$cox_pp
mod_itt <- cox_res$cox_itt
#record log result remeber to exp
exp(confint(mod_pp))
exp(confint(mod_itt))
# TTE method 
quantile(res_lr$RR_t59_itt, probs = 0.025) 
quantile(res_lr$RR_t59_itt, probs = 0.975)
quantile(res_lr$RR_t59_pp, probs = 0.025) 
quantile(res_lr$RR_t59_pp, probs = 0.975)
lr_res_pp$RR[60]
lr_res_itt$RR[60]
