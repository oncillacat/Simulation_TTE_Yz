rm(list =ls())
library(data.table)
library(dplyr)
library(FNN)
library(fastglm)
library(survival)
library(broom)
library(tidyr)
library(ggplot2)

logit_ps <- function(p) log(p / (1 - p))   # used in PSM caliper
MAX_TIME <- 60    # end of study (month 59 = last time point)
N_BOOT   <- 100   # number of bootstrap iterations
#=========================================================
# pipeline for scene 1 
#=========================================================
#method traditional
#----method traditional----
#dt <- generate_data()
load("dt.RData")
#data saved in 03232026 used in all presentations 
#pipeline create all 5 run and mark max_60 = 60 + N_boot =100 
#Note to self: delete max_60 = 60 + N_boot =100 if we don't use them in the end 
#* @param dt this is the only input
source("pipeline.R") 
ps_res<- run_psm(dt)
ps_matched <- ps_res$ps_matched
df_ps <- ps_res$df
cox_res <- run_cox(ps_matched,df_ps) # a list of hr, km curve, include both ITT & AT treatment 

res_lr <- data.table(bsample = NULL,
                     RD_t59  = NULL,
                     RR_t59  = NULL)
pb = txtProgressBar(min = 0, max = 100, initial = 0)
for (i in 1:100){
  #method TTE 
  unique_ids <- unique(dt$id)
  sampled_ids <- sample(unique_ids, size = length(unique_ids), replace = TRUE)
  id_map <- data.table(id = sampled_ids, boot_id = seq_along(sampled_ids))
  boot_dt <- dt[id_map, on = "id", allow.cartesian = TRUE]
  
  df_IPW <- run_ipw_boot(boot_dt)
  lr_res <- run_pooled_lr(df_IPW)
  
  scalar <- data.table(
    bsample = i,
    RD_t59  = lr_res[cal_time == 59, RD],
    RR_t59  = lr_res[cal_time == 59, RR]
  )
  res_lr <- rbind(scalar,res_lr)
  setTxtProgressBar(pb,i)
}

df_ipw = run_ipw(dt)
lr_res <- run_pooled_lr(df_ipw)

#---- result extract -----
mod_pp <- cox_res$cox_pp
mod_itt <- cox_res$cox_itt
#record log result remeber to exp
exp(confint(mod_pp))
exp(confint(mod_itt))
# TTE method 
quantile(res_lr$RR_t59, probs = 0.025) 
quantile(res_lr$RR_t59, probs = 0.975)
lr_res$
