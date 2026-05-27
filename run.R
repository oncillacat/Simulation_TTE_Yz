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
dt <- read.csv('Senario 2 t rr exp(0.7) simulated data 050526.csv')
dt<-dt[,-1]
setDT(dt)
#=========================================================
# pipeline for scene 1 
#=========================================================
#method traditional
#----method traditional----
dt <- generate_data( nsim     = 100000,
                     K        = 60,
                     A_params = c(-1.5, 0,5, -4.8, 0),       # beta_0_A, beta_L_A, beta_t_A, beta_t2_A
                     Y_params = c(-6, 0.01, 0.5, 0, 0, 0,0)) # beta_0_Y , beta_t_Y  ,beta_L_Y,
                                                             # beta_t2_Y ,beta_A_Y,beta_tA_Y, beta_t2A_Y
#load("dt.RData") 
#data saved in 03232026 used in all presentations with RR = 1
#pipeline create all 5 run and mark max_60 = 60 + N_boot =100 
#Note to self: delete max_60 = 60 + N_boot =100 if we don't use them in the end 
#* @param dt this is the only input
source("pipeline.R") 
ps_res<- run_psm(dt)
ps_matched <- ps_res$ps_matched
df_ps <- ps_res$df
cox_res <- run_cox(ps_matched,df_ps) # a list of hr, km curve, include both ITT & AT treatment 

res_lr <- data.table(bsample = NULL,
                     RD_t59_itt  = NULL,
                     RR_t59_itt  = NULL,
                     RD_t59_pp  = NULL,
                     RD_t59_pp  = NULL)

pb = txtProgressBar(min = 0, max = 100, initial = 0)
for (i in 1:100){
  #method TTE 
  unique_ids <- unique(dt$id)
  sampled_ids <- sample(unique_ids, size = length(unique_ids), replace = TRUE)
  id_map <- data.table(id = sampled_ids, boot_id = seq_along(sampled_ids))
  boot_dt <- dt[id_map, on = "id",allow.cartesian = T]
  
  res_IPW_boot <- run_ipw_boot(boot_dt)
  df_pp <- res_IPW_boot$df_ipw_pp
  df_itt <- res_IPW_boot$df_ipw_itt
  lr_res_itt <- run_pooled_lr(df_itt)
  lr_res_pp <- run_pooled_lr(df_pp)
  
  scalar <- data.table(
    bsample = i,
    RD_t59_itt  = lr_res_itt[cal_time == 59, RD],
    RR_t59_itt  = lr_res_itt[cal_time == 59, RR],
    RR_t59_pp = lr_res_pp[cal_time == 59, RR],
    RD_t59_pp  =lr_res_itt[cal_time == 59, RD]
  )
  res_lr <- rbind(scalar,res_lr)
  setTxtProgressBar(pb,i)
}

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
