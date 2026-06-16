'follow-up day extraction 060726'
'scen 2 TRR =0.4 '
rm(list = ls())
load('T040530A.RData')
library(dplyr)
library(data.table)
library(broom)
library(tidyr)
library(ggplot2)

'trad method'
ps_matched <- cox_res$ps_matched

'ITT'
hist(ps_matched$follow_ITT, main = "Scen2 TRR = 0.4 follow up trad ITT")
mean(ps_matched$follow_ITT)
median(ps_matched$follow_ITT)
c(quantile(ps_matched$follow_ITT, probs = 0.25),quantile(ps_matched$follow_ITT, probs = 0.75))

'PP'
hist(ps_matched$follow_up, main = "Scen2 TRR = 0.4 follow up trad PP")
mean(ps_matched$follow_up)
median(ps_matched$follow_up)
c(quantile(ps_matched$follow_up, probs = 0.25),quantile(ps_matched$follow_up, probs = 0.75))


'method TTE'
'ITT'
dat <- df_ipw_itt%>%dplyr::group_by(id)%>% summarise(follow_up = max(tstop))
hist(dat$follow_up,main = 'Scen2 TRR =0.4 follow up ITT TTE')
mean(dat$follow_up)
median(dat$follow_up)
c(quantile(dat$follow_up, probs = 0.25),quantile(dat$follow_up, probs = 0.75))

'pp'
dat <- df_ipw_pp%>%dplyr::group_by(id)%>% summarise(follow_up = max(tstop))
hist(dat$follow_up,main = 'Scen2 TRR =0.4 follow up PP TTE')
mean(dat$follow_up)
median(dat$follow_up)
c(quantile(dat$follow_up, probs = 0.25),quantile(dat$follow_up, probs = 0.75))

'end scen 2 TRR =0.4'
