'KP curve extraction'
rm(list = ls())
load('T100528A.RData')
library(dplyr)
library(data.table)
library(tidyr)
library(ggplot2)
library(survminer)
ps_matched <- cox_res$ps_matched

ps_plot <- ps_matched[, c("follow_up", "outcome_flag", "trt")]

km_fit <- survfit(
  Surv(follow_up, outcome_flag) ~ trt,
  data = ps_plot)

plot(
  km_fit,
  col = c("#7AC5CD","#7AC5CD"),
  lty     = c(2, 1),
  lwd     = 2,
  xlab    = "Months since index date",
  ylab    = "Probability of remaining event-free",
  main    = "Kaplan-Meier curves Tradition (PP) TRR =0.4"
)

legend("bottomleft",
       legend = c("Treated", "Control"),
       lty    = c(2, 1),
       lwd    = 2)

dat <- df_ipw_pp%>%
  dplyr::group_by(id)%>% 
  summarise(follow_up = max(tstop), trt = unique(A0),outcome  = any(outcome == 1))%>%
  mutate(index = 0) 
dat$outcome <- as.integer(dat$outcome)
km_ipw <- survfit(
  Surv(index, follow_up, outcome) ~ trt,
  data    = dat
)

plot(
  km_ipw,
  col = c("#7AC5CD","#7AC5CD"),
  lty     = c(2, 1),
  lwd     = 2,
  xlab    = "Months since index date",
  ylab    = "Probability of remaining event-free",
  main    = "Kaplan-Meier curves TTE (PP) TRR =0.4"
)

legend("bottomleft",
       legend = c("Treated", "Control"),
       lty    = c(2, 1),
       lwd    = 2)
#ITT 
ps_plot <- ps_matched[, c("follow_ITT", "outcome_ITT", "trt")]

km_fit <- survfit(
  Surv(ffollow_ITT, outcome_ITT) ~ trt,
  data = ps_plot)

plot(
  km_fit,
  lty     = c(2, 1),
  lwd     = 2,
  xlab    = "Months since index date",
  ylab    = "Probability of remaining event-free",
  main    = "Kaplan-Meier curves Tradition (ITT) TRR =1"
)

legend("bottomleft",
       legend = c("Control", "Treated"),
       lty    = c(2, 1),
       lwd    = 2)

'TTE ITT'
dat <- df_ipw_itt%>%
  dplyr::group_by(id)%>% 
  summarise(follow_up = max(tstop), trt = unique(A0),outcome  = any(outcome == 1))%>%
  mutate(index = 0) 
dat$outcome <- as.integer(dat$outcome)
km_ipw <- survfit(
  Surv(index, follow_up, outcome) ~ trt,
  data    = dat
)

plot(
  km_ipw,
  lty     = c(2, 1),
  lwd     = 2,
  xlab    = "Months since index date",
  ylab    = "Probability of remaining event-free",
  main    = "Kaplan-Meier curves TTE (ITT) TRR =1"
)

legend("bottomleft",
       legend = c("Control", "Treated"),
       lty    = c(2, 1),
       lwd    = 2)





col_09 <- "#08306B"  # darkest
col_07 <- "#2171B5"  # mid
col_04 <- "#9ECAE1"  # lightest