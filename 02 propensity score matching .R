#02 propensity score matching draft 
# Load libraries
library(dplyr)
s
# Estimate Propensity Scores using logistic regression
ps_model <- glm(Y ~ L + A + time + t2,
                data = dt,
                family = binomial(link = "logit"))

# Predicted propensity sores
data$propensity_score <- predict(ps_model, type = "response")

# View first rows
head(data)

# Split treated and control groups
treated <- data %>% filter(treatment == 1)
control <- data %>% filter(treatment == 0)

# For each treated unit, find closest control by propensity score
match_index <- sapply(treated$propensity_score, function(ps) {
  which.min(abs(control$propensity_score - ps))
})

matched_control <- control[match_index, ]

# Combine matched data
matched_data <- rbind(treated, matched_control)

# Estimate Treatment Effect
treated_mean <- mean(matched_data$outcome[matched_data$treatment == 1])
control_mean <- mean(matched_data$outcome[matched_data$treatment == 0])
effect <- treated_mean - control_mean

cat(sprintf("Estimated Treatment Effect: %.2f\n", effect))

# with replacement 
# 1 to 1 matching 

