perform_iptw_matching <- function(data) {
  
  # Identify the treatment period
  treatment_period <- min(data$Time[data$treated == 1])
  
  # Compute pre-treatment means 
  pre_treatment_means <- data %>%
    filter(Time < treatment_period) %>%
    group_by(Observation) %>%
    summarize(pre_mean = mean(value))
  
  # Merge pre-treatment means with the original data
  merged_data <- merge(data, pre_treatment_means, by = "Observation")
  
  # Fit logistic regression to calculate propensity scores based on pre-treatment means
  model <- glm(treated ~ pre_mean, data = merged_data, family = binomial())
  merged_data$propensity_score <- predict(model, type = "response")
  
  # Compute IPTW weights
  merged_data <- merged_data %>%
    mutate(weight = ifelse(grepl("Treated", Observation), 1 / propensity_score, 1 / (1 - propensity_score)))
  merged_data <- merged_data %>%
    mutate(weight = ifelse(treated == 1, 1 / propensity_score, 1 / (1 - propensity_score)))
  
  
  # Compute weighted observations
  merged_data <- merged_data %>%
    mutate(weighted_value = value * weight)
  
  final_df <- merged_data %>% select(Observation, Time, weighted_value, treated) %>%
    rename(value = weighted_value) %>% arrange(Time)
  
  # Return data with IPTW weights applied
  return(final_df)
}