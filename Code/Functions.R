library(synthdid)
library(rngtools)
library(xtable)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(reshape2)
library(patchwork)
library(fixest)

#######################################################
# Fix y-axis limit for all plotting functions








######################################################

#rm(list=ls())

########## Function to apply TWFE Regression ##########

twfe_formula <- as.formula("value ~ treated | Observation + Time")

twfe <- function(data, formula) {
  # Define the model
  model <- feols(formula, data = data)
  # Return only treated coefficient
  return(coef(summary(model)))
}


########## Simulation functions ##########

##### Function to simulate data with mean-zero innovations, different/identical 
### pre-/post trends and static/dynamic homogeneous treatment: 

### Inputs:
# N: Number of control units
# T: Number of time periods
# n_treated: Number of treated units
# treatment_period: Period in which treatment is applied
# pre_control_trend: Pre-treatment trend for control units
# pre_treated_trend: Pre-treatment trend for treated units
# post_control_trend: Post-treatment trend for control units
# post_treated_trend: Post-treatment trend for treated units
# treatment_effect: Treatment effect
# effect: Type of treatment effect (static or dynamic) as string input
# treated_mean: Mean value of initial random draw for treated units -> IPTW

treatment_simulation <- function(N, T, n_treated, treatment_period, pre_control_trend, pre_treated_trend,
                                        post_control_trend, post_treated_trend, treatment_effect, effect = "static", treated_mean = 0) {
  # Set up data frame with N rows and T columns
  control <- data.frame(matrix(ncol = T, nrow = N))
  # Set up separate data frame for treated units
  treated <- data.frame(matrix(ncol = T, nrow = n_treated))
  # Name dfs correctly
  colnames(control) <- paste(1:T)
  rownames(control) <- paste("Observation", 1:N, sep = "")
  colnames(treated) <- paste(1:T)
  rownames(treated) <- paste("Treated", 1:n_treated, sep = "")
  # Draw N random draws as starting points
  control[,1] <- rnorm(N, 0, 1)
  treated[,1] <- rnorm(n_treated, treated_mean, 1)
  # Fill rest of columns as previous value + upwards time-trend added and random draw and if treated, add treatment effect
  for (t in 2:T) {
    # Control units: Always add mean zero innovation
    control[,t] <- control[,t-1] + rnorm(N, 0, 1) +
      # If treatment is not yet applied, add pre-treatment control trend
      ifelse(t < treatment_period, pre_control_trend, 0) +
      # If treatment has been applied, add post-treatment control trend
      ifelse(t >= treatment_period, post_control_trend, 0)
    # Treated units: Always add mean zero innovation
    treated[,t] <- treated[,t-1] + rnorm(n_treated, 0, 1) +
      # If treatment is not yet applied, add pre-treatment treated trend
      ifelse(t < treatment_period, pre_treated_trend, 0) +
      # if we are in treatment period, add treatment effect depending on effect type
      ifelse(effect == "static" && t == treatment_period, treatment_effect, 0) +
      ifelse(effect == "dynamic" && t >= treatment_period, treatment_effect, 0) +
      # If we are in or after treatment period, add post-treatment treated trend
      ifelse(t >= treatment_period, post_treated_trend, 0) 
      
  }
  # Combine dfs
  df <- rbind(control, treated)
  # Pivot to long format
  df <- df %>% mutate(Observation = rownames(df)) %>% melt(id.vars = "Observation") %>%
    rename(Time = variable)  %>% arrange(Time)
  df$Time <- as.integer(df$Time)
  # Add treatment indicator period
  df$treated <- ifelse(df$Observation %in% rownames(treated) & as.numeric(df$Time) >= treatment_period, 1, 0) 
  df$treated <- as.integer(df$treated)
  
  return(df)
}



##### Function to simulate data with mean-zero innovations, different/identical 
### pre-/post trends and static/dynamic UNIT-SPECIFIC, time-constant treatment effects

### Additional inputs:
# deviation: Standard deviation of the treatment effect distribution
# can be understood as degree of heterogeneity in treatment effects

heterogeneous_treatment_simulation <- function(N, T, n_treated, treatment_period, pre_control_trend, pre_treated_trend,
                                 post_control_trend, post_treated_trend, treatment_effect, effect = "static", deviation = 1, treated_mean = 0) {
  # Set up data frame with N rows and T columns
  control <- data.frame(matrix(ncol = T, nrow = N))
  # Set up separate data frame for treated units
  treated <- data.frame(matrix(ncol = T, nrow = n_treated))
  # Name dfs correctly
  colnames(control) <- paste(1:T)
  rownames(control) <- paste("Observation", 1:N, sep = "")
  colnames(treated) <- paste(1:T)
  rownames(treated) <- paste("Treated", 1:n_treated, sep = "")
  # Draw N random draws as starting points
  control[,1] <- rnorm(N, 0, 1)
  treated[,1] <- rnorm(n_treated, treated_mean, 1)
  # Draw vector of heterogeneous treatment effects
  heterogeneity <- rnorm(n_treated, treatment_effect, deviation)
  # Fill rest of columns as previous value + upwards time-trend added and random draw and if treated, add treatment effect
  for (t in 2:T) {
    # Control units: Always add mean zero innovation
    control[,t] <- control[,t-1] + rnorm(N, 0, 1) +
      # If treatment is not yet applied, add pre-treatment control trend
      ifelse(t < treatment_period, pre_control_trend, 0) +
      # If treatment has been applied, add post-treatment control trend
      ifelse(t >= treatment_period, post_control_trend, 0)
    # Treated units: Always add mean zero innovation
    treated[,t] <- treated[,t-1] + rnorm(n_treated, 0, 1) +
      ifelse(t < treatment_period, pre_treated_trend, 0) +
      # if we are in treatment period, add treatment effect depending on effect type
      # Careful to disentangle if statement from vectorized operation!!!!!
      if (effect == "static" & t == treatment_period) heterogeneity else 0 +
      if (effect == "dynamic" & t >= treatment_period) heterogeneity else 0 +
      ifelse(t >= treatment_period, post_treated_trend, 0)
    
    
  }
  # Combine dfs
  df <- rbind(control, treated)
  # Pivot to long format
  df <- df %>% mutate(Observation = rownames(df)) %>% melt(id.vars = "Observation") %>%
    rename(Time = variable)  %>% arrange(Time)
  df$Time <- as.integer(df$Time)
  # Add treatment indicator period
  df$treated <- ifelse(df$Observation %in% rownames(treated) & as.numeric(df$Time) >= treatment_period, 1, 0) 
  df$treated <- as.integer(df$treated)
  
  return(df)
}


##### Function to simulate data with mean-zero innovations, different/identical 
### pre-/post trends and static/dynamic treatment effects thar are heterogeneous
### on unit AND time level

overall_heterogeneous_treatment_simulation <- function(N, T, n_treated, treatment_period, pre_control_trend, pre_treated_trend,
                                 post_control_trend, post_treated_trend, treatment_effect, effect = "static", deviation, treated_mean = 0) {
  # Set up data frame with N rows and T columns
  control <- data.frame(matrix(ncol = T, nrow = N))
  # Set up separate data frame for treated units
  treated <- data.frame(matrix(ncol = T, nrow = n_treated))
  # Name dfs correctly
  colnames(control) <- paste(1:T)
  rownames(control) <- paste("Observation", 1:N, sep = "")
  colnames(treated) <- paste(1:T)
  rownames(treated) <- paste("Treated", 1:n_treated, sep = "")
  # Draw N random draws as starting points
  control[,1] <- rnorm(N, 0, 1)
  treated[,1] <- rnorm(n_treated, treated_mean, 1)
  # Fill rest of columns as previous value + upwards time-trend added and random draw and if treated, add treatment effect
  for (i in 2:T) {
    # Control units: Always add mean zero innovation
    control[,i] <- control[,i-1] + rnorm(N, 0, 1) +
      # If treatment is not yet applied, add pre-treatment control trend
      ifelse(i < treatment_period, pre_control_trend, 0) +
      # If treatment has been applied, add post-treatment control trend
      ifelse(i >= treatment_period, post_control_trend, 0)
    # Treated units: Always add mean zero innovation
    treated[,i] <- treated[,i-1] + rnorm(n_treated, 0, 1) +
      ifelse(i < treatment_period, pre_treated_trend, 0) +
      # if we are in treatment period, add treatment effect depending on effect type
      if (effect == "static" & i == treatment_period) rnorm(n_treated, treatment_effect, deviation) else 0 +
      if (effect == "dynamic" & i >= treatment_period) rnorm(n_treated, treatment_effect, deviation) else 0 +
      # If we are in or after treatment period, add post-treatment treated trend
      ifelse(i >= treatment_period, post_treated_trend, 0)
  }
  # Combine dfs
  df <- rbind(control, treated)
  # Pivot to long format
  df <- df %>% mutate(Observation = rownames(df)) %>% melt(id.vars = "Observation") %>%
    rename(Time = variable)  %>% arrange(Time)
  df$Time <- as.integer(df$Time)
  # Add treatment indicator period
  df$treated <- ifelse(df$Observation %in% rownames(treated) & as.numeric(df$Time) >= treatment_period, 1, 0) 
  df$treated <- as.integer(df$treated)
  
  return(df)
}


########## Helper functions ##########

##### Function to perform ITWP matching ##########
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


##### Function to compute dynamic (accumulated) treatment effect #####
# Only works for homogeneous treatment effects
dynamic_treatment_effect <- function(treatment_effect, treatment_period, T) {
  # Define the number of periods after treatment
  num_periods <- T - treatment_period + 1
  # Compute the dynamic treatment effect using vectorized operations
  total_effect <- sum(treatment_effect * (1:num_periods))
  return(total_effect/num_periods)
}


##### Repeated simulation and estimation functions #####

### Function to simulate and estimate repeatedly
repeated_simulation <- function(simulation_function, simulation_params, iterations = 1000, iptw = FALSE) {
  # Initialize empty list to store dfs created in each iteration
  dfs = list()
  estimates_df = data.frame(matrix(ncol = 3, nrow = iterations))
  colnames(estimates_df) <- c("did", "sc", "sdid")
  
  for (i in 1:iterations) {
    # Simulate data
    data <- do.call(simulation_function, simulation_params)
    # If IPTW is True, apply respective matching function
    if (iptw == TRUE) {
      data <- perform_iptw_matching(data)
    }
    # Store data
    dfs[[i]] = data
    setup <- panel.matrices(data)
    estimates = lapply(estimators, function(estimator) { estimator(setup$Y,
                                                                   setup$N0, setup$T0) } )
    # Store estimates
    estimates_df[i, 1] = unlist(estimates)["did"]
    estimates_df[i, 2] = unlist(estimates)["sc"]
    estimates_df[i, 3] = unlist(estimates)["sdid"]
  }
  return(list("dfs" = dfs, "estimates" = estimates_df))
}


### Function to contrast time consumption of SDiD, SC and TWFE 
contrast_simulation <- function(simulation_function, simulation_params, iterations = 1000) {

  # Setup df for results
  times_df = data.frame(matrix(ncol = 4, nrow = iterations))
  colnames(times_df) <- c("N", "SDiD", "SC", "TWFE")
  
  # Perform estimation
  for (i in 1:iterations) {
   
     # Simulate data
    data <- do.call(simulation_function, simulation_params)
    
    # Start SDiD time
    sdid_start <- Sys.time()
    setup_sdid <- panel.matrices(data)
    # SDiD estimation
    sdid <- synthdid_estimate(setup_sdid$Y, setup_sdid$N0, setup_sdid$T0)
    # Stop time 
    sdid_end <- Sys.time()
    times_df[i, 2] <- sdid_end - sdid_start
    

    # Start SC time
    sc_start <- Sys.time()
    setup_sc <- panel.matrices(data)
    # SC estimation
    sc <- sc_estimate(setup_sc$Y, setup_sc$N0, setup_sc$T0)
    sc_end <- Sys.time()
    times_df[i, 3] <- sc_end - sc_start

    
    # Start TWFE time
    twfe_start <- Sys.time()
    # TWFE estimation
    twfe_est <- twfe(data, twfe_formula)
    twfe_end <- Sys.time()
    times_df[i, 4] <- twfe_end - twfe_start

  }
  # Store number of observations
  N <- as.numeric(simulation_params["n_treated"]) + as.numeric(simulation_params["N"])
  times_df[, 1] = N
  
  return(times_df)
}


### Function to simulate and estimate repeatedly including IPTW matching
repeated_simulation_iptw <- function(simulation_function, simulation_params, iterations = 1000) {
  dfs = list()
  estimates_df = data.frame(matrix(ncol = 3, nrow = iterations))
  colnames(estimates_df) <- c("did", "sc", "sdid")
  
  for (i in 1:iterations) {
    # Simulate data
    data <- do.call(simulation_function, simulation_params)
    
    # Apply IPTW matching
    data <- perform_iptw_matching(data)
   
     # Store data
    dfs[[i]] = data
    setup <- panel.matrices(data)
    estimates = lapply(estimators, function(estimator) { estimator(setup$Y,
                                                                   setup$N0, setup$T0) } )
    # Store estimates
    estimates_df[i, 1] = unlist(estimates)["did"]
    estimates_df[i, 2] = unlist(estimates)["sc"]
    estimates_df[i, 3] = unlist(estimates)["sdid"]
    
  }
  return(list("dfs" = dfs, "estimates" = estimates_df))
}


########## Plotting functions ##########

### All individual observations over time
plot_individual <- function(df) {
  # Create column indicating treated and control units
  df <- df %>% mutate(Type = ifelse(grepl("Treated", Observation), "Treated", "Control"))
  
  # Plot them separately to highlight treated observations
  ggplot(df, aes(x = Time, y = value, group = Observation, color = Type)) +
    geom_line() +
    geom_vline(xintercept = treatment_period, linetype = "dashed", linewidth = 1.5) +
    labs(title = "Development of Simulated Data Over Time",
         x = "Year",
         y = "Dependent Variable") +
    scale_linewidth_manual(values = c("Treated" = 1.5, "Control" = 1)) +
    scale_color_manual(values = c("Treated" = "red", "Control" = "black")) +
    theme_minimal() +
    theme(legend.position = "right",
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"), 
          axis.title = element_text(size = 12, face = "bold"),  
          axis.text = element_text(size = 10),
          panel.border = element_blank(),
          legend.title = element_text(size = 12, face = "bold", hjust = 0.5),
          legend.text = element_text(size = 10)) +
    scale_x_continuous(breaks = seq(min(df$Time), max(df$Time), by = 1))
}


### Means of treated and control group over time
plot_grouped <- function(df) {
  # Compute mean of treated units
  treated_means <- df %>% filter(grepl("Treated", Observation)) %>% 
    group_by(Time) %>% summarize(treated_avg = mean(value))
  # Compute mean of control group
  control_means <- df %>% filter(!grepl("Treated", Observation)) %>% 
    group_by(Time) %>% summarize(control_avg = mean(value))
  # Combine the two dfs
  average_data <- merge(treated_means, control_means)
  ggplot(data = average_data) +
    geom_line(aes(x = Time, y = control_avg, color = "Control mean"), linewidth = 1) +
    geom_line(aes(x = Time, y = treated_avg, color = "Treated mean"), linewidth = 1) +
    geom_vline(xintercept = treatment_period, linetype = "dashed", linewidth = 1.5) +
    labs(title = "Development of grouped means over time",
         x = "Year",
         y = "Dependent variable") +
    theme_minimal() +
    theme(legend.position = "right",
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"), 
          axis.title = element_text(size = 12, face = "bold"),  
          axis.text = element_text(size = 10)) +
    scale_color_manual(values = c("Treated mean" = "red", "Control mean" = "black")) +
    scale_x_continuous(breaks = seq(min(df$Time), max(df$Time), by = 1)) #+
  # Set ylim
  #ylim(95, 110)
}


### Means of treated and control group over time, over all simulations
plot_grouped_simulations <- function(dfs) {
  dfs <- do.call(rbind, dfs)
  # Compute mean of treated units
  treated_means <- dfs %>% filter(grepl("Treated", Observation)) %>% 
    group_by(Time) %>% summarize(treated_avg = mean(value))
  # Compute mean of control group
  control_means <- dfs %>% filter(!grepl("Treated", Observation)) %>% 
    group_by(Time) %>% summarize(control_avg = mean(value))
  # Combine the two dfs
  average_data <- merge(treated_means, control_means)
  ggplot(data = average_data) +
    geom_line(aes(x = Time, y = control_avg, color = "Control mean"), linewidth = 1) +
    geom_line(aes(x = Time, y = treated_avg, color = "Treated mean"), linewidth = 1) +
    geom_vline(xintercept = treatment_period, linetype = "dashed", linewidth = 1.5) +
    labs(title = paste0("Grouped means using ", iterations, " simulations with ", n_treated, " treated and ", N, " control units"),
         x = "Year",
         y = "Dependent variable",
         color = "Group",
         size = 14) +
    theme_minimal() +
    theme(legend.position = "right",
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"), 
          axis.title = element_text(size = 12, face = "bold"),  
          axis.text = element_text(size = 10),
          legend.text = element_text(size = 10), 
          legend.title = element_text(size = 12, face = "bold"),
          legend.title.align = 0.5) +
    scale_color_manual(values = c("Treated mean" = "red", "Control mean" = "black")) +
    scale_x_continuous(breaks = seq(min(df$Time), max(df$Time), by = 1)) #+
  # Set ylim
  #ylim(95, 110)
}

### Function to plot distribution of estimates
plot_estimates <- function(estimates_df, effect = 'static') {
  # If effect is dynamic, compute accumulated effect over time
  if (effect == 'dynamic') {
    treatment_effect <- dynamic_treatment_effect(treatment_effect, treatment_period, T)
  }
  
  ggplot(data = estimates_df) +
    geom_density(aes(x = did, color = "DiD"), linewidth = 1) +
    geom_density(aes(x = sc, color = "SC"), linewidth = 1) +
    geom_density(aes(x = sdid, color = "SDiD"), linewidth = 1) +
    geom_vline(xintercept = treatment_effect, linetype = "dashed", linewidth = 1.5) + 
    labs(title = paste0("Estimate distribution using ", iterations, " simulations ", n_treated, " treated and ", N, " control units"),
         x = "Estimate",
         y = "Density",
         color = "Method") +
    theme_minimal() +
    theme(legend.position = "right",
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"), 
          axis.title = element_text(size = 12, face = "bold"),  
          axis.text = element_text(size = 10),
          legend.title = element_text(face = "bold"),
          legend.text = element_text(size = 12)) +
    scale_color_manual(values = c("DiD" = "#1a80bb", "SC" = "#b8b8b8", "SDiD" = "#ea801c"),
                       labels = c("DiD" = "DiD", "SC" = "SC", "SDiD" = "SDiD")) +
    guides(color = guide_legend(title = "Estimator", title.position = "top", title.hjust = 0.5))
}

### Function to plot bias of estimates
plot_bias <- function(df, effect = 'static', y_limit = ylim(0, 0.5)) {
  
  # If effect is dynamic, compute accumulated effect over time
  if (effect == 'dynamic') {
    treatment_effect <- dynamic_treatment_effect(treatment_effect, treatment_period, T)
  }
    # Compute bias per method
  df <- df %>% mutate(did_bias = abs(did - treatment_effect),
                      sc_bias = abs(sc - treatment_effect),
                      sdid_bias = abs(sdid - treatment_effect))

  # Plot biases against smaple size
  bias_plot <- ggplot(data = df, aes(x = N)) +
    geom_line(aes(y = did_bias, color = "DiD"), linewidth = 1.5) +
    geom_line(aes(y = sc_bias, color = "SC"), linewidth = 1.5) +
    geom_line(aes(y = sdid_bias, color = "SDiD"), linewidth = 1.5) +
    labs(title = paste0("Distribution of bias per method depending on population size"),
         x = "Control units",
         y = "Bias",
         color = "Method") +
    theme_minimal() +
    theme(legend.position = "right",
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"), 
          axis.title = element_text(size = 12, face = "bold"),  
          axis.text = element_text(size = 10),
          legend.title = element_text(face = "bold"),
          legend.text = element_text(size = 12)) +
    scale_color_manual(values = c("DiD" = "#1a80bb", "SC" = "#b8b8b8", "SDiD" = "#ea801c"),
                       labels = c("DiD" = "DiD", "SC" = "SC", "SDiD" = "SDiD")) +
    guides(color = guide_legend(title = "Estimator", title.position = "top", title.hjust = 0.5)) +
    y_limit
  
  return(list("bias_df" = df, "plot" = bias_plot))
}

