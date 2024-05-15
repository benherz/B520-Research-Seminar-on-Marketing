library(synthdid)
library(rngtools)
library(future)
library(doFuture)
library(future.batchtools)
library(xtable)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(reshape2)
library(patchwork)

#######################################################
# Fix y-axis limit for all plotting functions








######################################################

#rm(list=ls())

########## Simulation functions ##########

### Function with mean zero innovations and linear time trends
trend_simulation <- function(N, T, n_treated, treatment_period, control_trend, treated_trend) {
  # Set up dataframe with N rows and T columns
  df <- data.frame(matrix(ncol = T, nrow = N))
  # Set up separate dataframe for treated units
  treated <- data.frame(matrix(ncol = T, nrow = n_treated))
  # Name dfs correctly
  colnames(df) <- paste(1:T)
  rownames(df) <- paste("Observation", 1:N, sep = "")
  colnames(treated) <- paste(1:T)
  rownames(treated) <- paste("Treated", 1:n_treated, sep = "")
  # Draw N random draws as starting points
  df[,1] <- rnorm(N, 0, 1)
  treated[,1] <- rnorm(n_treated, 0, 1)
  # Fill rest of columns as previous value + upwards time-trend added (%) and random draw
  for (i in 2:T) {
    df[,i] <- df[,i-1] + control_trend + rnorm(N, 0, 1)
    treated[,i] <- treated[,i-1] + treated_trend + rnorm(n_treated, 0, 1)
  }
  
  # Combine dfs
  df <- rbind(df, treated)
  # Pivot to long format
  df <- df %>% mutate(Observation = rownames(df)) %>% melt(id.vars = "Observation") %>%
    rename(Time = variable)  %>% arrange(Time)
  df$Time <- as.integer(df$Time)
  # Add treatment indicator period
  df$treated <- ifelse(df$Observation %in% rownames(treated) & as.numeric(df$Time) >= treatment_period, 1, 0) 
  df$treated <- as.integer(df$treated)
  
  return(df)
}

### Function with mean zero innovations, treatment and violation of parallel-trends-assumption
diff_trend_simulation <- function(N, T, n_treated, treatment_period, control_trend, treated_trend, treatment_effect) {
  # Set up dataframe with N rows and T columns
  df <- data.frame(matrix(ncol = T, nrow = N))
  # Set up separate dataframe for treated units
  treated <- data.frame(matrix(ncol = T, nrow = n_treated))
  # Name dfs correctly
  colnames(df) <- paste(1:T)
  rownames(df) <- paste("Observation", 1:N, sep = "")
  colnames(treated) <- paste(1:T)
  rownames(treated) <- paste("Treated", 1:n_treated, sep = "")
  # Draw N random draws as starting points
  df[,1] <- rnorm(N, 0, 1)
  treated[,1] <- rnorm(n_treated, 0, 1)
  # Fill rest of columns as previous value + time-trend added (%) and random draw
  for (i in 2:T) {
    df[,i] <- df[,i-1] + rnorm(N, 0, 1) +
      ifelse(i >= treatment_period, control_trend, 0)
    treated[,i] <- treated[,i-1] + rnorm(n_treated, 0, 1) +
      ifelse(i >= treatment_period, treated_trend, 0) +
      ifelse(i == treatment_period, treatment_effect, 0)
  }
  
  # Combine dfs
  df <- rbind(df, treated)
  # Pivot to long format
  df <- df %>% mutate(Observation = rownames(df)) %>% melt(id.vars = "Observation") %>%
    rename(Time = variable)  %>% arrange(Time)
  df$Time <- as.integer(df$Time)
  # Add treatment indicator period
  df$treated <- ifelse(df$Observation %in% rownames(treated) & as.numeric(df$Time) >= treatment_period, 1, 0) 
  df$treated <- as.integer(df$treated)
  
  return(df)
}

### Function with mean-zero innovations, linear time trend and static (one-time) treatment effect
static_treatment_simulation <- function(N, T, n_treated, treatment_period, control_trend,
                                 treated_trend, treatment_effect) {
  # Set up dataframe with N rows and T columns
  df <- data.frame(matrix(ncol = T, nrow = N))
  # Set up separate dataframe for treated units
  treated <- data.frame(matrix(ncol = T, nrow = n_treated))
  # Name dfs correctly
  colnames(df) <- paste(1:T)
  rownames(df) <- paste("Observation", 1:N, sep = "")
  colnames(treated) <- paste(1:T)
  rownames(treated) <- paste("Treated", 1:n_treated, sep = "")
  # Draw N random draws as starting points
  df[,1] <- rnorm(N, 0, 1)
  treated[,1] <- rnorm(n_treated, 0, 1)
  # Fill rest of columns as previous value + upwards time-trend added and random draw and if treated, add treatment effect
  for (i in 2:T) {
    df[,i] <- df[,i-1] + control_trend + rnorm(N, 0, 1)
    treated[,i] <- treated[,i-1] + treated_trend + rnorm(n_treated, 0, 1) +
      ifelse(i == treatment_period, treatment_effect, 0)
  }
  # Combine dfs
  df <- rbind(df, treated)
  # Pivot to long format
  df <- df %>% mutate(Observation = rownames(df)) %>% melt(id.vars = "Observation") %>%
    rename(Time = variable)  %>% arrange(Time)
  df$Time <- as.integer(df$Time)
  # Add treatment indicator period
  df$treated <- ifelse(df$Observation %in% rownames(treated) & as.numeric(df$Time) >= treatment_period, 1, 0) 
  df$treated <- as.integer(df$treated)
  
  return(df)
}

### Function with mean-zero innovations, linear time trend and dynamic treatment effect
dynamic_treatment_simulation <- function(N, T, n_treated, treatment_period, control_trend,
                                        treated_trend, treatment_effect) {
  # Set up dataframe with N rows and T columns
  df <- data.frame(matrix(ncol = T, nrow = N))
  # Set up separate dataframe for treated units
  treated <- data.frame(matrix(ncol = T, nrow = n_treated))
  # Name dfs correctly
  colnames(df) <- paste(1:T)
  rownames(df) <- paste("Observation", 1:N, sep = "")
  colnames(treated) <- paste(1:T)
  rownames(treated) <- paste("Treated", 1:n_treated, sep = "")
  # Draw N random draws as starting points
  df[,1] <- rnorm(N, 0, 1)
  treated[,1] <- rnorm(n_treated, 0, 1)
  # Fill rest of columns as previous value + upwards time-trend added and random draw and if treated, add treatment effect
  for (i in 2:T) {
    df[,i] <- df[,i-1] + control_trend + rnorm(N, 0, 1)
    treated[,i] <- treated[,i-1] + treated_trend + rnorm(n_treated, 0, 1) +
      ifelse(i >= treatment_period, treatment_effect, 0)
  }
  # Combine dfs
  df <- rbind(df, treated)
  # Pivot to long format
  df <- df %>% mutate(Observation = rownames(df)) %>% melt(id.vars = "Observation") %>%
    rename(Time = variable)  %>% arrange(Time)
  df$Time <- as.integer(df$Time)
  # Add treatment indicator period
  df$treated <- ifelse(df$Observation %in% rownames(treated) & as.numeric(df$Time) >= treatment_period, 1, 0) 
  df$treated <- as.integer(df$treated)
  
  return(df)
}


### Function to simulate data with mean-zero innovations and heterogenuous, static treatment effects
heterogeneous_static_simulation <- function(N, T, n_treated, treatment_period, control_trend,
                                 treated_trend, treatment_effect) {
  # Set up dataframe with N rows and T columns
  df <- data.frame(matrix(ncol = T, nrow = N))
  # Set up separate dataframe for treated units
  treated <- data.frame(matrix(ncol = T, nrow = n_treated))
  # Name dfs correctly
  colnames(df) <- paste(1:T)
  rownames(df) <- paste("Observation", 1:N, sep = "")
  colnames(treated) <- paste(1:T)
  rownames(treated) <- paste("Treated", 1:n_treated, sep = "")
  # Draw N random draws as starting points
  df[,1] <- rnorm(N, 0, 1)
  treated[,1] <- rnorm(n_treated, 0, 1)
  # Fill rest of columns as previous value + time-trend added + random draw and if treated, add heterogenuous treatment effect
  for (i in 2:T) {
    df[,i] <- df[,i-1] + control_trend + rnorm(N, 0, 1)
    treated[,i] <- treated[,i-1] + treated_trend + rnorm(n_treated, 0, 1) +
      ifelse(i == treatment_period, rnorm(n_treated, treatment_effect, 1), 0)
  }

  # Combine dfs
  df <- rbind(df, treated)
  # Pivot to long format
  df <- df %>% mutate(Observation = rownames(df)) %>% melt(id.vars = "Observation") %>%
    rename(Time = variable)  %>% arrange(Time)
  df$Time <- as.integer(df$Time)
  # Add treatment indicator period
  df$treated <- ifelse(df$Observation %in% rownames(treated) & as.numeric(df$Time) >= treatment_period, 1, 0) 
  df$treated <- as.integer(df$treated)
  
  return(df)
}



### Function to simulate data with mean-zero innovations and heterogenuous, dynamic treatment effects
heterogeneous_dynamic_simulation <- function(N, T, n_treated, treatment_period, control_trend,
                                            treated_trend, treatment_effect) {
  # Set up dataframe with N rows and T columns
  df <- data.frame(matrix(ncol = T, nrow = N))
  # Set up separate dataframe for treated units
  treated <- data.frame(matrix(ncol = T, nrow = n_treated))
  # Name dfs correctly
  colnames(df) <- paste(1:T)
  rownames(df) <- paste("Observation", 1:N, sep = "")
  colnames(treated) <- paste(1:T)
  rownames(treated) <- paste("Treated", 1:n_treated, sep = "")
  # Draw N random draws as starting points
  df[,1] <- rnorm(N, 0, 1)
  treated[,1] <- rnorm(n_treated, 0, 1)
  # Fill rest of columns as previous value + time-trend added + random draw and if treated, add heterogenuous treatment effect
  for (i in 2:T) {
    df[,i] <- df[,i-1] + control_trend + rnorm(N, 0, 1)
    treated[,i] <- treated[,i-1] + treated_trend + rnorm(n_treated, 0, 1) +
      ifelse(i >= treatment_period, rnorm(n_treated, treatment_effect, 1), 0)
  }
  
  # Combine dfs
  df <- rbind(df, treated)
  # Pivot to long format
  df <- df %>% mutate(Observation = rownames(df)) %>% melt(id.vars = "Observation") %>%
    rename(Time = variable)  %>% arrange(Time)
  df$Time <- as.integer(df$Time)
  # Add treatment indicator period
  df$treated <- ifelse(df$Observation %in% rownames(treated) & as.numeric(df$Time) >= treatment_period, 1, 0) 
  df$treated <- as.integer(df$treated)
  
  return(df)
}


### Function to simulate data with mean-zero innovations and heterogenuous, dynamic unit-specific treatment effects
heterogeneous_dynamic_unit_simulation <- function(N, T, n_treated, treatment_period, control_trend,
                                             treated_trend, treatment_effect) {
  # Set up dataframe with N rows and T columns
  df <- data.frame(matrix(ncol = T, nrow = N))
  # Set up separate dataframe for treated units
  treated <- data.frame(matrix(ncol = T, nrow = n_treated))
  # Name dfs correctly
  colnames(df) <- paste(1:T)
  rownames(df) <- paste("Observation", 1:N, sep = "")
  colnames(treated) <- paste(1:T)
  rownames(treated) <- paste("Treated", 1:n_treated, sep = "")
  # Draw N random draws as starting points
  df[,1] <- rnorm(N, 0, 1)
  treated[,1] <- rnorm(n_treated, 0, 1)
  # Draw vector of hetrogeneous treatment effects
  heterogeneity <- rnorm(n_treated, treatment_effect, 1)
  # Fill rest of columns as previous value + time-trend added + random draw and if treated, add heterogenuous treatment effect
  for (i in 2:T) {
    df[,i] <- df[,i-1] + control_trend + rnorm(N, 0, 1)
    treated[,i] <- treated[,i-1] + treated_trend + rnorm(n_treated, 0, 1) +
      ifelse(i >= treatment_period, heterogeneity, 0)
  }
  
  # Combine dfs
  df <- rbind(df, treated)
  # Pivot to long format
  df <- df %>% mutate(Observation = rownames(df)) %>% melt(id.vars = "Observation") %>%
    rename(Time = variable)  %>% arrange(Time)
  df$Time <- as.integer(df$Time)
  # Add treatment indicator period
  df$treated <- ifelse(df$Observation %in% rownames(treated) & as.numeric(df$Time) >= treatment_period, 1, 0) 
  df$treated <- as.integer(df$treated)
  
  return(df)
}

### Function to simulate and estimate repeatedly
repeated_simulation <- function(simulation_function, simulation_params, iterations = 1000) {
  dfs = list()
  estimates_df = data.frame(matrix(ncol = 3, nrow = iterations))
  colnames(estimates_df) <- c("did", "sc", "sdid")
  for (i in 1:iterations) {
    # Simulate data
    data <- do.call(simulation_function, simulation_params)
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
  # Subset treated units
  treated_units <- df %>% filter(grepl("Treated", Observation))
  # Subset control units
  control_units <- df %>% filter(!grepl("Treated", Observation))
  # Plot them separately to highlight treated observations
  ggplot() +
    geom_line(data = treated_units, aes(x = Time, y = value, color = Observation), linewidth = 2) +
    geom_line(data = control_units, aes(x = Time, y = value, color = Observation), linewidth = 1) +
    geom_vline(xintercept = treatment_period, linetype = "dashed", size = 1.5) +
    labs(title = "Development of simulated data over time",
         x = "Year",
         y = "Dependent variable") +
    theme_minimal() +
    theme(legend.position = "right",
          plot.title = element_text(hjust = 0.5, size = 20), 
          axis.title = element_text(size = 16),  
          axis.text = element_text(size = 14)) #+  
    #scale_color_manual(values = c("Treated" = "red", "Control" = "black"))
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
    geom_vline(xintercept = treatment_period, linetype = "dashed", size = 1.5) +
    labs(title = "Development of grouped means over time",
         x = "Year",
         y = "Dependent variable") +
    theme_minimal() +
    theme(legend.position = "right",
          plot.title = element_text(hjust = 0.5, size = 20, face = "bold"), 
          axis.title = element_text(size = 16, face = "bold"),  
          axis.text = element_text(size = 14)) +
    scale_color_manual(values = c("Treated mean" = "red", "Control mean" = "black")) #+
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
    geom_vline(xintercept = treatment_period, linetype = "dashed", size = 1.5) +
    labs(title = "Development of grouped means over all simulations",
         x = "Year",
         y = "Dependent variable",
         color = "Group") +
    theme_minimal() +
    theme(legend.position = "right",
          plot.title = element_text(hjust = 0.5, size = 20, face = "bold"), 
          axis.title = element_text(size = 16, face = "bold"),  
          axis.text = element_text(size = 14)) +
    scale_color_manual(values = c("Treated mean" = "red", "Control mean" = "black")) #+
  # Set ylim
  #ylim(95, 110)
}

### Function to plot distribution of estimates
plot_estimates <- function(estimates_df) {
  ggplot(data = estimates_df) +
    geom_density(aes(x = did, color = "DiD"), size = 1) +
    geom_density(aes(x = sc, color = "SC"), size = 1) +
    geom_density(aes(x = sdid, color = "SDiD"), size = 1) +
    geom_vline(xintercept = treatment_effect, linetype = "dashed", size = 1.5) +
    labs(title = paste0("Distribution of estimates using ", iterations, " simulations"),
         x = "Estimate",
         y = "Density",
         color = "Method") +
    theme_minimal() +
    theme(legend.position = "right",
          plot.title = element_text(hjust = 0.5, size = 20, face = "bold"), 
          axis.title = element_text(size = 16, face = "bold"),  
          axis.text = element_text(size = 14),
          legend.title = element_text(face = "bold"),
          legend.text = element_text(size = 12)) +
    scale_color_manual(values = c("DiD" = "#E41A1C", "SC" = "#377EB8", "SDiD" = "#4DAF4A"),
                       labels = c("DiD" = "DiD", "SC" = "SC", "SDiD" = "SDiD")) +
    guides(color = guide_legend(title = "Estimator", title.position = "top", title.hjust = 0.5))
}


########## Other functions ##########

### Function to compute dynamic treatment effect
dynamic_treatment_effect <- function(treatment_effect, treatment_period, T) {
  # Define the number of periods after treatment
  num_periods <- T - treatment_period + 1
  # Compute the dynamic treatment effect using vectorized operations
  total_effect <- sum(treatment_effect * (1:num_periods))
  return(total_effect/num_periods)
}

########## Test-wise estimation ##########

# Set parameter values
#N <- 30
#T <- 30
#n_treated <- 5
#treatment_period <- 30
#trend <- 0.05
#control_trend <- 0.1
#treated_trend <- 0.1
#mean = 0
#sd = 1
#treatment_effect = 2

# Simulate and visualize data
#df = treatment_simulation(N, T, n_treated, treatment_period, control_trend, treated_trend, treatment_effect)
#plot_grouped(df)

# Define estimators to be used
#estimators = list(did=did_estimate,
                  #sc=sc_estimate,
                  #sdid=synthdid_estimate)

# Convert to required format
#setup = panel.matrices(df)

# Check structure
#head(df)

# Compute estimates
#estimates = lapply(estimators, function(estimator) { estimator(setup$Y,
                                                               #setup$N0, setup$T0) } )
#estimates

########## Application ##########

### Function to run each simulation and estimation multiple times

#results <- repeated_simulation(treatment_simulation, list(N, T, n_treated, treatment_period, control_trend, treated_trend, treatment_effect), 10)
#dfs <- results$dfs
#estimates_df <- results$estimates
#plot_estimates(estimates_df)
#plot_grouped_simulations(dfs)

# Compute column-wise mean of estimates
#colMeans(estimates_df)
