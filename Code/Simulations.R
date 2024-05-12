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

rm(list=ls())

########## Simulation functions ##########

### Function for stagnant dependent variable
stagnant_simulation <- function(N, T, mean = 0, sd = 1, n_treated = 1, treatment_period) {
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
  df[,1] <- rnorm(N, mean, sd)
  treated[,1] <- rnorm(n_treated, mean, sd)
  ##### Fill rest of columns with starting values without any change #####
  for (i in 2:T) {
    df[,i] <- df[,i-1]
    treated[,i] <- treated[,i-1]
  }
  ##### Code in between is where most difference lies #####
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

### Function with mean zero innovations at each time point
noise_simulation <- function(N, T, mean = 0, sd = 1, n_treated, treatment_period) {
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
  df[,1] <- rnorm(N, mean, sd)
  treated[,1] <- rnorm(n_treated, mean, sd)
  # Fill rest of columns as previous value + random draw
  for (i in 2:T) {
    df[,i] <- df[,i-1] + rnorm(N, mean, sd)
    treated[,i] <- treated[,i-1] + rnorm(n_treated, mean, sd)
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


### Function with mean zero innovations and identical linear time trend
same_trend_simulation <- function(N, T, mean = 0, sd = 1, n_treated, treatment_period, trend) {
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
  df[,1] <- rnorm(N, mean, sd)
  treated[,1] <- rnorm(n_treated, mean, sd)
  # Fill rest of columns as previous value + upwards time-trend added (%) and random draw
  for (i in 2:T) {
    df[,i] <- df[,i-1] + trend + rnorm(N, mean, sd)
    treated[,i] <- treated[,i-1] + trend + rnorm(n_treated, mean, sd)
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

### Function with mean. zero innovations and different linear time trend
different_trend_simulation <- function(N, T, mean = 0, sd = 1, n_treated, treatment_period, control_trend, treated_trend) {
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
  df[,1] <- rnorm(N, mean, sd)
  treated[,1] <- rnorm(n_treated, mean, sd)
  # Fill rest of columns as previous value + upwards time-trend added (%) and random draw
  for (i in 2:T) {
    df[,i] <- df[,i-1] + control_trend + rnorm(N, mean, sd)
    treated[,i] <- treated[,i-1] + treated_trend + rnorm(n_treated, mean, sd)
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

### Function with mean-zero innovations, identical linear time trend and treatment effect
treatment_simulation <- function(N, T, mean = 0, sd = 1, n_treated, treatment_period, trend, treatment_effect) {
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
  df[,1] <- rnorm(N, 100, sd) ###############################################################
  treated[,1] <- rnorm(n_treated, 100, sd)
  # Fill rest of columns as previous value + upwards time-trend added (%) and random draw and if treated, add treatment effect
  for (i in 2:T) {
    df[,i] <- df[,i-1] + trend + rnorm(N, mean, sd)
    treated[,i] <- treated[,i-1] + trend + rnorm(n_treated, mean, sd) + ifelse(i >= treatment_period, treatment_effect, 0)
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



# Test-wise simulation

# Set parameter values
N <- 300
T <- 15
n_treated <- 100
treatment_period <- 10


# Set values for distribution of initial draws
mean <- 0
sd <- 1

########## Plotting functions ##########

### All individual observations over time
plot_individual <- function(df) {
  # Subset treated units
  treated_units <- df %>% filter(grepl("Treated", Observation))
  # Subset control units
  control_units <- df %>% filter(!grepl("Treated", Observation))
  # Plot them separately to highlight treated observations
  ggplot() +
    geom_line(data = treated_units, aes(x = Time, y = value, color = Observation), size = 3) +
    geom_line(data = control_units, aes(x = Time, y = value, color = Observation), size = 1) +
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
    geom_line(aes(x = Time, y = control_avg, color = "Control mean"), size = 1) +
    geom_line(aes(x = Time, y = treated_avg, color = "Treated mean"), size = 1) +
    geom_vline(xintercept = treatment_period, linetype = "dashed", size = 1.5) +
    labs(title = "Development of means over time",
         x = "Year",
         y = "Dependent variable") +
    theme_minimal() +
    theme(legend.position = "right",
          plot.title = element_text(hjust = 0.5, size = 20), 
          axis.title = element_text(size = 16),  
          axis.text = element_text(size = 14)) +
    scale_color_manual(values = c("Treated mean" = "red", "Control mean" = "black")) #+
  # Set ylim
  #ylim(95, 110)
}

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
    geom_line(aes(x = Time, y = control_avg, color = "Control mean"), size = 1) +
    geom_line(aes(x = Time, y = treated_avg, color = "Treated mean"), size = 1) +
    geom_vline(xintercept = treatment_period, linetype = "dashed", size = 1.5) +
    labs(title = "Development of means over time",
         x = "Year",
         y = "Dependent variable") +
    theme_minimal() +
    theme(legend.position = "right",
          plot.title = element_text(hjust = 0.5, size = 20), 
          axis.title = element_text(size = 16),  
          axis.text = element_text(size = 14)) +
    scale_color_manual(values = c("Treated mean" = "red", "Control mean" = "black")) #+
  # Set ylim
  #ylim(95, 110)
}

########## Test-wise estimation ##########
# Real-life dataset for comparison
og_data <- read.csv2("Data/california_prop99.csv") 
og_data$PacksPerCapita <- as.numeric(og_data$PacksPerCapita)
# Rename California to treated, everything else remains the same
og_data$State <- ifelse(og_data$State == "California", "Treated", og_data$State)
og_data <- og_data %>% rename(value = PacksPerCapita, Time = Year, Observation = State)

treatment_period = 1989
plot_grouped(og_data)

# Set parameter values
N <- 30
T <- 30
n_treated <- 5
treatment_period <- 30
trend <- 0.05
control_trend <- -0.1
treated_trend <- 0.2
mean = 0
sd = 1
treatment_effect = 2

# Define estimators to be used
estimators = list(did=did_estimate,
                  sc=sc_estimate,
                  sdid=synthdid_estimate)

# Simulate data
df_stagnant = stagnant_simulation(N, T, mean, sd, n_treated, treatment_period)
df_noise = noise_simulation(N, T, mean, sd, n_treated, treatment_period)
df_same_trend = same_trend_simulation(N, T, mean, sd, n_treated, treatment_period, trend)
df_different_trend = different_trend_simulation(N, T, mean, sd, n_treated, treatment_period, control_trend, treated_trend)
df_treatment = treatment_simulation(N, T, mean, sd, n_treated, treatment_period, trend, treatment_effect)



# Convert to required format
setup_og = panel.matrices(og_data)
setup_2 = panel.matrices(df_stagnant)
setup_3 = panel.matrices(df_noise)
setup_4 = panel.matrices(df_same_trend)
setup_5 = panel.matrices(df_different_trend)
setup_6 = panel.matrices(df_treatment)

# Check structure
head(setup_og)
head(setup_2)
head(setup_3)
head(setup_4)
head(setup_5)
head(setup_6)


plot_grouped(og_data)
head(df_stagnant)
head(og_data)

# Compute estimates
estimates_og = lapply(estimators, function(estimator) { estimator(setup_og$Y,
                                                                  setup_og$N0, setup_og$T0) } )

#estimates2 = lapply(estimators, function(estimator) { estimator(setup_2$Y,
                                                               #setup_2$N0, setup_2$T0) } )

estimates3 = lapply(estimators, function(estimator) { estimator(setup_3$Y,
                                                               setup_3$N0, setup_3$T0) } )

estimates4 = lapply(estimators, function(estimator) { estimator(setup_4$Y,
                                                               setup_4$N0, setup_4$T0) } )

estimates5 = lapply(estimators, function(estimator) { estimator(setup_5$Y,
                                                               setup_5$N0, setup_5$T0) } )

estimates6 = lapply(estimators, function(estimator) { estimator(setup_6$Y,
                                                               setup_6$N0, setup_6$T0) } )

########## Application ##########

iterations = 100
estimates_df = data.frame(matrix(ncol = 3, nrow = iterations))
dfs = list()
for (i in 1:iterations) {
  # Simulate a dataframe
  data <- treatment_simulation(N,T, mean, sd, n_treated, treatment_period, trend, treatment_effect)
  # Conduct estimation
  setup = panel.matrices(data)
  estimates = lapply(estimators, function(estimator) { estimator(setup$Y,
                                                                      setup$N0, setup$T0) } )
  # Store estimates
  estimates_df[i, 1] = unlist(estimates)["did"]
  estimates_df[i, 2] = unlist(estimates)["sc"]
  estimates_df[i, 3] = unlist(estimates)["sdid"]
  # Rename columns
  colnames(estimates_df) <- c("did", "sc", "sdid")
  dfs[[i]] = data
}

# Function to plot distribution of estimates
plot_estimates <- function(estimates_df) {
  ggplot(data = estimates_df) +
    geom_density(aes(x = did, color = "DiD"), size = 1) +
    geom_density(aes(x = sc, color = "SC"), size = 1) +
    geom_density(aes(x = sdid, color = "SDiD"), size = 1) +
    geom_vline(xintercept = treatment_effect, linetype = "dashed", size = 1.5) +
    labs(title = "Distribution of estimates for different estimators using 1000 simulations",
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

plot_estimates(estimates_df)


# Compute column-wise mean of estimates
colMeans(estimates_df)








########################################
dfs[1]
test <- dfs[1]
test

plot_grouped(test)
class(test)
test2 <- as.data.frame(test)
test2

plot_grouped(test2)

test2
# Combine all simulations in one dataframe with meaned observations
all_data <- do.call(rbind, dfs)
all_data
plot_grouped(all_data)


plot_grouped_simulations(dfs)
