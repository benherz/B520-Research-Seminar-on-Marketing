library(future)
library(doFuture)
library(future.batchtools)


##### Code to plot development of individual observations over time
individual_plot <- ggplot(data = df, aes(x = Time, y = value, color = Observation)) +
  geom_line() +
  geom_vline(xintercept = treatment_period, linetype = "dashed", size = 1.5) +
  labs(title = "Development of simulated data over time",
       x = "Year",
       y = "Size") +
  theme_minimal() +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5, size = 20), 
        axis.title = element_text(size = 16),  
        axis.text = element_text(size = 14)) +  
  scale_color_viridis_d() 

individual_plot



##### Code to plot means of control and treated
plot_grouped <- ggplot(data = average_data) +
  geom_line(aes(x = Time, y = control_avg, color = "Control mean")) +
  geom_line(aes(x = Time, y = treated_avg, color = "Treated mean")) +
  geom_vline(xintercept = treatment_period, linetype = "dashed", size = 1.5) +
  labs(title = "Development of means over time",
       x = "Year",
       y = "Size") +
  theme_minimal() +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5, size = 20), 
        axis.title = element_text(size = 16),  
        axis.text = element_text(size = 14)) +
  scale_color_manual(values = c("Control mean" = "red", "Treated mean" = "blue"))



##### Code to plot density distribution of estimates
test <- melt(estimates_df, id.vars = c("did", "sc", "sdid"))
estimates_df2 <- melt(estimates_df)
ggplot(data = estimates_df2) +
  geom_density(aes(x = value, color = variable), size = 1) +
  geom_vline(xintercept = treatment_effect, linetype = "dashed", size = 1.5) +
  labs(title = "Distribution of estimates for different estimators using 1000 noise simulations",
       x = "Estimate",
       y = "Density",
       color = "Method") +
  theme_minimal() +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5, size = 20), 
        axis.title = element_text(size = 16),  
        axis.text = element_text(size = 14)) +
  scale_color_manual(values = c("did" = "red", "sc" = "blue", "sdid" = "green"))


##### Function for stagnant dependent variable
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


##### Code to run a simulation function multiple times
iterations = 10
estimates_df = data.frame(matrix(ncol = 3, nrow = iterations))
dfs = list()
for (i in 1:iterations) {
  # Simulate a dataframe
  data <- treatment_simulation(N, T, n_treated, treatment_period, control_trend, treated_trend, treatment_effect)
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







##### Weird function for iptw matching 
perform_iptw_matching <- function(data) {
  
  # Conclude treatment period
  treatment_period <- min(data$Time[data$treated == 1])
  
  # Compute overall mean
  overall_mean <- mean(data$value)
  
  # Compute pre-treatment means for treated units
  treated_means <- data %>% filter(grepl("Treated", Observation) & Time < treatment_period) %>%
    group_by(Observation) %>% summarize(pre_mean = mean(value))
  
  # Compute pre-treatment means for control units
  control_means <- data %>% filter(!grepl("Treated", Observation) & Time < treatment_period) %>%
    group_by(Observation) %>% summarize(pre_mean = mean(value))
  
  # Compute inverse probability of receiving treatment as distance from overall mean
  treated_means$iptw <- 1 / abs(treated_means$pre_mean - overall_mean)
  control_means$iptw <- 1 / (1- abs(control_means$pre_mean - overall_mean))
  
  # Perform logistic activation (probability should be between 0 and 1)
  treated_means$activation <- 1 / (1 + exp(-treated_means$iptw))
  control_means$activation <- 1 / (1 + exp(-control_means$iptw))
  
  # Merge data
  mean_data <- rbind(treated_means, control_means)
  merged_data <- merge(data, mean_data, by = "Observation") #%>% 
  #select(-pre_mean, -iptw) 
  
  # Apply IPTW weights
  merged_data$value <- merged_data$value * merged_data$activation
  
  # Select only weighted data
  weighted_data <- merged_data# %>% select(-activation)
  
  # Return data with IPTW weights applied
  return(weighted_data)
}

test <- perform_iptw_matching(data)











###################################
# Code to test how population sizes affect estimates
##### N = 10, n_Treated = 2 #####
N <- 10
n_treated <- 2
simulation_params <- list(N = N, T = T, n_treated = n_treated, treatment_period = treatment_period,
                          pre_control_trend = pre_control_trend, pre_treated_trend = pre_treated_trend,
                          post_control_trend = post_control_trend, post_treated_trend = post_treated_trend,
                          treatment_effect = treatment_effect)
# Simulation & estimation
results_10_2 <- repeated_simulation(simulation_function = treatment_simulation,
                                    simulation_params = simulation_params, iterations = iterations)
estimates_10_2 <- results_10_2$estimates
df_10_2 <- results_10_2$dfs[[index]]

write.csv(df_10_2, paste0(path, "/Output/Data/df_10_2.csv"))
write.csv(estimates_10_2, paste0(path, "/Output/Data/estimates_10_2.csv"))

plot_estimates(estimates_10_2) %>%
  ggsave(filename = paste0(path, "/Output/Plots/estimates_plot_10_2.png"), width = 10, height = 5, units = "in", dpi = 300)
plot_grouped_simulations(results_10_2$dfs) %>%
  ggsave(filename = paste0(path, "/Output/Plots/data_plot_10_2.png"), width = 10, height = 5, units = "in", dpi = 300)


##### N = 50, n_Treated = 10 #####
N <- 50
n_treated <- 10
simulation_params <- list(N = N, T = T, n_treated = n_treated, treatment_period = treatment_period,
                          pre_control_trend = pre_control_trend, pre_treated_trend = pre_treated_trend,
                          post_control_trend = post_control_trend, post_treated_trend = post_treated_trend,
                          treatment_effect = treatment_effect)
# Simulation & estimation
results_50_10 <- repeated_simulation(simulation_function = treatment_simulation,
                                     simulation_params = simulation_params, iterations = iterations)
estimates_50_10 <- results_50_10$estimates
df_50_10 <- results_50_10$dfs[[index]]

write.csv(df_50_10, paste0(path, "/Output/Data/df_50_10.csv"))
write.csv(estimates_50_10, paste0(path, "/Output/Data/estimates_50_10.csv"))

plot_estimates(estimates_50_10) %>%
  ggsave(filename = paste0(path, "/Output/Plots/estimates_plot_50_10.png"), width = 10, height = 5, units = "in", dpi = 300)
plot_grouped_simulations(results_50_10$dfs) %>%
  ggsave(filename = paste0(path, "/Output/Plots/data_plot_50_10.png"), width = 10, height = 5, units = "in", dpi = 300)



### N = 100, n_Treated = 20 ###
N <- 100
n_treated <- 20
simulation_params <- list(N = N, T = T, n_treated = n_treated, treatment_period = treatment_period,
                          pre_control_trend = pre_control_trend, pre_treated_trend = pre_treated_trend,
                          post_control_trend = post_control_trend, post_treated_trend = post_treated_trend,
                          treatment_effect = treatment_effect)
# Simulation & estimation
results_100_20 <- repeated_simulation(simulation_function = treatment_simulation,
                                      simulation_params = simulation_params, iterations = iterations)
estimates_100_20 <- results_100_20$estimates
df_100_20 <- results_100_20$dfs[[index]]

write.csv(df_100_20, paste0(path, "/Output/Data/df_100_20.csv"))
write.csv(estimates_100_20, paste0(path, "/Output/Data/estimates_100_20.csv"))

plot_estimates(estimates_100_20) %>%
  ggsave(filename = paste0(path, "/Output/Plots/estimates_plot_100_20.png"), width = 10, height = 5, units = "in", dpi = 300)
plot_grouped_simulations(results_100_20$dfs) %>%
  ggsave(filename = paste0(path, "/Output/Plots/data_plot_100_20.png"), width = 10, height = 5, units = "in", dpi = 300)



### N = 500, n_Treated = 100 ###
N <- 500
n_treated <- 100
simulation_params <- list(N = N, T = T, n_treated = n_treated, treatment_period = treatment_period,
                          pre_control_trend = pre_control_trend, pre_treated_trend = pre_treated_trend,
                          post_control_trend = post_control_trend, post_treated_trend = post_treated_trend,
                          treatment_effect = treatment_effect)
# Simulation & estimation
results_500_100 <- repeated_simulation(simulation_function = treatment_simulation,
                                       simulation_params = simulation_params, iterations = iterations)
estimates_500_100 <- results_500_100$estimates
df_500_100 <- results_500_100$dfs[[index]]

write.csv(df_500_100, paste0(path, "/Output/Data/df_500_100.csv"))
write.csv(estimates_500_100, paste0(path, "/Output/Data/estimates_500_100.csv"))

plot_estimates(estimates_500_100) %>% 
  ggsave(filename = paste0(path, "/Output/Plots/estimates_plot_500_100.png"), width = 10, height = 5, units = "in", dpi = 300)
plot_grouped_simulations(results_500_100$dfs) %>%
  ggsave(filename = paste0(path, "/Output/Plots/data_plot_500_100.png"), width = 10, height = 5, units = "in", dpi = 300)

### N = 1000, n_Treated = 200 ###
N <- 1000
n_treated <- 200
simulation_params <- list(N = N, T = T, n_treated = n_treated, treatment_period = treatment_period,
                          pre_control_trend = pre_control_trend, pre_treated_trend = pre_treated_trend,
                          post_control_trend = post_control_trend, post_treated_trend = post_treated_trend,
                          treatment_effect = treatment_effect)
# Simulation & estimation
results_1000_200 <- repeated_simulation(simulation_function = treatment_simulation,
                                        simulation_params = simulation_params, iterations = iterations)
estimates_1000_200 <- results_1000_200$estimates
df_1000_200 <- results_1000_200$dfs[[index]]

write.csv(df_1000_200, paste0(path, "/Output/Data/df_1000_200.csv"))
write.csv(estimates_1000_200, paste0(path, "/Output/Data/estimates_1000_200.csv"))

plot_estimates(estimates_1000_200) %>%
  ggsave(filename = paste0(path, "/Output/Plots/estimates_plot_1000_200.png"), width = 10, height = 5, units = "in", dpi = 300)
plot_grouped_simulations(results_1000_200$dfs) %>%
  ggsave(filename = paste0(path, "/Output/Plots/data_plot_1000_200.png"), width = 10, height = 5, units = "in", dpi = 300)


### N = 5000, n_Treated = 1000 ###
N <- 5000
n_treated <- 1000

simulation_params <- list(N = N, T = T, n_treated = n_treated, treatment_period = treatment_period,
                          pre_control_trend = pre_control_trend, pre_treated_trend = pre_treated_trend,
                          post_control_trend = post_control_trend, post_treated_trend = post_treated_trend,
                          treatment_effect = treatment_effect)
# Simulation & estimation
results_5000_1000 <- repeated_simulation(simulation_function = treatment_simulation,
                                         simulation_params = simulation_params, iterations = iterations)
estimates_5000_1000 <- results_5000_1000$estimates
df_5000_1000 <- results_5000_1000$dfs[[index]]

write.csv(df_5000_1000, paste0(path, "/Output/Data/df_5000_1000.csv"))
write.csv(estimates_5000_1000, paste0(path, "/Output/Data/estimates_5000_1000.csv"))

plot_estimates(estimates_5000_1000) %>%
  ggsave(filename = paste0(path, "/Output/Plots/estimates_plot_5000_1000.png"), width = 10, height = 5, units = "in", dpi = 300)
plot_grouped_simulations(results_5000_1000$dfs) %>%
  ggsave(filename = paste0(path, "/Output/Plots/data_plot_5000_1000.png"), width = 10, height = 5, units = "in", dpi = 300)


#######################################
# Dynamic treatment effect random treatment
### 1.2 Dynamic treatment 

#### 1.2.1 Effect of sample size on estimates

```{r}
# Set seed anew, because chunks will be run separately
set.seed(100)

# Set parameters for simulation
T <- 10
n_treated <- 1
treatment_effect <- 1
control_sizes <- c(5, 10, 50, 100, 500, 1000)

# Set treatment effect to be dynamic (applied repeatedly in each period >= treatment_period)
effect = "dynamic"

# Treatment happens in T = 8
treatment_period <- T - 2

# Initialize df to store info about estimates per iteration
means <- data.frame(matrix(ncol = 5, nrow = length(control_sizes)))
sds <- data.frame(matrix(ncol = 5, nrow = length(control_sizes)))
colnames(means) <- c("Effect", "N", "did", "sc", "sdid")
colnames(sds) <- c("Effect", "N","did", "sc", "sdid")

for(i in 1:length(control_sizes)){
  
  # Grab N and n_treated 
  N <- control_sizes[i]
  
  # Set up list of simulation parameters
  simulation_params <- list(N = control_sizes[i], T = T, n_treated = n_treated, treatment_period = treatment_period,
                            pre_control_trend = pre_control_trend, pre_treated_trend = pre_treated_trend,
                            post_control_trend = post_control_trend, post_treated_trend = post_treated_trend,
                            treatment_effect = treatment_effect, effect = effect)
  
  # Simulation & estimation
  results <- repeated_simulation(simulation_function = treatment_simulation,
                                 simulation_params = simulation_params, iterations = iterations)
  estimates <- results$estimates
  df <- results$dfs[[index]]
  
  # Compute mean and sd of estimates and add average treatment effect
  effect <- dynamic_effect(treatment_effect, treatment_period, T)
  means[i, ] <- c(effect, N, colMeans(estimates))
  sds[i, ] <- c(effect, N, apply(estimates, 2, sd))
  
  write.csv(df, paste0(path, "/Output/Data/df_N_", N, "_ho_dyn_random.csv"))
  write.csv(estimates, paste0(path, "/Output/Data/estimates_N_", N, "_ho_dyn_random.csv"))
  
  plot_estimates(estimates) %>%
    ggsave(filename = paste0(path, "/Output/Plots/estimates_plot_N_", N, "_ho_dyn_random.png"), width = 10, height = 5, units = "in", dpi = 300)
  plot_grouped_simulations(results$dfs) %>%
    ggsave(filename = paste0(path, "/Output/Plots/data_plot_N_", N, "_ho_dyn_random.png"), width = 10, height = 5, units = "in", dpi = 300)
  # Print progress info
  print(paste0("Finished simulation for N = ", N, " and n_treated = ", n_treated, " at ", Sys.time()))
}

write.csv(means, paste0(path, "/Output/Data/means_per_N_ho_dyn_random.csv"))
write.csv(sds, paste0(path, "/Output/Data/sds_per_N_ho_dyn_random.csv"))
```


#### 1.2.2 Effect of pre-treatment periods on estimates

```{r}
set.seed(100)
N <- 8
T <- 10
n_treated <- 2
treatment_effect <- 1
iterations <- 1000

# Setting a lower value than 3 is not possible 
treatment_periods = c(3, 4, 5, 6, 7, 8, 9, 10)

# Initialize df to store info about estimates per iteration
means <- data.frame(matrix(ncol = 5, nrow = length(treatment_periods)))
sds <- data.frame(matrix(ncol = 5, nrow = length(treatment_periods)))
colnames(means) <- c("Effect", "treatment_period", "did", "sc", "sdid")
colnames(sds) <- c("EFfect", "treatment_period", "did", "sc", "sdid")

for(treatment_period in treatment_periods){
  simulation_params <- list(N = N, T = T, n_treated = n_treated, treatment_period = treatment_period,
                            pre_control_trend = pre_control_trend, pre_treated_trend = pre_treated_trend,
                            post_control_trend = post_control_trend, post_treated_trend = post_treated_trend,
                            treatment_effect = treatment_effect, effect = effect)
  
  # Extract N and n_treated for file naming
  N <- simulation_params$N
  n_treated <- simulation_params$n_treated
  
  # Simulation & estimation
  results <- repeated_simulation(simulation_function = treatment_simulation,
                                 simulation_params = simulation_params, iterations = iterations)
  # Assign results
  estimates <- results$estimates
  df <- results$dfs[[index]]
  
  # Compute mean and sd of estimates
  effect <- dynamic_effect(treatment_effect, treatment_period, T)
  means[treatment_period, ] <- c(effect, treatment_period, colMeans(estimates))
  sds[treatment_period, ] <- c(effect, treatment_period, apply(estimates, 2, sd))
  
  write.csv(df, paste0(path, "/Output/Data/df_treatment_period_", treatment_period, "_ho_dyn_random.csv"), )
  write.csv(estimates, paste0(path, "/Output/Data/estimates_treatment_period_", treatment_period,"_ho_dyn_random.csv"))
  
  plot_estimates(estimates) %>%
    ggsave(filename = paste0(path, "/Output/Plots/estimates_plot_treatment_period_",
                             treatment_period, "_ho_dyn_random.png"), width = 10, height = 5, units = "in", dpi = 300)
  plot_grouped_simulations(results$dfs) %>%
    ggsave(filename = paste0(path, "/Output/Plots/data_plot_treatment_period_",
                             treatment_period, "_ho_dyn_random.png"), width = 10, height = 5, units = "in", dpi = 300)
}

# Collect means and sd data
write.csv(means, paste0(path, "/Output/Data/means_treatment_period_", treatment_period, "_ho_dyn_random.csv"))
write.csv(sds, paste0(path, "/Output/Data/sds_treatment_period_", treatment_period, "_ho_dyn_random.csv"))
```


#### 1.2.3 Effect of heterogeneity (in treatment effect) on estimates


```{r}
set.seed(100)

# Set up new parameter configurations
heterogeneity <- c(1, 2, 3, 4, 5)
# Set overall sample size to 100, to allow for 2% and 5% share of treated units
N <- 20
treated_sizes <- c(1, 2)

# Set overall sample size to 100, to allow for 5% and 10% share of treated units
N <- 20
treated_sizes <- c(1, 2)

# Initialize df to store info about estimates per iteration
means <- data.frame(matrix(ncol = 6, nrow = combinations))
sds <- data.frame(matrix(ncol = 6, nrow = combinations))
colnames(means) <- c("Effect", "Heterogeneity","treated", "did", "sc", "sdid")
colnames(sds) <- c("Effect", "Heterogeneity", "treated", "did", "sc", "sdid")

for(treated_size in treated_sizes){
  
  n_treated <- treated_size
  
  for (deviation in heterogeneity){
    
    deviation <- deviation
    
    simulation_params <- list(N = N, T = T, n_treated = n_treated, treatment_period = treatment_period,
                              pre_control_trend = pre_control_trend, pre_treated_trend = pre_treated_trend,
                              post_control_trend = post_control_trend, post_treated_trend = post_treated_trend,
                              treatment_effect = treatment_effect, deviation = deviation, effect = effect)
    
    # Simulation & estimation
    results <- repeated_simulation(simulation_function = heterogeneous_treatment_simulation,
                                   simulation_params = simulation_params, iterations = iterations)
    # Assign results
    estimates <- results$estimates
    df <- results$dfs[[index]]
    
    # Compute mean and sd of estimates
    effect <- dynamic_effect(treatment_effect, treatment_period, T)
    means[treatment_period, ] <- c(effect, heterogeneity, n_treated, colMeans(estimates))
    sds[treatment_period, ] <- c(effect, heterogeneity, n_treated, apply(estimates, 2, sd))
    
    write.csv(df, paste0(path, "/Output/Data/df_dev_", heterogeneity, "_n_treated_", n_treated ,"dyn_random.csv"), )
    write.csv(estimates, paste0(path, "/Output/Data/estimates_dev_", heterogeneity, "_n_treated_", n_treated, "dyn_random.csv"))
    
    plot_estimates(estimates) %>%
      ggsave(filename = paste0(path, "/Output/Plots/estimates_plot_dev_",
                               heterogeneity, "_n_treated_", n_treated, "dyn_random.png"), width = 10, height = 5, units = "in", dpi = 300)
    plot_grouped_simulations(results$dfs) %>%
      ggsave(filename = paste0(path, "/Output/Plots/data_plot_dev_",
                               heterogeneity, "_n_treated_", n_treated, "dyn_random.png"), width = 10, height = 5, units = "in", dpi = 300)
    # Print progress info
    print(paste0("Finished simulation for deviation = ", deviation, " and n_treated = ", n_treated, " at ", Sys.time()))
  }
}

# Collect means and sd data
write.csv(means, paste0(path, "/Output/Data/means_per_N_he_dyn_random.csv"))
write.csv(sds, paste0(path, "/Output/Data/sds_per_N_he_dyn_random.csv"))
```
