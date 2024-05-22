source("Functions.R")

# Load in real-life dataset as a reference
og_data <- read.csv2("Data/california_prop99.csv") 
og_data$PacksPerCapita <- as.numeric(og_data$PacksPerCapita)
og_data$State <- ifelse(og_data$State == "California", "Treated", og_data$State)
og_data <- og_data %>% rename(value = PacksPerCapita, Time = Year, Observation = State)

# Format for synthdid-package
setup_og = panel.matrices(og_data)

# Define estimators to be used
estimators = list(did=did_estimate,
                  sc=sc_estimate,
                  sdid=synthdid_estimate)

# Estimation process
estimates = lapply(estimators, function(estimator) { estimator(setup_og$Y,
                                                               setup_og$N0, setup_og$T0) } )

unlist(estimates)

# Run TWFE as a check
twfe(og_data, twfe_formula)

a = synthdid_estimate(setup_og$Y, setup_og$N0, setup_og$T0)
b = sc_estimate(setup_og$Y, setup_og$N0, setup_og$T0)
c = did_estimate(setup_og$Y, setup_og$N0, setup_og$T0)
a
b
c

###################################
# Set number of simulation-and-estimation iterations
iterations = 100

# Set parameter values
N <- 50
T <- 10
n_treated <- 10
treatment_period <- T - 2
treated_trend <- 0
control_trend <- 0
treatment_effect <- 1


contrast_params <- list(N = N, T = T, n_treated = n_treated, treatment_period = treatment_period, 
                     control_trend = control_trend, treated_trend = treated_trend)

# Apply simulation
contrast_results <- contrast_simulation(simulation_function = trend_simulation,
                                     simulation_params = contrast_params, iterations = iterations)

# Set up population sizes to loop over
population_sizes <- list(10, 50, 100)#, 250, 500, 1000, 2500, 5000, 10000, 25000)
contrast_results = tibble()
all_means = tibble()

for(population in population_sizes) {
  # 90% of population will serve as control, 10% will be treated
   contrast_params <- list(N = 0.9 * population, T = T, n_treated = 0.1 * population, treatment_period = treatment_period, 
                     control_trend = control_trend, treated_trend = treated_trend, treatment_effect = treatment_effect)
   results <- contrast_simulation(simulation_function = static_treatment_simulation,
                                     simulation_params = contrast_params, iterations = iterations)
   all_means <- rbind(all_means, colMeans(results))
   contrast_results <- rbind(contrast_results, results)

}
colnames(contrast_results) <- c("N", "SDiD", "SC", "TWFE")
colnames(all_means) <- c("N", "SDiD", "SC", "TWFE")

# Save to csv
write.csv(contrast_results, "contrast_results2.csv")
write.csv(all_means, "all_means.csv")

a <- contrast_results %>% group_by(N) %>% summarise_all(mean) %>% as.data.frame()



# Plot a
ggplot(a, aes(x = N)) + 
  geom_line(aes(y = SDiD, color = "SDiD"), size = 1.2) + 
  geom_line(aes(y = SC, color = "SC"), size = 1.2) + 
  geom_line(aes(y = TWFE, color = "TWFE"), size = 1.2) + 
  labs(
    title = "Duration of Estimation Processes Depending on Population Size",
    x = "Population Size",
    y = "Time (s)",
    color = "Method"
  ) + 
  scale_color_manual(values = c("SDiD" = "red", "SC" = "blue", "TWFE" = "green")) + 
  theme_minimal() + 
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    legend.position = "right"
  )


# 101 check 
df <- static_treatment_simulation(N = 900, T = T, n_treated = 100, treatment_period = treatment_period, 
                       control_trend = control_trend, treated_trend = treated_trend, treatment_effect = treatment_effect)

estimators = list(sdid=synthdid_estimate)

time <- Sys.time()
setup <- panel.matrices(df)
# estimate
lapply(estimators, function(estimator) { estimator(setup$Y, setup$N0, setup$T0) } )
time2 <- Sys.time()
time2 - time




##### Simulation without trend or treatment #####

# Set number of simulation-and-estimation iterations
iterations = 100

# Set parameter values
N <- 30
T <- 30
n_treated <- 5
treatment_period <- T 
treated_trend <- 0
control_trend <- 0

# Collect parameters in list
noise_params <- list(N = T, T = T, n_treated = n_treated, treatment_period = treatment_period, 
                     control_trend = control_trend, treated_trend = treated_trend)

# Apply simulation
noise_results <- repeated_simulation(simulation_function = trend_simulation,
                                     simulation_params = noise_params, iterations = iterations)
# Visualize data
plot_grouped_simulations(noise_results$dfs)

# Visualize results
treatment_effect <- 0
plot_estimates(noise_results$estimates)


##### Simulation with trend but no treatment #####

# Set additional parameter values
control_trend <- 0.1
treated_trend <- 0.1

# Collect parameters in list
trend_params <- list(N = T, T = T, n_treated = n_treated, treatment_period = treatment_period,
                     control_trend = control_trend, treated_trend = treated_trend)

# Apply simulation
trend_results <- repeated_simulation(simulation_function = trend_simulation,
                                     simulation_params = trend_params, iterations = iterations)

# Visualize data
plot_grouped_simulations(trend_results$dfs)

# Visualize results
treatment_effect <- 0
plot_estimates(trend_results$estimates)

# Compute mean and sd of estimates
colMeans(trend_results$estimates)
apply(trend_results$estimates, 2, sd)


##### Simulation with same trend and static treatment #####

# Dynamic treatment effect results as (1 + 2 + 3)/3 for treatment = 1 and 
# treatment_period = T - 2

# Dynamic treatment effect results as (2 + 4 + 6 + 8)/4 for treatment = 2 and
# treatment_period = T - 3


# Set additional parameter values
control_trend <- 0.1
treated_trend <- 0.1
treatment_effect <- 1
treatment_period <- T - 2
iterations <- 1000

# Collect parameters in list
treatment_params <- list(N = T, T = T, n_treated = n_treated, treatment_period = treatment_period,
                         control_trend = control_trend, treated_trend = treated_trend, treatment_effect = treatment_effect)

# Apply simulation
static_treatment_results <- repeated_simulation(simulation_function = static_treatment_simulation,
                                     simulation_params = treatment_params, iterations = iterations)

# Visualize data
a <- plot_grouped_simulations(static_treatment_results$dfs)

# Visualize results
b <- plot_estimates(static_treatment_results$estimates)

# Compute mean and sd of estimates
colMeans(static_treatment_results$estimates)
apply(static_treatment_results$estimates, 2, sd)

# Save plots
ggsave(filename = "static_treatment_results.png", plot = b, width = 10, height = 5, units = "in", dpi = 300)
ggsave(filename = "static_treatment_data.png", plot = a, width = 10, height = 5, units = "in", dpi = 300)

##### Simulation with different trend and static treatment #####

# Set additional parameter values
control_trend <- 0.1
treated_trend <- 0.2

# Collect parameters in list
diff_trend_params <- list(N = T, T = T, n_treated = n_treated, treatment_period = treatment_period,
                          control_trend = control_trend, treated_trend = treated_trend, treatment_effect = treatment_effect)

# Apply simulation
diff_trend_results <- repeated_simulation(simulation_function = static_treatment_simulation,
                                     simulation_params = diff_trend_params, iterations = iterations)

# Visualize data
plot_grouped_simulations(diff_trend_results$dfs)

# Visualize results
plot_estimates(diff_trend_results$estimates)

# Compute mean and sd of estimates
colMeans(diff_trend_results$estimates)
apply(diff_trend_results$estimates, 2, sd)



##### Simulation with same trend and dynamic treatment #####

# Set additional parameter values
control_trend <- 0.1
treated_trend <- 0.1
treatment_effect <- 1
treatment_period <- T - 2

# Collect parameters in list
diff_params <- list(N = T, T = T, n_treated = n_treated, treatment_period = treatment_period,
                    control_trend = control_trend, treated_trend = treated_trend, treatment_effect = treatment_effect)

# Apply simulation
dynamic_treatment_results <- repeated_simulation(simulation_function = dynamic_treatment_simulation,
                                     simulation_params = diff_params, iterations = iterations)

# Visualize data
plot_grouped_simulations(dynamic_treatment_results$dfs)

# Visualize results
plot_estimates(dynamic_treatment_results$estimates)

# Compute mean and sd of estimates
colMeans(dynamic_treatment_results$estimates)
apply(dynamic_treatment_results$estimates, 2, sd)

##### Simulation with different pre-trend and static treatment #####

# Set additional parameter values
control_trend <- 0.1
treated_trend <- - 0.1
treatment_effect <- 1
treatment_period <- T - 2

# Collect parameters in list
diff_pre_trend_params <- list(N = T, T = T, n_treated = n_treated, treatment_period = treatment_period,
                          control_trend = control_trend, treated_trend = treated_trend, treatment_effect = treatment_effect)

# Apply s
diff_trend_results <- repeated_simulation(simulation_function = static_treatment_simulation,
                                     simulation_params = diff_trend_params, iterations = iterations)

# Visualize data
c <- plot_grouped_simulations(diff_trend_results$dfs)

# Visualize results
d <- plot_estimates(diff_trend_results$estimates)

# Compute mean and sd of estimates
colMeans(diff_trend_results$estimates)
apply(diff_trend_results$estimates, 2, sd)

ggsave(filename = "diff_pre_trend_results.png", plot = d, width = 10, height = 5, units = "in", dpi = 300)
ggsave(filename = "diff_pre_trend_data.png", plot = c, width = 10, height = 5, units = "in", dpi = 300)

c
d
##### Simulation with violation of parallel-trends assumption #####

# Set parameter values
T <- 30
N <- 30
n_treated <- 5
treatment_period <- T - 10
control_trend <- - 0.1
treated_trend <- 0.1
treatment_effect <- 1
iterations <- 1000

# Collect parameters in list
diff_trend_params <- list(N = T, T = T, n_treated = n_treated, treatment_period = treatment_period,
                          control_trend = control_trend, treated_trend = treated_trend, treatment_effect = treatment_effect)

# Apply simulation
diff_trend_results <- repeated_simulation(simulation_function = diff_trend_simulation,
                                     simulation_params = diff_trend_params, iterations = iterations)

# Visualize data
plot_grouped_simulations(diff_trend_results$dfs)

# Visualize results
plot_estimates(diff_trend_results$estimates)

# Compute mean and sd of estimates
colMeans(diff_trend_results$estimates)
apply(diff_trend_results$estimates, 2, sd)



##### Simulation with heterogeneous static treatment effect #####

# Set parameter values
T <- 30
N <- 30
n_treated <- 5
treatment_period <- T - 10
control_trend <- 0.1
treated_trend <- 0.1
treatment_effect <- 1
iterations <- 100


range <- 1:10
dynamic_treatment_effect(range, 1, 10)

# Collect parameters in list
heterogeneous_static_params <- list(N = T, T = T, n_treated = n_treated, treatment_period = treatment_period,
                          control_trend = control_trend, treated_trend = treated_trend, treatment_effect = treatment_effect)

# Apply simulation
heterogeneous_static_results <- repeated_simulation(simulation_function = heterogeneous_static_simulation,
                                     simulation_params = heterogeneous_static_params, iterations = iterations)

# Visualize data
plot_grouped_simulations(heterogeneous_static_results$dfs)

# Visualize results
plot_estimates(heterogeneous_static_results$estimates)

# Compute mean and sd of estimates
colMeans(heterogeneous_dynamic_results$estimates)
apply(heterogeneous_dynamic_results$estimates, 2, sd)


##### Simulation with heterogeneous dynamic treatment effect #####

# Set parameter values
T <- 30
N <- 30
n_treated <- 5
treatment_period <- T - 10
control_trend <- 0.1
treated_trend <- 0.1

dynamic_treatment_effect(treatment_effect, 1, 10)
  
# Collect parameters in list
heterogeneous_dynamic_params <- list(N = T, T = T, n_treated = n_treated, treatment_period = treatment_period,
                          control_trend = control_trend, treated_trend = treated_trend, treatment_effect = treatment_effect)

# Apply simulation
start_time <- Sys.time()
heterogeneous_dynamic_results <- repeated_simulation(simulation_function = heterogeneous_dynamic_unit_simulation,
                                     simulation_params = heterogeneous_dynamic_params, iterations = iterations)
end_time <- Sys.time()
time_elapsed <- end_time - start_time
time_elapsed

# Visualize data
plot_grouped_simulations(heterogeneous_dynamic_results$dfs)

# Visualize results
plot_estimates(heterogeneous_dynamic_results$estimates)

# Compute mean and sd of estimates
colMeans(heterogeneous_dynamic_results$estimates)
apply(heterogeneous_dynamic_results$estimates, 2, sd)
df <- heterogeneous_dynamic_simulation(N = 1, T = 30, n_treated = 20, 
      treatment_period = 10, control_trend = 0.1, treated_trend = 0.1, treatment_effect = 1)


plot_individual(df)


df <- heterogeneous_static_simulation(N = 1, T = 30, n_treated = 20, 
      treatment_period = 10, control_trend = 0.1, treated_trend = 0.1, treatment_effect = 1)
df2 <- heterogeneous_dynamic_simulation(N = 1, T = 30, n_treated = 20, 
      treatment_period = 10, control_trend = 0.1, treated_trend = 0.1, treatment_effect = 1)
df3 <- heterogeneous_dynamic_unit_simulation(N = 1, T = 30, n_treated = 20, 
      treatment_period = 10, control_trend = 0.1, treated_trend = 0.1, treatment_effect = 1)
plot_individual(df)
plot2 <- plot_individual(df2)
plot3 <- plot_individual(df3)
plot_individual(df3)

plot2
plot3

treatment_period <- 10
df <- diff_trend_simulation(N = 100, T = 50, n_treated = 30, 
      treatment_period = 10, control_trend = 0.1, treated_trend = 0.1, treatment_effect = 1)
plot_grouped(df)
