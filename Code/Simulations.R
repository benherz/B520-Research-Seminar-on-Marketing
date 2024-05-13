rm(list=ls())
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
iterations <- 200

# Collect parameters in list
treatment_params <- list(N = T, T = T, n_treated = n_treated, treatment_period = treatment_period,
                         control_trend = control_trend, treated_trend = treated_trend, treatment_effect = treatment_effect)

# Apply simulation
static_treatment_results <- repeated_simulation(simulation_function = static_treatment_simulation,
                                     simulation_params = treatment_params, iterations = iterations)

# Visualize data
plot_grouped_simulations(static_treatment_results$dfs)

# Visualize results
plot_estimates(static_treatment_results$estimates)

# Compute mean and sd of estimates
colMeans(static_treatment_results$estimates)
apply(static_treatment_results$estimates, 2, sd)


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



treatment_period <- 10





df <- heterogeneous_static_simulation(N = 1, T = 30, n_treated = 20, 
      treatment_period = 10, control_trend = 0.1, treated_trend = 0.1, treatment_effect = 1)
df2 <- heterogeneous_dynamic_simulation(N = 1, T = 30, n_treated = 20, 
      treatment_period = 10, control_trend = 0.1, treated_trend = 0.1, treatment_effect = 1)
df3 <- heterogeneous_dynamic_unit_simulation(N = 1, T = 30, n_treated = 20, 
      treatment_period = 10, control_trend = 0.1, treated_trend = 0.1, treatment_effect = 1)
plot_individual(df)
plot2 <- plot_individual(df2)
plot3 <- plot_individual(df3)


plot2
plot3
