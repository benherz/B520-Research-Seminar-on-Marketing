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
iterations = 1000

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


noise_params <- list(N = T, T = T, n_treated = n_treated, treatment_period = treatment_period)

noise_results <- repeated_simulation(simulation_function = noise_simulation,
                                     simulation_params = noise_params, iterations = iterations)




