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