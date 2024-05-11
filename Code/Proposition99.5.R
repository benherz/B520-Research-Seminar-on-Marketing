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

# Clear working space
rm(list=ls())

# Load in data and get variables of interest
dataset <- read.csv2("Data/california_prop99.csv") 
dataset$PacksPerCapita <- as.numeric(dataset$PacksPerCapita)
california <- dataset %>% filter(State == "California")
data <- dataset %>% filter(State != "California")

States <- unique(data$State)
data$PacksPerCapita <- as.numeric(data$PacksPerCapita)
california$PacksPerCapita <- as.numeric(california$PacksPerCapita)

# Min and max consumption as bounds for random draws of initial consumption
min_consumption <- min(data$PacksPerCapita) %>% as.numeric()
max_consumption <- max(data$PacksPerCapita) %>% as.numeric()

# Compute rate of consumption increase per year per state (in %)
diffs <- data.frame()
for (state in States){
  df = data %>% filter(State == state)
  for (i in 2:length(df$PacksPerCapita) -1){
    diff <- (df[i, "PacksPerCapita"] / df[i-1, "PacksPerCapita"])
    diffs[i -1, state] <- diff
  }
  }

# Compute range of increases as bounds for random draws of increment per period
min_increase_total <- min(diffs)
max_increase_total <- max(diffs)

# Compute min and max increase per state for more realistic depiction in data
min_increase_per_state <- apply(diffs, 2, min)
max_increase_per_state <- apply(diffs, 2, max)
increase_df <- rbind(min_increase_per_state, max_increase_per_state) %>% 
  as.data.frame()


# Create simulated data frame
T <- 31

df <- data.frame()
for (state in States) {
  df[state,1] <- runif(1, min = min_consumption, max = max_consumption) 
    for (t in 2:T){
      df[state,t] <- runif(1, min = increase_df["min_increase_per_state", state],
                           max = increase_df["max_increase_per_state", state]) *
        df[state, t-1]
    }
  }













###################################
simulation <- function(N) {
  datasets <- vector("list", length = N)
  for (n in 1:N) {
    df <- data.frame(row.names = States)
    for (state in States) {
      df[state,1] <- runif(1, min = min_consumption, max = max_consumption) 
      for (t in 2:T){
        df[state,t] <- runif(1, min = increase_df["min_increase_per_state", state],
                             max = increase_df["max_increase_per_state", state]) *
          df[state, t-1]
      }
    }
    datasets[n] <- df
  }
  return(datasets)
}
###################################

# Rename columns 
colnames(df) <- paste(min(data$Year):max(data$Year))
colnames(df)

# Add row names as column for id 
df$State <- rownames(df)

# Melt to long format
df_long <- melt(df, id.vars = "State")

# Rename columns
df_long <- df_long %>%
  rename(Year = variable, PacksPerCapita = value)

# Add treatment column
df_long$treated <- 0

# Change variables to correct format
df_long$Year <- as.integer(as.character(df_long$Year))
df_long$treated <- as.integer(df_long$treated)

# Check data structure
head(df_long)
head(california)

str(df_long)
str(california)

# Combine dataframes
#final_df <- rbind(df_long, california) %>% order(df_long$Year)


############################
# Combine the dataframes
combined_df <- rbind(df_long, california)

# Order the dataframe by 'Year'
combined_df <- combined_df[order(combined_df$Year), ]
head(combined_df)

treated_unit <- "California"
dependent_variable <- "PacksPerCapita"


### Plot development of dependent variable over time for each observation
synth_plot <- ggplot(filter(combined_df, State != treated_unit), aes(x = Year, y = PacksPerCapita, color = State)) +
  geom_line() +
  # Add California separately to highlight line
  geom_vline(xintercept = 1989, linetype = "dashed", size = 1.5) + 
  #geom_text(data = NULL, aes(x = 1988, y = 500, label = "Proposition 99", color = "black"), angle = 90, vjust = 1, size = 4) +
  geom_line(data = filter(combined_df, State == treated_unit), size = 2) +
  labs(title = "Development of PacksPerCapita over Time in simulated dataset",
       x = "Year",
       y = "Packs Per Capita") +
  theme_minimal() +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5, size = 20),  # Adjust title font size
        axis.title = element_text(size = 16),  # Adjust axis labels font size
        axis.text = element_text(size = 14)) +  
  scale_color_viridis_d() 
head(combined_df)

###  Original df 
og_plot <- ggplot(filter(dataset, State != treated_unit), aes(x = Year, y = PacksPerCapita, color = State)) +
  geom_line() +
  # Add California separately to highlight line
  geom_vline(xintercept = 1989, linetype = "dashed", size = 1.5) + 
  #geom_text(data = NULL, aes(x = 1988, y = 500, label = "Proposition 99", color = "black"), angle = 90, vjust = 1, size = 4) +
  geom_line(data = filter(dataset, State == treated_unit), size = 2) +
  labs(title = "Development of PacksPerCapita over Time in actual dataset",
       x = "Year",
       y = "Packs Per Capita") +
  theme_minimal() +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5, size = 20),  # Adjust title font size
        axis.title = element_text(size = 16),  # Adjust axis labels font size
        axis.text = element_text(size = 14)) +  
  scale_color_viridis_d() 



# Display both plots 
synth_plot + og_plot


### Plot California vs. equally weighted control group 
equal_control <- filter(dataset, State != "California") %>% group_by(Year) %>%
  summarise(Control = mean(PacksPerCapita))

equal_control_plot_data <- cbind(equal_control, filter(dataset, State == "California")["PacksPerCapita"]) %>%
  rename(California =  PacksPerCapita)


ggplot(data = equal_control_plot_data) +
  geom_line(aes(x = Year, y = Control, color = "Control"), size = 1) +  # Assign color label for legend
  geom_line(aes(x = Year, y = California, color = "California"), size = 1) +  # Assign color label for legend
  geom_vline(xintercept = 1989, linetype = "dashed", size = 1) +
  labs(title = "Development of PacksPerCapita in equally weighted control group vs. California",
       x = "Year",
       y = "Packs Per Capita",
       color = "Line") +  # Adjust legend title
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5, size = 16),  # Adjust title font size
        axis.title = element_text(size = 12),  # Adjust axis labels font size
        axis.text = element_text(size = 10))

###################################
# Define estimators to be used
estimators = list(did=did_estimate,
                  sc=sc_estimate,
                  sdid=synthdid_estimate)

# Convert from panel data to matrix
setup = panel.matrices(combined_df)
head(setup)

# Compute estimates
estimates = lapply(estimators, function(estimator) { estimator(setup$Y,
                                                               setup$N0, setup$T0) } )

head(estimates)



##############################
# Compute standard errors for Inference
standard.errors = mapply(function(estimate, name) {
  set.seed(12345)
  if(name == 'mc') { mc_placebo_se(setup$Y, setup$N0, setup$T0) }
  else {             sqrt(vcov(estimate, method='placebo'))     }
}, estimates, names(estimators))

head(standard.errors)

# Creating output table
california.table = rbind(unlist(estimates), unlist(standard.errors))
rownames(california.table) = c('estimate', 'standard error')
colnames(california.table) = toupper(names(estimators))
round(california.table, digits=1)

# DiD plot
synthdid_plot(estimates[1:3], facet.vertical=FALSE,
              control.name='control', treated.name='california',
              lambda.comparable=TRUE, se.method = 'none',
              trajectory.linetype = 1, line.width=.75, effect.curvature=-.4,
              trajectory.alpha=.7, effect.alpha=.7,
              diagram.alpha=1, onset.alpha=.7) +
  theme(legend.position=c(.26,.07), legend.direction='horizontal',
        legend.key=element_blank(), legend.background=element_blank(),
        strip.background=element_blank(), strip.text.x = element_blank())

# T-test for significance, H0: Effect is 0.
hypothesized_parameter = 0
t_did = (california.table[1,1] - hypothesized_parameter) / california.table[2,1]
t_sc = (california.table[1,2] - hypothesized_parameter) / california.table[2,2]
t_sdid = (california.table[1,3] - hypothesized_parameter) / california.table[2,3]
