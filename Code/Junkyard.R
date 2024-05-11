# Code to plot development of individual observations over time
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

# Code to plot means of control and treated
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