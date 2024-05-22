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
library(fixest)




# Read in dataset and get an overview
data <- read.csv2("Data/california_prop99.csv")
head(data)

# Define estimators to be used
estimators = list(did=did_estimate,
                  sc=sc_estimate,
                  sdid=synthdid_estimate)

# Convert from panel data to matrix
setup = panel.matrices(california_prop99)
head(setup)

# Compute estimates
estimates = lapply(estimators, function(estimator) { estimator(setup$Y,
                                                               setup$N0, setup$T0) } )

head(estimates)


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
plot <- synthdid_plot(estimates[1:3], facet.vertical=FALSE,
                      control.name='control', treated.name='california',
                      lambda.comparable=TRUE, se.method = 'none',
                      trajectory.linetype = 1, line.width=.75, effect.curvature=-.4,
                      trajectory.alpha=.7, effect.alpha=.7,
                      diagram.alpha=1, onset.alpha=.7) +
  theme(legend.position=c(.26,.07), legend.direction='horizontal',
        legend.text = element_text(size = 18),
        axis.title = element_text(size = 14), axis.text = element_text(size = 16),
        legend.key=element_blank(), legend.background=element_blank(),
        strip.background=element_blank(), strip.text.x = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 24)) +
  ggtitle('Proposition 99 Effects on Smoking Prevalence')

plot


ggsave("Output/Proposition99_plot.png", plot, width = 20, height=9, dpi=1000)

# T-test for significance, H0: Effect is 0.
hypothesized_parameter = 0
t_did = (california.table[1,1] - hypothesized_parameter) / california.table[2,1]
t_sc = (california.table[1,2] - hypothesized_parameter) / california.table[2,2]
t_sdid = (california.table[1,3] - hypothesized_parameter) / california.table[2,3]

ggsave("Output/Proposition99_effects.png", plot, width=8, height=6, dpi=300)
