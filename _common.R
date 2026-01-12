# load required packages
library(rstan)
library(tidyverse)
library(ggplot2)

# source core functions
# source("R/model_management.R")
source("R/data_prep.R")
source("R/diagnostics.R")
source("R/bayes_factor.R")
source("R/stan_fitting.R")

# # initialization of model registry (one-time)
# initialize_model_registry()

# uniform ggplot theme
theme_set(theme_minimal() + 
  theme(text = element_text(size = 11)))

# rstan options setting
rstan::rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())