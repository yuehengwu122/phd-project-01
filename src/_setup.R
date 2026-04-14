# ============================================================
# Project Setup: Load packages and source all modules
#
# Usage: source("src/_setup.R") at the top of any analysis script.
# ============================================================

# --- Required packages ---
library(rstan)
library(ggplot2)
library(logspline)

# --- Source project modules (order matters) ---
source("src/data_prep.R")
source("src/moment_prior.R")
source("src/model_fitting.R")
source("src/bayes_factor.R")
source("src/diagnostics.R")

# --- Global options ---
rstan::rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

theme_set(
  theme_minimal() +
    theme(text = element_text(size = 11))
)
