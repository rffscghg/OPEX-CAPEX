# Main script to run OPEX-CAPEX model

library(tidyverse)
library(hms)
library(parallel)
library(plotly)
library(RColorBrewer)
library(orca)

source("src/monte.r")
source("src/vfi.r")
source("src/utils.r")
source("src/parallelize.r")
source("src/surfaces.r")

# source("test/test_monte.r")
# source("test/test_vfi.r")
# source("test/test_utils.r")

# Parameters to adjust computational load

k_g_multiples  = seq(0.2, 4, by = .2)     # Relate state space to central value
c_f_multiples  = seq(0.2, 4, by = .2)     # Relate state space to central value
t = 10                                      # Lifespan/number of assets
n_mc = 1000                                  # Number of reps for Monte Carlo
t_mc = 100                                  # Number of timesteps for Monte Carlo
save_value_funcs = TRUE

# Parameters for plots

plot_multiples = c(0.2, 1.8)

source("src/compute.r")
source("src/figures.r")
