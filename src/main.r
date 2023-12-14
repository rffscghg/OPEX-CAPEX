# Main script to run OPEX-CAPEX model

library(tidyverse)
library(hms)
library(parallel)
library(plotly)
library(orca)
library(RColorBrewer)

source("src/monte.r")
source("src/vfi.r")
source("src/utils.r")
source("src/parallelize.r")

# source("test/test_monte.r")
# source("test/test_vfi.r")
# source("test/test_utils.r")

# Parameters to adjust computational intensity

k_g_multiples  = seq(.1, 3.1, length = 4) # Relate state space to central value
c_f_multiples  = seq(.1, 3.1, length = 4) # Relate state space to central value
t = 7                                     # Lifespan/number of assets

# Set up parallel compute

cl <- makeCluster(3, outfile = "") # number of worker processes

clusterEvalQ(cl, library(hms))

clusterEvalQ(cl, source("src/monte.r"))
clusterEvalQ(cl, source("src/vfi.r"))
clusterEvalQ(cl, source("src/utils.r"))
clusterEvalQ(cl, source("src/parallelize.r"))

clusterExport(cl, c("k_g_multiples", "c_f_multiples", "t"))

# Create grid of parameterizations

scenarios <- read_csv("data/scenarios.csv", col_select = -1) # See doc/README

options = tibble(
        opt_name = c("fossil-only", "green-only", "both-begin-fossil", "both-begin-green"),
        option = c("f", "g", "all", "all"),
        start_assets = c("f", "g", "f", "g")
    ) %>% 
    mutate(start_assets = strrep(start_assets, t - 1))

grid <- expand_grid(
    scenario = c("neutral", "power-plant", "vehicle"),
    opt_name = c("fossil-only", "green-only", "both-begin-fossil", "both-begin-green")
    ) %>%
    left_join(scenarios) %>%
    left_join(options) %>%
    expand_grid(k_g_multiples, c_f_multiples)

# Solve VFI for each scenario

V_func_params <- distinct(grid, across(k_g:option)) %>%
    mutate(V_id = row_number()) # So we can key the Monte Carlo to the right V

V_funcs <- parApply(cl, V_func_params, MARGIN = 1, FUN = function(x){
    parallel_vfi(
        params = x,
        c_f_multiples = c_f_multiples,
        k_g_multiples = k_g_multiples,
        t = t
    )
})

# Add V_id back into grid and pass V_funcs to worker processes

grid_w_id <- left_join(grid, V_func_params)

clusterExport(cl, "V_funcs")

# Run Monte Carlo simulations

mc_stats <- parApply(cl, grid_w_id, MARGIN = 1, FUN = function(x){
    monte_carlo_npv_stats(
        params = x,
        V_list = V_funcs,
        t = t
    )
})

# Shut down parallel compute

stopCluster(cl)
