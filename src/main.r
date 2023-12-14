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

# Parameters to adjust computational load

k_g_multiples  = seq(0.1, 3.1, length = 11) # Relate state space to central value
c_f_multiples  = seq(0.1, 3.1, length = 11) # Relate state space to central value
t = 7                                      # Lifespan/number of assets
n_mc = 50                                 # Number of reps for Monte Carlo
t_mc = 50                                  # Number of timesteps for Monte Carlo
save_value_funcs = TRUE

# Set up parallel compute

cl <- makeCluster(5, outfile = "") # number of worker processes

clusterEvalQ(cl, library(hms))

clusterEvalQ(cl, source("src/monte.r"))
clusterEvalQ(cl, source("src/vfi.r"))
clusterEvalQ(cl, source("src/utils.r"))
clusterEvalQ(cl, source("src/parallelize.r"))

clusterExport(cl, c("k_g_multiples", "c_f_multiples", "t", "n_mc", "t_mc"))

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
        t = t, 
        verbose = TRUE
    )
})

if (save_value_funcs) {

    write_V_time <- gsub("[:]", "" , round(Sys.time()), perl=TRUE)
    save(V_funcs, file = paste0("output/V-funcs ", write_V_time,".Rdata"))

}

# Add V_id back into grid and pass V_funcs to worker processes

grid_w_id <- left_join(grid, V_func_params)

clusterExport(cl, "V_funcs")

# Run Monte Carlo simulations

mc_stats <- parApply(cl, grid_w_id, MARGIN = 1, FUN = function(x){
    monte_carlo_npv_stats(
        params = x,
        c_f_multiples = c_f_multiples,
        k_g_multiples = k_g_multiples,
        V_list = V_funcs,
        t = t, 
        n_mc = n_mc, 
        t_mc = t_mc
    )
})

# Shut down parallel compute

stopCluster(cl)

# Add mc_stats to grid

results <- bind_cols(
    grid, 
    tibble(
        SD_PV =         unlist(lapply(mc_stats, function(x) x$SD.PV)),
        SD_PV_near =    unlist(lapply(mc_stats, function(x) x$SD.PV.near)),
        E_PV =          unlist(lapply(mc_stats, function(x) x$E.PV)),
        E_PV_near =     unlist(lapply(mc_stats, function(x) x$E.PV.near))
    )
)

write_time <- gsub("[:]", "" , round(Sys.time()), perl=TRUE)
write_csv(results, paste0("output/tidy-results ", write_time,".csv"))

# Make plots

p1 <- results %>%
    ggplot(aes(x = k_g_multiples*k_g, y = c_f_multiples*c_f, fill = SD_PV/E_PV)) +
    geom_raster() +
    facet_wrap(scenario~opt_name, scales = "free")

ggsave("figures/p1.png", p1)
