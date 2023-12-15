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

# Set up parallel compute

cl <- makeCluster(detectCores() - 1, outfile = "") # number of worker processes

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

plot_data <- results %>%
    filter(
        k_g_multiples >= plot_multiples[1],
        k_g_multiples <= plot_multiples[2],
        c_f_multiples >= plot_multiples[1],
        c_f_multiples <= plot_multiples[2],
    )

plot_scenarios <- plot_data %>%
    group_by(scenario, opt_name, k_g, c_f) %>%
    summarise(zmax = max(SD_PV)) %>% # Decided to set z limits manually
    mutate(
        xtitle = if_else(
            scenario == "vehicle",
            'Green CAPEX<br>     ($)',
            'Green CAPEX<br>     ($M)'
        ),
        ytitle = if_else(
            scenario == "vehicle",
            'Fossil OPEX<br> ($/year)',
            'Fossil OPEX<br> ($M/year)'
        ),
        ztitle = if_else(
            scenario == "vehicle",
            '$',
            '$M'
        ),
        ztop = if_else(
            scenario == "vehicle",
            80000,
            1200
        ),
    )

# SD_PV

for (i in 1:nrow(plot_scenarios)) {

    save_surface_plot(
        coords = results_to_SD_PV_xyz(plot_data, plot_scenarios$scenario[i], plot_scenarios$opt_name[i]),
        title = paste("Std. Dev. of NPV costs", plot_scenarios$scenario[i], plot_scenarios$opt_name[i], sep = ", "),
        scene = list(
            xaxis = list(
                title = plot_scenarios$xtitle[i],
                range = plot_scenarios$k_g[i] * plot_multiples
            ),
            yaxis = list(
                title = plot_scenarios$ytitle[i],
                range = plot_scenarios$c_f[i] * plot_multiples
            ),
            zaxis = list(
                title = plot_scenarios$ztitle[i],
                range = c(0, plot_scenarios$ztop[i])
            ),
            camera=list(eye=list(x=1.25*-1*1.5, y=1.25*-1*1.5, z=1.25*0.75*1.5))
            # default angles for x, y, and z are 1.25. Multiply by proportions to adjust
        ),
        file = paste0("figures/", plot_scenarios$scenario[i], "---", plot_scenarios$opt_name[i],".png")
    )

}

# Delta plots

# Change in SD from option, in dollars
scene_sd_delta = list(xaxis=list(title='Green CAPEX<br>     ($M)', range = c(45, 855)),
                      yaxis=list(title='Fossil OPEX<br> ($M/year)', range = c(2,48)),
                      zaxis=list(title='$M', range=c(-1250,250)),
                      camera=list(eye=list(x=1.25*-1*1.5, y=1.25*-1*1.5, z=1.25*0.75*1.5))) # default angles for x, y, and z are 1.25. Multiply by proportions to adjust

scene_ev_delta = list(xaxis=list(title='Green CAPEX<br>     ($)', range = c(2000, 78000)),
                yaxis=list(title='Fossil OPEX<br> ($/year)', range = c(100,1900)),
                zaxis=list(title='$', range=c(-1e5,2e4)),
                camera=list(eye=list(x=1.25*-1*1.5, y=1.25*-1*1.5, z=1.25*0.75*1.5))) # default angles for x, y, and z are 1.25. Multiply by proportions to adjust

delta_plot_scenarios <- expand_grid(plot_scenarios, opt_b = plot_scenarios$opt_name) %>%
    filter(
        opt_name != opt_b, 
        !str_detect(opt_name, "only"),
        !(str_detect(opt_name, "fossil")&str_detect(opt_b, "green")),
        !(str_detect(opt_name, "begin-green")&str_detect(opt_b, "fossil-only"))
    ) %>% 
    distinct_all()

for (i in 1:nrow(delta_plot_scenarios)) {
    if (delta_plot_scenarios$scenario[i] == "vehicle") {scene <- scene_ev_delta} else {scene <- scene_sd_delta}
    save_surface_plot(
        coords = a_minus_b_SD_PV_xyz(
            results, 
            delta_plot_scenarios$scenario[i], 
            delta_plot_scenarios$opt_name[i],
            delta_plot_scenarios$opt_b[i]
        ),
        title = "",
        scene = scene,
        file = paste0(
            "figures/", 
            delta_plot_scenarios$scenario[i], 
            "---delta---", 
            delta_plot_scenarios$opt_name[i],
            "---",
            delta_plot_scenarios$opt_b[i],
            ".png"
        ),
        color_scale = list(c(0, 1), c("blue", "#c0b3a2"))
    )
}
