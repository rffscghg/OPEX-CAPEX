# Main script to run OPEX-CAPEX model

library(tidyverse)
library(hms)

source("src/monte.r")
source("src/vfi.r")
source("src/utils.r")

source("test/test_vfi.r")
source("test/test_utils.r")

# Fossil exposure heatmaps

load("data/v_init_f_exposure_26_26_8.RData") # Loads results_1. Workaround.

results_1 <- vfi(
    c_f_vals = seq(50, 100, by = 2),
    k_g_vals = seq(500, 1000, by = 20),
    k_f = 400,
    c_g = 3,
    sigma_cf = .05,
    sigma_kg = .05,
    t = 4,
    const_scrap = TRUE,
    max_iter = 1000,
    threshold = 1e-5,
    verbose = TRUE,
    V_init = results_1
)

p1 <- tidy_V(results_1) %>%
    group_by(f_exposure, c_f, k_g) %>%
    summarise(value = mean(value)) %>%
    ggplot(aes(x = c_f, y = k_g, fill = value/max(value))) +
    geom_raster() +
    facet_wrap(~paste0("fossil-fuel exposure: ", f_exposure)) +
    scale_fill_viridis_c()

ggsave("figures/fossil_exposure.png", p1)

# Monte Carlo model

test <- monte_carlo(
    c_f_vals = seq(10, 60, by = 2),
    k_g_vals = seq(100, 600, by = 20),
    k_f = 500,
    c_g = 1,
    sigma_cf = .05,
    sigma_kg = .05,
    t = 4,
    const_scrap = TRUE,
    max_iter = 1000,
    threshold = 1e-3,
    verbose = TRUE,
    V_init = if (exists("test")) test$value_func,
    start_assets = "fgf"
)
