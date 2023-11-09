# Main script to run OPEX-CAPEX model
library(tidyverse)
library(hms)

source("src/vfi.r")
source("src/utils.r")

source("test/test_vfi.r")
source("test/test_utils.r")

# Fossil exposure heatmaps

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
    filter(option == "V_min") %>%
    group_by(f_exposure, c_f, k_g) %>%
    summarise(value = mean(value)) %>%
    ggplot(aes(x = c_f, y = k_g, fill = value/max(value))) +
    geom_raster() +
    facet_wrap(~paste0("fossil-fuel exposure: ", f_exposure)) +
    scale_fill_viridis_c()

ggsave("figures/fossil_exposure.png", p1)
