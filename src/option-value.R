# Main script to run OPEX-CAPEX model

library(tidyverse)
library(hms)

source("src/monte.r")
source("src/monte_skip.r")
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

# Monte Carlo model ("f" begins as more attractive, "g" improves over time)
c_f_vals = seq(4, 80, by = 10)
k_g_vals = seq(40, 800, by = 100)
length(c_f_vals)
mu_cf = .01
mu_kg = -.01

# Both portfolios
test_t4 <- monte_carlo(
  c_f_vals = c_f_vals,
  k_g_vals = k_g_vals,
  k_f = 250,
  c_g = 2,
  mu_cf = mu_cf,
  mu_kg = mu_kg,
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
testV_t4 = test_t4$value_func

test_t10 <- monte_carlo(
  c_f_vals = c_f_vals,
  k_g_vals = k_g_vals,
  k_f = 250,
  c_g = 2,
  mu_cf = mu_cf,
  mu_kg = mu_kg,
  sigma_cf = .05,
  sigma_kg = .05,
  t = 10,
  const_scrap = TRUE,
  max_iter = 1000,
  threshold = 1e-3,
  verbose = TRUE,
  V_init = if (exists("testV_t10_f")) testV_t10_f$value_func,
  start_assets = "fgf"
)
testV_t10 = test_t10$value_func

# Only-fossil portfolios
testV_t4_f <- vfi(
  c_f_vals = c_f_vals,
  k_g_vals = k_g_vals,
  k_f = 250,
  c_g = 2,
  mu_cf = mu_cf,
  mu_kg = mu_kg,
  sigma_cf = .05,
  sigma_kg = .05,
  t = 4,
  option='f',
  const_scrap = TRUE,
  max_iter = 1000,
  threshold = 1e-3,
  verbose = TRUE,
  V_init = if (exists("test_t4")) test_t4$value_func
  # start_assets = "fgf"
)
str(testV_f)

testV_t10_f <- vfi(
  c_f_vals = c_f_vals,
  k_g_vals = k_g_vals,
  k_f = 250,
  c_g = 2,
  mu_cf = mu_cf,
  mu_kg = mu_kg,
  sigma_cf = .05,
  sigma_kg = .05,
  t = 10,
  option='f',
  const_scrap = TRUE,
  max_iter = 1000,
  threshold = 1e-3,
  verbose = TRUE,
  V_init = if (exists("testV_t10_f")) testV_t10_f$value_func
  # start_assets = "fgf"
)

green_t4_OV = testV_t4_f$V_min - testV_t4$V_min
green_t10_OV = testV_t10_f$V_min - testV_t10$V_min 

library(plotly)

green_t4_OV[1:5,1:5,string2bin('ggg')]
green_t4_OV[1:5,1:5,string2bin('fff')]

# layout(add_surface(plot_ly(z=kernel2dsmooth(x=100*SD.pct.diff.near, kernel.type='disk', r=5)[-(1:3),-c(1:5, 49:51)], y=c_f_vals[-(1:3)], x=k_g_vals[-c(1:5, 49:51)])),
layout(add_surface(plot_ly(z=green_t4_OV[,,string2bin('ggg')], y=c_f_vals, x=k_g_vals)),
       title = 'Green Option Value, L=10',
       legend = list(title=list(text='$M Difference')),
       scene = list(xaxis=list(title='Green CAPEX'),
                    yaxis=list(title='Fossil OPEX'),
                    zaxis=list(title='$M')))

layout(add_surface(plot_ly(z=green_t10_OV[,,string2bin('fff')]-green_t4_OV[,,string2bin('fff')] , y=c_f_vals, x=k_g_vals)),
       title = 'Difference in green Option Value, L=10 vs. L=4',
       legend = list(title=list(text='$ Difference')),
       scene = list(xaxis=list(title='Green CAPEX'),
                    yaxis=list(title='Fossil OPEX'),
                    zaxis=list(title='$')))


