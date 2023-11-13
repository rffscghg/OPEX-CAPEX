rm(list=ls())
# Main script to run OPEX-CAPEX model
library(tidyverse)
library(hms)

source("src/vfi.r")
source("src/utils.r")

source("test/test_vfi.r")
source("test/test_utils.r")

# Fossil exposure heatmaps

# load("data/v_init_f_exposure_26_26_8.RData") # loads results_1. Workaround.

c_f_vals = seq(1, 50, by = 2)
k_g_vals = seq(10, 600, by = 20)

results_1 <- vfi(
    c_f_vals = c_f_vals,
    k_g_vals = k_g_vals,
    k_f = 275,
    c_g = 1,
    sigma_cf = .05,
    sigma_kg = .05,
    t = 2,
    const_scrap = TRUE,
    max_iter = 1000,
    threshold = 1e-5,
    verbose = TRUE
    # V_init = results_1
)

p1 <- tidy_V(results_1) %>%
    group_by(f_exposure, c_f, k_g) %>%
    summarise(value = mean(value)) %>%
    ggplot(aes(x = c_f, y = k_g, fill = value/max(value))) +
    geom_raster() +
    facet_wrap(~paste0("fossil-fuel exposure: ", f_exposure)) +
    scale_fill_viridis_c() + 
    ggtitle('Risk Neutral, Cost Setup')

p1

# ggsave("figures/fossil_exposure.png", p1)

# Find Y such that Y*L > max[max{k_g, k_f} + L*max{c_g, c_f}]
Y = max(k_g_vals, 275)/2 + max(c_f_vals, 1)
Y

results_risk_neutral <- vfi(
  c_f_vals = c_f_vals,
  k_g_vals = k_g_vals,
  k_f = 275,
  c_g = 1,
  sigma_cf = .05,
  sigma_kg = .05,
  t = 2,
  const_scrap = TRUE,
  max_iter = 1000,
  threshold = 1e-5,
  verbose = TRUE,
  V_init = results_1,
  risk_averse = TRUE,
  eta = 0,
  output_value = Y
)

tidy_V(results_risk_neutral) %>%
  group_by(f_exposure, c_f, k_g) %>%
  summarise(value = mean(value)) %>%
  ggplot(aes(x = c_f, y = k_g, fill = value/max(value))) +
  geom_raster() +
  facet_wrap(~paste0("fossil-fuel exposure: ", f_exposure)) +
  scale_fill_viridis_c() +
  ggtitle('Risk-averse setup, eta=0')

# results_risk_neutral$V_min - results_1$V_min

results_risk_averse <- vfi(
  c_f_vals = c_f_vals,
  k_g_vals = k_g_vals,
  k_f = 275,
  c_g = 1,
  sigma_cf = .05,
  sigma_kg = .05,
  t = 2,
  const_scrap = TRUE,
  max_iter = 1000,
  threshold = 1e-5,
  V_init = results_risk_neutral,
  # V_init = results_risk_averse,
  verbose = TRUE,
  risk_averse = TRUE,
  eta = 1,
  output_value = Y
)

tidy_V(results_risk_averse) %>%
  group_by(f_exposure, c_f, k_g) %>%
  summarise(value = mean(value)) %>%
  ggplot(aes(x = c_f, y = k_g, fill = value/max(value))) +
  geom_raster() +
  facet_wrap(~paste0("fossil-fuel exposure: ", f_exposure)) +
  scale_fill_viridis_c() +
  ggtitle('Risk-averse setup, eta=1')


# tidy_V(results_risk_averse1.01) %>%
#   group_by(f_exposure, c_f, k_g) %>%
#   summarise(value = mean(value)) %>%
#   ggplot(aes(x = c_f, y = k_g, fill = value/max(value))) +
#   geom_raster() +
#   facet_wrap(~paste0("fossil-fuel exposure: ", f_exposure)) +
#   scale_fill_viridis_c() +
#   ggtitle('Risk-averse setup, eta=1.01')

results_more_risk_averse <- vfi(
  c_f_vals = c_f_vals,
  k_g_vals = k_g_vals,
  k_f = 275,
  c_g = 1,
  sigma_cf = .05,
  sigma_kg = .05,
  t = 2,
  const_scrap = TRUE,
  max_iter = 1000,
  threshold = 1e-5,
  verbose = TRUE,
  V_init = results_risk_averse,
  risk_averse = TRUE,
  eta = 2,
  output_value = Y
)

tidy_V(results_more_risk_averse) %>%
  group_by(f_exposure, c_f, k_g) %>%
  summarise(value = mean(value)) %>%
  ggplot(aes(x = c_f, y = k_g, fill = value/max(value))) +
  geom_raster() +
  facet_wrap(~paste0("fossil-fuel exposure: ", f_exposure)) +
  scale_fill_viridis_c() +
  ggtitle('Risk-averse setup, eta=2')


opt_basic = 1*(results_1$V_f > results_1$V_g) # note: lower costs, reverse interp of the others
opt_risk_neutral = 1*(results_risk_neutral$V_f < results_risk_neutral$V_g)
opt_risk_averse = 1*(results_risk_averse$V_f < results_risk_averse$V_g)
opt_more_risk_averse = 1*(results_more_risk_averse$V_f < results_more_risk_averse$V_g)

# Overall choices not strongly affected. But what about when you're heavy on one
# or the other? 'fff'
mean(opt_basic)
mean(opt_risk_neutral)
mean(opt_risk_averse)
mean(opt_more_risk_averse)

mean(opt_risk_neutral[,,string2bin('f')])
mean(opt_risk_averse[,,string2bin('f')])
mean(opt_more_risk_averse[,,string2bin('f')])

mean(opt_risk_neutral[,,string2bin('g')])
mean(opt_risk_averse[,,string2bin('g')])
mean(opt_more_risk_averse[,,string2bin('g')])


table(opt_basic - opt_risk_neutral)

# all fossil
image(z=opt_basic[,,string2bin('f')], y=k_g_vals, x=c_f_vals, main='f')
image(z=opt_risk_neutral[,,string2bin('f')], y=k_g_vals, x=c_f_vals, main='f')
image(z=opt_more_risk_averse[,,string2bin('f')], y=k_g_vals, x=c_f_vals, main='f')
# when all fossil, risk aversion pushes for more fossil

# all green
image(z=opt_basic[,,string2bin('g')], y=k_g_vals, x=c_f_vals, main='g')
image(z=opt_risk_neutral[,,string2bin('g')], y=k_g_vals, x=c_f_vals, main='g')
image(z=opt_more_risk_averse[,,string2bin('g')], y=k_g_vals, x=c_f_vals, main='g')
# when all green, risk aversion pushes for more green

# wtf?

# to do: 
# I think the odd result owes to the fact that with sigma_f = sigma_g, 
# from a LCOE standpoint the uncertainty in g is greater because
# ~all of g's costs are subject to 5% uncertainty, whereas only the 
# OPEX component of f is subject to 5% uncertainty.
# Reparameterize from original model and see if the result stands






