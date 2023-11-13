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
k_g_vals = seq(10, 800, by = 20)
t = 4
thresh = 1e-3

results_1 <- vfi(
    c_f_vals = c_f_vals,
    k_g_vals = k_g_vals,
    k_f = 275,
    c_g = 1,
    sigma_cf = .05,
    sigma_kg = .05,
    t = t,
    const_scrap = TRUE,
    max_iter = 1000,
    threshold = thresh,
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
Y = max(k_g_vals, 275)/t + max(c_f_vals, 1)
Y

results_risk_neutral <- vfi(
  c_f_vals = c_f_vals,
  k_g_vals = k_g_vals,
  k_f = 275,
  c_g = 1,
  sigma_cf = .05,
  sigma_kg = .05,
  t = t,
  const_scrap = TRUE,
  max_iter = 1000,
  threshold = thresh,
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
  t = t,
  const_scrap = TRUE,
  max_iter = 1000,
  threshold = thresh,
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
  t = t,
  const_scrap = TRUE,
  max_iter = 1000,
  threshold = thresh,
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

results_more_risk_averse$V_min[25,30,1]

opt_basic = 1*(results_1$V_f > results_1$V_g) # note: lower costs, reverse interp of the others
opt_risk_neutral = 1*(results_risk_neutral$V_f < results_risk_neutral$V_g)
opt_risk_averse = 1*(results_risk_averse$V_f < results_risk_averse$V_g)
opt_more_risk_averse = 1*(results_more_risk_averse$V_f < results_more_risk_averse$V_g)

opex_f <- c_f_vals * exp(0) * 1 / (1 + .1)
opex_g <- 1 * 1 / (1 + .1)

# Compute Brownian motion density matrix and total expenses
phimat <- phi(c_f_vals, k_g_vals, 0, 0, 0.05, 0.05, t)
sum_f_vals <- 275 + opex_f
sum_g_vals <- k_g_vals + opex_g
const_scrap = T
n_states <- ifelse(const_scrap, 2^(t-1), 1)

# debugonce(.value)
.value(opex_f, sum_f_vals, opex_g, sum_g_vals, t, n_states, r=0.1, V=results_more_risk_averse$V_min, phimat, option='all', 
       const_scrap, risk_averse=T, eta=2, output_value=Y)
# Overall choices not strongly affected. But what about when you're heavy on one
# or the other? 'fff'
mean(opt_basic)
mean(opt_risk_neutral)
mean(opt_risk_averse)
mean(opt_more_risk_averse)

mean(opt_risk_neutral[,,string2bin(paste(rep('f', t-1), collapse=''))])
mean(opt_risk_averse[,,string2bin(paste(rep('f', t-1), collapse=''))])
mean(opt_more_risk_averse[,,string2bin(paste(rep('f', t-1), collapse=''))])

mean(opt_risk_neutral[,,string2bin(paste(rep('g', t-1), collapse=''))])
mean(opt_risk_averse[,,string2bin(paste(rep('g', t-1), collapse=''))])
mean(opt_more_risk_averse[,,string2bin(paste(rep('g', t-1), collapse=''))])

table(opt_basic - opt_risk_neutral)

# Always better in g state than in f state
# Note: first figure is flipped bc it's in terms of costs, second and third are in terms of utility
hist(results_1$V_min[,,string2bin(paste(rep('f', t-1), collapse=''))]/results_1$V_min[,,string2bin(paste(rep('g', t-1), collapse=''))] - 1)
hist(results_risk_averse$V_min[,,string2bin(paste(rep('g', t-1), collapse=''))]/results_risk_averse$V_min[,,string2bin(paste(rep('f', t-1), collapse=''))] - 1)
hist(results_more_risk_averse$V_min[,,string2bin(paste(rep('g', t-1), collapse=''))]/results_more_risk_averse$V_min[,,string2bin(paste(rep('f', t-1), collapse=''))] - 1)
# note: the effects on utility are much smaller than the effects on costs
# because of diminishing marginal utility
# Further, utility effects depend on the state space, since we set
# the output value based on the most expensive point in the grid.
# This means that the larger the grid, the larger the endowment, and hence
# smaller the effects of the decision on utility


# all fossil
image(z=opt_basic[,,string2bin(paste(rep('f', t-1), collapse=''))], y=k_g_vals, x=c_f_vals, main='f')
image(z=opt_risk_neutral[,,string2bin(paste(rep('f', t-1), collapse=''))], y=k_g_vals, x=c_f_vals, main='f')
image(z=opt_more_risk_averse[,,string2bin(paste(rep('f', t-1), collapse=''))], y=k_g_vals, x=c_f_vals, main='f')
# when all fossil, risk aversion pushes for more fossil

# all green
image(z=opt_basic[,,string2bin(paste(rep('g', t-1), collapse=''))], y=k_g_vals, x=c_f_vals, main='g')
image(z=opt_risk_neutral[,,string2bin(paste(rep('g', t-1), collapse=''))], y=k_g_vals, x=c_f_vals, main='g')
image(z=opt_more_risk_averse[,,string2bin(paste(rep('g', t-1), collapse=''))], y=k_g_vals, x=c_f_vals, main='g')
# when all green, risk aversion pushes for more green

# to do: 
# 0) Think about why we're getting the result that risk aversion
#  pushes for less diversification....
#
# 1) Could the odd result owe to the fact that with sigma_f = sigma_g, 
# from a LCOE standpoint the uncertainty in g is greater because
# ~all of g's costs are subject to 5% uncertainty, whereas only the 
# OPEX component of f is subject to 5% uncertainty.
# Reparameterize from original model and see if the result stands
#
# 2) A key result is that all else equal, in NPV terms it is better to be in a
# greener portfolio (although this came as a now-sunk cost in past periods).
# Show how the time horizon affects this result.
#
# 3) is it an interesting result that risk aversion doesn't change optimal 
# strategy very much? Why doesn't risk aversion push for more diversification?





