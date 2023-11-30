# Main script to run OPEX-CAPEX model
rm(list=ls())
setwd('/Users/Owner/Documents/GitHub/OPEX-CAPEX/')
setwd('/Users/prest/GitHub/OPEX-CAPEX/')
root = '/Users/prest/'
library(tidyverse)
library(hms)

source("src/monte.r")
# source("src/monte_skip.r")
source("src/vfi.r")
source("src/utils.r")

# source("test/test_vfi.r")
# source("test/test_utils.r")

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

# ggsave("figures/fossil_exposure.png", p1)

# Monte Carlo model ("f" begins as more attractive, "g" improves over time)
# c_f_vals = seq(5, 85, by = 10)
# k_g_vals = seq(50, 850, by = 100)
c_f_vals = seq(5, 45, by = 5)
k_g_vals = seq(50, 850, by = 50)
length(c_f_vals)
# mu_cf = .01
# mu_kg = -.01
mu_cf = 0
mu_kg = 0
k_f = 173.5 # for grid going up to 85
k_f = 296 # for grid going up to 45
c_g = 0
sigma_kg = 0.05
sigma_cf = 0.0703 # for grid going up to 85
sigma_cf = 0.1223 # for grid going up to 45
t10 = 10
t4 = 4
q = 1
r = 0.1
DF4 = 1/(1+r)^(1:t4)
DF10 = 1/(1+r)^(1:t10)

## Mean
# t=4
median(k_g_vals)*exp(mu_kg*t4) + sum(c_g*q*DF4) # green
k_f + median(c_f_vals)*exp(mu_cf*t4)*sum(exp(mu_cf*(1:t4))*q*DF4) # fossil

# t=10- aligned
median(k_g_vals)*exp(mu_kg*t10) + sum(c_g*q*DF10) # green
k_f + median(c_f_vals)*exp(mu_cf*t10)*sum(exp(mu_cf*(1:t10))*q*DF10) # fossil

## Variance
c_f_sim = random_walk_gbm(n=(2*t10)*100000, mu=mu_cf, sigma=sigma_cf, t=2*t10, x0=median(c_f_vals))

# t=4
median(k_g_vals)^2*exp(2*mu_kg*t4) * (exp(sigma_cf^2 * t4) - 1) # green
var(k_f + apply(c_f_sim, 2, function(x) sum(x[(t4+1):(2*t4)]*q*DF4)))

# t=10 - aligned
c_f_sim = random_walk_gbm(n=(2*t10)*100000, mu=mu_cf, sigma=sigma_cf, t=2*t10, x0=median(c_f_vals))
median(k_g_vals)^2*exp(2*mu_kg*t10) * (exp(sigma_kg^2 * t10) - 1) # green
var(k_f + apply(c_f_sim, 2, function(x) sum(x[(t10+1):(2*t10)]*q*DF10))) # fossil
rm(c_f_sim)

# Both portfolios
test_t4 <- monte_carlo(
  c_f_vals = c_f_vals,
  k_g_vals = k_g_vals,
  k_f = k_f,
  c_g = c_g,
  mu_cf = mu_cf,
  mu_kg = mu_kg,
  sigma_cf = sigma_cf,
  sigma_kg = sigma_kg,
  t = 4,
  t_mc = 100,
  const_scrap = TRUE,
  max_iter = 1000,
  threshold = 1e-3,
  verbose = TRUE,
  V_init = if (exists("test_t4")) test_t4$value_func,
  start_assets = "fff"
)
# test_t4_big = test_t4

test_t10 <- monte_carlo(
  c_f_vals = c_f_vals,
  k_g_vals = k_g_vals,
  k_f = k_f,
  c_g = c_g,
  mu_cf = mu_cf,
  mu_kg = mu_kg,
  sigma_cf = sigma_cf,
  sigma_kg = sigma_kg,
  t = 10,
  t_mc = 100,
  const_scrap = TRUE,
  max_iter = 1000,
  threshold = 1e-3,
  verbose = TRUE,
  V_init = if (exists("test_t10")) test_t10$value_func,
  start_assets = paste(rep("f", 9), collapse='')
)

# Only-fossil portfolios
test_t4_f <- monte_carlo(
  c_f_vals = c_f_vals,
  k_g_vals = k_g_vals,
  k_f = k_f,
  c_g = c_g,
  mu_cf = mu_cf,
  mu_kg = mu_kg,
  sigma_cf = sigma_cf,
  sigma_kg = sigma_kg,
  t = 4,
  t_mc = 100,
  option = 'f',
  const_scrap = TRUE,
  max_iter = 1000,
  threshold = 1e-3,
  verbose = TRUE,
  V_init = if (exists("test_t4_f")) test_t4_f$value_func,
  start_assets = "fff"
)

test_t10_f <- monte_carlo(
  c_f_vals = c_f_vals,
  k_g_vals = k_g_vals,
  k_f = k_f,
  c_g = c_g,
  mu_cf = mu_cf,
  mu_kg = mu_kg,
  sigma_cf = sigma_cf,
  sigma_kg = sigma_kg,
  t = 10,
  t_mc = 100,
  option = 'f',
  const_scrap = TRUE,
  max_iter = 1000,
  threshold = 1e-3,
  verbose = TRUE,
  V_init = if (exists("test_t10_f")) test_t10_f$value_func,
  start_assets = paste(rep("f", 9), collapse='')
)

# Only-green portfolios
test_t4_g <- monte_carlo(
  c_f_vals = c_f_vals,
  k_g_vals = k_g_vals,
  k_f = k_f,
  c_g = c_g,
  mu_cf = mu_cf,
  mu_kg = mu_kg,
  sigma_cf = sigma_cf,
  sigma_kg = sigma_kg,
  t = 4,
  t_mc = 100,
  option = 'g',
  const_scrap = TRUE,
  max_iter = 1000,
  threshold = 1e-3,
  verbose = TRUE,
  V_init = if (exists("test_t4_g")) test_t4_g$value_func,
  start_assets = "ggg"
)

test_t10_g <- monte_carlo(
  c_f_vals = c_f_vals,
  k_g_vals = k_g_vals,
  k_f = k_f,
  c_g = c_g,
  mu_cf = mu_cf,
  mu_kg = mu_kg,
  sigma_cf = sigma_cf,
  sigma_kg = sigma_kg,
  t = 10,
  t_mc = 100,
  option = 'g',
  const_scrap = TRUE,
  max_iter = 1000,
  threshold = 1e-3,
  verbose = TRUE,
  V_init = if (exists("test_t10_g")) test_t10_g$value_func,
  start_assets = paste(rep("g", 9), collapse='')
)

agg_cost = function(mat, horizon=NULL) {
  if (is.null(horizon)) horizon = nrow(mat)
  apply(mat, 2, function(x) sum(x[1:horizon]/(1+0.1)^c(seq_along(x[1:horizon]) - 1)))
}

## Iterate over the full grid 
E.PV.near <- SD.PV.near <- E.PV <- SD.PV <- array(NA, dim=c(length(c_f_vals), length(k_g_vals), 3, 2),
                                                  dimnames=list(c_f=c_f_vals, 
                                                                k_g=k_g_vals, 
                                                                option=c('all','g','f'),
                                                                portfolio=c('all-f', 'all-g')))

st = Sys.time()
for (i in 1:length(c_f_vals)) {          
  for (j in 1:length(k_g_vals)) {
    # for (k in 1:(2^(9-1))) {
    # for (k in string2bin('fffffffff')) { 
      for (o in c('all','g','f')) {
        for (k in 1:2) {
        if (o=='all') {
          vf_temp = test_t10$value_func
          # k = string2bin('fffffffff')
        }
        if (o=='g') {
          vf_temp = test_t10_g$value_func
          # k = string2bin('ggggggggg')
        }
        if (o=='f') {
          vf_temp = test_t10_f$value_func
          # k = string2bin('fffffffff')
        }
        mc_temp <- monte_carlo(
          c_f_vals = c_f_vals,
          k_g_vals = k_g_vals,
          k_f = k_f,
          c_g = c_g,
          sigma_cf = sigma_cf,
          sigma_kg = sigma_kg,
          t = 10,
          t_mc = 100,
          r = r,
          option = o,
          const_scrap = TRUE,
          max_iter = 1000,
          threshold = 1e-3,
          verbose = TRUE,
          V_init = if (exists("vf_temp")) vf_temp,
          start_cf = c_f_vals[i],            # Starting value for c_f
          start_kg = k_g_vals[j],            # Starting value for k_g
          start_assets = bin2string(ifelse(k==1, string2bin('fffffffff'), string2bin('ggggggggg')), t=10)
        )
        
        rel.cost.temp = mc_temp$realized_costs
        PV.Costs.near = apply(rel.cost.temp[1:10,], 2, function(x) sum(x*(1+r)^-(1:10)))
        PV.Costs = apply(rel.cost.temp, 2, function(x) sum(x*(1+r)^-(1:nrow(rel.cost.temp))))
        
        E.PV.near[i,j,o,k] = mean(PV.Costs.near)
        SD.PV.near[i,j,o,k] = sd(PV.Costs.near)
        E.PV[i,j,o,k] = mean(PV.Costs)
        SD.PV[i,j,o,k] = sd(PV.Costs)
      }
    }
    message('k_g = ', k_g_vals[j],' done at ', Sys.time())
  }
  message('c_f = ', c_f_vals[i],' done at ', Sys.time())
}
ed = Sys.time()
difftime(ed, st)
# save.image(paste0(root,'/OneDrive - rff/Documents - RPE-Electric Power/OPEX CAPEX Price Risk/output/working_opex_capex_data_',Sys.Date(),'.RData'))
load(paste0(root,'/OneDrive - rff/Documents - RPE-Electric Power/OPEX CAPEX Price Risk/output/working_opex_capex_data_2023-11-29.RData'))
# check:
hist(E.PV[,,'all','all-f']/test_t10$value_func$V_min[,,string2bin('fffffffff')])
hist(E.PV[,,'f','all-f']/test_t10$value_func$V_min[,,string2bin('fffffffff')])
hist(E.PV[,,'g','all-g']/test_t10$value_func$V_min[,,string2bin('ggggggggg')])
# can be different along the edges of the parameter space

# figures for paper:
# 1) show value functions for green only, fossil only, both, and delta (green vs fossil option value) (1x3 panel)
# 2) show monte carlo std deviation grids for the same combination. this shows option value for std dev. Another 1x3 panel. 
# This shows the effect of adding green to the choice set on the std dev (distinct from actually buying it)
# 3) Finally,  for the case with both options, show the std dev of the all-f portfolio, 
#   the all-g portfolio, and their difference. 
#   This shows the effect of actually buying the green portfolio on std dev (difference in std dev btwn all-g and all-f)
#    for both near-term change in std dev and long-term (2x3 panels)
# This yields a 4x3 panel grid.
# Do this same figure under 3 scenarios: aligned (above), and empirical analogues to ICE/EV and NG/Wind
# This will result in 3 figures, each with 4x3 panels.
#   in general, show all-f starting portfolio for fossil only, all-g for green, and all-f for both

ncols = 100
bks.V = seq(min(test_t10_g$value_func$V_min), max(test_t10_g$value_func$V_min), length.out=ncols+1)
bks.near = seq(0.5*min(SD.PV.near), max(SD.PV.near), length.out=ncols+1)
bks.far = seq(0.5*min(SD.PV), max(SD.PV), length.out=ncols+1)
bks.delta = seq(min(SD.PV[,,'all','all-f']-SD.PV[,,'all','all-g']), 
                max(SD.PV[,,'all','all-f']-SD.PV[,,'all','all-g']), length.out=ncols+1)
bks.delta = c(seq(-50, 0, length.out=ncols/2), 0, seq(0, 200, length.out=ncols/2))
# This figure:
par(mfrow=c(3,4))
image(test_t10_f$value_func$V_min[,,string2bin('fffffffff')], x=c_f_vals, y=k_g_vals, main='Fossil Only, Long-run', breaks = bks.V, col = hcl.colors(ncols, 'Reds', rev=T))
image(test_t10_g$value_func$V_min[,,string2bin('ggggggggg')], x=c_f_vals, y=k_g_vals, main='Green Only, Long-run', breaks = bks.V, col = hcl.colors(ncols, 'Reds', rev=T))
image(test_t10$value_func$V_min[,,string2bin('fffffffff')], x=c_f_vals, y=k_g_vals, main='Both Options, Long-run\nStarting with All-Fossil Portfolio', breaks = bks.V, col = hcl.colors(ncols, 'Reds', rev=T))
image(test_t10$value_func$V_min[,,string2bin('ggggggggg')], x=c_f_vals, y=k_g_vals, main='Both Options, Long-run\nStarting with All-Green Portfolio', breaks = bks.V, col = hcl.colors(ncols, 'Reds', rev=T))
image(SD.PV[,,'f','all-f'], x=c_f_vals, y=k_g_vals, main='SD, Fossil Only, Long-run\nStarting with All-Fossil Portfolio', breaks = bks.far, col = hcl.colors(ncols, 'Reds', rev=T))
image(SD.PV[,,'g','all-g'], x=c_f_vals, y=k_g_vals, main='SD, Green Only, Long-run\nStarting with All-Green Portfolio', breaks = bks.far, col = hcl.colors(ncols, 'Reds', rev=T))
image(SD.PV[,,'all','all-f'], x=c_f_vals, y=k_g_vals, main='SD, Both, Long-run\nStarting with All-Fossil Portfolio', breaks = bks.far, col = hcl.colors(ncols, 'Reds', rev=T))
image(SD.PV[,,'all','all-g'], x=c_f_vals, y=k_g_vals, main='SD, Both, Long-run\nStarting with All-Green Portfolio', breaks = bks.far, col = hcl.colors(ncols, 'Reds', rev=T))
image(SD.PV[,,'all','all-f']-SD.PV[,,'all','all-g'], x=c_f_vals, y=k_g_vals, main='SD, Both, Long-run\nGreen Procurement Benefit', breaks = bks.delta, col = hcl.colors(ncols, 'Red-Green', rev=F))
image(SD.PV.near[,,'all','all-f']-SD.PV.near[,,'all','all-g'], x=c_f_vals, y=k_g_vals, main='SD, Both, Short-run\nGreen Procurement Benefit', breaks = bks.delta, col = hcl.colors(ncols, 'Red-Green', rev=F))

library(plotly)
library(orca)
# 3D Surface version
# PV Value function
out_dir = paste0(root,'/OneDrive - rff/Documents - RPE-Electric Power/OPEX CAPEX Price Risk/output/figures/')
scene_vf = list(xaxis=list(title='Green CAPEX'),
               yaxis=list(title='Fossil OPEX'),
               zaxis=list(title='$B', range=c(0,10)),
               camera=list(eye=list(x=1.25*-1, y=1.25*-1, z=1.25*0.75))) # default angles for x, y, and z are 1.25. Multiply by proportions to adjust
layout(add_surface(plot_ly(z=test_t10_f$value_func$V_min[,,string2bin('fffffffff')]/1e3, y=c_f_vals, x=k_g_vals)),
       title = 'NPV Costs, Fossil-only', 
       scene = scene_vf,
       legend = list(title=list(text='$B')))
orca(file = paste0(out_dir,"my_plot.png")) # need to fix this error
layout(add_surface(plot_ly(z=test_t10_g$value_func$V_min[,,string2bin('ggggggggg')]/1e3, y=c_f_vals, x=k_g_vals)),
       title = 'NPV Costs, Green-only',
       scene = scene_vf,
       legend = list(title=list(text='$B')))
layout(add_surface(plot_ly(z=test_t10$value_func$V_min[,,string2bin('fffffffff')]/1e3, y=c_f_vals, x=k_g_vals)),
       title = 'NPV Costs, Both Options, starting with All-Fossil Portfolio',
       scene = scene_vf,
       legend = list(title=list(text='$B')))
layout(add_surface(plot_ly(z=test_t10$value_func$V_min[,,string2bin('ggggggggg')]/1e3, y=c_f_vals, x=k_g_vals)),
       title = 'NPV Costs, Both Options, starting with All-Green Portfolio',
       scene = scene_vf,
       legend = list(title=list(text='$B')))

# Std Dev under each scenario, long-run
scene_sd = list(xaxis=list(title='Green CAPEX'),
               yaxis=list(title='Fossil OPEX'),
               zaxis=list(title='$M', range=c(0,1000)),
               camera=list(eye=list(x=1.25*-1, y=1.25*-1, z=1.25*0.75))) # default angles for x, y, and z are 1.25. Multiply by proportions to adjust
layout(add_surface(plot_ly(z=SD.PV[,,'f','all-f'], y=c_f_vals, x=k_g_vals)),
       title = 'Std. Dev. of NPV Costs, Fossil-only',
       scene = scene_sd,
       legend = list(title=list(text='$M')))
layout(add_surface(plot_ly(z=SD.PV[,,'g','all-g'], y=c_f_vals, x=k_g_vals)),
       title = 'Std. Dev. of NPV Costs, Green-only',
       scene = scene_sd,
       legend = list(title=list(text='$M')))
layout(add_surface(plot_ly(z=SD.PV[,,'all','all-f'], y=c_f_vals, x=k_g_vals)),
       title = 'Std. Dev. of NPV Costs, Both Options, starting with All-Fossil Portfolio',
       scene = scene_sd,
       legend = list(title=list(text='$M')))
layout(add_surface(plot_ly(z=SD.PV[,,'all','all-g'], y=c_f_vals, x=k_g_vals)),
       title = 'Std. Dev. of NPV Costs, Both Options, starting with All-Green Portfolio',
       scene = scene_sd,
       legend = list(title=list(text='$M')))

# Std Dev under each scenario, short-run
layout(add_surface(plot_ly(z=SD.PV.near[,,'f','all-f'], y=c_f_vals, x=k_g_vals)),
       title = 'Std. Dev. of NPV Costs, Fossil-only\nNear Term Only',
       scene = scene_sd,
       legend = list(title=list(text='$M')))
layout(add_surface(plot_ly(z=SD.PV.near[,,'g','all-g'], y=c_f_vals, x=k_g_vals)),
       title = 'Std. Dev. of NPV Costs, Green-only\nNear Term Only',
       scene = scene_sd,
       legend = list(title=list(text='$M')))
layout(add_surface(plot_ly(z=SD.PV.near[,,'all','all-f'], y=c_f_vals, x=k_g_vals)),
       title = 'Std. Dev. of NPV Costs, Both Options, starting with All-Fossil Portfolio\nNear Term Only',
       scene = scene_sd,
       legend = list(title=list(text='$M')))
layout(add_surface(plot_ly(z=SD.PV.near[,,'all','all-g'], y=c_f_vals, x=k_g_vals)),
       title = 'Std. Dev. of NPV Costs, Both Options, starting with All-Green Portfolio\nNear Term Only',
       scene = scene_sd,
       legend = list(title=list(text='$M')))


# To do: changes (in % and $) in MC due to 1) each added option and 2) actual change in portfolio
# Here is a start:
scene_pct = list(xaxis=list(title='Green CAPEX'),
                yaxis=list(title='Fossil OPEX'),
                zaxis=list(title='%', range=c(-100,5)),
                camera=list(eye=list(x=1.25*-1, y=1.25*-1, z=1.25*0.75)))
# Option value
layout(add_surface(plot_ly(z=(SD.PV.near[,,'all','all-f']/SD.PV.near[,,'f','all-f']-1)*100, y=c_f_vals, x=k_g_vals)),
       title = 'Value of Adding Green Option\nNear Term',
       scene = scene_pct,
       legend = list(title=list(text='%')))
layout(add_surface(plot_ly(z=(SD.PV.near[,,'all','all-f']/SD.PV.near[,,'g','all-f']-1)*100, y=c_f_vals, x=k_g_vals)),
       title = 'Value of Adding Fossil Option\nNear Term',
       scene = scene_pct,
       legend = list(title=list(text='%')))
# Procurement value ($)
scene_sd2 = list(xaxis=list(title='Green CAPEX'),
                yaxis=list(title='Fossil OPEX'),
                zaxis=list(title='$M', range=c(-25,200)),
                camera=list(eye=list(x=1.25*-1, y=1.25*-1, z=1.25*0.75))) # default angles for x, y, and z are 1.25. Multiply by proportions to adjust
layout(add_surface(plot_ly(z=SD.PV.near[,,'all','all-f']-SD.PV.near[,,'all','all-g'], y=c_f_vals, x=k_g_vals)),
       title = 'Value of Procuring Green Portfolio\nNear Term',
       scene = scene_sd2,
       legend = list(title=list(text='%')))
layout(add_surface(plot_ly(z=SD.PV[,,'all','all-f']-SD.PV[,,'all','all-g'], y=c_f_vals, x=k_g_vals)),
       title = 'Value of Procuring Green Portfolio\nLong Term',
       scene = scene_sd2,
       legend = list(title=list(text='%')))
# bks.near <- bks.near <- seq(min(E.PV.near), max(E.PV), length.out=ncols+1)
par(mfrow=c(1,3))
# image(E.PV.near[,,'f'], x=c_f_vals, y=k_g_vals, main='Fossil Only, Short-run', breaks = bks.near, col = hcl.colors(ncols, 'Reds', rev=T))
# image(E.PV.near[,,'g'], x=c_f_vals, y=k_g_vals, main='Green Only, Short-run', breaks = bks.near, col = hcl.colors(ncols, 'Reds', rev=T))
# image(E.PV.near[,,'all'], x=c_f_vals, y=k_g_vals, main='Both, Short-run', breaks = bks.far, col = hcl.colors(ncols, 'Reds', rev=T))
image(E.PV[,,'f'], x=c_f_vals, y=k_g_vals, main='Fossil Only, Long-run', breaks = bks.far, col = hcl.colors(ncols, 'Reds', rev=T))
image(E.PV[,,'g'], x=c_f_vals, y=k_g_vals, main='Green Only, Long-run', breaks = bks.far, col = hcl.colors(ncols, 'Reds', rev=T))
image(E.PV[,,'all'], x=c_f_vals, y=k_g_vals, main='Both, Long-run', breaks = bks.far, col = hcl.colors(ncols, 'Reds', rev=T))


# bks.near <- bks.near <- seq(min(SD.PV.near), max(SD.PV), length.out=ncols+1)
par(mfrow=c(2,3))
image(SD.PV.near[,,'f'], x=c_f_vals, y=k_g_vals, main='Fossil Only, Short-run', breaks = bks.near, col = hcl.colors(ncols, 'Reds', rev=T))
image(SD.PV.near[,,'g'], x=c_f_vals, y=k_g_vals, main='Green Only, Short-run', breaks = bks.near, col = hcl.colors(ncols, 'Reds', rev=T))
image(SD.PV.near[,,'all'], x=c_f_vals, y=k_g_vals, main='Both, Short-run', breaks = bks.far, col = hcl.colors(ncols, 'Reds', rev=T))
image(SD.PV[,,'f'], x=c_f_vals, y=k_g_vals, main='Fossil Only, Long-run', breaks = bks.far, col = hcl.colors(ncols, 'Reds', rev=T))
image(SD.PV[,,'g'], x=c_f_vals, y=k_g_vals, main='Green Only, Long-run', breaks = bks.far, col = hcl.colors(ncols, 'Reds', rev=T))
image(SD.PV[,,'all'], x=c_f_vals, y=k_g_vals, main='Both, Long-run', breaks = bks.far, col = hcl.colors(ncols, 'Reds', rev=T))

# run_mc = function(c_f, k_g, start_asset) {
#   return(monte_carlo(
#     c_f_vals = c_f_vals,
#     k_g_vals = k_g_vals,
#     k_f = k_f,
#     c_g = c_g,
#     mu_cf = mu_cf,
#     mu_kg = mu_kg,
#     sigma_cf = sigma_cf,
#     sigma_kg = sigma_kg,
#     t = 10,
#     t_mc = 100,
#     const_scrap = TRUE,
#     max_iter = 1000,
#     threshold = 1e-3,
#     verbose = TRUE,
#     V_init = if (exists("test_t10")) test_t10$value_func,
#     start_cf = c_f,
#     start_kg = k_g,
#     start_assets = start_asset
#   ))
# }
# # debugonce(monte_carlo)
# library(parallel)
# 
# cl <- makeCluster(detectCores())
# run_mc(c_f=c_f_vals[2], k_g=k_g_vals[2], start_asset=paste(rep("f", 9), collapse=''))
# ?meshgrid
# data_list = list('c_f_vals' = c_f_vals,
#                  'k_g_vals' = k_g_vals,
#                  'start_asset'= bin2string(1:2^(10-1), t=10))
# # parLapply(cl, )
# stopCluster(cl)
# 
# median(test_t4_g$V_g[1,])
# test_t4_g$value_func$V_min[index_nearest(median(c_f_vals), c_f_vals),
#                            index_nearest(median(k_g_vals), k_g_vals),
#                           string2bin('ggg')]
(median(k_g_vals) + 4*2)*(1+.1)/.1 
mean(agg_cost(test_t4_g$realized_costs))

range(test_t4_g$realized_costs/test_t4_g$k_g)

plot(apply(test_t4_g$realized_costs, 1, mean), type='l')
lines(apply(test_t4_g$k_g + 4*2, 1, mean), type='l', col=2)
lines(apply(test_t4_g$random_kg + 4*2, 1, mean), type='l', col=3)

plot(apply(test_t4_g$realized_costs, 1, sd), type='l')
lines(apply(test_t4_g$k_g + 4*2, 1, sd), type='l', col=2)
lines(apply(test_t4_g$random_kg + 4*2, 1, sd), type='l', col=3)
# true sd understated due to gridding

# debugonce(monte_carlo)
# option value
green_t4_OV = test_t4_f$value_func$V_min - test_t4$value_func$V_min
green_t10_OV = test_t10_f$value_func$V_min - test_t10$value_func$V_min 

opt_strat = 1*(test_t10$value_func$V_f>test_t10$value_func$V_g)
range(apply(opt_strat, 3, mean))
image(opt_strat[,,1], x=c_f_vals, y=k_g_vals)

h4 = 100; h10 = 100
# h4 = 4; h10 = 10
t4_rc = agg_cost(test_t4$realized_costs, horizon = h4)
t10_rc = agg_cost(test_t10$realized_costs, horizon = h10)
t4_f_rc = agg_cost(test_t4_f$realized_costs, horizon = h4)
t10_f_rc = agg_cost(test_t10_f$realized_costs, horizon = h10)
t4_g_rc = agg_cost(test_t4_g$realized_costs, horizon = h4)
t10_g_rc = agg_cost(test_t10_g$realized_costs, horizon = h10)

# compare value function to simulated analogue (i.e., the mean)
mean(t4_rc)
test_t4$value_func$V_min[index_nearest(median(c_f_vals), c_f_vals),
                         index_nearest(median(k_g_vals), k_g_vals),
                         string2bin('fff')]

mean(t4_g_rc)
test_t4_g$value_func$V_min[index_nearest(median(c_f_vals), c_f_vals),
                         index_nearest(median(k_g_vals), k_g_vals),
                         string2bin('fff')]

# On a 4-year horizon, expected costs are smaller under fossil than green (and of course smallest with both)
mean(t4_f_rc)
mean(t4_g_rc)
mean(t4_rc)

# On a 10-year horizon, expected costs are larger under fossil than green
mean(t10_f_rc)
mean(t10_g_rc)
mean(t10_rc)

# On a 4-year horizon, std dev of costs are largest with green only
sd(t4_f_rc)
sd(t4_g_rc)
sd(t4_rc)
message('\n')
# On a 10-year horizon, std dev of costs are largest with fossil only options
sd(t10_f_rc)
sd(t10_g_rc)
sd(t10_rc)

green_t4_OV[1:5,1:5,string2bin('ggg')]
green_t4_OV[1:5,1:5,string2bin('fff')]
green_t10_OV[1:5,1:5,string2bin('ggg')]
green_t10_OV[1:5,1:5,string2bin('fff')]

# layout(add_surface(plot_ly(z=kernel2dsmooth(x=100*SD.pct.diff.near, kernel.type='disk', r=5)[-(1:3),-c(1:5, 49:51)], y=c_f_vals[-(1:3)], x=k_g_vals[-c(1:5, 49:51)])),
layout(add_surface(plot_ly(z=green_t10_OV[,,string2bin('ggg')], y=c_f_vals, x=k_g_vals)),
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

# What is the key question:
# How does moving towards a greener portfolio affect not just the mean
# but also the variance in costs?
# So I think this is about the variance in the Monte Carlo, but what compared to what?
# I think it is the difference in variance between having the option and not.
# This is not just the option value, which is based on expected values, but
# rather based on the Monte Carlo uncertainty. 
# That is, this is the Std Dev analogue of option value:
# What is bigger: the reduction in uncertainty from adding the green
# option or from adding the fossil option? And over what time period?
# Question: for what starting point? ggg or fff? The natural one is today's
# mostly fossil starting point.

p1 <- tidy_V(test_t10$value_func) %>%
  filter(f_exposure %in% c(0,15,30,45)) %>% 
  group_by(f_exposure, c_f, k_g) %>%
  summarise(value = mean(value)) %>%
  ggplot(aes(x = c_f, y = k_g, fill = value)) +
  geom_raster() +
  facet_wrap(~paste0("fossil-fuel exposure: ", f_exposure)) +
  scale_fill_viridis_c()
p1

ggsave("figures/fossil_exposure.png", p1)
