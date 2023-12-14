# Main script to run OPEX-CAPEX model
rm(list=ls())
palette(c(rgb(4, 39, 60, maxColorValue=255), # RFF black
          rgb(116, 100, 94, maxColorValue=255), # RFF brown
          rgb(80, 177, 97, maxColorValue=255), # RFF green
          rgb(136, 196, 244, maxColorValue=255), # RFF blue
          rgb(255, 102, 99, maxColorValue=255), # RFF coral
          rgb(118, 94, 165, maxColorValue=255), # RFF purple
          rgb(244, 162, 95, maxColorValue=255), # RFF orange
          rgb(235, 211, 103, maxColorValue=255), # RFF yellow
          rgb(224, 228, 231, maxColorValue=255))) # RFF light gray

script.start = Sys.time()
setwd('/Users/Owner/Documents/GitHub/OPEX-CAPEX/')
root = '/Users/Owner/'
# setwd('/Users/prest/GitHub/OPEX-CAPEX/')
# root = '/Users/prest/'
library(tidyverse)
library(hms)

source("src/monte.r")
# source("src/monte_skip.r")
source("src/vfi.r")
source("src/utils.r")

# source("test/test_vfi.r")
# source("test/test_utils.r")

# Fossil exposure heatmaps

# load("data/v_init_f_exposure_26_26_8.RData") # Loads results_1. Workaround.
# 
# results_1 <- vfi(
#   c_f_vals = seq(50, 100, by = 2),
#   k_g_vals = seq(500, 1000, by = 20),
#   k_f = 400,
#   c_g = 3,
#   sigma_cf = .05,
#   sigma_kg = .05,
#   t = 4,
#   const_scrap = TRUE,
#   max_iter = 1000,
#   threshold = 1e-3,
#   verbose = TRUE,
#   V_init = results_1
# )
# 
# p1 <- tidy_V(results_1) %>%
#   group_by(f_exposure, c_f, k_g) %>%
#   summarise(value = mean(value)) %>%
#   ggplot(aes(x = c_f, y = k_g, fill = value/max(value))) +
#   geom_raster() +
#   facet_wrap(~paste0("fossil-fuel exposure: ", f_exposure)) +
#   scale_fill_viridis_c()
#
# ggsave("figures/fossil_exposure.png", p1)

# Monte Carlo model ("f" begins as more attractive, "g" improves over time)
# c_f_vals = seq(5, 85, by = 10)
# k_g_vals = seq(50, 850, by = 100)
# c_f_vals = seq(5, 45, by = 5)
# k_g_vals = seq(50, 850, by = 50)
# Regular grid
{
  c_f_vals = seq(5, 45, by = 2.5)
  k_g_vals = seq(50, 850, by = 25)
  c_f_vals_core = c_f_vals
  k_g_vals_core = k_g_vals
}
# Fuller grid
{
  c_f_vals = seq(5, 45+20, by = 2.5)
  k_g_vals = seq(50, 850+200, by = 25)
  idx_cf = (1:length(c_f_vals))[-which(c_f_vals>45)]
  idx_kg = (1:length(k_g_vals))[-which(k_g_vals>850)]
  c_f_vals_core = c_f_vals[idx_cf]
  k_g_vals_core = k_g_vals[idx_kg]
}
length(c_f_vals)
length(k_g_vals)
median(c_f_vals_core) # $25
median(k_g_vals_core) # $450
# mu_cf = .01
# mu_kg = -.01
mu_cf = 0
mu_kg = 0
# k_f = 173.5 # for grid going up to 85
# k_f = 296 # for grid going up to 45
k_f = median(k_g_vals_core)
c_g = median(c_f_vals_core)
# c_g = 0
sigma_kg = 0.05
# sigma_cf = 0.0703 # for grid going up to 85
sigma_cf = 0.1223 # for grid going up to 45
t10 = 10
t4 = 4
q = 1
r = 0.1
DF4 = 1/(1+r)^(1:t4)
DF10 = 1/(1+r)^(1:t10)

## Mean
# t=4
median(k_g_vals_core)*exp(mu_kg*t4) + sum(c_g*q*DF4) # green
k_f + median(c_f_vals_core)*exp(mu_cf*t4)*sum(exp(mu_cf*(1:t4))*q*DF4) # fossil

# t=10- aligned
median(k_g_vals_core)*exp(mu_kg*t10) + sum(c_g*q*DF10) # green
k_f + median(c_f_vals_core)*exp(mu_cf*t10)*sum(exp(mu_cf*(1:t10))*q*DF10) # fossil

## Variance
c_f_sim = random_walk_gbm(n=(2*t10)*100000, mu=mu_cf, sigma=sigma_cf, t=2*t10, x0=median(c_f_vals_core))

# t=4
median(k_g_vals_core)^2*exp(2*mu_kg*t4) * (exp(sigma_cf^2 * t4) - 1) # green
var(k_f + apply(c_f_sim, 2, function(x) sum(x[(t4+1):(2*t4)]*q*DF4)))

# t=10 - aligned
# c_f_sim = random_walk_gbm(n=(2*t10)*100000, mu=mu_cf, sigma=sigma_cf, t=2*t10, x0=median(c_f_vals_core))
median(k_g_vals_core)^2*exp(2*mu_kg*t10) * (exp(sigma_kg^2 * t10) - 1) # green
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
# to do: fix the monte carlo so it can be run with only f as an option. what happened??
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
          skipVFI = TRUE,
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
script.end = Sys.time()
difftime(script.end, script.start) 
# 1.7 hours with state space of 17*33*512 = 287k
# 3.8 hours with state space of 25*41*512 = 525k
dim(test_t10$value_func$V_min)
prod(dim(test_t10$value_func$V_min))
# save.image(paste0(root,'/OneDrive - rff/Documents - RPE-Electric Power/OPEX CAPEX Price Risk/output/working_opex_capex_data_matched_costs_',Sys.Date(),'.RData'))
load(paste0(root,'/OneDrive - rff/Documents - RPE-Electric Power/OPEX CAPEX Price Risk/output/working_opex_capex_data_matched_costs_2023-12-05.RData'))

# check:
hist(E.PV[idx_cf, idx_kg,'all','all-f']/test_t10$value_func$V_min[idx_cf, idx_kg,string2bin('fffffffff')])
hist(E.PV[idx_cf, idx_kg,'f','all-f']/test_t10$value_func$V_min[idx_cf, idx_kg,string2bin('fffffffff')])
hist(E.PV[idx_cf, idx_kg,'g','all-g']/test_t10$value_func$V_min[idx_cf, idx_kg,string2bin('ggggggggg')])
# can be different along the edges of the parameter space
opt_strat = 1*(test_t10$value_func$V_f[idx_cf, idx_kg,]<test_t10$value_func$V_g[idx_cf, idx_kg,])
apply(opt_strat, 3, mean)
table(opt_strat[idx_cf, idx_kg,1]==opt_strat[idx_cf, idx_kg,512])
image(opt_strat[idx_cf,idx_kg,1], x=c_f_vals_core, y=k_g_vals_core)
points(x=c_g, y=k_f, pch=19, cex=2, col=4)

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
bks.delta = seq(min(SD.PV[idx_cf, idx_kg,'all','all-f']-SD.PV[idx_cf, idx_kg,'all','all-g']), 
                max(SD.PV[idx_cf, idx_kg,'all','all-f']-SD.PV[idx_cf, idx_kg,'all','all-g']), length.out=ncols+1)
bks.delta = c(seq(-10, 0, length.out=ncols/2), 0, seq(0, 250, length.out=ncols/2))
# This figure:
par(mfrow=c(3,4))
image(test_t10_f$value_func$V_min[idx_cf, idx_kg,string2bin('fffffffff')], x=c_f_vals_core, y=k_g_vals_core, main='Fossil Only, Long-run', breaks = bks.V, col = hcl.colors(ncols, 'Reds', rev=T))
image(test_t10_g$value_func$V_min[idx_cf, idx_kg,string2bin('ggggggggg')], x=c_f_vals_core, y=k_g_vals_core, main='Green Only, Long-run', breaks = bks.V, col = hcl.colors(ncols, 'Reds', rev=T))
image(test_t10$value_func$V_min[idx_cf, idx_kg,string2bin('fffffffff')], x=c_f_vals_core, y=k_g_vals_core, main='Both Options, Long-run\nStarting with All-Fossil Portfolio', breaks = bks.V, col = hcl.colors(ncols, 'Reds', rev=T))
image(test_t10$value_func$V_min[idx_cf, idx_kg,string2bin('ggggggggg')], x=c_f_vals_core, y=k_g_vals_core, main='Both Options, Long-run\nStarting with All-Green Portfolio', breaks = bks.V, col = hcl.colors(ncols, 'Reds', rev=T))
image(SD.PV[idx_cf, idx_kg,'f','all-f'], x=c_f_vals_core, y=k_g_vals_core, main='SD, Fossil Only, Long-run\nStarting with All-Fossil Portfolio', breaks = bks.far, col = hcl.colors(ncols, 'Reds', rev=T))
image(SD.PV[idx_cf, idx_kg,'g','all-g'], x=c_f_vals_core, y=k_g_vals_core, main='SD, Green Only, Long-run\nStarting with All-Green Portfolio', breaks = bks.far, col = hcl.colors(ncols, 'Reds', rev=T))
image(SD.PV[idx_cf, idx_kg,'all','all-f'], x=c_f_vals_core, y=k_g_vals_core, main='SD, Both, Long-run\nStarting with All-Fossil Portfolio', breaks = bks.far, col = hcl.colors(ncols, 'Reds', rev=T))
image(SD.PV[idx_cf, idx_kg,'all','all-g'], x=c_f_vals_core, y=k_g_vals_core, main='SD, Both, Long-run\nStarting with All-Green Portfolio', breaks = bks.far, col = hcl.colors(ncols, 'Reds', rev=T))
image(SD.PV[idx_cf, idx_kg,'all','all-f']-SD.PV[idx_cf, idx_kg,'all','all-g'], x=c_f_vals_core, y=k_g_vals_core, main='SD, Both, Long-run\nGreen Procurement Benefit', breaks = bks.delta, col = hcl.colors(ncols, 'Red-Green', rev=F))
image(SD.PV.near[idx_cf, idx_kg,'all','all-f']-SD.PV.near[idx_cf, idx_kg,'all','all-g'], x=c_f_vals_core, y=k_g_vals_core, main='SD, Both, Short-run\nGreen Procurement Benefit', breaks = bks.delta, col = hcl.colors(ncols, 'Red-Green', rev=F))
image(opt_strat[idx_cf, idx_kg,1], x=c_f_vals_core, y=k_g_vals_core, main='Optimal Strategy', col = hcl.colors(ncols, 'Reds', rev=T))

# At point where costs are aligned, which has greater option value? procurement value?
# Check approximate alignment:
test_t10_f$value_func$V_min[median(idx_cf), median(idx_kg),c(1,512)]
test_t10_g$value_func$V_min[median(idx_cf), median(idx_kg),c(1,512)]
test_t10$value_func$V_min[median(idx_cf), median(idx_kg),c(1,512)]

E.PV[median(idx_cf), median(idx_kg),,]
SD.PV[median(idx_cf), median(idx_kg),,]

green_ov = test_t10$value_func$V_min - test_t10_f$value_func$V_min
fossil_ov = test_t10$value_func$V_min - test_t10_g$value_func$V_min

green_ov[median(idx_cf), median(idx_kg),]
fossil_ov[median(idx_cf), median(idx_kg),]

E.PV.highn <- SD.PV.highn <-  array(NA, dim=c(3, 2),
                                    dimnames=list(option=c('all','g','f'),
                                                  portfolio=c('all-f', 'all-g')))
for (o in c('all','g','f')) {
  for (k in 1:2) {
  if (o=='all') {
    vf_temp = test_t10$value_func
  }
  if (o=='g') {
    vf_temp = test_t10_g$value_func
  }
  if (o=='f') {
    vf_temp = test_t10_f$value_func
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
    n_mc=1e4,
    r = r,
    option = o,
    const_scrap = TRUE,
    max_iter = 1000,
    threshold = 1e-3,
    verbose = TRUE,
    V_init = if (exists("vf_temp")) vf_temp,
    skipVFI = TRUE,
    start_cf = median(c_f_vals_core),            # Starting value for c_f
    start_kg = median(k_g_vals_core),            # Starting value for k_g
    start_assets = bin2string(ifelse(k==1, string2bin('fffffffff'), string2bin('ggggggggg')), t=10)
  )
  rel.cost.temp = mc_temp$realized_costs
  PV.Costs = apply(rel.cost.temp, 2, function(x) sum(x*(1+r)^-(1:nrow(rel.cost.temp))))
  
  E.PV.highn[o,k] = mean(PV.Costs)
  SD.PV.highn[o,k] = sd(PV.Costs)
  }
}
SD.PV[median(idx_cf), median(idx_kg),,]
SD.PV.highn

SD.PV[median(idx_cf), median(idx_kg),'all',] - SD.PV[median(idx_cf), median(idx_kg),'f',]
SD.PV.highn['all',]-SD.PV.highn['f',] # fossil SD option value
SD.PV[median(idx_cf), median(idx_kg),'all',] - SD.PV[median(idx_cf), median(idx_kg),'g',]
SD.PV.highn['all',]-SD.PV.highn['g',] # fossil SD option value



green_ov_sd = SD.PV[c(1, median(idx_cf), max(idx_cf)), c(1, median(idx_kg), max(idx_kg)), 'all', 'all-f'] - 
  SD.PV[c(1, median(idx_cf), max(idx_cf)), c(1, median(idx_kg), max(idx_kg)), 'f', 'all-f']
# fossil_ov_sd = SD.PV[c(1, median(idx_cf), max(idx_cf)), c(1, median(idx_kg), max(idx_kg)), 'all', 'all-f'] - 
#   SD.PV[c(1, median(idx_cf), max(idx_cf)), c(1, median(idx_kg), max(idx_kg)), 'g', 'all-f']
fossil_ov_sd = SD.PV[c(1, median(idx_cf), max(idx_cf)), c(1, median(idx_kg), max(idx_kg)), 'all', 'all-g'] - 
  SD.PV[c(1, median(idx_cf), max(idx_cf)), c(1, median(idx_kg), max(idx_kg)), 'g', 'all-g']

green_ov_sd
fossil_ov_sd
fossil_ov_sd<green_ov_sd
par(mfrow=c(1,2))
barplot(c(fossil_ov_sd[2,2], green_ov_sd[2,2]), names.arg=c('Fossil OV','Green OV'))
barplot(c(SD.PV.highn['all','all-g'] - SD.PV.highn['g','all-g'],
          SD.PV.highn['all','all-f'] - SD.PV.highn['f','all-f']), names.arg=c('Fossil OV','Green OV'))
# green_ov_sd_near = SD.PV.near[median(idx_cf), median(idx_kg), 'all', ] - SD.PV.near[median(idx_cf), median(idx_kg), 'f', ]
# fossil_ov_sd_near = SD.PV.near[median(idx_cf), median(idx_kg), 'all', ] - SD.PV.near[median(idx_cf), median(idx_kg), 'g', ]

library(plotly)
library(orca)
library(RColorBrewer)

# display.brewer.all()
# colscale = brewer.pal(3, 'YlOrRd')
colscale = 'RdOrYl'
# colscale = 'Plasma'
# colscale = NULL

# 3D Surface version
# PV Value function
out_dir = paste0(root,'/OneDrive - rff/Documents - RPE-Electric Power/OPEX CAPEX Price Risk/output/figures/matched_costs/')
# dir.create(out_dir)
setwd(out_dir)
# out_dir=='C:/Users/prest//OneDrive - rff/Documents - RPE-Electric Power/OPEX CAPEX Price Risk/output/figures/'
scene_vf = list(xaxis=list(title='Green CAPEX<br>     ($M)'),
                yaxis=list(title='Fossil OPEX<br> ($M/year)'),
                zaxis=list(title='$B', range=c(0,10)),
                camera=list(eye=list(x=1.25*-1*1.5, y=1.25*-1*1.5, z=1.25*0.75*1.5))) # default angles for x, y, and z are 1.25. Multiply by proportions to adjust
p1 = layout(add_surface(plot_ly(z=test_t10_f$value_func$V_min[idx_cf, idx_kg,string2bin('fffffffff')]/1e3, y=c_f_vals_core, x=k_g_vals_core), colorscale = colscale),
       title = list(text='NPV Costs, Fossil-only', y = 0.9),
       margin = list(t = 0),
       scene = scene_vf,
       legend = list(title='$B', font=list(size=12)))
orca(p1, file = "value_function_fossil_only.png", scale=0.75, width=800*0.8, height=800*0.8) 
p2 = layout(add_surface(plot_ly(z=test_t10_g$value_func$V_min[idx_cf, idx_kg,string2bin('ggggggggg')]/1e3, y=c_f_vals_core, x=k_g_vals_core), colorscale = colscale),
       title = list(text='NPV Costs, Green-only', y = 0.9),
       scene = scene_vf,
       legend = list(title=list(text='$B')))
orca(p2, file = "value_function_green_only.png", scale=0.75, width=800*0.8, height=800*0.8) 
p3 = layout(add_surface(plot_ly(z=test_t10$value_func$V_min[idx_cf, idx_kg,string2bin('fffffffff')]/1e3, y=c_f_vals_core, x=k_g_vals_core), colorscale = colscale),
       title = list(text='NPV Costs, Both Options<br>Starting with All-Fossil Portfolio', y = 0.9),
       scene = scene_vf,
       legend = list(title=list(text='$B')))
orca(p3, file = "value_function_both_fossil_portfolio.png", scale=0.75, width=800*0.8, height=800*0.8) 
p4 = layout(add_surface(plot_ly(z=test_t10$value_func$V_min[idx_cf, idx_kg,string2bin('ggggggggg')]/1e3, y=c_f_vals_core, x=k_g_vals_core), colorscale = colscale),
       title = list(text='NPV Costs, Both Options<br>Starting with All-Green Portfolio', y = 0.9),
       scene = scene_vf,
       legend = list(title=list(text='$B')))
orca(p4, file = "value_function_both_green_portfolio.png", scale=0.75, width=800*0.8, height=800*0.8) 

# Std Dev under each scenario, long-run
library(smoothie)
library(ks)
# Function to apply moving average smoothing to a matrix
smooth_matrix <- function(mat, window_size = 3) {
  n <- nrow(mat)
  m <- ncol(mat)
  smoothed_mat <- matrix(NA, nrow = n, ncol = m)
  
  for (i in 1:n) {
    for (j in 1:m) {
      # Extract the window around each element (handling edge cases)
      row_start <- max(1, i - floor(window_size / 2))
      row_end <- min(n, i + floor(window_size / 2))
      col_start <- max(1, j - floor(window_size / 2))
      col_end <- min(m, j + floor(window_size / 2))
      
      # Calculate the average of the window
      smoothed_mat[i, j] <- mean(mat[row_start:row_end, col_start:col_end])
    }
  }
  return(smoothed_mat)
}

scene_sd = list(xaxis=list(title='Green CAPEX<br>     ($M)'),
                yaxis=list(title='Fossil OPEX<br> ($M/year)'),
                zaxis=list(title='$M', range=c(0,1000)),
                camera=list(eye=list(x=1.25*-1*1.5, y=1.25*-1*1.5, z=1.25*0.75*1.5))) # default angles for x, y, and z are 1.25. Multiply by proportions to adjust
p5 = layout(add_surface(plot_ly(z=smooth_matrix(SD.PV[,,'f','all-f']),
                   y=c_f_vals_core, x=k_g_vals_core), colorscale = colscale),
       title = list(text='Std. Dev. of NPV Costs, Fossil-only', y = 0.9),
       scene = scene_sd,
       legend = list(title=list(text='$M')))
orca(p5, file = "std_dev_fossil_only.png", scale=0.75, width=800*0.8, height=800*0.8)
p6 = layout(add_surface(plot_ly(z=smooth_matrix(SD.PV[, ,'g','all-g'])[idx_cf, idx_kg], 
                           y=c_f_vals_core, x=k_g_vals_core), colorscale = colscale),
       title = list(text='Std. Dev. of NPV Costs, Green-only', y = 0.9),
       scene = scene_sd,
       legend = list(title=list(text='$M')))
orca(p6, file = "std_dev_green_only.png", scale=0.75, width=800*0.8, height=800*0.8)
p7 = layout(add_surface(plot_ly(z=smooth_matrix(SD.PV[, ,'all','all-f'])[idx_cf, idx_kg], 
                                y=c_f_vals_core, x=k_g_vals_core), colorscale = colscale),
       title = list(text='Std. Dev. of NPV Costs, Both Options<br>Starting with All-Fossil Portfolio', y = 0.9),
       scene = scene_sd,
       legend = list(title=list(text='$M')))
orca(p7, file = "std_dev_both_fossil_portfolio.png", scale=0.75, width=800*0.8, height=800*0.8)
p8 = layout(add_surface(plot_ly(z=smooth_matrix(SD.PV[, ,'all','all-g'])[idx_cf, idx_kg], 
                                y=c_f_vals_core, x=k_g_vals_core), colorscale = colscale),
       title = list(text='Std. Dev. of NPV Costs, Both Options<br>Starting with All-Green Portfolio', y = 0.9),
       scene = scene_sd,
       legend = list(title=list(text='$M')))
orca(p8, file = "std_dev_both_green_portfolio.png", scale=0.75, width=800*0.8, height=800*0.8)

# Std Dev under each scenario, short-run
p9 = layout(add_surface(plot_ly(z=smooth_matrix(SD.PV.near[,,'f','all-f']),
                                y=c_f_vals_core, x=k_g_vals_core), colorscale = colscale),
            title = list(text='Std. Dev. of NPV Costs, Fossil-only<br>Near-term Only', y = 0.9),
            scene = scene_sd,
            legend = list(title=list(text='$M')))
orca(p9, file = "std_dev_fossil_only_nearterm.png", scale=0.75, width=800*0.8, height=800*0.8)
p10 = layout(add_surface(plot_ly(z=smooth_matrix(SD.PV.near[, ,'g','all-g'])[idx_cf, idx_kg], 
                                y=c_f_vals_core, x=k_g_vals_core), colorscale = colscale),
            title = list(text='Std. Dev. of NPV Costs, Green-only<br>Near-term Only', y = 0.9),
            scene = scene_sd,
            legend = list(title=list(text='$M')))
orca(p10, file = "std_dev_green_only_nearterm.png", scale=0.75, width=800*0.8, height=800*0.8)
p11 = layout(add_surface(plot_ly(z=smooth_matrix(SD.PV.near[, ,'all','all-f'])[idx_cf, idx_kg], 
                                y=c_f_vals_core, x=k_g_vals_core), colorscale = colscale),
            title = list(text='Std. Dev. of NPV Costs, Both Options<br>Starting with All-Fossil Portfolio<br>Near-term Only', y = 0.9),
            scene = scene_sd,
            legend = list(title=list(text='$M')))
orca(p11, file = "std_dev_both_fossil_portfolio_nearterm.png", scale=0.75, width=800*0.8, height=800*0.8)
p12 = layout(add_surface(plot_ly(z=smooth_matrix(SD.PV.near[, ,'all','all-g'])[idx_cf, idx_kg], 
                                y=c_f_vals_core, x=k_g_vals_core), colorscale = colscale),
            title = list(text='Std. Dev. of NPV Costs, Both Options<br>Starting with All-Green Portfolio<br>Near-term Only', y = 0.9),
            scene = scene_sd,
            legend = list(title=list(text='$M')))
orca(p12, file = "std_dev_both_green_portfolio_nearterm.png", scale=0.75, width=800*0.8, height=800*0.8)


# Changes in std dev from adding 1) green option, 2) green portfolio
scene_pct = list(xaxis=list(title='Green CAPEX<br>     ($M)'),
                yaxis=list(title='Fossil OPEX<br> ($M/year)'),
                zaxis=list(title='%', range=c(-100,100)),
                camera=list(eye=list(x=1.25*-1*1.5, y=1.25*-1*1.5, z=1.25*0.75*1.5))) # default angles for x, y, and z are 1.25. Multiply by proportions to adjust
# to do here: 
# 1) fix the color scale so that reductions are always blue (or green), and not red
# 2) deal with the outlier in p13 at the lowest value of c_g and low-ish values of k_g
# 3) do the same thing for procurement value (SD.PV[,,'all','all-g'] - SD.PV[,,'all','all-f']) in both % and $

# Chang in SD from option, in percent (%)
green_option_pct_reduc = (SD.PV[,,'all','all-f']/SD.PV[,,'f','all-f']-1)*100
p13 = layout(add_surface(plot_ly(z=smooth_matrix(green_option_pct_reduc, window_size = 6)[idx_cf, idx_kg],
                                y=c_f_vals_core, x=k_g_vals_core), colorscale = colscale),
            title = list(text='Percent Change in Std. Dev. of NPV Costs<br>from Adding Green Option', y = 0.9),
            scene = scene_pct,
            legend = list(title=list(text='%')))
orca(p13, file = "pct_change_std_dev_green_option.png", scale=0.75, width=800*0.8, height=800*0.8)

fossil_option_pct_reduc = (SD.PV[,,'all','all-g']/SD.PV[,,'g','all-g']-1)*100
p14 = layout(add_surface(plot_ly(z=smooth_matrix(fossil_option_pct_reduc, window_size = 6)[idx_cf, idx_kg],
                                 y=c_f_vals_core, x=k_g_vals_core), colorscale = colscale),
             title = list(text='Percent Change in Std. Dev. of NPV Costs<br>from Adding Fossil Option', y = 0.9),
             scene = scene_pct,
             legend = list(title=list(text='%')))
orca(p14, file = "pct_change_std_dev_fossil_option.png", scale=0.75, width=800*0.8, height=800*0.8)

image(green_option_pct_reduc<fossil_option_pct_reduc)
image(green_option_dol_reduc<fossil_option_dol_reduc)

# Change in SD from option, in dollars
scene_sd_delta = list(xaxis=list(title='Green CAPEX<br>     ($M)'),
                      yaxis=list(title='Fossil OPEX<br> ($M/year)'),
                      zaxis=list(title='$M', range=c(-900,250)),
                      camera=list(eye=list(x=1.25*-1*1.5, y=1.25*-1*1.5, z=1.25*0.75*1.5))) # default angles for x, y, and z are 1.25. Multiply by proportions to adjust

green_option_dol_reduc = SD.PV[,,'all','all-f']-SD.PV[,,'f','all-f']
p15 = layout(add_surface(plot_ly(z=smooth_matrix(green_option_dol_reduc, window_size = 6)[idx_cf, idx_kg],
                                 y=c_f_vals_core, x=k_g_vals_core), colorscale = colscale),
             title = list(text='Dollar Change in Std. Dev. of NPV Costs<br>from Adding Green Option', y = 0.9),
             scene = scene_sd_delta,
             legend = list(title=list(text='$')))
orca(p15, file = "dol_change_std_dev_green_option.png", scale=0.75, width=800*0.8, height=800*0.8)

fossil_option_dol_reduc = SD.PV[,,'all','all-g']-SD.PV[,,'g','all-g']
p16 = layout(add_surface(plot_ly(z=smooth_matrix(fossil_option_dol_reduc, window_size = 8)[idx_cf, idx_kg],
                                 y=c_f_vals_core, x=k_g_vals_core), colorscale = colscale),
             title = list(text='Dollar Change in Std. Dev. of NPV Costs<br>from Adding Fossil Option', y = 0.9),
             scene = scene_sd_delta,
             legend = list(title=list(text="$M")))
orca(p16, file = "dol_change_std_dev_fossil_option.png", scale=0.75, width=800*0.8, height=800*0.8)

# Procurement value
green_own_pct_reduc = (SD.PV[,,'all','all-g']/SD.PV[,,'all','all-f']-1)*100
p17 = layout(add_surface(plot_ly(z=smooth_matrix(green_own_pct_reduc)[idx_cf, idx_kg],
                                 y=c_f_vals_core, x=k_g_vals_core), colorscale = colscale),
             title = list(text='Percent Change in Std. Dev. of NPV Costs<br>from Procuring Green Asset', y = 0.9),
             scene = scene_pct,
             legend = list(title=list(text='%')))
orca(p17, file = "pct_change_std_dev_green_procurement.png", scale=0.75, width=800*0.8, height=800*0.8)

green_own_dol_reduc = SD.PV[,,'all','all-g']-SD.PV[,,'all','all-f']
p18 = layout(add_surface(plot_ly(z=smooth_matrix(green_own_dol_reduc)[idx_cf, idx_kg],
                                 y=c_f_vals_core, x=k_g_vals_core), colorscale = colscale),
             title = list(text='Dollar Change in Std. Dev. of NPV Costs<br>from Procuring Green Asset', y = 0.9),
             scene = scene_sd_delta,
             legend = list(title=list(text='$M')))
orca(p18, file = "dol_change_std_dev_green_procurement.png", scale=0.75, width=800*0.8, height=800*0.8)

# Near term
green_own_pct_reduc_near = (SD.PV.near[,,'all','all-g']/SD.PV.near[,,'all','all-f']-1)*100
p17 = layout(add_surface(plot_ly(z=smooth_matrix(green_own_pct_reduc_near)[idx_cf, idx_kg],
                                 y=c_f_vals_core, x=k_g_vals_core), colorscale = colscale),
             title = list(text='Percent Change in Std. Dev. of NPV Costs<br>from Procuring Green Asset<br>Near-term Only', y = 0.9),
             scene = scene_pct,
             legend = list(title=list(text='%')))
orca(p17, file = "pct_change_std_dev_green_procurement_near.png", scale=0.75, width=800*0.8, height=800*0.8)

green_own_dol_reduc_near = SD.PV.near[,,'all','all-g']-SD.PV.near[,,'all','all-f']
p18 = layout(add_surface(plot_ly(z=smooth_matrix(green_own_dol_reduc_near)[idx_cf, idx_kg],
                                 y=c_f_vals_core, x=k_g_vals_core), colorscale = colscale),
             title = list(text='Dollar Change in Std. Dev. of NPV Costs<br>from Procuring Green Asset<br>Near-term Only', y = 0.9),
             scene = scene_sd_delta,
             legend = list(title=list(text='$M')))
orca(p18, file = "dol_change_std_dev_green_procurement_near.png", scale=0.75, width=800*0.8, height=800*0.8)


i=median(idx_cf); j=median(idx_kg)
# i=length(idx_cf); j=1
# i=1; j=length(idx_kg)

c_f_vals[i]; k_g_vals[j]
# In general, option value is larger than procurement value.
# Option value is the value from adding optionality
# Call the sd in an all-g eqm state A, the SD in a all-f fossil-only strategy B,
# and the SD in the all-f eqm state C. B>C by optionality. So -B<-C and hence A-B<A-C.

# Expectation (in theory should replace this with the value function)
E.PV[i,j,'all','all-g']-E.PV[i,j,'f','all-f'] # option value plus owning it
E.PV[i,j,'all','all-f']-E.PV[i,j,'f','all-f'] # "pure" option value only
E.PV[i,j,'all','all-g']-E.PV[i,j,'all','all-f'] # just procurement value
E.PV[i,j,'f','all-f']; E.PV[i,j,'all','all-f']; E.PV[i,j,'all','all-g']
Es_plot = c(E.PV[i,j,'f','all-f'], E.PV[i,j,'g','all-g'], E.PV[i,j,'all','all-f'], E.PV[i,j,'all','all-g'])
names(Es_plot) = c('Fossil-only\nStrategy', 'Green-only\nStrategy',
                    'Both Options,\nAll-Fossil Chosen', 'Both Options,\nAll-Green Chosen')

# Standard Deviations 
SD.PV[i,j,'all','all-g']-SD.PV[i,j,'f','all-f'] # option value plus owning it
SD.PV[i,j,'all','all-f']-SD.PV[i,j,'f','all-f'] # "pure" option value only
SD.PV[i,j,'all','all-g']-SD.PV[i,j,'all','all-f'] # just procurement value
SD.PV[i,j,'f','all-f']; SD.PV[i,j,'all','all-f']; SD.PV[i,j,'all','all-g']
SDs_plot = c(SD.PV[i,j,'f','all-f'], SD.PV[i,j,'g','all-g'], SD.PV[i,j,'all','all-f'], SD.PV[i,j,'all','all-g'])
names(SDs_plot) = c('Fossil-only\nStrategy', 'Green-only\nStrategy',
                    'Both Options,\nAll-Fossil Chosen', 'Both Options,\nAll-Green Chosen')

dev.off()
svg(filename=paste0(out_dir, 'std_dev_barchats.svg'), width=8)
barplot(SDs_plot, ylab='Standard Deviation in NPV Costs ($M)', col=4)
dev.off()

# Near-term only
SD.PV.near[i,j,'all','all-g']-SD.PV.near[i,j,'f','all-f'] # option value plus owning it
SD.PV.near[i,j,'all','all-f']-SD.PV.near[i,j,'f','all-f'] # "pure" option value only
SD.PV.near[i,j,'all','all-g']-SD.PV.near[i,j,'all','all-f'] # just procurement value
SD.PV.near[i,j,'f','all-f']; SD.PV.near[i,j,'all','all-f']; SD.PV.near[i,j,'all','all-g']
SDs_plot_near = c(SD.PV.near[i,j,'f','all-f'], SD.PV.near[i,j,'g','all-g'], SD.PV.near[i,j,'all','all-f'], SD.PV.near[i,j,'all','all-g'])
names(SDs_plot_near) = c('Fossil-only\nStrategy', 'Green-only\nStrategy',
                    'Both Options,\nAll-Fossil Chosen', 'Both Options,\nAll-Green Chosen')

dev.off()

svg(filename=paste0(out_dir, 'std_dev_barchats_nearterm.svg'), width=8)
par(mfrow=c(1,2))
barplot(SDs_plot, ylab='Standard Deviation in NPV Costs ($M)', col=4)
barplot(SDs_plot_near, ylab='Standard Deviation in NPV Costs ($M)\nNear Term Only', col=4)
dev.off()

svg(filename=paste0(out_dir, 'std_dev_barchats_near_vs_long_term.svg'), width=8)
barplot(rbind(SDs_plot, SDs_plot_near), beside=T, ylab='Standard Deviation in NPV Costs ($M)', col=c(4,5))
grid(nx=NA, ny=NULL, lty=1, col='gray90')
barplot(rbind(SDs_plot, SDs_plot_near), beside=T, col=c(4,5), add=T)
legend(y=700, x=5, horiz=F, legend=c('Long-run', 'First 10 Years Only'), fill=4:5, x.intersp = 0.25, bty='n')
dev.off()

# Very similar expected values at the median/aligned values
test_t10_f$value_func$V_min[i,j,string2bin('fffffffff')]
test_t10_g$value_func$V_min[i,j,string2bin('ggggggggg')]
test_t10$value_func$V_min[i,j,string2bin('ggggggggg')]
test_t10$value_func$V_min[i,j,string2bin('fffffffff')]
plot(test_t10$value_func$V_min[9,21,], ylim=c(0,7500), type='l')

# Owning fossil vs all-g
SD.PV[i,j,'all','all-f']-SD.PV[i,j,'g','all-g'] # option value plus owning it
SD.PV[i,j,'all','all-g']-SD.PV[i,j,'g','all-g'] # "pure" option value only
SD.PV[i,j,'all','all-f']-SD.PV[i,j,'all','all-g'] # just procurement value
SD.PV[i,j,'g','all-g']; SD.PV[i,j,'all','all-g']; SD.PV[i,j,'all','all-f']


# Option value vs procurement value. When is option value bigger (i.e., more negative)?
# It looks like in the regions when green
image(z=(SD.PV[,,'all','all-f']-SD.PV[,,'f','all-f'])<(SD.PV[,,'all','all-g']-SD.PV[,,'all','all-f']),
      x=c_f_vals, y=k_g_vals)

# # Option value (%)
# # Green
# layout(add_surface(plot_ly(z=(SD.PV[idx_cf, idx_kg,'all','all-f']/SD.PV[idx_cf, idx_kg,'f','all-f']-1)*100, y=c_f_vals_core, x=k_g_vals_core), colorscale = colscale),
#        title = 'Value of Adding Green Option\nLong Term',
#        scene = scene_pct,
#        legend = list(title=list(text='%')))
# # Fossil
# layout(add_surface(plot_ly(z=(SD.PV[idx_cf, idx_kg,'all','all-g']/SD.PV[idx_cf, idx_kg,'g','all-g']-1)*100, y=c_f_vals_core, x=k_g_vals_core), colorscale = colscale),
#        title = 'Value of Adding Fossil Option\nLong Term',
#        scene = scene_pct,
#        legend = list(title=list(text='%')))


# Option value ($)
layout(add_surface(plot_ly(z=(SD.PV[idx_cf, idx_kg,'all','all-f']-SD.PV[idx_cf, idx_kg,'f','all-f']), y=c_f_vals_core, x=k_g_vals_core), colorscale = colscale),
       title = 'Value of Adding Green Option\nLong Term',
       scene = scene_sd_delta,
       legend = list(title=list(text='$M')))
layout(add_surface(plot_ly(z=(SD.PV[idx_cf, idx_kg,'all','all-f']-SD.PV[idx_cf, idx_kg,'g','all-f']), y=c_f_vals_core, x=k_g_vals_core), colorscale = colscale),
       title = 'Value of Adding Fossil Option\nLong Term',
       scene = scene_sd_delta,
       legend = list(title=list(text='$M')))

# Procurement value ($)
scene_sd2 = list(xaxis=list(title='Green CAPEX'),
                 yaxis=list(title='Fossil OPEX'),
                 zaxis=list(title='$M', range=c(-5,250)),
                 camera=list(eye=list(x=1.25*-1, y=1.25*-1, z=1.25*0.75))) # default angles for x, y, and z are 1.25. Multiply by proportions to adjust
scene_pct2 = list(xaxis=list(title='Green CAPEX'),
                 yaxis=list(title='Fossil OPEX'),
                 zaxis=list(title='%', range=c(-100,5)),
                 camera=list(eye=list(x=1.25*-1, y=1.25*-1, z=1.25*0.75)))

layout(add_surface(plot_ly(z=SD.PV[idx_cf, idx_kg,'all','all-f']-SD.PV[idx_cf, idx_kg,'all','all-g'], y=c_f_vals_core, x=k_g_vals_core), colorscale = colscale),
       title = 'Reduction in SD ($) from Procuring Green Portfolio\nLong Term',
       scene = scene_sd2,
       legend = list(title=list(text='$M')))
layout(add_surface(plot_ly(z=(SD.PV[idx_cf, idx_kg,'all','all-g']/SD.PV[idx_cf, idx_kg,'all','all-f']-1)*100, y=c_f_vals_core, x=k_g_vals_core), colorscale = colscale),
       title = 'Reduction in SD (%) from Procuring Green Portfolio\nLong Term',
       scene = scene_pct,
       legend = list(title=list(text='%')))
# once we have the set of figures we want, run it 4 times:
# 1) the above "normalized" parameterization
# 2) a similar one but where the central starting values for k_g and k_f are the same, and same for c_f and c_g,
# so the *only* difference is uncertainty
# 3) with "realistic" parameter values for EVs/ICEs
# 4) with "realistic" parameter vavlues for NG/wind-solar

# bks.near <- bks.near <- seq(min(E.PV.near), max(E.PV), length.out=ncols+1)
par(mfrow=c(1,3))
# image(E.PV.near[idx_cf, idx_kg,'f'], x=c_f_vals_core, y=k_g_vals_core, main='Fossil Only, Short-run', breaks = bks.near, col = hcl.colors(ncols, 'Reds', rev=T))
# image(E.PV.near[idx_cf, idx_kg,'g'], x=c_f_vals_core, y=k_g_vals_core, main='Green Only, Short-run', breaks = bks.near, col = hcl.colors(ncols, 'Reds', rev=T))
# image(E.PV.near[idx_cf, idx_kg,'all'], x=c_f_vals_core, y=k_g_vals_core, main='Both, Short-run', breaks = bks.far, col = hcl.colors(ncols, 'Reds', rev=T))
image(E.PV[idx_cf, idx_kg,'f'], x=c_f_vals_core, y=k_g_vals_core, main='Fossil Only, Long-run', breaks = bks.far, col = hcl.colors(ncols, 'Reds', rev=T))
image(E.PV[idx_cf, idx_kg,'g'], x=c_f_vals_core, y=k_g_vals_core, main='Green Only, Long-run', breaks = bks.far, col = hcl.colors(ncols, 'Reds', rev=T))
image(E.PV[idx_cf, idx_kg,'all'], x=c_f_vals_core, y=k_g_vals_core, main='Both, Long-run', breaks = bks.far, col = hcl.colors(ncols, 'Reds', rev=T))


# bks.near <- bks.near <- seq(min(SD.PV.near), max(SD.PV), length.out=ncols+1)
par(mfrow=c(2,3))
image(SD.PV.near[idx_cf, idx_kg,'f'], x=c_f_vals_core, y=k_g_vals_core, main='Fossil Only, Short-run', breaks = bks.near, col = hcl.colors(ncols, 'Reds', rev=T))
image(SD.PV.near[idx_cf, idx_kg,'g'], x=c_f_vals_core, y=k_g_vals_core, main='Green Only, Short-run', breaks = bks.near, col = hcl.colors(ncols, 'Reds', rev=T))
image(SD.PV.near[idx_cf, idx_kg,'all'], x=c_f_vals_core, y=k_g_vals_core, main='Both, Short-run', breaks = bks.far, col = hcl.colors(ncols, 'Reds', rev=T))
image(SD.PV[idx_cf, idx_kg,'f'], x=c_f_vals_core, y=k_g_vals_core, main='Fossil Only, Long-run', breaks = bks.far, col = hcl.colors(ncols, 'Reds', rev=T))
image(SD.PV[idx_cf, idx_kg,'g'], x=c_f_vals_core, y=k_g_vals_core, main='Green Only, Long-run', breaks = bks.far, col = hcl.colors(ncols, 'Reds', rev=T))
image(SD.PV[idx_cf, idx_kg,'all'], x=c_f_vals_core, y=k_g_vals_core, main='Both, Long-run', breaks = bks.far, col = hcl.colors(ncols, 'Reds', rev=T))

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
(median(k_g_vals_core) + 4*2)*(1+.1)/.1 
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
image(opt_strat[idx_cf, idx_kg,1], x=c_f_vals_core, y=k_g_vals_core)

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
test_t4$value_func$V_min[index_nearest(median(c_f_vals_core), c_f_vals_core),
                         index_nearest(median(k_g_vals_core), k_g_vals_core),
                         string2bin('fff')]

mean(t4_g_rc)
test_t4_g$value_func$V_min[index_nearest(median(c_f_vals_core), c_f_vals_core),
                           index_nearest(median(k_g_vals_core), k_g_vals_core),
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

# layout(add_surface(plot_ly(z=kernel2dsmooth(x=100*SD.pct.diff.near, kernel.type='disk', r=5)[-(1:3),-c(1:5, 49:51)], y=c_f_vals_core[-(1:3)], x=k_g_vals_core[-c(1:5, 49:51)])),
layout(add_surface(plot_ly(z=green_t10_OV[idx_cf, idx_kg,string2bin('ggg')], y=c_f_vals_core, x=k_g_vals_core), colorscale = colscale),
       title = 'Green Option Value, L=10',
       legend = list(title=list(text='$M Difference')),
       scene = list(xaxis=list(title='Green CAPEX'),
                    yaxis=list(title='Fossil OPEX'),
                    zaxis=list(title='$M')))

layout(add_surface(plot_ly(z=green_t10_OV[idx_cf, idx_kg,string2bin('fff')]-green_t4_OV[idx_cf, idx_kg,string2bin('fff')] , y=c_f_vals_core, x=k_g_vals_core), colorscale = colscale),
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
