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

############# Make figures
out_dir = paste0(root,'/OneDrive - rff/Documents - RPE-Electric Power/OPEX CAPEX Price Risk/output/figures/matched_costs/')
out_dir = paste0(root,'/OneDrive - rff/Documents - RPE-Electric Power/OPEX CAPEX Price Risk/output/figures/matched_costs_test/')
dir.create(out_dir)
setwd(out_dir)
# out_dir=='C:/Users/prest//OneDrive - rff/Documents - RPE-Electric Power/OPEX CAPEX Price Risk/output/figures/'

### Optimal strategy matrix
library(colorspace)
opt_strat = 1*(test_t10$value_func$V_f[idx_cf, idx_kg,]<test_t10$value_func$V_g[idx_cf, idx_kg,])

# Note: conceptually, it can vary by portfolio, but in practice it doesn't
apply(opt_strat, 3, mean)
table(opt_strat[idx_cf, idx_kg,1]==opt_strat[idx_cf, idx_kg, dim(opt_strat)[3]])

svg(filename=paste0(out_dir, 'optimal_strategy.svg'), width=8, height=8)
image(opt_strat[idx_cf,idx_kg,1], x=c_f_vals_core, y=k_g_vals_core, col=lighten(c(3,2),0.5),
      xlab='Fossil OPEX ($M/year)', ylab='Green CAPEX ($M)')
text(x=30, y=200, labels='Green preferred', col=darken(3,0.5), font=2, cex=1.5)
text(x=15, y=700, labels='Fossil preferred', col=darken(2,0.5), font=2, cex=1.5)
points(x=median(c_f_vals_core), y=median(k_g_vals_core), pch=19, cex=2, col=4)
dev.off()

### Barcharts of SD of costs
i=median(idx_cf); j=median(idx_kg)
c_f_vals[i]; k_g_vals[j]

# Standard Deviations 
SDs_plot = c(SD.PV[i,j,'f','all-f'], SD.PV[i,j,'g','all-g'], SD.PV[i,j,'all','all-f'], SD.PV[i,j,'all','all-g'])
names(SDs_plot) = c('Fossil-only\nStrategy', 'Green-only\nStrategy',
                    'Both Options,\nAll-Fossil Chosen', 'Both Options,\nAll-Green Chosen')

# Barchart, long-term SD
svg(filename=paste0(out_dir, 'std_dev_barchats.svg'), width=8)
barplot(SDs_plot, ylab='Standard Deviation in NPV Costs ($M)', col=4)
dev.off()

# Near-term only
SDs_plot_near = c(SD.PV.near[i,j,'f','all-f'], SD.PV.near[i,j,'g','all-g'], SD.PV.near[i,j,'all','all-f'], SD.PV.near[i,j,'all','all-g'])
names(SDs_plot_near) = c('Fossil-only\nStrategy', 'Green-only\nStrategy',
                         'Both Options,\nAll-Fossil Chosen', 'Both Options,\nAll-Green Chosen')

dev.off()

# "the uncertainty reduction from replacing the fossil portfolio with a green 
# one is larger in the near-term (red bars) than in the long-term (blue bars),
# both in terms of total dollars and as a percent."
diff(SDs_plot_near[3:4]); SDs_plot_near[4]/SDs_plot_near[3]-1
diff(SDs_plot[3:4]); SDs_plot[4]/SDs_plot[3]-1

# Barchart, near-term SD
svg(filename=paste0(out_dir, 'std_dev_barchats_nearterm.svg'), width=8)
par(mfrow=c(1,2))
barplot(SDs_plot, ylab='Standard Deviation in NPV Costs ($M)', col=4)
barplot(SDs_plot_near, ylab='Standard Deviation in NPV Costs ($M)\nNear Term Only', col=4)
dev.off()

# Long-term and near-term side-by-side
svg(filename=paste0(out_dir, 'std_dev_barchats_near_vs_long_term.svg'), width=8)
barplot(rbind(SDs_plot, SDs_plot_near), beside=T, ylab='Standard Deviation in NPV Costs ($M)', col=c(4,5))
grid(nx=NA, ny=NULL, lty=1, col='gray90')
barplot(rbind(SDs_plot, SDs_plot_near), beside=T, col=c(4,5), add=T)
legend(y=700, x=5, horiz=F, legend=c('Long-run', 'First 10 Years Only'), fill=4:5, x.intersp = 0.25, bty='n')
dev.off()

#### 3D Surface plots
library(plotly)
library(orca)
library(RColorBrewer)

# display.brewer.all()
colscale = 'RdOrYl'

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

## Std Dev under each scenario
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
# To do here: 
# 1) fix the color scale so that reductions are always blue (or green), and not red

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

image(green_option_pct_reduc<fossil_option_pct_reduc)
image(green_option_dol_reduc<fossil_option_dol_reduc)

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

# Fossil exposure plot
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
