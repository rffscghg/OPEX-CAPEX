# Main script to run OPEX-CAPEX model

out.dir = 'figures'

library(tidyverse)
library(hms)

source("src/monte.r")
source("src/monte_skipVFI.r")
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

# load(paste0(out.dir, 'solved_t4_cf10-60-by-1_kg100-600-by-10.RData')) # 3 minutes
# load(paste0(out.dir, 'solved_t10_cf10-60-by-1_kg100-600-by-10.RData')) # 168 mins (60x)

test_t4 <- monte_carlo(
  c_f_vals = seq(10, 60, by = 1),
  k_g_vals = seq(100, 600, by = 10),
  k_f = 300,
  c_g = 1,
  sigma_cf = .05,
  sigma_kg = .05,
  t = 4,
  r = 0.1,
  const_scrap = TRUE,
  max_iter = 1000,
  threshold = 1e-3,
  verbose = TRUE,
  V_init = if (exists("test_t4")) test_t4$value_func,
  start_assets = "fgf"
)

opt_strat_t4 = 1*(test_t4$value_func$V_f[,,string2bin('fgf')] > test_t4$value_func$V_g[,,string2bin('fgf')])
image(opt_strat_t4, x=seq(10, 60, by = 1), y=seq(100, 600, by = 10), 
      xlab=expression('c'[f]), ylab=expression('k'[g]))
diff(apply( 1*(test_t4$value_func$V_f > test_t4$value_func$V_g), 3, mean))

# test_t10 <- monte_carlo(
#   c_f_vals = seq(10, 60, by = 1),
#   k_g_vals = seq(100, 600, by = 10),
#   k_f = 300,
#   c_g = 1,
#   sigma_cf = .05,
#   sigma_kg = .05,
#   t = 10,
#   r = 0.1,
#   const_scrap = TRUE,
#   max_iter = 1000,
#   threshold = 1e-3,
#   verbose = TRUE,
#   V_init = if (exists("test_t10")) test_t10$value_func,
#   # skipVFI=TRUE,
#   start_assets = "fgf"
# )
#
# save(test_t4,   file=paste0(out.dir, 'solved_t4_cf10-60-by-1_kg100-600-by-10.RData')) # 3 minutes
# save(test_t10, file=paste0(out.dir, 'solved_t10_cf10-60-by-1_kg100-600-by-10.RData')) # 208 mins (70x)

# key insight from the monte carlo:
#    when you're starting deep in the optimal-g (or optimal-f) region and running the scenario
#    for >L periods, the starting asset portfolio doesn't matter too much.
#    Why? Because after t years you have turned over//flushed out the original 
#    portfolio and replaced it with an optimal one.
#    So perhaps we want to just separately consider the near-term and long-term
#    variance, where the long-run variance is likely smaller (as a %)
#    Further, more-g portfolios should be more insulated in the near term than more-f
#    portfolios, but both equally insulated in the long term. THIS IS THE KEY RESULT!
#    In other words, optimization ensures that the long-run doesn't matter,
#    but in the near-term "g" portfolios are more protected from risk

opt_strat_t10 = 1*(test_t10$value_func$V_f[,,string2bin('fgf')] > test_t10$value_func$V_g[,,string2bin('fgf')])
image(opt_strat_t10, x=seq(10, 60, by = 1), y=seq(100, 600, by = 10), 
      xlab=expression('c'[f]), ylab=expression('k'[g]))
diff(apply( 1*(test_t10$value_func$V_f > test_t10$value_func$V_g), 3, mean))

opt2 = opt_strat_t4 + opt_strat_t10

c_f_vals = seq(10, 60, by = 1)
k_g_vals = seq(100, 600, by = 10)

png(paste0(out.dir,"figures/optimal_strategies_v2.png"),
    width = 480*2, height = 480*2, pointsize=12*1.5)
image(opt2, x=c_f_vals, y=k_g_vals, col=c(colorspace::lighten('#533600', 0.2),
                                          colorspace::lighten('#90B88E', 0),
                                          colorspace::lighten('darkgreen', 0.75)),
      ylab='Green Capital Cost', xlab='Fossil Operating Cost')
text(x=20, y=500, labels='Fossil optimal\nunder L=4 & L=10')
text(x=30, y=250, labels='Green optimal under L=4 & L=10')
text(x=45, y=500, labels='Green optimal\nunder L=10 only')
dev.off()

str(test_t4$realized_costs)
str(test_t4$pick_f)

png(paste0(out.dir,"figures/illustrative_brownian_motion.png"),
    width = 480*2, height = 480*2, pointsize=12*1.5)
image(opt_strat_t4, x=seq(10, 60, by = 1), y=seq(100, 600, by = 10), 
      col=c(colorspace::lighten('#533600', 0.2),
            colorspace::lighten('darkgreen', 0.75)),
      ylab='Green Capital Cost', xlab='Fossil Operating Cost')
lines(x=c(median(c_f_vals), test_t4$c_f[,174]), y=c(median(k_g_vals),test_t4$k_g[,174]), col=2, type='b', pch=19)
lines(x=c(median(c_f_vals), test_t4$c_f[,161]), y=c(median(k_g_vals),test_t4$k_g[,161]), col=4, type='b', pch=19)
points(x=median(c_f_vals), y=median(k_g_vals), pch=19, col=7, cex=2.5)
dev.off()


# Run Monte Carlo across full grid
r=0.1
E.PV.near <- SD.PV.near <- test_t4$value_func$V_min
E.PV <- SD.PV <- test_t4$value_func$V_min
st = Sys.time()
for (i in 1:length(c_f_vals)) {          
  for (j in 1:length(k_g_vals)) {
    for (k in 1:(2^(4-1))) { 
      
      test_t4.temp <- monte_carlo(
        c_f_vals = seq(10, 60, by = 1),
        k_g_vals = seq(100, 600, by = 10),
        k_f = 300,
        c_g = 1,
        sigma_cf = .05,
        sigma_kg = .05,
        t = 4,
        r = r,
        const_scrap = TRUE,
        max_iter = 1000,
        threshold = 1e-3,
        verbose = TRUE,
        V_init = if (exists("test_t4")) test_t4$value_func,
        skipVFI=TRUE,
        start_cf = c_f_vals[i],            # Starting value for c_f
        start_kg = k_g_vals[j],            # Starting value for k_g
        start_assets = bin2string(k, t=4)
      )
      
      rel.cost.temp = test_t4.temp$realized_costs
      PV.Costs.near = apply(rel.cost.temp[1:4,], 2, function(x) sum(x*(1+r)^-(1:4)))
      PV.Costs = apply(rel.cost.temp, 2, function(x) sum(x*(1+r)^-(1:nrow(rel.cost.temp))))
      
      E.PV.near[i,j,k] = mean(PV.Costs.near)
      SD.PV.near[i,j,k] = sd(PV.Costs.near)
      E.PV[i,j,k] = mean(PV.Costs)
      SD.PV[i,j,k] = sd(PV.Costs)
    }
    message('k_g = ', k_g_vals[j],' done at ', Sys.time())
  }
  message('c_f = ', c_f_vals[i],' done at ', Sys.time())
}
ed = Sys.time()
difftime(ed, st)

find.boundary = function(opt_strat) {
  nk = ncol(opt_strat)
  nf = nrow(opt_strat)
  y0.idx = which.min(opt_strat[1,]==1)
  y0 = mean(as.numeric(dimnames(opt_strat)[[2]][c(y0.idx, y0.idx+1)]))
  y1.idx = which.min(opt_strat[nk,]==1)
  y1 = mean(as.numeric(dimnames(opt_strat)[[2]][c(y1.idx, y1.idx+1)]))
  y1 = as.numeric(dimnames(opt_strat)[[2]][which.min(opt_strat[nk,]==1)])
  x0 = as.numeric(dimnames(opt_strat)[[1]][1])
  x1 = as.numeric(dimnames(opt_strat)[[1]][nf])
  
  lmfit = lm(c(y1,y0) ~ c(x1, x0))
  # predict(lmfit, newdata=as.numeric(dimnames(opt_strat)[[1]]))
  return(lmfit)
}

# save(E.PV.near,   file=paste0(out.dir, 'montecarlo_EPVnear_t4_cf10-60-by-1_kg100-600-by-10.RData')) # 2.7 hours
# save(SD.PV.near,   file=paste0(out.dir, 'montecarlo_SDPVnear_t4_cf10-60-by-1_kg100-600-by-10.RData'))
# save(E.PV,   file=paste0(out.dir, 'montecarlo_EPV_t4_cf10-60-by-1_kg100-600-by-10.RData')) # 2.7 hours
# save(SD.PV,   file=paste0(out.dir, 'montecarlo_SDPV_t4_cf10-60-by-1_kg100-600-by-10.RData'))

library(plotly)
add_surface(plot_ly(z=SD.PV.near[,,string2bin('ggg')], x=c_f_vals, y=k_g_vals))
add_surface(plot_ly(z=SD.PV.near[,,string2bin('ggg')]-SD.PV.near[,,string2bin('fff')], x=c_f_vals, y=k_g_vals))
add_surface(plot_ly(z=SD.PV.near[,,string2bin('ggg')]/SD.PV.near[,,string2bin('fff')] - 1, x=c_f_vals, y=k_g_vals))

SD.pct.diff.near = SD.PV.near[,,string2bin('ggg')]/SD.PV.near[,,string2bin('fff')] - 1
SD.pct.diff.long = SD.PV[,,string2bin('ggg')]/SD.PV[,,string2bin('fff')] - 1
mean(SD.pct.diff.near<0)
mean(SD.pct.diff.long<0)

# Unsmoothed
add_surface(plot_ly(z=SD.pct.diff.near, y=c_f_vals, x=k_g_vals))
library(smoothie)

# % Diff
# Near term
# Note: the smoothing produces odd effects near the boundary. Drop some of these for visualization
layout(add_surface(plot_ly(z=kernel2dsmooth(x=100*SD.pct.diff.near, kernel.type='disk', r=5)[-(1:3),-c(1:5, 49:51)], y=c_f_vals[-(1:3)], x=k_g_vals[-c(1:5, 49:51)])),
       title = 'Near-term % Difference in Std. Dev. Costs',
       legend = list(title=list(text='% Difference')),
       scene = list(xaxis=list(title='Green CAPEX'),
                    yaxis=list(title='Fossil OPEX'),
                    zaxis=list(title='% Difference')))
# Long term
layout(add_surface(plot_ly(z=kernel2dsmooth(x=100*SD.pct.diff.long, kernel.type='disk', r=9)[-(1:3),-c(1:5, 49:51)], y=c_f_vals[-(1:3)], x=k_g_vals[-c(1:5, 49:51)])),
       title = 'Long-term % Difference in Std. Dev. Costs',
       legend = list(title=list(text='% Difference')),
       scene = list(xaxis=list(title='Green CAPEX'),
                    yaxis=list(title='Fossil OPEX'),
                    zaxis=list(title='% Difference', range=c(-30,5))))

SD.diff.near = SD.PV.near[,,string2bin('ggg')]-SD.PV.near[,,string2bin('fff')] 
SD.diff.long = SD.PV[,,string2bin('ggg')]-SD.PV[,,string2bin('fff')] 

# $ Diff
# Near term
layout(add_surface(plot_ly(z=kernel2dsmooth(x=SD.diff.near, kernel.type='disk', r=5)[,-c(1:5, 49:51)], y=c_f_vals[-(1:3)], x=k_g_vals[-c(1:5, 49:51)])),
       title = 'Near-term $M Difference in Std. Dev. Costs',
       legend = list(title=list(text='$M Difference')),
       scene = list(xaxis=list(title='Green CAPEX'),
                    yaxis=list(title='Fossil OPEX'),
                    zaxis=list(title='$M Difference', range=c(-15, 5))))
# Long term
layout(add_surface(plot_ly(z=kernel2dsmooth(x=SD.diff.long, kernel.type='disk', r=9)[,-c(1:5, 49:51)], y=c_f_vals[-(1:3)], x=k_g_vals[-c(1:5, 49:51)])),
       title = 'Long-term $M Difference in Std. Dev. Costs',
       legend = list(title=list(text='$M Difference')),
       scene = list(xaxis=list(title='Green CAPEX'),
                    yaxis=list(title='Fossil OPEX'),
                    zaxis=list(title='$M Difference', range=c(-15, 5))))
