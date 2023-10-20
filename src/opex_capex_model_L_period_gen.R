palette(c(rgb(4, 39, 60, maxColorValue=255), # RFF black
          rgb(116,100,94, maxColorValue=255), # RFF brown
          rgb(80,177,97, maxColorValue=255), # RFF green
          rgb(136,196,244, maxColorValue=255), # RFF blue
          rgb(255, 102, 99, maxColorValue=255), # RFF coral
          rgb(118, 94, 165, maxColorValue=255), # RFF purple
          rgb(244, 162, 95, maxColorValue=255), # RFF orange
          rgb(235, 211, 103, maxColorValue=255), # RFF yellow
          rgb(224, 228, 231, maxColorValue=255))) # RFF light gray

rm(list=ls())
set.seed(1)
# Discount rate
r = 0.1

# Plot parameter
n_colors = 200

#### Capital costs
## Green cost parameters
# Capital
k.g.0 = 450
mu.g = 0
# mu.g = -0.10
sigma.g = 0.05
# sigma.g = 0.05*10
# Operating
c.g = 0

## Fossil cost parameters
# Capital
# k.f = 277.5
k.f = 278 # for L=10

# Operating
# c.f.0 = 19 # for L=1
c.f.0 = 28 # for L=10
mu.f = 0
# mu.f = 0.10
sigma.f = 0.1118 # for L=10. To do: update
# sigma.f = 0.09133 # for L=1
# sigma.f = 0.09133/1000 # for L=1

# Asset life
L = 10 # To do: update
# L = 1
DF = 1/(1+r)^(1:L)

# Output/year
q = 1 # for L=10 (To do)
# q = 10

# Try this: what if the second asset has the same cost profile as fossil, but
# just was independent?
# k.g.0 = k.f
# c.g = c.f.0
# sigma.g = sigma.f


# Simulate cost process
# GBM function
generate.GBM = function(mu, sigma, T=100, B=1000, S0=1) {
  t = 1:T
  
  # Compute geometric brownian motion
  eps <- matrix(rnorm(T*B), ncol = B, nrow = T)  
  
  # S0 = 100 # starting value
  # mu = 0.02
  # sigma = 0.01
  
  x <- exp((mu - sigma * sigma / 2) + sigma * eps)
  x[1:5,1:5]
  
  S <- apply(rbind(rep(S0, B), x), 2, cumprod)
  return(S)
}
set.seed(1)
c.f = generate.GBM(mu=mu.f, sigma=sigma.f, S0=c.f.0, T=20, B=10000)

# Expected costs in L years
k.g.0*exp(mu.g*L) + sum(c.g*q*DF) # green
k.f + c.f.0*sum(exp(mu.f*L)*exp(mu.f*(1:L))*q*DF)  #  fossil
c.f.0*exp(mu.f*L); mean(c.f[L+1,]) # mean matches
c.f.0^2*exp(2*mu.f*L)*(exp(sigma.f^2*L) - 1); var(c.f[L+1,]) # var ~matches

# Variance in costs in L years
k.g.0^2*exp(2*mu.g*L)*(exp(sigma.g^2 * L) - 1) # green
c.f.replace.sum = k.f + apply(c.f, 2, function(x) sum(x[(L+2):(2*L+1)]*q*DF))
var(c.f.replace.sum) # fossil

# Expected LCOE as of today
k.g.0/sum((q*DF))
(k.f + sum(c.f.0*exp(mu.f*(1:L))*q*DF))/sum(q*DF)

# Set up matrices with conditional probability densities for k.g[t+1] given k.g[t],
# and similarly for c.f
# Discretize into a given number states of each, with the probability of going from 
# price level i to level j, Phi.g[i,j] for green and Phi.f[i,j] for fossil.
k.g.grid = seq(100, 800, by=10)
c.f.grid = seq(1, 40, by=1)

grid.length.g = length(k.g.grid)
grid.length.f = length(c.f.grid)
grid.length.g*grid.length.f # state space size: 2,840 points

# Create a function giving probability distribution of price
dgbm = function(x, mu, sigma, t, x0) {
  return(dlnorm(x, meanlog = log(x0) + (mu-1/2*sigma^2)*t, sdlog = sigma*sqrt(t)))
  # return(1/sqrt(2*pi)*1/(x*sigma*sqrt(t))*exp(-(log(x)-log(x0)-(mu-0.5*sigma^2)*t)^2/(2*sigma^2*t))) # same
  # https://en.wikipedia.org/wiki/Geometric_Brownian_motion
}
t=10
# Check:
dgbm(k.g.grid, mu=mu.g, sigma=sigma.g, t=t, x0=k.g.0)/
  dlnorm(k.g.grid, meanlog=log(k.g.0) + (mu.g-1/2*sigma.g^2)*t, sdlog=sigma.g*sqrt(t))

# Compute diffusion matrices
# (To do: change to also calculate these for 1 through L?)
Phi.g <- Phi.g.10 <- matrix(0, nrow=grid.length.g, ncol=grid.length.g)
rownames(Phi.g) <- colnames(Phi.g) <- k.g.grid
rownames(Phi.g.10) <- colnames(Phi.g.10) <- k.g.grid

Phi.f <- Phi.f.L <- matrix(0, nrow=grid.length.f, ncol=grid.length.f)
rownames(Phi.f) <- colnames(Phi.f) <- c.f.grid
rownames(Phi.f.L) <- colnames(Phi.f.L) <- c.f.grid

# for c.f
for (i in 1:grid.length.f) {
  Phi.f[i,] =    dgbm(c.f.grid, mu=mu.f, sigma=sigma.f, x0=c.f.grid[i], t=1)/sum(dgbm(c.f.grid, mu=mu.f, sigma=sigma.f, x0=c.f.grid[i], t=1))
  Phi.f.L[i,] = dgbm(c.f.grid, mu=mu.f, sigma=sigma.f, x0=c.f.grid[i], t=10)/sum(dgbm(c.f.grid, mu=mu.f, sigma=sigma.f, x0=c.f.grid[i], t=L))
}
# for k.g
for (i in 1:grid.length.g) {
  Phi.g[i,] =    dgbm(k.g.grid, mu=mu.g, sigma=sigma.g, x0=k.g.grid[i], t=1)/sum(dgbm(k.g.grid, mu=mu.g, sigma=sigma.g, x0=k.g.grid[i], t=1))
  Phi.g.10[i,] = dgbm(k.g.grid, mu=mu.g, sigma=sigma.g, x0=k.g.grid[i], t=10)/sum(dgbm(k.g.grid, mu=mu.g, sigma=sigma.g, x0=k.g.grid[i], t=L))
}
range(apply(Phi.g, 1, sum))
range(apply(Phi.g.10, 1, sum))

## Compute joint conditional probability distributions
# This will give a matrix of probabilities of jumping to a (k.g[t+h], c.f[t+h]) pair,
# given the current (k.g[t], c.f[t]) pair. We need such a matrix for every (k.g[t], c.f[t]) pair
# Hence we will have an array of matrices. Both the array and each matrix will have 
# rows corresponding to c.f.grid and columns matching k.g.grid.
# This should match the dimensions of the optimal strategy vector.
dummy.mat = matrix(0, nrow=grid.length.f, ncol=grid.length.g)
Phi.f.g <- Phi.f.g.10 <-  array(list(dummy.mat), dim=c(grid.length.f, grid.length.g), dimnames=list(c.f.grid, k.g.grid))

rep.row<-function(x,n) matrix(rep(x,each=n),nrow=n)
rep.col<-function(x,n) matrix(rep(x,each=n), ncol=n, byrow=TRUE)

for (i in 1:grid.length.f) for (j in 1:grid.length.g) {
  # note: this for t=1 step ahead
  # Phi.f.g[[i,j]] = rep.row(Phi.f[i,], grid.length.g) * rep.col(Phi.g[j,], grid.length.f) # old
  Phi.f.g[[i,j]] = rep.row(Phi.g[j,], grid.length.f) * rep.col(Phi.f[i,], grid.length.g) 
  Phi.f.g[[i,j]]  = Phi.f.g[[i,j]] / sum(Phi.f.g[[i,j]] ) # normalize to a probability matrix
  
  # for t=10 steps ahead (To do: change to calculate these for 1 through L?)
  # Phi.f.g.10[[i,j]] = rep.row(Phi.f.L[i,], grid.length.g) * rep.col(Phi.g.10[j,], grid.length.f) # old
  Phi.f.g.10[[i,j]] = rep.row(Phi.g.10[j,], grid.length.f) * rep.col(Phi.f.L[i,], grid.length.g) 
  Phi.f.g.10[[i,j]]  = Phi.f.g.10[[i,j]] / sum(Phi.f.g.10[[i,j]] ) # normalize to a probability matrix
  
  rownames(Phi.f.g[[i,j]]) <- rownames(Phi.f.g.10[[i,j]]) <-  c.f.grid
  colnames(Phi.f.g[[i,j]]) <- colnames(Phi.f.g.10[[i,j]]) <-  k.g.grid
}

# Visualize
image(Phi.f.g[['35','760']], x=c.f.grid, y=k.g.grid) # probability of going from (c.f=$35, k.g=$760) to every other pair
image(Phi.f.g.10[['35','760']], x=c.f.grid, y=k.g.grid)
image(Phi.f.g[['10','320']], x=c.f.grid, y=k.g.grid)
image(Phi.f.g.10[['10','320']], x=c.f.grid, y=k.g.grid)
range(apply(Phi.f.g, 1:2, function(x) sum(x[[1]]))) # check that every matrix's rows sum to 1
range(apply(Phi.f.g, 1:2, function(x) sum(x[[1]]))) # each matrix's probabilities should sum to 1

# But if you're sufficiently away from the edges, its essentially fine
plot(Phi.g['680',] ~ k.g.grid, type='l', ylim=c(0,0.50), col=5)
lines(Phi.g['120',] ~ k.g.grid, type='l', col=5)
lines(Phi.g['500',] ~ k.g.grid, type='l', col=5)
lines(Phi.g.10['680',] ~ k.g.grid, type='l', col=4)
lines(Phi.g.10['120',] ~ k.g.grid, type='l', col=4)
lines(Phi.g.10['500',] ~ k.g.grid, type='l', col=4)

# Create optimal strategy matrix for each (k.g, c.f)
# (initial values are arbitrary. They will be overwritten on each iteration)
opt.strat = array(NA, dim=c(grid.length.f, grid.length.g),
                  dimnames=list(c.f.grid, k.g.grid)) # optimal strategies: 1=fossil, 2=green
dim(opt.strat)
# initial guess: above-average c.f values and below-average k.g values will 
# encourage green investments that is, the lower triangle of opt.strat (high fossil cost/low green)
# is 2, and the rest is 1.
opt.strat[upper.tri(opt.strat)] = 1
opt.strat[lower.tri(opt.strat, diag=T)] = 2
table(opt.strat)

grid.size = grid.length.f*grid.length.g
# Initiate value function. One value for each state. initial value is arbitrary
V = array(runif(grid.size), dim=c(grid.length.f, grid.length.g),
          dimnames=list(c.f.grid, k.g.grid))

# check that dimensions of V, opt.strat, and Phi.f.g[[i,j]] match
if (!identical(dimnames(V), dimnames(opt.strat))) stop('Error on opt.strat')
if (!identical(dimnames(V), dimnames(Phi.f.g[[1,1]]))) stop('Error on Phi.f.g')

# Note: Phi.g.10 is not currently used in the model because we have assumed
# an asset life of 1 period. When we extend it to L>1, this will change. (To do)
V.calc = function(V) {
  Vnew = array(runif(grid.size), dim=c(grid.length.f, grid.length.g),
               dimnames=list(c.f.grid, k.g.grid))
  
  C.array = array(NA, dim=c(2, length(c.f.grid), length(k.g.grid)), 
                  dimnames = list(c('fossil','green'),
                                  c.f.grid, k.g.grid))
  
  # For each potential c.f and k.g value, calculate npv costs and the optimal investment
  for (i in 1:grid.length.f) for (j in 1:grid.length.g) {
    # To do: change C.f to include expected path of c.f through h=1,..,L, rather than just 1 year
    C.f = k.f + sum(c.f.grid[i]*exp(mu.f*(1:L))*q*DF) + sum(Phi.f.g.10[[i, j]]*V)/(1+r)^L
    C.g = k.g.grid[j] + sum(c.g*q*DF) + sum(Phi.f.g.10[[i, j]]*V)/(1+r)^L
    
    C.vec = c(C.f, C.g)
    opt.strat[i,j] = which.min(C.vec)
    C.array[,i,j] = C.vec
    
    Vnew[i,j] = min(C.vec)
  }
  
  return(list(Vnew=Vnew, opt.strat=opt.strat, npvcost=C.array))
}
delta = 1
iter = 0

st = Sys.time()
while (delta>1e-6) {
  Vnew = V.calc(V)$Vnew
  delta = max(Vnew-V)
  V = Vnew
  iter = iter+1
}
ed = Sys.time()
difftime(ed, st) 
# on Brian's machine as of 10/12 with a state space size of 2,840, this took 10 seconds to converge.
iter

V.out = V.calc(V)

delta_V = V.out$npvcost['green',,]-V.out$npvcost['fossil',,]
delta_V_breaks = c(seq(min(delta_V),0,length = n_colors/2), 0, seq(0,max(delta_V),length=n_colors/2))

# Examine solved model
png(filename = paste0('output/figures/npv_value_L',L,'_',Sys.Date(),'.png'), width = 480*3, height=480, pointsize = 24)
par(mfrow=c(1,3))
image(z=-V.out$Vnew, x=c.f.grid, y=k.g.grid, main='NPV of costs (optimal)', col=hcl.colors(n_colors, palette='Heat'), xlab='Fossil Operating Cost', ylab='Green Capital Cost')

image(delta_V, x=c.f.grid, y=k.g.grid,
      main='Green advantage\n(difference in $ NPV cost)',
      xlab='Fossil Operating Cost', ylab='Green Capital Cost',
      col=rev(hcl.colors(n_colors, palette='Green-Orange')), breaks = delta_V_breaks)

image(-V.out$opt.strat, x=c.f.grid, y=k.g.grid, main='Optimal current-period strategy', col=hcl.colors(n_colors, palette='Green-Orange'), xlab='Fossil Operating Cost', ylab='Green Capital Cost')
dev.off()

# because our c.f.0 and k.g.0 values equated the expectation and variance of costs
# between the two assets, we are on an edge case where very small changes in either
# in the level of either cost parameter can flip the decision at these prices
mean(V)
mean(V.out$npvcost['green',,]-V.out$npvcost['fossil',,])
# Optimal strategy 
V.out$opt.strat[as.character(c.f.0), as.character(k.g.0)]
V.out$opt.strat[as.character((c.f.0-2):(c.f.0+2)), as.character(seq(k.g.0-20,k.g.0+20, by=10))]

V.out$Vnew[as.character(c.f.0), as.character(k.g.0)]
V.out$Vnew[as.character((c.f.0-2):(c.f.0+2)), as.character(seq(k.g.0-20,k.g.0+20, by=10))]

### Run with only green, and again with only fossil
# Do the optimization with green as the only option. 
# The easiest but brute-forcey way to do this is to simply replace C.f with C.g's value, so the options are equivalent.
V.calc.onlyg = function(V) {
  Vnew = array(runif(grid.size), dim=c(grid.length.f, grid.length.g),
               dimnames=list(c.f.grid, k.g.grid))
  
  C.array = array(NA, dim=c(2, length(c.f.grid), length(k.g.grid)), 
                  dimnames = list(c('fossil','green'),
                                  c.f.grid, k.g.grid))
  
  # For each potential c.f and k.g value, calculate npv costs and the optimal investment
  for (i in 1:grid.length.f) for (j in 1:grid.length.g) {
    # C.f = k.f + sum(c.f.grid[i]*exp(mu.f*(1:L))*q*DF) + sum(Phi.f.g.10[[i, j]]*V)/(1+r)^L
    C.g = k.g.grid[j] + sum(c.g*q*DF) + sum(Phi.f.g.10[[i, j]]*V)/(1+r)^L
    C.f = C.g
    
    C.vec = c(C.f, C.g)
    opt.strat[i,j] = which.min(C.vec)
    C.array[,i,j] = C.vec
    
    Vnew[i,j] = min(C.vec)
  }
  
  return(list(Vnew=Vnew, opt.strat=opt.strat, npvcost=C.array))
}

iter = 0
delta = 1
st = Sys.time()
while (delta>1e-6) {
  Vnew.onlyg = V.calc.onlyg(V)$Vnew
  delta = max(Vnew.onlyg-V)
  V = Vnew.onlyg
  iter = iter+1
}
ed = Sys.time()
difftime(ed, st) 

V.g = V.calc.onlyg(V)

par(mfrow=c(2,2))
image(z=V.g$Vnew, x=c.f.grid, y=k.g.grid, main='NPV of costs (optimal)',xlab='Fossil Operating Cost', ylab='Green Capital Cost')
points(x=c.f.0, y=k.g.0, pch=19, cex=3, col=4)

image(V.g$npvcost['green',,]-V.g$npvcost['fossil',,], x=c.f.grid, y=k.g.grid,
      main='Green advantage (difference in $ NPV cost)',
      xlab='Fossil Operating Cost', ylab='Green Capital Cost',
      col=rev(hcl.colors(200, palette='Red-Green')))

image(V.g$opt.strat, x=c.f.grid, y=k.g.grid, main='Optimal strategy', col=hcl.colors(200, palette='Green-Orange'), xlab='Fossil Operating Cost', ylab='Green Capital Cost')
points(x=c.f.0, y=k.g.0, pch=19, cex=3, col=4)

# Do the optimization with fossil as the only option. 
# The easiest but brute-forcey way to do this is to simply replace C.g with C.f's value, so the options are equivalent.
V.calc.onlyf = function(V) {
  Vnew = array(runif(grid.size), dim=c(grid.length.f, grid.length.g),
               dimnames=list(c.f.grid, k.g.grid))
  
  C.array = array(NA, dim=c(2, length(c.f.grid), length(k.g.grid)), 
                  dimnames = list(c('fossil','green'),
                                  c.f.grid, k.g.grid))
  
  # For each potential c.f and k.g value, calculate npv costs and the optimal investment
  for (i in 1:grid.length.f) for (j in 1:grid.length.g) {
    C.f = k.f + sum(c.f.grid[i]*exp(mu.f*(1:L))*q*DF) + sum(Phi.f.g.10[[i, j]]*V)/(1+r)^L
    # C.g = k.g.grid[j] + sum(c.g*q*DF) + sum(Phi.f.g.10[[i, j]]*V)/(1+r)^L
    C.g = C.f
    
    C.vec = c(C.f, C.g)
    opt.strat[i,j] = which.min(C.vec)
    C.array[,i,j] = C.vec
    
    Vnew[i,j] = min(C.vec)
  }
  
  return(list(Vnew=Vnew, opt.strat=opt.strat, npvcost=C.array))
}

iter = 0
delta = 1
st = Sys.time()
while (delta>1e-6) {
  Vnew.onlyg = V.calc.onlyf(V)$Vnew
  delta = max(Vnew.onlyg-V)
  V = Vnew.onlyg
  iter = iter+1
}
ed = Sys.time()
difftime(ed, st) 

V.f = V.calc.onlyf(V)

par(mfrow=c(2,2))
image(z=V.f$Vnew, x=c.f.grid, y=k.g.grid, main='NPV of costs (optimal)',xlab='Fossil Operating Cost', ylab='Green Capital Cost')
points(x=c.f.0, y=k.g.0, pch=19, cex=3, col=4)

image(V.f$npvcost['green',,]-V.f$npvcost['fossil',,], x=c.f.grid, y=k.g.grid,
      main='Green advantage (difference in $ NPV cost)',
      xlab='Fossil Operating Cost', ylab='Green Capital Cost',
      col=rev(hcl.colors(200, palette='Red-Green')))

image(V.f$opt.strat, x=c.f.grid, y=k.g.grid,  main='Optimal strategy', col=hcl.colors(200, palette='Green-Orange'), xlab='Fossil Operating Cost', ylab='Green Capital Cost')
points(x=c.f.0, y=k.g.0, pch=19, cex=3, col=4)

green.value = V.out$Vnew - V.f$Vnew
fossil.value = V.out$Vnew - V.g$Vnew
delta_value = green.value - fossil.value # negative means green is better

delta_breaks = c(seq(min(delta_value),0,length = n_colors/2), 0, seq(0,max(delta_value),length=n_colors/2))

png(filename = paste0('output/figures/option_value_L',L,'_',Sys.Date(),'.png'), width = 480*3, height=480, pointsize = 24)
# png(filename = paste0('output/figures/option_value_L',L,'_drift_',Sys.Date(),'.png'), width = 480*3, height=480, pointsize = 24)
par(mfrow=c(1,3))
image(z=fossil.value, x=c.f.grid, y=k.g.grid, main='Fossil option value', col=hcl.colors(n_colors, palette='BluGrn'), xlab='Fossil Operating Cost', ylab='Green Capital Cost')
image(z=green.value, x=c.f.grid, y=k.g.grid, main='Green option value', col=hcl.colors(n_colors, palette='BluGrn'), xlab='Fossil Operating Cost', ylab='Green Capital Cost')
image(z=delta_value, x=c.f.grid, y=k.g.grid, main='Green option value over fossil', col=hcl.colors(n_colors, palette='Green-Orange'), xlab='Fossil Operating Cost', ylab='Green Capital Cost', breaks = delta_breaks)
dev.off()

green.value[as.character((c.f.0-2):(c.f.0+2)), as.character(seq(k.g.0-20,k.g.0+20, by=10))]
fossil.value[as.character((c.f.0-2):(c.f.0+2)), as.character(seq(k.g.0-20,k.g.0+20, by=10))]
