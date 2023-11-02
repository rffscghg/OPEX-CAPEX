# Value function iteration
vfi <- function(
    c_f_vals,                   # State space for fossil-fuel operating costs
    k_g_vals,                   # State space for green-energy capital costs
    k_f = 1,                    # Fossil-fuel capital costs
    c_g = 0,                    # Green-energy operating costs
    mu_cf = 0,                  # Drift for fossil-fuel operating costs (geometric Brownian motion)
    mu_kg = 0,                  # Drift for green-energy capital costs (geometric Brownian motion)
    sigma_cf = 1,               # Volatility for fossil-fuel operating costs (geometric Brownian motion)
    sigma_kg = 1,               # Volatility for green-energy capital costs (geometric Brownian motion)
    q = 1,                      # Output
    t = 1,                      # Number of timesteps
    r = 0.1,                    # Discount rate
    option = "all",             # Which options to choose from
    threshold = 1e-6,           # Fit threshold (value function iteration)
    max_iter = 100,             # Maximum number of iterations (value function iteration)
    verbose = TRUE              # Print supplementary information to the console
    ) {

    phi <- phi(c_f_vals, k_g_vals, mu_cf, mu_kg, sigma_cf, sigma_kg, t)

    value_V <- function(V) value(c_f_vals, k_g_vals, k_f, c_g, mu_cf, mu_kg, sigma_cf, sigma_kg, q, t, r, V, phi, option)

    V <- value_V()

    delta <- 1
    iter <- 0

    t_start <- Sys.time()

    while ((delta > 1e-6) & (iter < max_iter)) {
        V_new <- value_V(V$V_min)
        delta <- max(abs(V_new$V_min - V$V_min))
        V <- V_new
        iter <- iter + 1
        if (verbose) cat("iteration: ",iter,"\n")
    }

    t_run <- Sys.time() - t_start
    
    if (verbose) cat(iter, " iterations yielded a fit to a precision of ", delta," in ", round(t_run), "seconds \n")

    V <- lapply(V, function(x) {
        rownames(x) <- c_f_vals
        colnames(x) <- k_g_vals
        return(x)
    })

    return(V)

}

# Value function
value <- function(
    c_f_vals,                   # State space for fossil-fuel operating costs
    k_g_vals,                   # State space for green-energy capital costs
    k_f = 1,                    # Fossil-fuel capital costs
    c_g = 0,                    # Green-energy operating costs
    mu_cf = 0,                  # Drift for fossil-fuel operating costs (geometric Brownian motion)
    mu_kg = 0,                  # Drift for green-energy capital costs (geometric Brownian motion)
    sigma_cf = 1,               # Volatility for fossil-fuel operating costs (geometric Brownian motion)
    sigma_kg = 1,               # Volatility for green-energy capital costs (geometric Brownian motion)
    q = 1,                      # Output
    t = 1,                      # Number of timesteps
    r = 0.1,                    # Discount rate
    V,                          # Value function matrix (dimensions determined by c_f_vals and k_g_vals)
    phi,                        # Two-dimensional Brownian motion density matrix
    option = "all"              # Which options to choose from
    ) {

    if(missing(V)) {
        V <- matrix(
            runif(length(c_f_vals)*length(k_g_vals)), 
            length(c_f_vals), 
            length(k_g_vals)
        )
    }

    V_min <- V # Copy dimensions, values will be overwritten
    V_f <- V
    V_g <- V

    Vleft_f <- k_f + rowSums(tcrossprod(c_f_vals, exp(mu_cf * (1:t))*q*(1+r)^-(1:t)))

    Vleft_g <- k_g_vals + sum(c_g*q*(1+r)^-(1:t))

    # TODO: replace with array multiplication between phi and V to speed up computation (probably by a lot)
    for (i in 1:length(c_f_vals)) {          
        for (j in 1:length(k_g_vals)) {
            V_right <- sum(phi[,,i,j]*V)*(1+r)^-t
            V_f[i,j]  <- Vleft_f[i] + V_right
            V_g[i,j]  <- Vleft_g[j] + V_right
        }
    }

    if (option == "all") {
        V_min <- pmin(V_f, V_g)
    } else if (option == "f") {
        V_min <- V_f
    } else if (option == "g") {
        V_min <- V_g
    } else stop("Invalid option argument. Choose 'all', 'f', or 'g'.")


    return(list(
        V_min = V_min,
        V_f = V_f, 
        V_g = V_g
    ))

}

# Two-dimensional Brownian motion density matrix
phi <- function(
    c_f_vals,                   # State space for fossil-fuel operating costs
    k_g_vals,                   # State space for green-energy capital costs
    mu_cf = 0,                  # Drift for fossil-fuel operating costs (geometric Brownian motion)
    mu_kg = 0,                  # Drift for green-energy capital costs (geometric Brownian motion)
    sigma_cf = 1,               # Volatility for fossil-fuel operating costs (geometric Brownian motion)
    sigma_kg = 1,               # Volatility for green-energy capital costs (geometric Brownian motion)
    t = 1                       # Number of timesteps
    ) {
    
    phi_array <- array(
        data = NA, # to fill in
        dim = c(
            length(c_f_vals),
            length(k_g_vals),
            length(c_f_vals),
            length(k_g_vals)            
        ),
        dimnames = list(
            c_f = c_f_vals, 
            k_g = k_g_vals, 
            c_f_0 = c_f_vals, 
            k_g_0 = k_g_vals
        )
    )
    
    for (i in 1:length(c_f_vals)) {
        for (j in 1:length(k_g_vals)) {
            phi_c_f <- dgbm(c_f_vals, mu_cf, sigma_cf, t, c_f_vals[i])
            phi_k_g <- dgbm(k_g_vals, mu_kg, sigma_kg, t, k_g_vals[j])
            gbm_2d <- tcrossprod(phi_c_f, phi_k_g)
            phi_array[,,i,j] <- gbm_2d/sum(gbm_2d)
        }
    }

    return(phi_array)

}

# Geometric Brownian motion density function
dgbm <- function(x, mu, sigma, t, x0) {

    dlnorm(x, meanlog = log(x0) + (mu-1/2*sigma^2)*t, sdlog = sigma*sqrt(t))

}
