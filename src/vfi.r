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
    const_scrap = FALSE,        # Constant scrappage (`t` assets held at once, oldest replaced each timestep)
    threshold = 1e-6,           # Fit threshold (value function iteration)
    max_iter = 100,             # Maximum number of iterations (value function iteration)
    verbose = FALSE             # Print supplementary information to the console
    ) {

    # Initialize value function array
    V <- array(
        runif(length(c_f_vals)*length(k_g_vals)), 
        c(
            length(c_f_vals), 
            length(k_g_vals),
            1
        )
    )

    # Compute Brownian motion density matrix
    phi <- phi(c_f_vals, k_g_vals, mu_cf, mu_kg, sigma_cf, sigma_kg, t)

    # Calculate fossil-fuel operating and total expenses
    opex_f <- c_f_vals * exp(mu_cf) * q / (1 + r) # Single timestep
    if (const_scrap) {
        sum_f_vals <- k_f + opex_f
    } else {
        sum_f_vals <- k_f + rowSums(tcrossprod(c_f_vals, exp(mu_cf * (1:t))*q*(1+r)^-(1:t)))
    }

    # Calculate green operating and total expenses
    opex_g <- c_g * q / (1 + r) # Single timestep
    if (const_scrap) {
        sum_g_vals <- k_g_vals + opex_g
    } else {
        sum_g_vals <- k_g_vals + sum(c_g*q*(1+r)^-(1:t))
    }

    # Create a shorthand version of the value function
    value_V <- function(V, option = "all") {
        .value(opex_f, sum_f_vals, opex_g, sum_g_vals, t, r, V, phi, option, const_scrap)
    }

    # Set up for value function iteration
    delta <- 1
    iter <- 0
    t_start <- Sys.time()

    # Carry out value function iteration
    while ((delta > 1e-6) & (iter < max_iter)) {
        V_new <- value_V(V, option)
        delta <- max(abs(V_new - V))
        V <- V_new
        iter <- iter + 1
        if (verbose) cat("iteration", iter, "complete.\n")
    }

    # Calculate and display runtime
    t_run <- Sys.time() - t_start
    cat(iter, "iterations yielded a fit to a precision of", delta, "in", t_run, "seconds\n")

    # Return solved value function matrix
    return(list(
        V_min = as.matrix(V[,,1]),
        V_f = ifelse(option == "all", as.matrix(value_V(V, "f")[,,1]), NA), # Otherwise V_min = V_f
        V_g = ifelse(option == "all", as.matrix(value_V(V, "g")[,,1]), NA)  # Otherwise V_min = V_g
    ))

}

# Value function
.value <- function(
    opex_f,                     # Vector of discounted fossil-fuel single-period operating expense over c_f range
    sum_f_vals,                 # Vector of fossil-fuel total costs (sans replacement) over c_f range
    opex_g,                     # Discounted green single-period operating expense
    sum_g_vals,                 # Vector of green total costs (sans replacement) over k_g range
    t,                          # Number of timesteps
    r,                          # Discount rate
    V,                          # Value function array
    phi,                        # Two-dimensional Brownian motion density matrix
    option,                     # Which options to choose from
    const_scrap                 # Constant scrappage (`t` assets held at once, oldest replaced each timestep)
    ) {
    
    # Match dimensions
    V_f <- V
    V_g <- V
    
    # Iterate value function
    for (i in 1:length(sum_f_vals)) {          
        for (j in 1:length(sum_g_vals)) {
            for (k in 1:1) { # Change to 1:(2^(t-1))
                V_right <- sum(phi[,,i,j]*V[,,1])*(1+r)^-t # Change to V[,,k]
                
                N_f <- binary_digit_sum(k - 1) # So t = 1 -> k = 1 - > {N_f, N_g} = 0, i.e., no legacy assets
                N_g <- t - N_f - 1

                legacy <- N_f * opex_f[i] + N_g * opex_g

                V_f[i,j,1]  <- sum_f_vals[i] + V_right + legacy*const_scrap
                V_g[i,j,1]  <- sum_g_vals[j] + V_right + legacy*const_scrap
            }
        }
    }

    # Choose output
    if (option == "all") {
        V_out <- pmin(V_f, V_g)
    } else if (option == "f") {
        V_out <- V_f
    } else if (option == "g") {
        V_out <- V_g
    } else stop("Invalid option argument. Choose 'all', 'f', or 'g'.")

    return(V_out)

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
    
    # Initialize phi array
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
    
    # Calculate phi likelihoods from Brownian motion density function
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

binary_digit_sum <- function(x) {
    x_2 <- x
    while (x_2 > 0) {
        x_2 <- floor(x_2 / 2)
        x <- x - x_2
    }
    return(x)
}
