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

    # Initialize value function matrix
    V <- matrix(
        runif(length(c_f_vals)*length(k_g_vals)), 
        length(c_f_vals), 
        length(k_g_vals)
    )

    # Compute Brownian motion density matrix
    phi <- phi(c_f_vals, k_g_vals, mu_cf, mu_kg, sigma_cf, sigma_kg, t)

    # Sum operating and capital expenses
    sum_f_vals <- k_f + rowSums(tcrossprod(c_f_vals, exp(mu_cf * (1:t))*q*(1+r)^-(1:t)))
    sum_g_vals <- k_g_vals + sum(c_g*q*(1+r)^-(1:t))

    # Create a shorthand version of the value function
    value_V <- function(V, option = "all") value(sum_f_vals, sum_g_vals, t, r, V, phi, option)

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
    }

    # Calculate and display runtime
    t_run <- Sys.time() - t_start
    if (verbose) cat(iter, " iterations yielded a fit to a precision of ", delta," in ", t_run, "seconds \n")

    # Return solved value function matrix
    return(list(
        V_min = V,
        V_f = value_V(V, "f"), # TODO: Add tests for V_f and V_g (comparing to Brian's code)
        V_g = value_V(V, "g")
    ))

}

# Value function
value <- function(
    sum_f_vals,                 # Vector of fossil-fuel total costs (sans replacement) over c_f range
    sum_g_vals,                 # Vector of green total costs (sans replacement) over k_g range
    t = 1,                      # Number of timesteps
    r = 0.1,                    # Discount rate
    V,                          # Value function matrix (dimensions determined by c_f_vals and k_g_vals)
    phi,                        # Two-dimensional Brownian motion density matrix
    option = "all"              # Which options to choose from
    ) {

    V_f <- V
    V_g <- V

    V_array <- array(V, dim(phi))
    V_x_phi <- V_array * phi
    V_replacement <- apply(V_x_phi, 3:4, sum) * (1+r)^-t

    sum_f_mat <- matrix(
        rep(sum_f_vals, length(sum_g_vals)), 
        nrow = length(sum_f_vals)
    )

    sum_g_mat <- matrix(
        rep(sum_g_vals, length(sum_f_vals)), 
        nrow = length(sum_g_vals),
        byrow = TRUE
    )

    sum_min_mat <- pmin(sum_f_mat, sum_g_mat)

    # TODO: replace with array multiplication between phi and V to speed up computation (not sure how much)
    for (i in 1:length(sum_f_vals)) {          
        for (j in 1:length(sum_g_vals)) {
            V_f[i,j]  <- sum_f_vals[i] + V_replacement[i,j]
            V_g[i,j]  <- sum_g_vals[j] + V_replacement[i,j]
        }
    }

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
