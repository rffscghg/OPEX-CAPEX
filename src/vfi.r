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

    # Calculate the number of possible tuples of legacy assets, i.e., those already in operation
    n_states <- ifelse(const_scrap, 2^(t-1), 1)
    
    # Initialize value function array
    V <- array(
        runif(length(c_f_vals)*length(k_g_vals)), 
        c(
            length(c_f_vals), 
            length(k_g_vals),
            n_states
        ),
        list(
            c_f = c_f_vals,
            k_g = k_g_vals,
            legacy_state = 1:n_states
        )
    )

    # Compute single-timestep operating expenses (accounting for output, drift, and discounting)
    opex_f <- c_f_vals * exp(mu_cf) * q / (1 + r)
    opex_g <- c_g * q / (1 + r)

    # Compute Brownian motion density matrix and total expenses
    if (const_scrap) {

        phi <- phi(c_f_vals, k_g_vals, mu_cf, mu_kg, sigma_cf, sigma_kg, t = 1)
        sum_f_vals <- k_f + opex_f
        sum_g_vals <- k_g_vals + opex_g

    } else {

        phi <- phi(c_f_vals, k_g_vals, mu_cf, mu_kg, sigma_cf, sigma_kg, t)
        sum_f_vals <- k_f + rowSums(tcrossprod(c_f_vals, exp(mu_cf * (1:t))*q*(1+r)^-(1:t)))
        sum_g_vals <- k_g_vals + sum(c_g*q*(1+r)^-(1:t))

    }

    # Create a shorthand version of the value function
    value_V <- function(V, option = "all") {
        .value(opex_f, sum_f_vals, opex_g, sum_g_vals, t, n_states, r, V, phi, option, const_scrap)
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
    message(iter, " iterations yielded a fit to a precision of ", delta, " in ", toString(as_hms(round(t_run,1))), " (hh:mm:ss)")

    # Calculate single-option value functions
    V_f = if (option == "all") {value_V(V, "f")} else {NA} # Otherwise V_min = V_f
    V_g = if (option == "all") {value_V(V, "g")} else {NA} # Otherwise V_min = V_g
    
    # Return solved value function
    return(list(
        V_min = V,
        V_f = V_f,
        V_g = V_g
    ))

}

# Value function
.value <- function(
    opex_f,                     # Vector of discounted fossil-fuel single-period operating expense over c_f range
    sum_f_vals,                 # Vector of fossil-fuel total costs (sans replacement) over c_f range
    opex_g,                     # Discounted green single-period operating expense
    sum_g_vals,                 # Vector of green total costs (sans replacement) over k_g range
    t,                          # Number of timesteps
    n_states,                   # Number of possible tuples of legacy assets
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
            for (k in 1:n_states) { # The binary representation of k-1 can be thought of as the series of
                                    # assets, with the newest asset as the rightmost digit.
                                    # For example: t = 3 -> n_states = 4, k = 4 -> k - 1 = 3 = 11(base 2).
                                    # The modular arithmetic below transforms k - 1 by shifting all assets
                                    # left, adding an asset (1 for fossil), and dropping the leftmost asset.
                
                # For constant scrappage: look one timestep ahead for c_f, k_g, and A (set of t assets)
                # For single asset: look t timesteps ahead for c_f and k_g
                V_right_f <- sum(phi[,,i,j]*V[,,(2*k - 1)%%n_states + 1])*(1+r)^-ifelse(const_scrap,1,t)
                V_right_g <- sum(phi[,,i,j]*V[,,(2*k - 2)%%n_states + 1])*(1+r)^-ifelse(const_scrap,1,t)
                
                N_f <- binary_digit_sum(k - 1) # t=1 -> k=1 -> {N_f, N_g}={0,0}, i.e., no legacy assets
                N_g <- t - N_f - 1
                legacy <- N_f * opex_f[i] + N_g * opex_g

                V_f[i,j,k]  <- sum_f_vals[i] + V_right_f + legacy*const_scrap
                V_g[i,j,k]  <- sum_g_vals[j] + V_right_g + legacy*const_scrap

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

# Calculate the sum of the digits of x as a binary number
binary_digit_sum <- function(x) {

    x_2 <- x
    while (x_2 > 0) {
        x_2 <- floor(x_2 / 2)
        x <- x - x_2
    }

    return(x)

}
