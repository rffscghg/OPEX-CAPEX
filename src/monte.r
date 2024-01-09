# Monte Carlo simulation
monte_carlo <- function(
    c_f_vals,                               # State space for fossil-fuel operating costs
    k_g_vals,                               # State space for green-energy capital costs
    k_f = 1,                                # Fossil-fuel capital costs
    c_g = 0,                                # Green-energy operating costs
    mu_cf = 0,                              # Drift for fossil-fuel operating costs (geometric Brownian motion)
    mu_kg = 0,                              # Drift for green-energy capital costs (geometric Brownian motion)
    sigma_cf = 1,                           # Volatility for fossil-fuel operating costs (geometric Brownian motion)
    sigma_kg = 1,                           # Volatility for green-energy capital costs (geometric Brownian motion)
    q = 1,                                  # Output
    t = 1,                                  # Number of timesteps
    r = 0.1,                                # Discount rate
    option = "all",                         # Which options to choose from
    const_scrap = FALSE,                    # Constant scrappage (`t` assets held at once, oldest replaced each timestep)
    threshold = 1e-6,                       # Fit threshold (value function iteration)
    max_iter = 100,                         # Maximum number of iterations (value function iteration)
    verbose = FALSE,                        # Print supplementary information to the console
    V_init = NULL,                          # Starting values for value function iteration (`monte_carlo(...)$value_func` format) 
    skipVFI = FALSE,                        # option to skip the VFI step
    deterministic_prices = NULL,            # Deterministic prices (tibble or dataframe with `c_f` and `k_g` columns)
    n_mc = 1000,                            # Monte Carlo sample size
    t_mc = 10*t,                            # Monte Carlo time horizon
    start_cf = median(c_f_vals),            # Starting value for c_f
    start_kg = median(k_g_vals),            # Starting value for k_g
    start_assets = NULL,                    # Starting assets (e.g., "ffg" means two older fossil-fuel assets and one new green asset)
    running_in_parallel = FALSE             # Suppresses a warning that is annoying when iterated at length
    ) {

    # Calculate the number of possible tuples of legacy assets, i.e., those already in operation
    n_states <- ifelse(const_scrap, 2^(t-1), 1)

    # Assign starting state of legacy assets
    if (const_scrap & (is.null(start_assets))) {
        warning(paste0(
            "Starting assets were set to all fossil-fuel by default.",
            " Use `start_assets` to provide a string",
            " of starting assets one character shorter",
            " than the value of the `t` argument."))
        start_state <- 2^(t-1)
    } else if (!is.null(start_assets)) {
        start_state <- string2bin(start_assets)
    } else {
        start_state <- NULL
    }

    if (skipVFI) {
        V = V_init
    } else {
        # Calculate value function
        V <- vfi(
            c_f_vals,
            k_g_vals,
            k_f,
            c_g,
            mu_cf,
            mu_kg,
            sigma_cf,
            sigma_kg,
            q,
            t,
            r,
            option = option,
            const_scrap,
            threshold,
            max_iter,
            verbose,
            V_init
        )
    }
    
    if (is.null(deterministic_prices)) {
        # Sample random walks
        random_cf <- random_walk_gbm(n_mc, mu_cf, sigma_cf, t_mc, start_cf)
        random_kg <- random_walk_gbm(n_mc, mu_kg, sigma_kg, t_mc, start_kg)
    } else {
        # Use deterministic prices
        random_cf <- matrix(deterministic_prices$c_f)
        random_kg <- matrix(deterministic_prices$k_g)
        n_mc <- 1
        t_mc <- nrow(deterministic_prices)
    }

    # Copy random-walk dimensions to other variables
    c_f_index <- random_cf
    k_g_index <- random_cf
    c_f <- random_cf
    k_g <- random_cf
    V_f_mc <- random_cf
    V_g_mc <- random_cf
    decision_mc <- random_cf
    realized_costs <- random_cf
    legacy_state_mc <- random_cf

    t_start <- Sys.time()

    # Run model
    for (i in 1:nrow(random_cf)) {
        for(j in 1:ncol(random_cf)) {
            
            # Initialize legacy assets
            prior_legacy <- ifelse(i == 1, start_state, legacy_state_mc[i-1,j]) # Decisions and realized costs are based
                                                                                # on legacy assets at timestep `i - 1`

            # Go from a floating point random value to the nearest c_f x k_g grid cell
            c_f_index[i,j] <- index_nearest(random_cf[i,j], c_f_vals)
            k_g_index[i,j] <- index_nearest(random_kg[i,j], k_g_vals)
            c_f[i,j] <- c_f_vals[c_f_index[i,j]]
            k_g[i,j] <- k_g_vals[k_g_index[i,j]]
            
            # Tally legacy assets
            N_f <- binary_digit_sum(prior_legacy - 1) # `prior_legacy - 1` is the binary representation of the legacy state
            N_g <- t - 1 - N_f
            legacy_opex <- N_f*c_f[i,j]*q + N_g*c_g*q

            # Retrieve values from V and choose best option
            if (option == "all") {
                V_f_mc[i,j] <- V$V_f[c_f_index[i,j], k_g_index[i,j], prior_legacy]
                V_g_mc[i,j] <- V$V_g[c_f_index[i,j], k_g_index[i,j], prior_legacy]
                decision_mc[i,j] <- V_f_mc[i,j] < V_g_mc[i,j] # TRUE = fossil-fuel, FALSE = green, minimize V
            } else if (option == "f") {
                V_f_mc[i,j] <- V$V_min[c_f_index[i,j], k_g_index[i,j], prior_legacy]
                V_g_mc[i,j] <- NA
                decision_mc [i,j] <- TRUE
            } else if (option == "g") {
                V_f_mc[i,j] <- NA
                V_g_mc[i,j] <- V$V_min[c_f_index[i,j], k_g_index[i,j], prior_legacy]
                decision_mc[i,j] <- FALSE
            } else stop("The option argument must be set to 'all', 'f', or 'g'.")

            # Calculate realized costs and update legacy assets
            if (decision_mc[i,j]) {
                realized_costs[i,j] <- k_f + c_f[i,j]*q + legacy_opex
                legacy_state_mc[i,j] <- (2*prior_legacy - 1)%%n_states + 1 # Add fossil-fuel
            } else {
                realized_costs[i,j] <- k_g[i,j] + c_g*q + legacy_opex
                legacy_state_mc[i,j] <- (2*prior_legacy - 2)%%n_states + 1 # Add green
            }
        }
    }

    # Calculate and display runtime
    t_run <- Sys.time() - t_start
    message(
        "Completed a Monte Carlo simulation for a sample size of ", n_mc,
        " and a time horizon of ", t_mc, 
        " in ", toString(as_hms(round(t_run,1))), 
        " (hh:mm:ss) at ", Sys.time()
    )

    if (((min(c_f_index) == 1) | (max(c_f_index) == length(c_f_vals))) & !running_in_parallel) {
        warning(paste0(
            "The simulation reached the edge of the c_f state space. ",
            "Consider expanding the range to maintain consistent Brownian motion."
        ))
    }

    if (((min(k_g_index) == 1) | (max(k_g_index) == length(k_g_vals))) & !running_in_parallel) {
        warning(paste0(
            "The simulation reached the edge of the k_g state space. ",
            "Consider expanding the range to maintain consistent Brownian motion."
        ))
    }

    return(list(
        V_f = V_f_mc, 
        V_g = V_g_mc,
        c_f = c_f,
        k_g = k_g,
        random_cf = random_cf,
        random_kg = random_kg,
        pick_f = decision_mc, 
        realized_costs = realized_costs,
        legacy_state = legacy_state_mc, 
        value_func = V
    ))

}

# Find the index of the vector element closest to a value
index_nearest <- function(value, vector) {

    which(abs(value - vector) == min(abs(value - vector)))

}

# Make a matrix where each column is a GBM random walk and each row is a timestep
random_walk_gbm <- function(n, mu, sigma, t, x0) {

    noise <- matrix(
        data = rlnorm(n*t, meanlog = mu-1/2*sigma^2, sdlog = sigma),  # Always one timestep and 
        nrow = t,                                                     # independent of x0
        ncol = n
    )

    rand_walk <- x0*apply(noise, 2, cumprod)

    return(rand_walk)

}
