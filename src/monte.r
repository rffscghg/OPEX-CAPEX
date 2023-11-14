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
    const_scrap = FALSE,                    # Constant scrappage (`t` assets held at once, oldest replaced each timestep)
    threshold = 1e-6,                       # Fit threshold (value function iteration)
    max_iter = 100,                         # Maximum number of iterations (value function iteration)
    verbose = FALSE,                        # Print supplementary information to the console
    V_init = NULL,                          # Starting values for iteration, in the same format as this function's output
    n_mc = 1000,                            # Monte Carlo sample size
    t_mc = 10*t,                            # Monte Carlo time horizon
    start_cf = median(c_f_vals),            # Starting value for c_f
    start_kg = median(k_g_vals),            # Starting value for k_g
    start_assets = NULL                     # Starting assets (e.g., "ffg" means two older fossil-fuel assets and one new green asset)
    ) {

    if (const_scrap & (is.null(start_assets))) {

        warning(paste0(
            "Starting assets set to all fossil-fuel.",
            " Use `start_assets` to provide a string",
            " of starting assets one character shorter",
            " than the number of timesteps."))
        start_state <- 2^(t-1)

    } else if (!is.null(start_assets)) {

        start_state <- string2bin(start_assets)

    } else {

        start_state <- NULL

    }

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
        option = "all", # `monte_carlo()` depends on V_f and V_g output
        const_scrap,
        threshold,
        max_iter,
        verbose,
        V_init
    )

    # Sample random walks
    random_cf <- random_walk_gbm(n_mc, mu_cf, sigma_cf, t_mc, start_cf)
    random_kg <- random_walk_gbm(n_mc, mu_kg, sigma_kg, t_mc, start_kg)

    # Copy random-walk dimensions to output variables
    decision_mc <- random_cf
    V_f_mc <- random_cf
    legacy_state_mc <- random_cf

    return(decision_mc)

}

# Make a matrix where each column is a GBM random walk and each row is a timestep
random_walk_gbm <- function(n, mu, sigma, t, x0) {

    noise <- matrix(
        data = rlnorm(n*t, meanlog = mu-1/2*sigma^2, sdlog = sigma),  # Always one timestep and 
        nrow = t,                                                   # independent of x0
        ncol = n
    )

    rand_walk <- x0*apply(noise, 2, cumprod)

    return(rand_walk)

}
