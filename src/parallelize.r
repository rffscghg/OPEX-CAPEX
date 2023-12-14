parallel_vfi <- function(
    params,                 # All c, k, mu, sigma, and option values
    c_f_multiples,          # State space relative to central value
    k_g_multiples,          # State space relative to central value
    t,
    q = 1,
    r = 0.1,
    const_scrap = TRUE,
    threshold = 1e-3,
    max_iter = 1000,
    verbose = FALSE
) {

    V <- vfi(
        c_f_vals = as.numeric(params["c_f"]) * c_f_multiples,
        k_g_vals = as.numeric(params["k_g"]) * k_g_multiples,
        k_f = as.numeric(params["k_f"]),
        c_g = as.numeric(params["c_g"]),
        mu_cf = as.numeric(params["mu_f"]),
        mu_kg = as.numeric(params["mu_g"]),
        sigma_cf = as.numeric(params["sigma_f"]),
        sigma_kg = as.numeric(params["sigma_g"]),
        q = q,
        t = t,
        r = r,
        option = params["option"],          
        const_scrap = const_scrap,        
        threshold = threshold,
        max_iter = max_iter,
        verbose = verbose,
        V_init = NULL   
    )

    return(list(V$V_min))

}