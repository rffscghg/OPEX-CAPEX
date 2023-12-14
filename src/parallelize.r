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

    message("Beginning VFI for:\n", paste(names(params), params, sep = " = ", collapse = "\n"),"\n")

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

    return(V)

}

monte_carlo_npv_stats <- function(
    params,
    V_list,
    t,
    q = 1,
    r = 0.1,
    const_scrap = TRUE,
    skipVFI = TRUE,
    n_mc = 1000,                            # Monte Carlo sample size
    t_mc = 100,                             # Monte Carlo time horizon
    t_nt = 10,                              # "Near-term" time horizon (for NPV realized costs)
    verbose = FALSE
    ) {
    
    mc <- monte_carlo(
        c_f_vals = as.numeric(params["c_f"]) * sort(unique(as.numeric(params["c_f_multiples"]))),
        k_g_vals = as.numeric(params["k_g"]) * sort(unique(as.numeric(params["k_g_multiples"]))),
        start_cf = as.numeric(params["c_f"]) * as.numeric(params["c_f_multiples"]),
        start_kg = as.numeric(params["k_g"]) * as.numeric(params["k_g_multiples"]),
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
        verbose = verbose,             
        V_init = V_list[[as.numeric(params["V_id"])]],               
        skipVFI = skipVFI,             
        n_mc = n_mc,                   
        t_mc = t_mc,                   
        start_assets = params["start_assets"],
        running_in_parallel = TRUE
    )

    rc = mc$realized_costs

    PV.Costs.near = apply(
        rc[1:t_nt,], 
        MARGIN = 2, 
        function(x) sum(x*(1+r)^-(1:t_nt))
    )

    PV.Costs = apply(
        rc, 
        MARGIN = 2, 
        function(x) sum(x*(1+r)^-(1:nrow(rc)))
    )
    
    return(list(
        E.PV.near = mean(PV.Costs.near),
        SD.PV.near = sd(PV.Costs.near),
        E.PV = mean(PV.Costs),
        SD.PV = sd(PV.Costs)
    ))

}
