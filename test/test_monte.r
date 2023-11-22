# Test Monte Carlo functionality

test_mc_f <- monte_carlo(
    c_f_vals = seq(10, 100, by = 10),
    k_g_vals = seq(100, 1000, by = 100),
    k_f = 250,
    c_g = 2,
    sigma_cf = .01,
    sigma_kg = .01,
    t = 3,
    option = "f",
    const_scrap = FALSE,
    max_iter = 1000,
    threshold = 1e-2,
    verbose = TRUE,
    V_init = if (exists("test")) test$value_func,
    start_assets = "gg"
)
