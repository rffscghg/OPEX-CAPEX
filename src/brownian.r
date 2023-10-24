
# Geometric Brownian motion density function
dgbm <- function(x, mu, sigma, t, x0) {
    dlnorm(x, meanlog = log(x0) + (mu-1/2*sigma^2)*t, sdlog = sigma*sqrt(t))
}

# Two-dimensional Brownian motion density matrix
phi <- function(
    c_f_vals, 
    k_g_vals, 
    mu_cf = 0, 
    mu_kg = 0, 
    sigma_cf = 1, 
    sigma_kg = 1, 
    t = 1
    ) {
    phi_array = array(
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
