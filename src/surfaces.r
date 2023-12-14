save_surface_plot <- function(
    coords, 
    title, 
    scene, 
    file,
    color_scale = "RdOrYl"

    ) {

    p = layout(
        add_surface(plot_ly(x = coords$x, y = coords$y, z = coords$z), colorscale = color_scale),
        title = list(text = title, y = 0.9),
        scene = scene)

    orca(p, file = file, scale=0.75, width=800*0.8, height=800*0.8)

}

# Only works on current format of `results` in main.jl, not ideal
results_to_SD_PV_xyz <- function(
    results, 
    scenario_name, 
    opt_name_, 
    smooth_w = 3
    ) { 

    tibble <- results %>%
        filter(scenario == scenario_name, opt_name == opt_name_) %>%
        mutate(
            k_g_val = k_g * k_g_multiples,
            c_f_val = c_f * c_f_multiples
        ) %>%
        select(k_g_val, c_f_val, SD_PV) %>%
        pivot_wider(names_from = k_g_val, values_from = SD_PV)

    cfs <- tibble$c_f_val
    kgs <- as.numeric(colnames(tibble)[-1])
    sds <- as.matrix(tibble[,-1])

    return(list(
        x = kgs,
        y = cfs,
        z = smooth_matrix(sds, window_size = smooth_w)
    ))

}

# Function to apply moving average smoothing to a matrix
smooth_matrix <- function(mat, window_size = 3) {
  n <- nrow(mat)
  m <- ncol(mat)
  smoothed_mat <- matrix(NA, nrow = n, ncol = m)
  
  for (i in 1:n) {
    for (j in 1:m) {
      # Extract the window around each element (handling edge cases)
      row_start <- max(1, i - floor(window_size / 2))
      row_end <- min(n, i + floor(window_size / 2))
      col_start <- max(1, j - floor(window_size / 2))
      col_end <- min(m, j + floor(window_size / 2))
      
      # Calculate the average of the window
      smoothed_mat[i, j] <- mean(mat[row_start:row_end, col_start:col_end])
    }
  }
  return(smoothed_mat)
}
