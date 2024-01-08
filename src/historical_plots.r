# Function to save plots based on historical data

save_historical_plots <- function(
    data,                       # Needs to have date, c_f, and k_g columns and no NAs
    V_funcs,                    # Value functions, listed in order: "f" option, "g" option, "all" option
    V_func_params,              # Parameters of value functions, as a 3-row tibble: c_f, k_g, k_f, c_g, mu_f, mu_g, sigma_f, sigma_g, option
    t,                          # Number of assets
    y_axis_title_k_g,
    y_axis_title_c_f,
    y_axis_title_N_f,
    y_max_k_g,
    y_max_c_f,
    y_max_annual_cost,
    multiplier = 1,
    plot_filename
) {

    iterations <- c(1,2,3,3) # Run once for "f" option, once for "g" option, once for "all, start f", and once for "all, start g"

    # Plots of k_g and c_f data

    k_g_plot <- data %>%
        ggplot(aes(x = date, y = k_g*multiplier)) +
        geom_line(color = "#50B161") +
        geom_point(color = "#50B161") +
        theme_bw() +
        scale_y_continuous(
            limits = c(0, y_max_k_g),
            expand = c(0,0),
            labels = scales::label_dollar(scale_cut = scales::cut_short_scale())
        ) +
        labs(x = "Year", y = y_axis_title_k_g)

    c_f_plot <- data %>%
        ggplot(aes(x = date, y = c_f*multiplier)) +
        geom_line(color = "#ff6663") +
        geom_point(color = "#ff6663") +
        theme_bw() +
        scale_y_continuous(
            limits = c(0, y_max_c_f),
            expand = c(0,0),
            labels = scales::label_dollar(scale_cut = scales::cut_short_scale())
        ) +
        labs(x = "Year", y = y_axis_title_c_f)

    # Deterministic model run using `monte_carlo()`

    historical_mc_params <- V_func_params[iterations,] %>% # Select power-plant data
        mutate(
            start_assets = rep(
                c(
                    paste(rep("f",t-1), collapse = ""), 
                    paste(rep("g",t-1), collapse = "")
                ),
                2
            ))

    historical_results <- list()
    historical_N_f <- list()
    historical_realized_costs <- list()

    for (i in 1:length(iterations)) {

        historical_results[[i]] <- monte_carlo(
            c_f_vals = as.numeric(historical_mc_params[i,"c_f"]) * c_f_multiples,
            k_g_vals = as.numeric(historical_mc_params[i,"k_g"]) * k_g_multiples,
            k_f = as.numeric(historical_mc_params[i,"k_f"]),
            c_g = as.numeric(historical_mc_params[i,"c_g"]),
            mu_cf = as.numeric(historical_mc_params[i,"mu_f"]),
            mu_kg = as.numeric(historical_mc_params[i,"mu_g"]),
            sigma_cf = as.numeric(historical_mc_params[i,"sigma_f"]),
            sigma_kg = as.numeric(historical_mc_params[i,"sigma_g"]),
            q = 1,
            r = 0.1,
            t = t,
            const_scrap = TRUE,
            skipVFI = TRUE,
            deterministic_prices = data,
            option = historical_mc_params[i,"option"],          
            V_init = V_funcs[[iterations[i]]],               
            start_assets = historical_mc_params[i,"start_assets"],
        )

        historical_N_f[[i]] <- sapply(as.vector(historical_results[[i]]$legacy_state) - 1, binary_digit_sum)
        historical_realized_costs[[i]] <- as.vector(historical_results[[i]]$realized_costs)

    }

    # N_f figure

    tidy_hist_N_f <- bind_cols(historical_N_f)

    colnames(tidy_hist_N_f) <- paste0(historical_mc_params$option, "-start-", str_sub(historical_mc_params$start_assets,,1))

    N_f_plot <- tidy_hist_N_f %>%
        mutate(date = data$date) %>%
        pivot_longer(-date) %>%
        ggplot(aes(x = date, y = value, color = name)) +
        geom_line(show.legend = FALSE) +
        geom_point(show.legend = FALSE) +
        theme_bw() +
        scale_y_continuous(breaks = 0:10, minor_breaks = NULL) +
        scale_color_manual(values = c("#755EA6","#74645E","#ff6663","#50B161")) +
        labs(x = "Year", y = y_axis_title_N_f)

    # Annual costs

    tidy_hist_cost <- bind_cols(historical_realized_costs)

    colnames(tidy_hist_cost) <- paste0(historical_mc_params$option, "-start-", str_sub(historical_mc_params$start_assets,,1))

    annual_cost_plot <- tidy_hist_cost %>%
        mutate(date = data$date) %>%
        pivot_longer(-date) %>%
        ggplot(aes(x = date, y = value*multiplier, color = name)) +
        geom_line() +
        geom_point() +
        theme_bw() +
        scale_color_manual(values = c("#755EA6","#74645E","#ff6663","#50B161")) +
        scale_y_continuous(
            limits = c(0, y_max_annual_cost),
            expand = c(0,0),
            labels = scales::label_dollar(scale_cut = scales::cut_short_scale())
        ) +
        labs(x = "Year", y = "Annual costs", color = "") +
        theme(legend.position = "bottom")

    # Save plot

    ggsave(
        plot_filename, 
        plot_grid(
            k_g_plot, 
            c_f_plot, 
            N_f_plot,
            annual_cost_plot,
            ncol = 1, 
            align = "hv"), 
        width = 7, 
        height = 9
    )

}
