# Make figures

results <- read_csv("output/tidy-results 2023-12-14 231102.csv") # temporary
load("output/V-funcs 2023-12-14 222631.Rdata") # temporary

### Surface plots ###

surface_data <- results %>%
    filter(
        k_g_multiples >= plot_multiples[1],
        k_g_multiples <= plot_multiples[2],
        c_f_multiples >= plot_multiples[1],
        c_f_multiples <= plot_multiples[2],
    )

surface_scenarios <- surface_data %>%
    group_by(scenario, opt_name, k_g, c_f) %>%
    summarise(zmax = max(SD_PV)) %>% # Decided to set z limits manually
    mutate(
        xtitle = if_else(
            scenario == "vehicle",
            'Green CAPEX<br>     ($)',
            'Green CAPEX<br>     ($M)'
        ),
        ytitle = if_else(
            scenario == "vehicle",
            'Fossil OPEX<br> ($/year)',
            'Fossil OPEX<br> ($M/year)'
        ),
        ztitle = if_else(
            scenario == "vehicle",
            '$',
            '$M'
        ),
        ztop = if_else(
            scenario == "vehicle",
            70000,
            1300
        ),
    )

# SD_PV

for (i in 1:nrow(surface_scenarios)) {

    save_surface_plot(
        coords = results_to_SD_PV_xyz(surface_data, surface_scenarios$scenario[i], surface_scenarios$opt_name[i]),
        title = paste("Std. Dev. of NPV costs", surface_scenarios$scenario[i], surface_scenarios$opt_name[i], sep = ", "),
        scene = list(
            xaxis = list(
                title = surface_scenarios$xtitle[i]
            ),
            yaxis = list(
                title = surface_scenarios$ytitle[i]
            ),
            zaxis = list(
                title = surface_scenarios$ztitle[i],
                range = c(0, surface_scenarios$ztop[i])
            ),
            camera=list(eye=list(x=1.25*-1*1.5, y=1.25*-1*1.5, z=1.25*0.75*1.5))
            # default angles for x, y, and z are 1.25. Multiply by proportions to adjust
        ),
        file = paste0("figures/surfaces/", surface_scenarios$scenario[i], "---", surface_scenarios$opt_name[i],".png")
    )

}

# Delta plots

# Change in SD from option, in dollars
scene_sd_delta = list(xaxis=list(title='Green CAPEX<br>     ($M)'),
                      yaxis=list(title='Fossil OPEX<br> ($M/year)'),
                      zaxis=list(title='$M', range=c(-1050,50)),
                      camera=list(eye=list(x=1.25*-1*1.5, y=1.25*-1*1.5, z=1.25*0.75*1.5))) # default angles for x, y, and z are 1.25. Multiply by proportions to adjust

scene_ev_delta = list(xaxis=list(title='Green CAPEX<br>     ($)'),
                yaxis=list(title='Fossil OPEX<br> ($/year)'),
                zaxis=list(title='$', range=c(-6e4,1e4)),
                camera=list(eye=list(x=1.25*-1*1.5, y=1.25*-1*1.5, z=1.25*0.75*1.5))) # default angles for x, y, and z are 1.25. Multiply by proportions to adjust

delta_surface_scenarios <- expand_grid(surface_scenarios, opt_b = surface_scenarios$opt_name) %>%
    filter(
        opt_name != opt_b, 
        !str_detect(opt_name, "only"),
        !(str_detect(opt_name, "fossil")&str_detect(opt_b, "green")),
        !(str_detect(opt_name, "begin-green")&str_detect(opt_b, "fossil-only"))
    ) %>% 
    distinct_all()

for (i in 1:nrow(delta_surface_scenarios)) {
    if (delta_surface_scenarios$scenario[i] == "vehicle") {scene <- scene_ev_delta} else {scene <- scene_sd_delta}
    save_surface_plot(
        coords = a_minus_b_SD_PV_xyz(
            surface_data, 
            delta_surface_scenarios$scenario[i], 
            delta_surface_scenarios$opt_name[i],
            delta_surface_scenarios$opt_b[i]
        ),
        title = "",
        scene = scene,
        file = paste0(
            "figures/surfaces/", 
            delta_surface_scenarios$scenario[i], 
            "---delta---", 
            delta_surface_scenarios$opt_name[i],
            "---",
            delta_surface_scenarios$opt_b[i],
            ".png"
        ),
        color_scale = list(c(0, 1), c("blue", "#dfd8d4"))
    )
}

### Bar graphs ###

extremes <- results %>%
    filter(scenario == "neutral") %>%
    filter(
        k_g_multiples %in% c(min(plot_multiples), max(plot_multiples)),
        c_f_multiples %in% c(min(plot_multiples), max(plot_multiples)),
    )

central <- results %>%
    filter(scenario == "neutral") %>%
    filter(k_g_multiples == 1, c_f_multiples == 1)

bar_graph <- bind_rows(extremes, central) %>%
    mutate(opt_name = factor(
        opt_name, 
        levels = c("fossil-only", "green-only", "both-begin-fossil", "both-begin-green")
    )) %>%
    select(
        opt_name, 
        k_g, 
        c_f, 
        k_g_multiples, 
        c_f_multiples, 
        SD_PV, 
        SD_PV_near
    ) %>%
    pivot_longer(SD_PV:SD_PV_near) %>%
    mutate(
        name = factor(
            name, 
            levels = c("SD_PV", "SD_PV_near"), 
            labels = c("Long-run", "First 10 years only")
        ),
        CAPEX = fct_reorder(factor(paste0("Green CAPEX = $", k_g_multiples * k_g, "M")), -k_g_multiples),
        OPEX = fct_reorder(factor(paste0("Fossil OPEX = $", c_f_multiples * c_f, "M/yr")), c_f_multiples)
    ) %>%
    filter(name != "Long-run") %>%
    ggplot(aes(x = opt_name, y = value*1e6)) +
    geom_col(position = "dodge", fill = "#04273C", color = "#04273C") +
    facet_grid(CAPEX~OPEX) +
    theme_bw() +
    scale_y_continuous(
        breaks = c(0, 2.5e8, 5e8, 7.5e8, 1e9, 1.25e9), 
        labels = scales::label_dollar(scale_cut = scales::cut_short_scale())
    ) +
    labs(x = "", y = "Standard Deviation of NPV Costs", fill = "") +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1), 
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank()
    )

ggsave("figures/bar_graph.png", bar_graph, width = 7, height = 6)

### Historical data graphs ###

historical <- read_csv("data/historical.csv")

# Plots of k_g and c_f data

k_g_plot <- historical %>%
    ggplot(aes(x = date, y = k_g*1e6)) +
    geom_line(color = "#50B161") +
    geom_point(color = "#50B161") +
    theme_bw() +
    scale_y_continuous(
        limits = c(0, 9e8),
        expand = c(0,0),
        labels = scales::label_dollar(scale_cut = scales::cut_short_scale())
    ) +
    labs(x = "Year", y = "Wind Power CAPEX")

c_f_plot <- historical %>%
    ggplot(aes(x = date, y = c_f*1e6)) +
    geom_line(color = "#ff6663") +
    geom_point(color = "#ff6663") +
    theme_bw() +
    scale_y_continuous(
        limits = c(0, 7.5e7),
        expand = c(0,0),
        labels = scales::label_dollar(scale_cut = scales::cut_short_scale())
    ) +
    labs(x = "Year", y = "Natural Gas Power OPEX")

# Deterministic model run using `monte_carlo()`

historical_mc_params <- V_func_params[c(4,5,6,6),] %>% # Select power-plant data
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

for (i in 1:nrow(historical_mc_params)) {

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
        deterministic_prices = historical,
        option = historical_mc_params[i,"option"],          
        V_init = V_funcs[[as.numeric(historical_mc_params[i,"V_id"])]],               
        start_assets = historical_mc_params[i,"start_assets"],
    )

    historical_N_f[[i]] <- sapply(as.vector(historical_results[[i]]$legacy_state) - 1, binary_digit_sum)
    historical_realized_costs[[i]] <- as.vector(historical_results[[i]]$realized_costs)

}

# N_f figure

tidy_hist_N_f <- bind_cols(historical_N_f)

colnames(tidy_hist_N_f) <- paste0(historical_mc_params$option, "-start-", str_sub(historical_mc_params$start_assets,,1))

N_f_plot <- tidy_hist_N_f %>%
    mutate(date = historical$date) %>%
    pivot_longer(-date) %>%
    ggplot(aes(x = date, y = value, color = name)) +
    geom_line(show.legend = FALSE) +
    geom_point(show.legend = FALSE) +
    theme_bw() +
    scale_y_continuous(breaks = 0:10, minor_breaks = NULL) +
    scale_color_manual(values = c("#755EA6","#74645E","#ff6663","#50B161")) +
    labs(x = "Year", y = "# of legacy fossil-fuel plants")

# Annual costs

tidy_hist_cost <- bind_cols(historical_realized_costs)

colnames(tidy_hist_cost) <- paste0(historical_mc_params$option, "-start-", str_sub(historical_mc_params$start_assets,,1))

annual_cost_plot <- tidy_hist_cost %>%
    mutate(date = historical$date) %>%
    pivot_longer(-date) %>%
    ggplot(aes(x = date, y = value*1e6, color = name)) +
    geom_line(show.legend = FALSE) +
    geom_point(show.legend = FALSE) +
    theme_bw() +
    scale_color_manual(values = c("#755EA6","#74645E","#ff6663","#50B161")) +
    scale_y_continuous(
        limits = c(0, 1e9),
        expand = c(0,0),
        labels = scales::label_dollar(scale_cut = scales::cut_short_scale())
    ) +
    labs(x = "Year", y = "Annual costs", color = "")

# Save plot

ggsave(
    "figures/temporal.png", 
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
