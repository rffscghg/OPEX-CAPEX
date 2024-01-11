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

historical <- read_csv("data/tidy/historical.csv")

h_power <- historical %>%
    transmute(date, c_f = c_f_power, k_g = k_g_power) %>%
    filter(!is.na(c_f), !is.na(k_g))

h_vehic <- historical %>%
    transmute(date, c_f = c_f_vehicle, k_g = k_g_vehicle) %>%
    filter(!is.na(c_f), !is.na(k_g))

save_historical_plots(
    data = filter(h_power, date > 2012),
    V_funcs = V_funcs[4:6],
    V_func_params = V_func_params[4:6,],
    y_axis_title_k_g = "Wind Power CAPEX",
    y_axis_title_c_f = "Natural Gas Power OPEX",
    y_axis_title_N_f = "# of legacy fossil-fuel plants",
    y_max_k_g = 1e9,
    y_max_c_f = 6e7,
    y_max_capex_cost = 1e9,
    y_max_opex_cost = 6e8,
    legend_labels = c(
        "Both options,\nstarting with\nnatural gas",
        "Both options,\nstarting with\nwind",
        "Only natural gas",
        "Only wind"
    ),
    multiplier = 1e6,
    t = t,
    plot_filename = "figures/temporal_power.png"
)

save_historical_plots(
    data = filter(h_vehic, date > 2013),
    V_funcs = V_funcs[7:9],
    V_func_params = V_func_params[7:9,],
    y_axis_title_k_g = "Electric vehicle CAPEX",
    y_axis_title_c_f = "Gas-powered vehicle OPEX",
    y_axis_title_N_f = "# of legacy gas vehicles",
    y_max_k_g = 8e4,
    y_max_c_f = 2000,
    y_max_capex_cost = 8e4,
    y_max_opex_cost = 2e4,
    legend_labels = c(
        "Both options,\nstarting with ICEVs",
        "Both options,\nstarting with EVs",
        "Only gas vehicles",
        "Only EVs"
    ),
    t = t,
    plot_filename = "figures/temporal_vehicle.png"
)

### Overlaying historical data and optimal-strategy heatmaps

optimal_g_power <- V_funcs[[6]]$V_g < V_funcs[[6]]$V_f
optimal_g_vehic <- V_funcs[[9]]$V_g < V_funcs[[9]]$V_f

optimal_g_power_all <- apply(optimal_g_power, c(1,2), mean)
optimal_g_vehic_all <- apply(optimal_g_vehic, c(1,2), mean)

p_opt_power <- as.tibble(optimal_g_power_all) %>%
    mutate(c_f = rownames(optimal_g_power_all)) %>%
    pivot_longer(1:length(k_g_multiples), names_to = "k_g") %>%
    mutate(across(c(c_f, k_g), as.numeric)) %>%
    ggplot(aes(x = k_g*1e6, y = c_f*1e6, fill = factor(value, labels = c("Natural gas plant", "Wind plant")))) +
    geom_tile(color = "white") +
    geom_text(aes(x = k_g*1e6, y = c_f*1e6, label = date), data = filter(h_power, date %in% c(2010, 2015, 2020, 2022)), inherit.aes = FALSE, nudge_y = c(0,1e6,0,0), nudge_x = c(7e7,7e7,-7e7,-7e7)) +
    geom_path(aes(x = k_g*1e6, y = c_f*1e6), data = filter(h_power, date >= 2010), inherit.aes = FALSE, color = "#04273C") +
    geom_point(aes(x = k_g*1e6, y = c_f*1e6), data = filter(h_power, date >= 2010), inherit.aes = FALSE, color = "#04273C") +
    theme_bw() +
    scale_y_continuous(expand = c(0,0), labels = scales::label_dollar(scale_cut = scales::cut_short_scale())) +
    scale_x_continuous(expand = c(0,0), labels = scales::label_dollar(scale_cut = scales::cut_short_scale())) +
    scale_fill_manual(values = c( "#ff6663","#50B161")) +
    labs(x = "Wind power CAPEX", y = "Natural gas OPEX", fill = "Optimal strategy") +
    theme(aspect.ratio = 1, legend.position = "bottom")

ggsave("figures/optimal_power.png", p_opt_power, width = 7, height = 7)

p_opt_vehic <- as.tibble(optimal_g_vehic_all) %>%
    mutate(c_f = rownames(optimal_g_vehic_all)) %>%
    pivot_longer(1:length(k_g_multiples), names_to = "k_g") %>%
    mutate(across(c(c_f, k_g), as.numeric)) %>%
    ggplot(aes(x = k_g, y = c_f, fill = factor(value, labels = c("Gas vehicle", "Electric vehicle")))) +
    geom_tile(color = "white") +
    geom_text(aes(x = k_g, y = c_f, label = date), data = filter(h_vehic, date %in% c(2010, 2015, 2020, 2023)), inherit.aes = FALSE, nudge_y = c(90,90,0,0), nudge_x = c(3000,-3000,-6000,-6000)) +
    geom_path(aes(x = k_g, y = c_f), data = h_vehic, inherit.aes = FALSE, color = "#04273C") +
    geom_point(aes(x = k_g, y = c_f), data = h_vehic, inherit.aes = FALSE, color = "#04273C") +
    theme_bw() +
    scale_y_continuous(expand = c(0,0), labels = scales::label_dollar(scale_cut = scales::cut_short_scale())) +
    scale_x_continuous(expand = c(0,0), labels = scales::label_dollar(scale_cut = scales::cut_short_scale())) +
    scale_fill_manual(values = c( "#ff6663","#50B161")) +
    labs(x = "Electric vehicle CAPEX", y = "Gas vehicle OPEX", fill = "Optimal strategy") +
    theme(aspect.ratio = 1, legend.position = "bottom")

ggsave("figures/optimal_vehicle.png", p_opt_vehic, width = 7, height = 7)
