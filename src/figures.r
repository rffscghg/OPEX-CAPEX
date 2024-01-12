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
        scenario == "neutral"
    )

surface_scenarios <- surface_data %>%
    mutate(opt_name_title = factor(
        opt_name, 
        levels = c("fossil-only", "green-only", "both-begin-fossil", "both-begin-green"),
        labels = c("(a) Fossil-only strategy","(b) Green-only strategy","(c) Optimal strategy, beginning with fossil","(d) Optimal strategy, beginning with green")
    )) %>%
    group_by(scenario, opt_name, opt_name_title, k_g, c_f) %>%
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

color_scales_sd <- list(
    both_begin_f = list(list(0, "#ffd1d0"), list(1, "#e65c59")),
    both_begin_g = list(list(0, "#dbedfc"), list(1, "#7ab0dc")),
    only_f = list(list(0, "#d5d1cf"), list(1, "#685a55")),
    only_g = list(list(0, "#cbe8d0"), list(1, "#489f57"))
)

for (i in 1:nrow(surface_scenarios)) {

    save_surface_plot(
        coords = results_to_SD_PV_xyz(surface_data, surface_scenarios$scenario[i], surface_scenarios$opt_name[i]),
        color_scale = color_scales_sd[[i]],
        title = surface_scenarios$opt_name_title[i],
        scene = list(
            xaxis = list(
                title = surface_scenarios$xtitle[i], 
                showline = TRUE, linewidth = 4, gridcolor = "black"
            ),
            yaxis = list(
                title = surface_scenarios$ytitle[i], 
                range = list(3,49), 
                showline = TRUE, linewidth = 4, gridcolor = "black"
            ),
            zaxis = list(
                title = surface_scenarios$ztitle[i],
                range = c(0, surface_scenarios$ztop[i]), 
                showline = TRUE, zerolinewidth = 4, linewidth = 4, gridcolor = "black", tickvals = list(200,400,600,800,1000,1200)
            ),
            camera=list(eye=list(x=1.25*-1*1.5, y=1.25*-1*1.5, z=1.25*0.75*1.5))
            # default angles for x, y, and z are 1.25. Multiply by proportions to adjust
        ),
        file = paste0("figures/surfaces/", surface_scenarios$scenario[i], "---", surface_scenarios$opt_name[i],".png")
    )

}

# Delta plots

color_scales_delta <- list(
    bbf_of = list(list(0, "#bbafd2"), list(1, "#6a5595")),
    bbg_bbf = list(list(0, "#fad1af"), list(1, "#dc9256")),
    bbg_go = list(list(0, "#f5e9b3"), list(1, "#d4be5d"))
)

titles_delta <- c(
    "(e) Reduced uncertainty when green option is<br>available versus unavailable, i.e., (c) - (a)",
    "(g) Reduced uncertainty after having<br>acquired green assets with both<br>options available, i.e., (d) - (c)",
    "(f) Reduced uncertainty when fossil option is<br>available versus unavailable, i.e., (d) - (b)"
)

# Change in SD from option, in dollars
scene_sd_delta = list(xaxis=list(title='Green CAPEX<br>     ($M)', 
                showline = TRUE, linewidth = 4, gridcolor = "black"),
                      yaxis=list(title='Fossil OPEX<br> ($M/year)', 
                showline = TRUE, linewidth = 4, gridcolor = "black", range = list(3,49)),
                      zaxis=list(title='$M', range=c(-1050,50), 
                showline = TRUE, zerolinewidth = 1, linewidth = 4, gridcolor = "black"),
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
        title = titles_delta[i],
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
        color_scale = color_scales_delta[[i]]
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

bar_graph <- central %>%
    mutate(opt_name = factor(
        opt_name, 
        levels = c("fossil-only", "green-only", "both-begin-fossil", "both-begin-green"),
        labels = c("Fossil-only\nstrategy","Green-only\nstrategy","Optimal strategy\nbeginning with\nfossil","Optimal strategy\nbeginning with\ngreen")
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
    filter(name == "Long-run") %>%
    ggplot(aes(x = opt_name, fill = opt_name, y = value*1e6)) +
    geom_col(position = "dodge", color = "black", show.legend = FALSE) +
    theme_bw() +
    scale_y_continuous(
        breaks = c(0, 2.5e8, 5e8, 7.5e8, 1e9, 1.25e9), minor_breaks = NULL,
        labels = scales::label_dollar(scale_cut = scales::cut_short_scale()),
        expand = c(0,0),
        limits = c(0, 8e8)
    ) +
    scale_fill_manual(values = c("#74645E","#50B161","#ff6663","#88c4f4")) +
    labs(x = "", y = "Standard Deviation of NPV Costs", fill = "") +
    theme(
        legend.position = c(0.85, 0.9),
        panel.grid.major.x = element_blank(),
        plot.background = element_rect(fill = "white", color = "white"), 
        axis.line = element_line(),
        legend.background = element_blank()
    )

ggsave("figures/central_bar_graph.svg", bar_graph, width = 5, height = 4)

extremes_bar_graph <- extremes %>%
    mutate(opt_name = factor(
        opt_name, 
        levels = c("fossil-only", "green-only", "both-begin-fossil", "both-begin-green"),
        labels = c("Fossil\nonly","Green\nonly","Optimal,\nstart\nfossil","Optimal,\nstart\ngreen")
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
        )
    ) %>%
    filter(name == "Long-run") %>%
    mutate(
        CAPEX = fct_reorder(factor(paste0(c("Low", "Low", "High", "High")," Green CAPEX ($", k_g_multiples * k_g, "M)")), -k_g_multiples),
        OPEX = fct_reorder(factor(paste0(c("Low", "High")," Fossil OPEX ($", c_f_multiples * c_f, "M/yr)")), c_f_multiples)
    ) %>%
    filter(name == "Long-run") %>%
    ggplot(aes(x = opt_name, fill = opt_name, y = value*1e6)) +
    geom_col(position = "dodge", show.legend = FALSE, color = "black") +
    facet_grid(CAPEX~OPEX) +
    theme_bw() +
    scale_y_continuous(
        breaks = c(0, 2.5e8, 5e8, 7.5e8, 1e9, 1.25e9), minor_breaks = NULL,
        labels = scales::label_dollar(scale_cut = scales::cut_short_scale()),
        expand = c(0,0),
        limits = c(0, 1.3e9)
    ) +
    scale_fill_manual(values = c("#74645E","#50B161","#ff6663","#88c4f4")) +
    labs(x = "", y = "Standard Deviation of NPV Costs", fill = "") +
    theme(
        plot.background = element_rect(fill = "white", color = "white"), 
        panel.grid.major.x = element_blank(),
        axis.line = element_line(),
        strip.background = element_blank(),
        panel.spacing = unit(.25, "inch")
    )

ggsave("figures/extremes_bar_graph.svg", extremes_bar_graph, width = 6, height = 5)

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
        "Natural gas only",
        "Wind only",
        "Optimal strategy\nstarting with\nnatural gas",
        "Optimal strategy\nstarting with\nwind"
    ),
    multiplier = 1e6,
    t = t,
    plot_filename = "figures/temporal_power.svg"
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
        "Gas vehicles only",
        "Electric vehicles only",
        "Optimal strategy\nstarting with\ngas vehicles",
        "Optimal strategy\nstarting with\nelectric vehicles"
    ),
    t = t,
    plot_filename = "figures/temporal_vehicle.svg"
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
    ggplot(aes(x = k_g*1e6, y = c_f*1e6, fill = factor(value, labels = c("Natural gas plant     ", "Wind turbine")))) +
    geom_tile(color = "black") +
    geom_text(aes(x = k_g*1e6, y = c_f*1e6, label = date), data = filter(h_power, date %in% c(2010, 2015, 2020, 2022)), inherit.aes = FALSE, nudge_y = 1.1*c(0,1e6,0,0), nudge_x = 1.1*c(7e7,6e7,-7e7,-7e7), color = "white", size = 2.75) +
    geom_path(aes(x = k_g*1e6, y = c_f*1e6), data = filter(h_power, date >= 2010), inherit.aes = FALSE, color = "white") +
    geom_point(aes(x = k_g*1e6, y = c_f*1e6), data = filter(h_power, date >= 2010), inherit.aes = FALSE, color = "white") +
    theme_bw() +
    scale_y_continuous(expand = c(0,0), labels = scales::label_dollar(scale_cut = scales::cut_short_scale())) +
    scale_x_continuous(expand = c(0,0), labels = scales::label_dollar(scale_cut = scales::cut_short_scale())) +
    scale_fill_manual(values = c( "#74645e","#50B161")) +
    labs(x = "Wind power CAPEX", y = "Natural gas OPEX", fill = "Optimal strategy:") +
    theme(aspect.ratio = 1, legend.position = "bottom")

ggsave("figures/optimal_power.svg", p_opt_power, width = 6, height = 6)

p_opt_vehic <- as.tibble(optimal_g_vehic_all) %>%
    mutate(c_f = rownames(optimal_g_vehic_all)) %>%
    pivot_longer(1:length(k_g_multiples), names_to = "k_g") %>%
    mutate(across(c(c_f, k_g), as.numeric)) %>%
    ggplot(aes(x = k_g, y = c_f, fill = factor(value, labels = c("Gas vehicle     ", "Electric vehicle")))) +
    geom_tile(color = "black") +
    geom_text(aes(x = k_g, y = c_f, label = date), data = filter(h_vehic, date %in% c(2010, 2015, 2020, 2023)), inherit.aes = FALSE, nudge_y = 1.1*c(90,90,0,0), nudge_x = 1.1*c(3000,-3000,-5000,-5000), color = "white", size = 2.75) +
    geom_path(aes(x = k_g, y = c_f), data = h_vehic, inherit.aes = FALSE, color = "white") +
    geom_point(aes(x = k_g, y = c_f), data = h_vehic, inherit.aes = FALSE, color = "white") +
    theme_bw() +
    scale_y_continuous(expand = c(0,0), labels = scales::label_dollar(scale_cut = scales::cut_short_scale())) +
    scale_x_continuous(expand = c(0,0), labels = scales::label_dollar(scale_cut = scales::cut_short_scale())) +
    scale_fill_manual(values = c( "#74645e","#50B161")) +
    labs(x = "Electric vehicle CAPEX", y = "Gas vehicle OPEX", fill = "Optimal strategy:") +
    theme(aspect.ratio = 1, legend.position = "bottom")

ggsave("figures/optimal_vehicle.svg", p_opt_vehic, width = 6, height = 6)
