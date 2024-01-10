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
    data = h_power,
    V_funcs = V_funcs[4:6],
    V_func_params = V_func_params[4:6,],
    t = t,
    y_axis_title_k_g = "Wind Power CAPEX",
    y_axis_title_c_f = "Natural Gas Power OPEX",
    y_axis_title_N_f = "# of legacy fossil-fuel plants",
    y_max_k_g = 1e9,
    y_max_c_f = 1.1e8,
    y_max_annual_cost = 1.25e9,
    multiplier = 1e6,
    plot_filename = "figures/temporal_power.png"
)

save_historical_plots(
    data = h_vehic,
    V_funcs = V_funcs[7:9],
    V_func_params = V_func_params[7:9,],
    t = t,
    y_axis_title_k_g = "Electric vehicle CAPEX",
    y_axis_title_c_f = "Gas-powered vehicle OPEX",
    y_axis_title_N_f = "# of legacy gas-powered vehicles",
    y_max_k_g = 1.5e5,
    y_max_c_f = 2500,
    y_max_annual_cost = 1.6e5,
    plot_filename = "figures/temporal_vehicle.png"
)

save_boxplot(
    data = h_power,
    V_funcs = V_funcs[4:6],
    V_func_params = V_func_params[4:6,],
    t = t,
    x_axis_title = "# of fossil-fuel plants in year 1",
    y_max = 1e10,
    multiplier = 1e6,
    plot_filename = "figures/boxplot_power.png"
)

save_boxplot(
    data = h_vehic,
    V_funcs = V_funcs[7:9],
    V_func_params = V_func_params[7:9,],
    t = t,
    x_axis_title = "# of gas-powered vehicles in year 1",
    y_max = 8.5e5,
    plot_filename = "figures/boxplot_vehicle.png"
)

### Individual asset graphs

indiv_tot_f <- tibble(1:10)
indiv_tot_g <- tibble(1:10)

# Loop over all years in historical data
for (i in 1:nrow(h_vehic)) {
    indiv_tot_f[,i] <- h_vehic$c_f[i:(i+9)]
    indiv_tot_g[,i] <- h_vehic$k_g[i:(i+9)]
}

# Add k_f to first year of c_f
indiv_tot_f[1,] <- indiv_tot_f[1,] + V_func_params$k_f[7]

# Replace k_g with c_g in all but first year
indiv_tot_g[2:10,] <- V_func_params$c_g[7] # 7th row has parameters for vehicle example

colnames(indiv_tot_f) <- h_vehic$date
colnames(indiv_tot_g) <- h_vehic$date

tidy_tot_g <- indiv_tot_g %>%
    mutate(year = 1:nrow(indiv_tot_g)) %>%
    pivot_longer(-year, names_to = "start_year", values_to = "tot_g") %>%
    group_by(year) %>%
    mutate(yearly_min_g = min(tot_g, na.rm = TRUE), yearly_max_g = max(tot_g, na.rm = TRUE))

tidy_tot_f <- indiv_tot_f %>%
    mutate(year = 1:nrow(indiv_tot_f)) %>%
    pivot_longer(-year, names_to = "start_year", values_to = "tot_f") %>%
    group_by(year) %>%
    mutate(yearly_min_f = min(tot_f, na.rm = TRUE), yearly_max_f = max(tot_f, na.rm = TRUE))

p_line_indiv <- left_join(tidy_tot_g, tidy_tot_f) %>%
    ggplot() +
    geom_ribbon(aes(x = year, ymin = yearly_min_g, ymax = yearly_max_g), fill = "#50B161", alpha = .1) +
    geom_line(aes(x = year, y = tot_g, group = start_year), col = "#50B161", alpha = .5) +
    geom_ribbon(aes(x = year, ymin = yearly_min_f, ymax = yearly_max_f), fill = "#ff6663", alpha = .1) +
    geom_line(aes(x = year, y = tot_f, group = start_year), col = "#ff6663", alpha = .5) +
    theme_bw() +
    scale_y_continuous(
        limits = c(0, 1.5e5),
        expand = c(0,0),
        labels = scales::label_dollar(scale_cut = scales::cut_short_scale())
    ) +
    scale_x_continuous(breaks = 1:10, minor_breaks = NULL) +
    labs(x = "Year in vehicle lifespan", y = "Total OPEX + CAPEX costs")

ggsave("figures/individual_assets_line.png", p_line_indiv, height = 7, width = 7)

p_box_indiv <- left_join(tidy_tot_g, tidy_tot_f) %>%
    filter(start_year %in% 2010:2014) %>% # Filter out assets with NAs (in future years)
    group_by(start_year) %>%
    summarise(
        `Total lifespan cost of fossil-fuel vehicles` = sum(tot_f), 
        `Total lifespan cost of EVs` = sum(tot_g)
    ) %>%
    pivot_longer(`Total lifespan cost of fossil-fuel vehicles`:`Total lifespan cost of EVs`) %>%
    ggplot(aes(y = value, x = name, color = name, fill = name)) +
    geom_boxplot(alpha = .2, show.legend = FALSE) +
    theme_bw() +
    scale_y_continuous(
    limits = c(0, 1.5e5),
    expand = c(0,0),
    labels = scales::label_dollar(scale_cut = scales::cut_short_scale())
    ) +
    scale_color_manual(values = c("#50B161", "#ff6663"), aesthetics = c("colour", "fill")) +
    labs(x = "", y = "", title = "Total costs of vehicles purchased in 2010-2014 (10-year lifespans)")

ggsave("figures/individual_assets_box.png", p_box_indiv, height = 7, width = 7)
