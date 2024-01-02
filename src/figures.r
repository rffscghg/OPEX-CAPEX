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
    ggplot(aes(x = opt_name, fill = name, y = value*1e6)) +
    geom_col(position = "dodge") +
    facet_grid(CAPEX~OPEX) +
    theme_bw() +
    scale_y_continuous(
        breaks = c(0, 2.5e8, 5e8, 7.5e8, 1e9, 1.25e9), 
        labels = scales::label_dollar(scale_cut = scales::cut_short_scale())
    ) +
    scale_fill_manual(values = c("#04273C", "#88C4F4")) +
    labs(x = "", y = "Standard Deviation of NPV Costs", fill = "") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("figures/bar_graph.png", bar_graph, width = 7, height = 6)

### Historical data graphs ###

historical <- read_csv("data/historical.csv")

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

ggsave("figures/temporal.png", plot_grid(k_g_plot, c_f_plot, ncol = 1, labels = c("a", "b")), width = 7, height = 9)

V_f_pp <- V_funcs[[4]]$V_min
V_g_pp <- V_funcs[[5]]$V_min
V_all_pp <- V_funcs[[6]]$V_min

