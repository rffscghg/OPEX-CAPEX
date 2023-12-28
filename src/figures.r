# Make figures

plot_data <- results %>%
    filter(
        k_g_multiples >= plot_multiples[1],
        k_g_multiples <= plot_multiples[2],
        c_f_multiples >= plot_multiples[1],
        c_f_multiples <= plot_multiples[2],
    )

plot_scenarios <- plot_data %>%
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

for (i in 1:nrow(plot_scenarios)) {

    save_surface_plot(
        coords = results_to_SD_PV_xyz(plot_data, plot_scenarios$scenario[i], plot_scenarios$opt_name[i]),
        title = paste("Std. Dev. of NPV costs", plot_scenarios$scenario[i], plot_scenarios$opt_name[i], sep = ", "),
        scene = list(
            xaxis = list(
                title = plot_scenarios$xtitle[i]
            ),
            yaxis = list(
                title = plot_scenarios$ytitle[i]
            ),
            zaxis = list(
                title = plot_scenarios$ztitle[i],
                range = c(0, plot_scenarios$ztop[i])
            ),
            camera=list(eye=list(x=1.25*-1*1.5, y=1.25*-1*1.5, z=1.25*0.75*1.5))
            # default angles for x, y, and z are 1.25. Multiply by proportions to adjust
        ),
        file = paste0("figures/", plot_scenarios$scenario[i], "---", plot_scenarios$opt_name[i],".png")
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

delta_plot_scenarios <- expand_grid(plot_scenarios, opt_b = plot_scenarios$opt_name) %>%
    filter(
        opt_name != opt_b, 
        !str_detect(opt_name, "only"),
        !(str_detect(opt_name, "fossil")&str_detect(opt_b, "green")),
        !(str_detect(opt_name, "begin-green")&str_detect(opt_b, "fossil-only"))
    ) %>% 
    distinct_all()

for (i in 1:nrow(delta_plot_scenarios)) {
    if (delta_plot_scenarios$scenario[i] == "vehicle") {scene <- scene_ev_delta} else {scene <- scene_sd_delta}
    save_surface_plot(
        coords = a_minus_b_SD_PV_xyz(
            plot_data, 
            delta_plot_scenarios$scenario[i], 
            delta_plot_scenarios$opt_name[i],
            delta_plot_scenarios$opt_b[i]
        ),
        title = "",
        scene = scene,
        file = paste0(
            "figures/", 
            delta_plot_scenarios$scenario[i], 
            "---delta---", 
            delta_plot_scenarios$opt_name[i],
            "---",
            delta_plot_scenarios$opt_b[i],
            ".png"
        ),
        color_scale = list(c(0, 1), c("blue", "#dfd8d4"))
    )
}
