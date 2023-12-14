# Main script to run OPEX-CAPEX model

library(tidyverse)
library(hms)
library(parallel)
library(plotly)
library(orca)
library(RColorBrewer)

source("src/monte.r")
source("src/vfi.r")
source("src/utils.r")

source("test/test_monte.r")
source("test/test_vfi.r")
source("test/test_utils.r")

scenarios <- read_csv("data/scenarios.csv", col_select = -1) # Calculated in .../README.rmd document

options = tibble(
    opt_name = c("fossil-only", "green-only", "both-begin-fossil", "both-begin-green"),
    option = c("f", "g", "all", "all"),
    start_assets = c("f", "g", "f", "g")
)

grid <- expand_grid(
    scenario = c("neutral", "power-plant", "vehicle"),
    opt_name = c("fossil-only", "green-only", "both-begin-fossil", "both-begin-green")
    ) %>%
    left_join(scenarios) %>%
    left_join(options)
