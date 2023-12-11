OPEX-CAPEX model
================
Last updated by Jordan Wingenroth on
12/11/23

## Load VFI and Monte Carlo functions

``` r
library(tidyverse)
library(hms)

source("../src/utils.r")
source("../src/vfi.r")
source("../src/monte.r")
```

## Determine appropriate parameters

We intend to run our OPEX-CAPEX model for realistic scenarios. The most
important parameters we need to determine are $\sigma_f$, $\sigma_g$,
$\mu_f$, and $\mu_g$, which represent uncertanty and drift values for
future “fossil” OPEX and “green” CAPEX. We will also want to tune our
state space for OPEX and CAPEX costs to cover the range of realistic
values for actual investments, including fossil CAPEX and green OPEX,
which are constants in the current version of our model. We are planning
to present results for two scenarios: power plants and cars.

Our main source is the Rhodium Climate Outlook, which was just released
a couple weeks ago (November 30, 2023). We got these values from the
authors:

|                                      | Year | Mean    | StdDv  | 2022 value | Min    | Max     | Unit      |
|--------------------------------------|------|---------|--------|------------|--------|---------|-----------|
| Brent Crude Oil Price                | 2030 | 88.08   | 43.46  | 100.90     | 29.62  | 155.21  | \$/barrel |
| Henry Hub Natural Gas Price          | 2030 | 3.28    | 1.30   | 6.45       | 1.31   | 6.00    | \$/mmbtu  |
| Solar Overnight Capital Cost         | 2050 | 632.00  | 198.71 | 1458.61    | 287.58 | 980.19  | \$/kW     |
| Wind Land Overnight Capital Cost     | 2050 | 928.00  | 291.76 | 1519.08    | 422.28 | 1439.25 | \$/kW     |
| Wind Offshore Overnight Capital Cost | 2050 | 1718.00 | 712.65 | 5029.64    | 667.27 | 2868.02 | \$/kW     |
| EV Batteries                         | 2050 | 112.00  | 86.88  | 240.83     | 16.67  | 234.67  | \$/kWh    |

### Power plants

We’re going to compare onshore wind (green) with natural gas (fossil).
The “Wind Land Overnight Capital Cost” row above provides present-day
and 2050 prices for capital expense per kW of green power generation. We
will use the “Henry Hub Natural Gas Price” for fossil OPEX. We also plan
to incorporate fixed O&M expenses in the CAPEX figures, and some
dimensional analysis is required to get to the correct units. The
$\sigma$ and $\mu$ parameters do not require dimensional, so let’s start
there:

``` r
mu_formula <- function(x, x0, t, t0) log(x/x0)/(t-t0) # Adapted from "Model Documentation.docx"
sigma_formula <- function(x, sd_x, t, t0) sqrt(log((sd_x/x)^2 + 1)/(t-t0))

c_f     = 3.27788226897 # 2030 value
c_f_0   = 6.45          # 2022 value
sd_c_f  = 1.30384615402 # StdDv in 2030
t_f       = 2030          
t_f_0     = 2022

mu_c_f = mu_formula(c_f, c_f_0, t_f, t_f_0)
sigma_c_f = sigma_formula(c_f, sd_c_f, t_f, t_f_0)
```

### Vehicles

## Model runs
