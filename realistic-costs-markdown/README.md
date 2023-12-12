OPEX-CAPEX model
================
Last updated by Jordan Wingenroth on
12/12/23

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
future “green” CAPEX ($k_g$) and “fossil” OPEX ($c_f$). We will also
want to tune other parameters and ranges for OPEX and CAPEX costs to be
as realistic as possible, including fossil CAPEX ($k_f$) and green OPEX
($c_g$), which are constants in the current version of our model. We are
planning to present results for two scenarios: power plants and cars.

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
and 2050 prices for capital expense per kW of green power generation,
which will inform our $k_g$ variable. We will use the “Henry Hub Natural
Gas Price” to calculate $c_f$, although because 2022 was a price spike,
we will use the average of all of the 2023 monthly data currently
available in the [FRED](https://fred.stlouisfed.org/series/MHHNGSP)
database. We also plan to incorporate fixed O&M expenses in the CAPEX
figures, and some dimensional analysis is required to get to units that
can be meaningfully compared between the two technologies.

#### $\mu$ and $\sigma$

These parameters are unitless and do not require dimensional analysis.
The [NREL ATB](https://atb.nrel.gov/electricity/2023/land-based_wind)
estimates fixed O&M values to decline over time, with pretty high
uncertainty, so it seems wrong to assume the same values for present and
future $k_g$, so I am going to ignore that for this part.

``` r
FRED_HH_values <- c(3.27, 2.38, 2.31, 2.16, 2.15, 2.18, 2.55, 2.58, 2.64, 2.98, 2.71) # Jan to Nov

mu_formula <- function(x, x0, t, t0) log(x/x0)/(t-t0) # Adapted from "Model Documentation.docx"
sigma_formula <- function(x, sd_x, t, t0) sqrt(log((sd_x/x)^2 + 1)/(t-t0))

k_g         = 928                   # 2050 value
k_g_0       = 1519.08               # 2022 value
sd_k_g      = 291.7610619           # StdDv in 2050
t_g         = 2050
t_g_0       = 2022

mu_k_g      = mu_formula(k_g, k_g_0, t_g, t_g_0)
sigma_k_g   = sigma_formula(k_g, sd_k_g, t_g, t_g_0)

c_f         = 3.27788226897         # 2030 value
c_f_0       = mean(FRED_HH_values)  # Average of 2023 monthly values
sd_c_f      = 1.30384615402         # StdDv in 2030
t_f         = 2030          
t_f_0       = 2023

mu_c_f      = mu_formula(c_f, c_f_0, t_f, t_f_0)
sigma_c_f   = sigma_formula(c_f, sd_c_f, t_f, t_f_0)
```

This results in the following estimates:

| Variable | Drift ($\mu$) | Volatility ($\sigma$) |
|:---------|:--------------|:----------------------|
| $k_g$    | -0.0176       | 0.05802               |
| $c_f$    | 0.03659       | 0.14486               |

#### $k$ values

To convert units for $k$ values from \$/kW to investment cost for a
plant, we will make assumptions about capacity factor and the generation
($q$) size of the plant, and we will add in fixed O&M costs. Capacity
factors come from [Waiting for
Clarity](https://www.rff.org/publications/reports/waiting-for-clarity-how-a-price-on-carbon-can-inspire-investment/),
the 2021 report authored by Brian and his colleagues. We’ll run the
Monte Carlo simulation for a wide range of starting values for $k_g$ but
it would be good to have one, probably the central one, match the actual
value from Rhodium. We want to keep everything in the same start year,
so we will adjust that $k_{g,0}$ starting value to 2023 to match the
$c_{f,0}$ year, using the drift calculated above from the Rhodium
values.

We also need our $k_f$ value, which we will get from the [NREL Annual
Technology
Baseline](https://atb.nrel.gov/electricity/2023/fossil_energy_technologies),
specifically using the values for “NG Combined Cycle” plants. They have
two specifications for those plants, “H-Frame” and “F-Frame”, and their
starting year is 2021, so their 2023 projections are slightly uncertain.
I will use estimates in the center of their range for now, for both
CAPEX and fixed O&M costs, as well as for green fixed O&M costs. We can
change this later if need be, but it seems like it will have a
negligible effect.

Reminder: code chunks in this document can use variables from earlier
chunks.

``` r
dr              = 0.1                        # Our model uses a 10% discount rate
T               = 10                         # Our model assumes plants last 10 years
q               = 1                          # Our model assumes generation of 1 million MWh/year

# Convert $/kW to total cost in millions of $ (where q is in millions of MWh/year)
k_unit_conversion <- function(k, onm, CF, dr, T, q) {

    k_plus_onm  = k + sum(onm*(1-dr)^(1:T))  # Add in fixed O&M
    k_millionMW = k_plus_onm * 1000          # $/kW to millions of dollars per million MW
    k_capacity  = k_millionMW/CF             # Need to build 1/CF units to provide 1 unit of power
    k_milMWh_yr = k_capacity / 8760          # $1 per million MW = $1/8760 per million MWh/year

    return(k_milMWh_yr)    
}

k_g_2023        = exp(log(k_g_0) + mu_k_g)   # Estimate 2023 value from 2022 value
k_f_2023        = mean(c(1234.2, 1274.6))

onm_g_2023      = mean(28.8, 29.9)           # In $/kW-year units
onm_f_2023      = mean(c(30.4, 30.9))

CF_g            = .409
CF_f            = .550

k_g_adj         = k_unit_conversion(k_g_2023, onm_g_2023, CF_g, dr, T, q)
k_f_adj         = k_unit_conversion(k_f_2023, onm_f_2023, CF_f, dr, T, q)
```

This lands us at $k_g$ = \$463.71 million and $k_f$ = \$297.65 million.

We address future decline in $k_g$ as part of our model, but we assume a
constant $k_f$, which seems like it might bias results a bit based on
the NREL ATB projections.

#### $c$ values

Converting operating expenses is easier in terms of units, because we
are leaving out variable O&M and assuming that price per mmbtu scales
linearly. We could add in variable O&M but we’d want to do so for both
options to stay consistent, and I wasn’t sure what values to use for
green power. The best source I found was an [LBNL
report](https://eta-publications.lbl.gov/sites/default/files/opex_paper_final.pdf),
in case we do want to add it in.

So all we need to do to get from the Rhodium Henry Hub value to millions
of dollars per million MWh/year is multiply by a [heat
rate](https://www.eia.gov/electricity/annual/html/epa_08_02.html) that
was most recently estimated at 7,596 Btu/kWh, converted to the correct
units: 7.596 mmbtu/MWh. \$/MWh equals millions of dollars per million
MWh. Thus, the 2023 $c_{f,0}$ estimate converts to \$19.273 million per
year.

With zero variable O&M and zero fuel costs, we effectively assume
$c_{g,0}$ = 0.

### Vehicles

## Model runs
