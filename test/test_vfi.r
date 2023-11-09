# Compare value function iteration to Brian's original code, using coarser resolution for speed

message("\nBeginning single-asset accuracy tests.")

### L = 10 NPV
test_L10 <- vfi(
    c_f_vals = seq(1, 40, by = 3),
    k_g_vals = seq(100, 800, by = 20),
    k_f = 278,
    c_g = 0,
    sigma_cf = .1118,
    sigma_kg = .05,
    t = 10
)

test_original_L10 <- read.csv("test/original_L10_npv.csv")

max_error_L10 <- max(abs(test_L10$V_min - test_original_L10))

if (max_error_L10 > 1e-6) {
    stop(paste0("Value function iteration yielded a poor fit to previous methods, with an error of ", max_error_L10))
}

### L = 10 green option value over fossil option value with drift
test_L10_option_w_drift <- function(option = "all", V_init = NULL) {
    vfi(
        c_f_vals = seq(1, 40, by = 3),
        k_g_vals = seq(100, 800, by = 20),
        k_f = 278,
        c_g = 0,
        mu_cf = 0.1,
        mu_kg = -0.1,
        sigma_cf = .1118,
        sigma_kg = .05,
        t = 10,
        option = option
    )
}

opt_f <- test_L10_option_w_drift("f")
opt_g <- test_L10_option_w_drift("g")

option_delta <- opt_g$V_min - opt_f$V_min # TODO: check the sign here. Value versus cost can be confusing.

test_original_L10_option_drift <- read.csv("test/original_L10_option_delta_w_drift.csv")

max_error_L10_option_drift <- max(abs(option_delta - test_original_L10_option_drift))

if (max_error_L10_option_drift > 1e-6) {
    stop(paste0("Value function iteration yielded a poor fit to previous methods, with an error of ", max_error_L10_option_drift))
}

message("Accuracy tests passed.")

message("\nBeginning miscellaneous tests")

# Test initial value functionality

test_v_init <- test_L10_option_w_drift("f", V_init = opt_g)

test_original_v_init <- read.csv("test/L10_V_f_with_V_g_as_V_init.csv")

if (isFALSE(test_v_init$V_min[,,1] == test_original_v_init)) {
    stop("V_init functionality is broken.")
}

message("Miscellaneous tests passed.")

# Test runtime

message("\nBeginning single-asset runtime test. The current benchmark is 0.7 seconds.")
vfi(
    c_f_vals = seq(1, 41, by = 2),
    k_g_vals = seq(100, 900, by = 20),
    k_f = 278,
    c_g = 0,
    sigma_cf = .1118,
    sigma_kg = .05,
    t = 10
)

message("\nBeginning constant-scrappage runtime test. The current benchmark is 1.2 seconds.")
vfi(
    c_f_vals = seq(1, 41, by = 5),
    k_g_vals = seq(100, 900, by = 80),
    k_f = 278,
    c_g = 0,
    sigma_cf = .1118,
    sigma_kg = .05,
    t = 3,
    const_scrap = TRUE,
    max_iter = 1000
)
