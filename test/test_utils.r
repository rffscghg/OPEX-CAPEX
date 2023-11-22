# Test functions in utils.r

message("\nBeginning utils.r accuracy tests.")

# Test that VFI data-tidying function works for const_scrap == FALSE

legacy_0 <- tidy_V(test_L10) %>%
    distinct(legacy_str, legacy_state) %>%
    unlist()

if (isTRUE(legacy_0 != c("", "1"))) stop("VFI data-tidying function failed a test.")

# Test that Monte Carlo data-tidying function works for option == "f" and const_scrap == FALSE

tidy_mc_f <- tidy_mc(test_mc_f)

if (
    isFALSE(tidy_mc_f$pick_f == 1) | 
    isFALSE(is.na(tidy_mc_f$V_g)) |
    isFALSE(tidy_mc_f$legacy_state == 1)
) stop("Monte Carlo data-tidying function failed a test.")

message("Accuracy tests passed.")
