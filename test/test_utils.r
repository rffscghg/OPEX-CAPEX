# Test functions in figures.r

message("\nBeginning figures accuracy tests.")

# Test that data-tidying function works for const_scrap == FALSE

legacy_0 <- tidy_V(test_L10) %>%
    distinct(legacy_str, legacy_state) %>%
    unlist()
if (isTRUE(legacy_0 != c("", "1"))) stop("Data-tidying function behaved unexpectedly for `const_scrap == FALSE`")

message("Accuracy tests passed.")
