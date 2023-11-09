# Export "tidy" data from VFI output
tidy_V <- function(vfi_obj) {

    # Get c_f vals from first dimension
    c_f_vals <- as.numeric(dimnames(vfi_obj[[1]])$c_f)

    vfi_tidy <- vfi_obj %>%
        lapply(as_tibble) %>%
        lapply(mutate, c_f = c_f_vals) %>%
        lapply(pivot_longer, -c_f) %>%
        lapply(separate, name, c("k_g", "legacy_state"), sep = "\\.") %>%
        lapply(mutate, k_g = as.numeric(k_g), legacy_state = as.numeric(legacy_state))

    for (i in 1:length(vfi_tidy)) {
        vfi_tidy[[i]]$option <- names(vfi_tidy)[i]
    }

    return(bind_rows(vfi_tidy))

}

