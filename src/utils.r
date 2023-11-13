# Export "tidy" data from VFI output
tidy_V <- function(vfi_obj) {

    # Get c_f_vals from first dimension
    c_f_vals <- as.numeric(dimnames(vfi_obj[[1]])$c_f)

    # Calculate t
    t <- log(dim(vfi_obj[[1]])[3], base = 2) + 1

    vfi_tidy <- vfi_obj %>%
        lapply(as_tibble) %>%
        lapply(mutate, c_f = c_f_vals) %>%
        lapply(pivot_longer, -c_f) %>%
        lapply(separate, name, c("k_g", "legacy_state"), sep = "\\.") %>%
        lapply(
            mutate, 
            k_g = as.numeric(k_g), 
            legacy_state = as.numeric(legacy_state),
            legacy_str = bin2string(legacy_state, t)
        )

    for (i in 1:length(vfi_tidy)) {

        vfi_tidy[[i]]$f_exposure <- 0

        for (j in 1:nrow(vfi_tidy[[i]])) {
            
            # I tried to do this using mutate but ran into vectorization issues
            vfi_tidy[[i]]$f_exposure[j] <- cumul_years_left(vfi_tidy[[i]]$legacy_str[j])

        }

        vfi_tidy[[i]]$option <- names(vfi_tidy)[i]
    }

    return(bind_rows(vfi_tidy))

}

# Calculate the number of years remaining for fossil assets from a legacy_state
cumul_years_left = function(x) {

  return(sum(1:nchar(x) * (strsplit(x, split = NULL)[[1]]=='f')))

}

# Calculate the sum of the digits of x as a binary number
binary_digit_sum <- function(x) {

    x_2 <- x
    while (x_2 > 0) {
        x_2 <- floor(x_2 / 2)
        x <- x - x_2
    }

    return(x)

}

# Convert a legacy_state argument (and the number of timesteps) to an "fg" string for legacy assets
bin2string = function(x, t) {

    x = x - 1
    x = sapply(x, function(z) paste(rev(as.integer(intToBits(z))), collapse=""))
    x = gsub('1', 'f', x)
    x = gsub('0', 'g', x)
    x = substr(x, start=nchar(x)-t+1+1, stop=nchar(x))

    return(x)

}

# Convert an "fg" string for legacy assets to the corresponding legacy_state argument
string2bin = function(x) {

    x = gsub('f','1', x)
    x = gsub('g','0', x)
    x = strtoi(x, base=2)
    x = x+1

    return(x)

}
