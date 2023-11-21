# Export "tidy" data from Monte Carlo output
tidy_mc <- function(mc) {

    # Drop value function (last item)
    mc <- mc[-length(mc)]

    # Get variable names
    mc_names <- names(mc)

    # Create individual tidy tables (tibbles)
    mc <- mc %>%
        lapply(as_tibble) %>%
        lapply(mutate, timestep = row_number()) %>%
        lapply(pivot_longer, -timestep, names_to = "trialnum")

    # Name the columns with the values
    for (i in 1:length(mc)) {
        mc[[i]]$var <- mc_names[i]
        mc[[i]] <- pivot_wider(mc[[i]], names_from = var)
    }

    # Join into one tibble
    mc <- mc %>%
        reduce(inner_join, by = c("timestep", "trialnum")) %>%
        mutate(trialnum = as.integer(str_sub(trialnum,2)))

    return(mc)

}

# Export "tidy" data from VFI output
tidy_V <- function(
    V,                  # Data in the format of `vfi()` output
    index = 1           # Which of the three (V_min, V_f, V_g) output arrays to use
    ) {

    # Get c_f_vals from first dimension
    c_f_vals <- as.numeric(dimnames(V[[index]])$c_f)

    # Calculate t
    t <- log(dim(V[[index]])[3], base = 2) + 1

    V <- V[[index]] %>%
        as_tibble() %>%
        mutate(c_f = c_f_vals) %>%
        pivot_longer(-c_f) %>%
        separate(name, c("k_g", "legacy_state"), sep = "\\.") %>%
        mutate(
            k_g = as.numeric(k_g), 
            legacy_state = as.numeric(legacy_state),
            legacy_str = bin2string(legacy_state, t)
        )

    V$f_exposure <- 0

    for (i in 1:nrow(V)) {
        
        # I tried to do this using mutate but ran into vectorization issues
        V$f_exposure[i] <- cumul_years_left(V$legacy_str[i])

    }

    return(bind_rows(V))

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
