elecdemand <- fpp2::elecdemand
class(elecdemand) <- c("mts", "ts")
elecdemand <- tsibble::as_tsibble(elecdemand, gather = FALSE)
usethis::use_data(elecdemand, overwrite=TRUE)
