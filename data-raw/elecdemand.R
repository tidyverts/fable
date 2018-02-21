elecdemand <- fpp2::elecdemand
class(elecdemand) <- c("mts", "ts")
elecdemand <- as_tsibble(elecdemand, gather = FALSE)
use_data(elecdemand, overwrite=TRUE)
