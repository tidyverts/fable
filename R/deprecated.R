# These functions manage deprecated functions in unreleased packages.
# They are due to be removed when these packages are released.

default_time_units <- function(...){
  if(utils::packageVersion("tsibble") >= "0.8.9.9000") getExportedValue("tsibble", "default_time_units")(...) else tsibble::time_unit(...)
}

units_since <- function(...){
  if(utils::packageVersion("tsibble") >= "0.8.9.9000") as.double(...) else tsibble::units_since(...)
}
