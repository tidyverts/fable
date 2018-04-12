#' @inherit tsibble::guess_frequency
#' @param x A tsibble
#' @importFrom tsibble guess_frequency index
#' @importFrom dplyr pull
#' @examples
#' library(tsibble)
#' guess_frequency(pedestrian) 
#' @export
guess_frequency.tbl_ts <- function(x){
  x %>%
    pull(!!index(.)) %>%
    guess_frequency
}

#' Extract frequencies for common seasonal periods
#' 
#' @param x An object containing temporal data (such as a tsibble, interval, datetime and others.)
#' 
#' @return A named vector of frequencies appropriate for the provided data.
#' 
#' @references <https://robjhyndman.com/hyndsight/seasonal-periods/>
#' 
#' @examples 
#' common_periods(tsibble::pedestrian)
#' 
#' @export
common_periods <- function(x){
  UseMethod("common_periods")
}

#' @importFrom tsibble pull_interval
#' @export
common_periods.default <- function(x){
  common_periods(pull_interval(x))
}

#' @export
common_periods.tbl_ts <- function(x){
  common_periods(tsibble::interval(x))
}

#' @export
common_periods.interval <- function(x){
  freq_sec <- c(year = 31557600, week = 604800, day = 86400, hour = 3600, minute = 60, second = 1)
  switch(paste(names(x), collapse = ""),
         "unit" = c("none" = 1),
         "year" = c("year" = 1),
         "quarter" = c("year" = 4/x[["quarter"]]),
         "month" = c("year" = 12/x[["month"]]),
         "week" = c("year" = 52/x[["week"]]),
         "day" = c("year" = 365.25, "week" = 7)/x[["day"]],
         with(list(secs = freq_sec/sum(as.numeric(x)*freq_sec[names(x)])), {
           if(any(is.na(secs))){
             abort("Irregular time series provided")
           }
           secs[secs>1]
         })
  )
}

get_frequencies <- function(period, data, ...){
  UseMethod("get_frequencies")
}

get_frequencies.numeric <- function(period, data, ...){
  period
}

get_frequencies.character <- function(period, data, ...){
  frequencies <- common_periods(data)
  if(period == "all"){
    return(frequencies)
  }
  else if(period == "smallest"){
    return(frequencies[which.min(frequencies)])
  }
  else if(period == "largest"){
    return(frequencies[which.max(frequencies)])
  }
  else{
    out <- frequencies[period]
    if(any(is.na(out))){
      bad_freq <- period[which(is.na(out))]
      warn(sprintf(
        "Could not find %s for `c(%s)`, possible frequencies for this data are: `c(%s)`.\nUnknown frequencies have been ignored.",
        ifelse(length(bad_freq)==1, "an appropriate frequency", "appropriate frequencies"),
        paste0(paste0(paste0('"', bad_freq, '"'), collapse = ", ")),
        paste0(paste0(paste0('"', names(frequencies), '"'), collapse = ", "))
        )
      )
    }
    return(out[!is.na(out)])
  }
}