#' @inherit tsibble::guess_frequency
#' @param x A tsibble
#' @importFrom tsibble guess_frequency index
#' @importFrom dplyr pull
#' @export
guess_frequency.tbl_ts <- function(x){
  x %>%
    pull(!!index(.)) %>%
    guess_frequency
}

#' Extract potential time-series frequencies
#' 
#' @param x An object containing temporal data (such as a tsibble, interval, datetime and others.)
#' 
#' @return A vector of possible frequencies appropriate for the time intervals provided.
#' 
#' @references <https://robjhyndman.com/hyndsight/seasonal-periods/>
#' 
#' @examples 
#' possible_frequencies(tsibble::pedestrian)
#' 
#' @export
possible_frequencies <- function(x){
  UseMethod("possible_frequencies")
}

#' @export
possible_frequencies.default <- function(x){
  guess_frequency(x)
}

#' @export
possible_frequencies.tbl_ts <- function(x){
  x <- tsibble::interval(x)
  possible_frequencies(x)
}

#' @export
possible_frequencies.interval <- function(x){
  freq_sec <- c(year = 31557600, week = 604800, day = 86400, hour = 3600, minute = 60, second = 1)
  switch(paste(names(x), collapse = ""),
         "unit" = 1,
         "year" = 1,
         "quarter" = 4/x[["quarter"]],
         "month" = 12/x[["month"]],
         "week" = 52/x[["week"]],
         "day" = c(365.25,7)/x[["day"]],
         with(list(secs = freq_sec/sum(as.numeric(x)*freq_sec[names(x)])), {
           if(any(is.na(secs))){
             abort("Irregular time series provided")
           }
           as.numeric(secs)[secs>1]
         })
  )
}