#' Is an object constant?
#'
#' Returns true if the object's numerical values do not vary.
#'
#' @param x object to be tested
is.constant <- function(x) {
  x <- as.numeric(x)
  y <- rep(x[1], length(x))
  return(isTRUE(all.equal(x, y)))
}

assignSpecials <- function(x, env = caller_env()){
  x %>% 
    imap(function(.x, nm){
      if(length(.x) > 1) warn(sprintf("Only one special for `%s` is allowed, defaulting to the first usage", nm))
      .x[[1]] %>% 
        imap(function(.x, .y) assign(.y, .x, envir = env))
  })
}