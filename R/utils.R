flatten_first_args <- function(specials){
  specials %>%
    map(~ if(length(.x) > 1){stop("Only one special of each type is allowed for this model")} else {.x[[1]]}) %>%
    set_names(NULL) %>%
    unlist(recursive = FALSE) %>%
    as.list # If no args are provided, unlist removes list structure
}

#' @importFrom purrr imap reduce
enclass <- function(x, subclass = NULL, ...){
  dots_list(...) %>%
    imap(function(value, name) set_names(list(value), name)) %>%
    reduce(.init = x, # Add attributes (from ...)
           function(x, attr) {
             if (!is.null(attr[[1]])) {
               attr(x, names(attr)) <- attr[[1]]
             }
             x
           }) %>%
    add_class(subclass)
}

add_class <- function(x, new_class){
  `class<-`(x, union(new_class, class(x)))
}

rm_class <- function(x, class){ 
  `class<-`(x, setdiff(class(x), class))
} 

#' @importFrom utils tail
fc_idx <- function(idx, h){
  seq(tail(idx, 1), length.out = h + 1, by = time_unit(idx)) %>% tail(-1)
}

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
        imap(~ assign(.y, .x, envir = env))
  })
}