flatten_first_args <- function(specials){
  specials %>%
    map(~ if(length(.x) > 1){stop("Only one special of each type is allowed for this model")} else {.x[[1]]}) %>%
    set_names(NULL) %>%
    unlist(recursive = FALSE) %>%
    as.list # If no args are provided, unlist removes list structure
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