#' Extend a fitted model with new data
#'
#' @param object An object (such as a model) which can be extended with additional data
#' @param ... Additional arguments passed on to other methods (usually streamed data and options)
#'
#' @export
stream <- function(object, ...){
  UseMethod("stream")
}

#' @export
stream.mable <- function(object, data, ...){
  newdata <- data %>% 
    group_by(!!!syms(key_vars(object))) %>%
    nest %>%
    rename(.newdata = !!sym("data"))
  
  object %>%
    left_join(newdata) %>%
    mutate(
      data = map2(data, .newdata, bind_rows),
      model = map2(model, .newdata, stream, ...) %>% enclass("lst_mdl")
    ) %>%
    select(exclude(".newdata"))
}