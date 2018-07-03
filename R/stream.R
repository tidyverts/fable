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
  obj_vars <- key_vars(object)
  newdata <- data %>% 
    group_by(!!!syms(obj_vars)) %>%
    nest %>%
    rename(.newdata = !!sym("data"))
  
  if(length(key_vars(object)) == 0){
    object <- object %>%
      mutate(.newdata = newdata$.newdata)
  }
  else{
    object <- object %>%
      left_join(newdata, by = obj_vars)
  }
  object %>%
    mutate(
      data = map2(data, .newdata, dplyr::bind_rows),
      model = map2(model, .newdata, stream, ...) %>% enclass("lst_mdl")
    ) %>%
    select(exclude(".newdata"))
}