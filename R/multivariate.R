#' @importFrom purrr map_dfr
#' @importFrom dplyr mutate
#' @importFrom tsibble nest
multi_univariate <- function(data, cl){
  data %>% 
    group_by(!!!syms(key_vars(data))) %>%
    nest %>%
    mutate(model = map(data,
                       function(x){
                         # Re-evaluate cl in environment with split data
                         eval_tidy(
                           get_expr(cl), 
                           env = child_env(caller_env(), !!expr_text(get_expr(cl)$data) := x)
                         )
                       })
           ) %>%
    mutate(data = map(!!sym("model"), ~.x$data[[1]]),
           model = map(!!sym("model"), ~.x$model[[1]])) %>%
    enclass("mable")
}