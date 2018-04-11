#' @importFrom tsibble split_by key unkey
#' @importFrom purrr map_dfr
multi_univariate <- function(data, cl){
  data %>%
    split_by(!!!key(.)) %>%
    map_dfr(function(x){
      # Re-evaluate cl in environment with split data
      eval_tidy(
        get_expr(cl), 
        env = child_env(get_env(cl), !!expr_text(get_expr(cl)$data) := x)
      )
    })
}