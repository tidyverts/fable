#' @importFrom purrr map_dfr
multi_univariate <- function(data, cl){
  data %>%
    split_by(key_vars(.)) %>%
    map_dfr(function(x){
      # Re-evaluate cl in environment with split data
      eval_tidy(
        get_expr(cl), 
        env = child_env(get_env(cl), !!expr_text(get_expr(cl)$data) := x)
      )
    })
}