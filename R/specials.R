#' Create evaluation environment for specials
#' 
#' Allows extension packages to make use of the formula parsing of specials.
#' 
#' @param ... A named set of functions which used to parse formula inputs
#' @param parent_env The parent environment of the specials (to find other user objects)
#' @param required_specials The names of specials which must be provided (and if not, are included used with no inputs).
#' 
#' @export
new_specials_env <- function(..., parent_env = caller_env(), required_specials = NULL){
  child_env(parent_env, !!!dots_list(...)) %>%
    enclass(NULL, required_specials = required_specials)
}

tbl_xreg <- function(...){
  tibble(...)
}

exprs_xreg <- function(...){
  exprs(...)
}

no_xreg <- function(...){
  abort("Exogenous regressors are not supported for this model type.")
}

trend <- function(data, knots = NULL){
  idx <- data %>% pull(!!index(data)) %>% as.numeric
  tibble(trend = idx,
         !!!set_names(as.numeric(knots) %>% map(~ pmax(0, idx-.x)), map_chr(knots, ~ paste0("trend_",format(.x)))))
}

season <- function(data, period){
  period <- get_frequencies(period, data)
  idx <- seq_len(NROW(data))
  tibble(!!!set_names(seq_len(period-1) %>% map(~ (idx - .x)%%period == 0), paste0("season_",seq_len(period-1))))
}

fourier <- function(data, period, K){
  idx_num <- data %>% pull(!!index(data)) %>% as.numeric
  period <- get_frequencies(period, data)
  if (length(period) != length(K)) {
    abort("Number of periods does not match number of orders")
  }
  if (any(2 * K > period)) {
    abort("K must be not be greater than period/2")
  }
  
  fourier_exprs <- map2(as.numeric(period), K,
   function(period, K){
     set_names(seq_len(K)/period, paste0(seq_len(K), "_", round(period)))
   }) %>%
    invoke(c, .) %>%
    .[!duplicated(.)] %>%
    imap(function(p, name){
      out <- exprs(C = cospi(2 * !!p * idx_num))
      browser()
      if(abs(2 * p - round(2 * p)) > .Machine$double.eps){
        out <- c(out, exprs(S = sinpi(2 * !!p * idx_num)))
      }
      names(out) <- paste0(names(out), name)
      out
    }) %>%
    set_names(NULL) %>%
    unlist(recursive = FALSE)
  
  tibble(!!!fourier_exprs)
}