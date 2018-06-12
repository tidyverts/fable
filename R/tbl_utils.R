big_mark <- function (x, ...){
  mark <- if (identical(getOption("OutDec"), ",")) 
    "."
  else ","
  formatC(x, big.mark = mark, ...)
}

dim_tbl <- function (x){
  dim_x <- dim(x)
  format_dim <- purrr::map_chr(dim_x, big_mark)
  paste(format_dim, collapse = " x ")
}