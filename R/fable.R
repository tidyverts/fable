key_vars.fable <- function(x){
  syms(setdiff(colnames(x), c("data", "model", "forecast")))
}