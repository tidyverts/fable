context("setup-data.R")

USAccDeaths_tbl <- as_tsibble(USAccDeaths)
UKLungDeaths <- as_tsibble(cbind(mdeaths, fdeaths), pivot_longer = FALSE)

fc_mean <- function(x){
  if(inherits(x, "distribution")) mean(x) else x
}
