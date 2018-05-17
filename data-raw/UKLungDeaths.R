UKLungDeaths <- tsibble::as_tsibble(cbind(mdeaths, fdeaths), gather = FALSE)
usethis::use_data(UKLungDeaths, overwrite=TRUE)
