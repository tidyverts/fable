# fable (development version)

## Bug fixes

* Fixed NNETAR & VAR handling of missing values (#215).
* Fix ETS forecasting with forecast horizons less than the seasonal period (#219).
* Fixed season() special for non-seasonally based time indices (#220)
* Fix issue with simulation forecasting from damped ETS models.

# fable 0.1.1

## Improvements

* Added interpolation method for `MEAN()` model (#203).
* Added rolling mean option for `MEAN()` model (#204).

## Bug fixes

* Corrected forecast standard error for drift models.

# fable 0.1.0

* First release.

## New features

* Support for 9 models and relevant methods: `ARIMA`, `ETS`, `TSLM`, `MEAN`, `RW`, `NAIVE`, `SNAIVE`, `NNETAR`, `VAR`.
