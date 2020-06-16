# fable 0.2.1

This release coincides with v0.2.0 of the fabletools package, which contains
some substantial changes to the output of `forecast()` methods.
These changes to fabletools emphasise the distribution in the fable 
object. The most noticeable is a change in column names of the fable, with the
distribution now stored in the column matching the response variable, and the
forecast mean now stored in the `.mean` column. 
For a complete summary of these changes, refer to the fabletools v0.2.0 release
news: https://fabletools.tidyverts.org/news/index.html

## New features

* Added the `THETA()` method.

## Improvements

* Forecasts distributions are now provided by the distributional package. They
  are now more space efficient and allows calculation of distributional 
  statistics including the `mean()`, `median()`, `variance()`, `quantile()`,
  `cdf()`, and  `density()`.
* The uncertainty of the drift parameter in random walk models (`RW()`, 
  `NAIVE()` and `SNAIVE()`) is now included in data generated with `generate()`.
* Added Syntetos-Boylan and Shale-Boylan-Johnston variants of `CROSTON()` method.
* Performance improvements.

## Bug fixes

* Fixed issue with approximation being used when refitting ARIMA models and when
  a specific model is requested.
* Fixed `glance()` for `TSLM()` models when the data contains missing values.
* Fixed typo in `glance()` output of `ETS()` models.

## Breaking changes

* The sample path means are now used instead of analytical means when forecasts 
  are produced from sample paths.

# fable 0.2.0

## Improvements

* Added autoregressive modelling with `AR()`.
* Better handling of rank deficiency in `ARIMA()`.
* Added `generate.ARIMA()` method.
* Added bootstrap forecast paths for `ARIMA()` models.
* `ARIMA()` specials now allow specifying fixed coefficients via the `fixed` argument.
* Documentation improvements.

# fable 0.1.2

## Improvements

* Added `CROSTON()` for Croston's method of intermittent demand forecasting.
* Documentation improvements

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
