## Test environments
* local ubuntu 18.04 install, R 3.5.3, R 3.6.0
* ubuntu 14.04 (on travis-ci), R-devel, R 3.6.0, R 3.5.3
* macOS 10.13 (on travis-ci), R 3.6.1
* windows server 2012 R2 (on AppVeyor), R 3.6.1
* win-builder, R-devel, R 3.6.1

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

## Re-submission

* Fixed URL in vignette.
* Updated tests and examples to conditionally depend on feasts.
* Updated description to use '' for package names.
* Improved detail in package description to provide examples of included models.
* Removed author field from docs.
* Added Gabriel Caceres as ctb.
* Added \value to documentation for model functions and methods.
* Added more examples to model method documentation.
* Clarified models with their acryonyms in the Description.

> If there are references describing the (theoretical background of)
methods in your package, please add these in the Description field of
your DESCRIPTION file in the form

There are no published references for the package's methods yet.