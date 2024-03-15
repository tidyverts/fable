## Submission

Attempted to fix issue regarding defining R_NO_REMAP before R headers from C++ 
as advised by email from CRAN

Note that it is not possible for us to test this, but I think we have done 
what is required.

## Test environments
* local ubuntu 20.04 install, R 4.1.2
* ubuntu-latest (on GitHub actions), R-devel, R-release, R-oldrel
* macOS (on GitHub actions), R-release
* windows (on GitHub actions), R-release
* win-builder, R-devel

## R CMD check results

0 errors | 0 warnings | 0 notes

## Revdep checks

All reverse dependencies have been checked, none have changed to worse.
