## Submission

This patch release resolves outstanding issues on the CRAN check page introduced
by the dplyr v1.0.0 and tsibble v0.9.0 releases.

## Test environments
* local ubuntu 18.04 install, R 3.6.3
* ubuntu 16.04 (on GitHub actions), R 4.0.0, R 3.6.3, R 3.5.3
* macOS (on GitHub actions), R-devel, R 4.0.0
* windows (on GitHub actions), R 3.6.3
* win-builder, R-devel, R 4.0.0, R 3.6.3

## R CMD check results

0 errors | 0 warnings | 0 notes

## revdep checks

All revdeps have been checked. All changes to worse are due
to upstream changes in tsibble v0.9.0, and have fixed versions ready for 
submission.

## Resubmission

Resubmitting as requested now that distributional and fabletools are on CRAN.