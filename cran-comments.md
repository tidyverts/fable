## Test environments
* local ubuntu 18.04 install, R 3.5.3, R 3.6.0
* ubuntu 14.04 (on travis-ci), R-devel, R 3.6.0, R 3.5.3
* macOS 10.13 (on travis-ci), R 3.6.1
* windows server 2012 R2 (on AppVeyor), R 3.6.1

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

## Comments

This is a joint submission with feasts and fabletools. These packages are
interdependent in their examples and functionality, so all packages will need
to be installed prior to checking the packages.

fabletools should be installed first, and then fable and feasts.

Unfortunately due to this dependency, the packages could not be tested using 
win-builder. Instead these packages have been checked on Windows using AppVeyor.

Apologies for the added complexity in reviewing this submission.