# pwtest

Generate unweighted and prognosis-weighted statistics for tests of as-if random and tests of continuity proposed by [Bicalho, Bouyamourn, and Dunning (2025)](https://arxiv.org/abs/2205.10478) and additional diagnostic tools.

## Overview

Standard covariate-by-covariate tests in natural experiments can lead to false rejections or failure to reject tests of as-if random and tests of continuity in natural experiments. That is because covariates that are uninformative about potential outcomes provide little information about whether potential outcomes are balanced in expectation (in tests of as-if random), or whether the potential outcomes are continuous at the cuttoff in regression discontinuity designs (tests of continuity). Building on prior research on the prognostic or informative value of covariates, we develop a global test of as-if random (and a test of continuity) that upweights covariates associated with the outcome.

## Installation

You can install the development version of **pwtest** from GitHub with:

```r
# install.packages("remotes")
remotes::install_github("clarabicalho/pwtest")
```
