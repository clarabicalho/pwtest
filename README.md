# Prognosis-Weighted Balance Tests (`pwtest` package)

Generate unweighted and prognosis-weighted statistics for tests of as-if random and tests of continuity proposed by [Bicalho, Bouyamourn, and Dunning (2025)](https://arxiv.org/abs/2205.10478) and additional diagnostic tools.

## Overview

Standard covariate-by-covariate tests in natural experiments can lead to false rejections or failure to reject tests of as-if random and tests of continuity in natural experiments. That is because covariates that are uninformative about potential outcomes provide little information about whether potential outcomes are balanced in expectation (in tests of as-if random), or whether the potential outcomes are continuous at the cuttoff in regression discontinuity designs (tests of continuity). Building on prior research on the prognostic or informative value of covariates, we develop a global test of as-if random (and a test of continuity) that upweights covariates associated with the outcome.

## Installation

You can install the development version of **pwtest** from GitHub with:

```r
# install.packages("remotes")
remotes::install_github("clarabicalho/pwtest")
```

## Main functions

* `pwtest()` produces unweighted (optional) and prognosis-weighted statistics with standard errors and p-values for the _test of as-if random_.
*  `pwtest_rdd()` produces unweighted (optional) and prognosis-weighted statistics with standard errors and p-values for the _test of continuity_ (RD designs).
*  `prog_bal()` generates a plot of standardized covariate difference-in-means in the y-axis and prognosis weights as standardized coefficients from regression of outcome on control units. The plot offers a visual diagnostic for covariate-by-covariate standardized difference in means and prognosis weights from standardized coefficients of prognostic regression of $Y^C(0)$ on set of covariates. The black dots show covariates with significant p-values ($\alpha = 0.05$) two-tailed t-tests of difference in means. The red triangle indicates the value of the R-squared from the prognosis regression on the x-axis and the balance regression on the y-axis.

## Usage Example

We offer an example of the code syntax using open-source data available from Caughey and Sekhon (2011). The code produces the results for the Caughey and Sekhon study reported in Figure 1 in Bicalho, Bouyamourn and Dunning (2026). The data were downloaded in .Rdata format from the Harvard Dataverse at [https://dataverse.harvard.edu/file.xhtml?persistentId=doi:10.7910/DVN/8EYYA2/DRWB57](https://dataverse.harvard.edu/file.xhtml?persistentId=doi:10.7910/DVN/8EYYA2/DRWB57). 

```
    library(pwtest)
     
    # load reproduction data
    load(file = "RDReplication.RData")
    
    # filter as per RDrdobReplication.do
    cs_asif <- filter(x, Use == 1 & abs(DifDPct) < .5 & !is.na(DifDPct))
    cs_rd <- filter(x, Use == 1 & !is.na(DifDPct))

    # as-if random test
    delta_asif <- pwtest(data = cs_asif, 
                         covariates = c("DWinPrv", "DPctPrv", "DifDPPrv",
                          "IncDWNOM1", "ElcSwing", "DemInc", "NonDInc",
                          "PrvTrmsD", "PrvTrmsO", "RExpAdv", "DExpAdv",
                          "SoSDem", "GovDem", "VtTotPct"),
                         treatment = "DemWin", 
                         outcome = "DPctNxt", 
                         nsims = 500)

    # as-if random test results
    delta_asif$estimates
    
    # testing continuity of potential outcome in RD design
    delta_rd <- pwtest_rdd(
        data = cs_rd, 
        covariates = c("DWinPrv", "DPctPrv", "DifDPPrv",
                       "IncDWNOM1", "ElcSwing", "DemInc", "NonDInc",
                       "PrvTrmsD", "PrvTrmsO", "RExpAdv", "DExpAdv",
                       "SoSDem", "GovDem", "VtTotPct"),
        treatment = "DemWin",
        running_var = "DifDPct",
        outcome = "DPctNxt", 
        nsims = 500,
        se_type = "bootstrap"
    )

    # continuity test results
    delta_rd

    # example graph
    plot_pbal(delta_asif, label_option = "minmax", 
              show_color_legend = TRUE)
```
