#' Produces unweighted delta statistic (sum of standardized difference in means for all covariates across treatment and control groups) with standard errors and p-values
#' @param data data.frame containing covariates, treatment assignment, and outcome variable
#' @param covariates character vector of names of placebo variables.
#' @param treatment name of variable indicating (binary) treatment assigned
#' @param outcome name of outcome variable
#' @param n_bootstraps numeric scalar indicating number of bootstraps to use for the resampling-based p-values
#' @param se_type character. Which type of standard error to return.
#' @returns data.frame containing unweighted delta statistic, its analytic standard error and bootstrap p-value as well as the number of bootstraps used.
#' @export

uwtest <- function(data,
                   covariates = c("X1", "X2", "X3"),
                   treatment = "Z",
                   outcome = "Y",
                   n_bootstraps = 500,
                   se_type = c("bootstrap", "analytic")) {

  se_type <- match.arg(se_type)

  # observed statistics-------------------------------------
  # list-wise deletion
  datac <- data[complete.cases(data[,c(outcome,treatment,covariates)]),]

  # standardize data relative to entire study group (the finite population)

  data <- datac %>%
    mutate_at(.vars = c(outcome, covariates),
              .funs = scale) %>% as.data.frame()

  treat_i <- which(data[[treatment]] == 1)
  control_i <- which(data[[treatment]] == 0)

  # check if all values constant in any covariate
  ux_values <- sapply(covariates, function(var) length(unique(data[control_i,var]))==1L)

  if(any(ux_values)) stop("Variable(s) ",
                          paste0(covariates[ux_values], collapse = ", "),
                          " are constant in the non-missing control data.")

  # difference in means -----------------------------------------------------
  DIM <- sapply(covariates, function(x) {
    tryCatch({
      tres <- t.test(data[treat_i, x], data[control_i, x], na.rm = TRUE)
      dim <- tres$estimate[1] - tres$estimate[2]
      return(unname(dim))
    }, error = function(e) {
      message(sprintf("Error in t.test for covariate '%s': %s", x, e$message))
      return(NA)
    })
  })


  # standard errors ---------------------------------------------------------

  # check if any of the DIM are missing or NaN
  cov_nan <- names(DIM)[is.nan(DIM)]
  if(any(is.nan(DIM))) stop(paste0("The following covariates have no observations in either treatment or control conditions or covariate values are constant across both groups (so cannot be standardized): ",
                                   paste(cov_nan, collapse = ", "), ". Consider removing these covariates."))

  #analytic standard error
  Nc <- nrow(datac) # changed from data_uw
  z0c <- datac[[treatment]] == 0
  n1c <- sum(!z0c, na.rm = TRUE)
  n0c <- sum(z0c, na.rm = TRUE)

  # analytic standard error of unweighted delta
  if(length(covariates)>1L){

    varcov <- cov(datac[,covariates], use = "everything")
    sum_sigma2 <- sum(diag(varcov)*(nrow(datac)-1)/nrow(datac))
    sum_covs <- sum(varcov[upper.tri(varcov)]*(nrow(datac)-1)/nrow(datac))
    multiplier <- ((Nc^2)/(Nc-1))/(n1c*n0c)

    # analytic SE
    uwdelta_se <- sqrt(multiplier*sum(sum_sigma2, 2*sum_covs))

  } else { # when only one covariate
    uwdelta_se <- sqrt(var(datac[!z0c,covariates])/n1c +
                         var(datac[z0c,covariates])/n0c)
  }

  # bootstrapping -------------------------------------

  # resample from control group with replacement
  samples <- get_bootstrap_samples(control_i, length(treat_i), n_bootstraps)

  # standardize bootstrap population relative to control group SD
  data_stdc <- std_data(data, c(outcome, covariates), treatment)

  # obtain delta distribution from bootstrap samples
  bstats <- apply(samples, 2, function(z) {
    bsample <- data_stdc[z, ]
    bsample[[treatment]][1:length(treat_i)] <- 1

    treat_i <- which(bsample[[treatment]] == 1)
    control_i <- which(bsample[[treatment]] == 0)

    tryCatch({
      DIM <- sapply(covariates, function(x){
        tres <- t.test(bsample[treat_i, x], bsample[control_i, x], na.rm = TRUE)
        dim <- tres$estimate[1] - tres$estimate[2]
        return(unname(dim))
      })
      return(sum(DIM))
    }, error = function(e) return(NA))
  })

  # check for completeness and draw additional bootstrap samples if necessary
  nc_bootstraps <- sum(!is.na(bstats))

  while (nc_bootstraps < as.integer(n_bootstraps)) {
    add_samples <- get_bootstrap_samples(control_i, length(treat_i), n_bootstraps-nc_bootstraps)

    # obtain delta distribution from bootstrap samples
    add_bstats <- apply(add_samples, 2, function(z) {
      bsample <- data_stdc[z, ]
      bsample[[treatment]][1:length(treat_i)] <- 1

      treat_i <- which(bsample[[treatment]] == 1)
      control_i <- which(bsample[[treatment]] == 0)

      tryCatch({
        DIM <- sapply(covariates, function(x){
          tres <- t.test(bsample[treat_i, x], bsample[control_i, x], na.rm = TRUE)
          dim <- tres$estimate[1] - tres$estimate[2]
          return(unname(dim))
        })
        return(sum(DIM))
      }, error = function(e) return(NA))
    })

    bstats <- c(bstats[!is.na(bstats)], add_bstats[!is.na(add_bstats)])
    nc_bootstraps <- sum(!is.na(bstats))
  }

  # SE bootstrap distribution
  uwdelta_se_b <- sd(bstats, na.rm = TRUE) * (n_bootstraps - 1) / n_bootstraps

  # output ------------------------------------------------------------------

  #pwdelta estimates
  est_data <- data.frame(
    n_bootstraps = n_bootstraps,
    uwdelta = sum(DIM),
    uwdelta_se = switch(se_type,
      "analytic" = uwdelta_se,
      "bootstrap" = uwdelta_se_b),
    uwdelta_p = get_pvalue_uw(bstats, sum(DIM), n_bootstraps)
  )

  return(est_data)

}
