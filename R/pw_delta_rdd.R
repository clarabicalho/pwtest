#' Calculates prognosis-weighted delta for RD designs
#' @param data data.frame containing covariates, treatment assignment, and outcome variable
#' @param covariates character vector of covariate names
#' @param running_var character string with running variable name
#' @param treatment name of variable indicating (binary) treatment assigned
#' @param outcome name of outcome variable
#' @param standardize logical. Whether to standardize data inside function
#' @param simulation logical. Whether running the function on bootstrap sample
#' @param rd_estimator character. Whether to use the conventional ("h") or the bias-corrected local-polynomial point estimator ("b"). See `rdrobust()` for more details. Defaults to conventional estimate ("h").
#' @param ... arguments passed on to `rdrobust` function. If `rdrobust()` arguments `y` and `covs` are not specified, they will take the values of the variables defined by `outcome` and `covariates`, respectively. All other arguments, if not specified, will take the default values in `rdrobust()`

pw_delta_rdd <- function(data, covariates, running_var, treatment, outcome,
                         standardize = TRUE, simulation = FALSE,
                         rd_estimator = "h", ...){

  argg <- as.list(match.call())
  c <- NULL

  if("c" %in% names(argg)){
    warning("Arguments `treatment` and `c` both specified, will use `treatment` var to define treatment condition, but `c` will be passed onto `rdrobust()`. Please ensure the values coded in `treatment` are consistent with value of `c`.")
  }

  # obtain prognostic weights/coefficients for each covariate calculated for
  # all control units in the full data
  pw_full <- pw_rdd(data = data, covariates = covariates, treatment = treatment,
                outcome = outcome, standardize = standardize, simulation = simulation)

  # if(any(is.na(pw_full))) stop()
  # standardize data relative to entire study group (the finite population)
  # uses same standardization procedure for data as in non-RD case
  # REVIEW: does not standardize running variable
  if(standardize){
    if(simulation & length(unique(data[[outcome]]))==1L){
      data <- data %>%
        dplyr::mutate_at(.vars = c(covariates),
                         .funs = stdr) %>% as.data.frame()
    } else {
      data <- data %>%
        dplyr::mutate_at(.vars = c(outcome, covariates),
                         .funs = stdr) %>% as.data.frame()
    }
  }

  # fitted values of Y0 with prognosis weights
  # (estimated coefs from control group regression of Y0 on covariates)
  Y0hat <- as.matrix(data[,covariates])%*%as.matrix(pw_full, nrow = length(pw_full))

  # code values for rdrobust arguments
  if(!"y" %in% names(argg)) argg$y <- Y0hat
  if(!"x" %in% names(argg)) argg$x <- data[[running_var]]

  # rdrobust() inherits arguments from pwtest()
  rd_argg <- intersect(names(argg), names(formals(rdrobust)))
  # REVIEW: rdrobust does not take variables in `covariates` for the argument `covs`
  # unless `covs` is specified (separately)
  rd_out <- do.call("rdrobust", args = argg[rd_argg])

  # use the rdrobust output to extract bandwidth and
  # recalculate weights within the bandwidth
  # w/in bw: (re)estimate weights, (re)estimate Y0hat, estimate dii
  argg$h <- rd_out$bws[1,1] # optimal bandwidth
  cutoff <- ifelse("c" %in% names(argg), argg$c, 0)

  data_bw <- subset(data, data[[running_var]] >= cutoff - argg$h & data[[running_var]] <= cutoff + argg$h)

  pw_bw <- pw_rdd(data = data_bw, covariates = covariates, treatment = treatment,
              outcome = outcome, standardize = standardize, simulation = simulation)

  # fitted values of Y0 with prognosis weights within the bandwidth
  # (estimated coefs from control group regression of Y0 on covariates)
  Y0hat <- as.matrix(data[,covariates])%*%as.matrix(pw_bw, nrow = length(pw_bw))
  argg$y <- Y0hat # overwrite outcome with within-bw prognostic weights
  rd_argg <- intersect(names(argg), names(formals(rdrobust)))

  # run dii estimation on reweighted (within-bandwidth) fitted Y0hat
  rd_out <- do.call("rdrobust", args = argg[rd_argg])

  # pw delta as difference in intercepts for Y0hat
  if(rd_estimator == "h") pwdelta <- unname(rd_out$Estimate[,"tau.us"])
  if(rd_estimator == "b") pwdelta <- unname(rd_out$Estimate[,"tau.bc"])

  # UW delta estimates using the optimal or user-set bandwidth
  if(!rd_estimator %in% names(argg)){
    # Note: if not specified, the conventional or bias-corrected bandwidth passed onto covariate-by-covariate difference in intercepts is taken from the same data used to calculate the prognosis weighted delta
    argg[[rd_estimator]] <- unname(rd_out$bws[rd_estimator,])
  }

  # Note: calculates the difference in intercepts for covariates using the same bandwidth set by user or defaulted in rdrobust() with the fitted Y0hat.
  # All other values are rdrobust defaults if not set by user.
  dii_covs <- sapply(covariates, function(covariate){
    argg_cov <- argg
    argg_cov$y <- data[[covariate]]
    argg_new <- intersect(names(argg), names(formals(rdrobust)))
    rd <- do.call("rdrobust", args = argg_cov[argg_new])
    if(rd_estimator == "h") uwd <- unname(rd$Estimate[,"tau.us"])
    if(rd_estimator == "b") uwd <- unname(rd$Estimate[,"tau.bc"])
    return(uwd)
  })

  # R-squared from prognosis regression
  prog_mod_f <- paste(c(outcome, paste(covariates, collapse = " + ")), collapse = " ~ ")
  prog_mod <- lm(formula = prog_mod_f, data = data_bw[data_bw[[substitute(treatment)]] == 0, ])
  prog_Rsq <- summary(prog_mod)$r.squared

  # R-squared from balance regression
  bal_mod_f <- paste(c(treatment, paste(covariates, collapse = " + ")), collapse = " ~ ")
  bal_mod <- lm(formula = bal_mod_f, data = data_bw)
  bal_Rsq <- summary(bal_mod)$r.squared

  return(list(dii = dii_covs,
              uwdelta = sum(dii_covs),
              pw = pw_bw,
              pwdelta = pwdelta,
              pwdelta_se = rd_out$se,
              rdrobust_output = rd_out,
              prog_Rsq = prog_Rsq,
              bal_Rsq = bal_Rsq))

}
