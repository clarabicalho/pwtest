# Clean and standardize data before calculating coefficients of prognosis
# regression on control units
# Returns vector of labelled covariate (prognosis) weights

#' @importFrom stats sd
stdr <- function(x){
  (x - mean(x, na.rm = TRUE))/(sd(x, na.rm = TRUE)*(length(x)-1)/length(x))
}

#'Calculates prognosis weights from observed control-group sample
#' @param data data.frame containing covariates, treatment assignment, and outcome variable
#' @param covariates character vector of covariate names
#' @param treatment name of variable indicating (binary) treatment assigned
#' @param outcome name of outcome variable
#' @param standardize whether to standardize data inside function
#' @param simulation logical. Whether running the function on bootstrap sample
#' @importFrom magrittr %>%
#' @importFrom tidyr drop_na
#' @importFrom tidyselect all_of
#' @importFrom dplyr mutate_at filter
#' @importFrom stats coef var lm
pw <- function(data, covariates, treatment, outcome, standardize, simulation){

  # listwise deletion of observations with missing values
  data_pw <- data %>% tidyr::drop_na(tidyselect::all_of(c(outcome, covariates, treatment)))
  z0 <- data_pw[[treatment]] == 0

  # standardize control data for prognosis regression
  if(standardize){
    if(simulation & length(unique(data[[outcome]]))==1L){
      # standardize -covariates only- in simulation runs with fixed POs
      data_c <- data_pw %>% dplyr::filter(z0) %>%
        dplyr::mutate_at(.vars = c(covariates),
                  .funs = stdr) %>% as.data.frame()
    } else {
      data_c <- data_pw %>% dplyr::filter(z0) %>%
        dplyr::mutate_at(.vars = c(outcome, covariates),
                  .funs = stdr) %>% as.data.frame()
    }

  } else {
    data_c <- data_pw %>% dplyr::filter(z0)
  }

  # check if variance = 0 and return error
  check_var <- data_c %>%
    dplyr::select(tidyselect::all_of(c(covariates, outcome))) %>%
    apply(., 2, var, na.rm = TRUE)
  var_na <- names(check_var[is.na(check_var)])
  if(length(var_na)>=1L) stop(paste0("The following variables are constant in the control group among complete cases, so cannot be standardized: ",
                                     paste0(var_na, collapse = ", "),
                                     ". Consider an alternative, for example, excluding the covariate(s)."))
  # calculate prognosis weights
  X <- as.matrix(data_c[,covariates], ncol = length(covariates))
  Y <- as.matrix(data_c[,outcome], ncol = 1)
  pw <- coef(lm(Y ~ X-1))

  names(pw) <- covariates
  return(pw)
}

#' Calculates unweighted delta
#' @param data data.frame containing covariates, treatment assignment, and outcome variable
#' @param covariates character vector of covariate names
#' @param treatment name of variable indicating (binary) treatment assigned
#' @param outcome name of outcome variable
#' @param standardize whether to standardize data inside function
#' @param simulation logical. Whether running the function on bootstrap sample
#' @importFrom magrittr %>%
#' @importFrom tidyr drop_na
#' @importFrom tidyselect all_of
#' @importFrom stats cov var complete.cases
uw_delta <- function(data, covariates, treatment, outcome, standardize = TRUE, simulation = FALSE){

  # drop if missing treatment value
  data_uw <- data %>% tidyr::drop_na(tidyselect::all_of(c(treatment)))
  z0 <- data_uw[[treatment]] == 0

  # standardize data relative to entire study group (the finite population)
  if(standardize){
    if(simulation & length(unique(data[[outcome]]))==1L){
      data_uw <- data_uw %>%
        dplyr::mutate_at(.vars = c(covariates),
                  .funs = stdr) %>% as.data.frame()
    } else {
      data_uw <- data_uw %>%
        dplyr::mutate_at(.vars = c(outcome, covariates),
                  .funs = stdr) %>% as.data.frame()
    }
  }

  # covariate-by-covariate difference in means
  DIM <- sapply(covariates, function(x){
    d <- mean(data_uw[!z0, x], na.rm = TRUE) -
      mean(data_uw[z0, x], na.rm = TRUE)
    return(DIM = d)
  })
  # names(DIM) <- paste0("dim_", names(DIM))
  # unweighted delta
  uwdelta <- sum(DIM)

  # check if any of the DIM are missing or NaN
  cov_nan <- names(DIM)[is.nan(DIM)]
  if(any(is.nan(DIM))) stop(paste0("The following covariates have no observations in either treatment or control conditions or covariate values are constant across both groups (so cannot be standardized): ",
                                   paste(cov_nan, collapse = ", "), ". Consider removing these covariates."))

  # NOTE: list-wise deletion to get analytic standard error
  data_uwc <- data_uw[complete.cases(data_uw[,covariates]),]
  Nc <- nrow(data_uwc) # changed from data_uw
  z0c <- data_uwc[[treatment]] == 0
  n1c <- sum(!z0c, na.rm = TRUE)
  n0c <- sum(z0c, na.rm = TRUE)

  # analytic standard error of unweighted delta
  if(length(covariates)>1L){

    varcov <- cov(data_uwc[,covariates], use = "everything")
    sum_sigma2 <- sum(diag(varcov)*(nrow(data_uwc)-1)/nrow(data_uwc))
    sum_covs <- sum(varcov[upper.tri(varcov)]*(nrow(data_uwc)-1)/nrow(data_uwc))
    multiplier <- ((Nc^2)/(Nc-1))/(n1c*n0c)

    # analytic SE
    uwdelta_se <- sqrt(multiplier*sum(sum_sigma2, 2*sum_covs))

  } else { # when only one covariate
    uwdelta_se <- sqrt(var(data_uwc[!z0c,covariates])/n1c +
                         var(data_uwc[z0c,covariates])/n0c)
  }

  return(list(dim = DIM,
              uwdelta = uwdelta,
              uwdelta_se = uwdelta_se))
}

#' Calculates prognosis-weighted delta
#' @param data data.frame containing covariates, treatment assignment, and outcome variable
#' @param covariates character vector of covariate names
#' @param treatment name of variable indicating (binary) treatment assigned
#' @param outcome name of outcome variable
#' @param standardize whether to standardize data inside function
#' @param DIM difference in covariate means across treatment and control from `uw_delta()`
#' @param simulation logical. Whether running the function on bootstrap sample

pw_delta <- function(data, covariates, treatment, outcome, standardize = TRUE,
                     DIM, simulation = FALSE){

  pweights <- pw(data, covariates, treatment, outcome, standardize, simulation)
  pwdelta_j <- pweights*DIM
  # names(pwdelta_j) <- paste0("pw_", covariates)

  # prognosis weighted delta
  pwdelta <- sum(pweights*DIM, na.rm = TRUE)

  # R-squared from prognosis regression
  prog_mod_f <- paste(c(outcome, paste(covariates, collapse = " + ")), collapse = " ~ ")
  prog_mod <- lm(formula = prog_mod_f, data = data[data[[treatment]] == 0, ])
  prog_Rsq <- summary(prog_mod)$r.squared

  # R-squared from balance regression
  bal_mod_f <- paste(c(treatment, paste(covariates, collapse = " + ")), collapse = " ~ ")
  bal_mod <- lm(formula = bal_mod_f, data = data)
  bal_Rsq <- summary(bal_mod)$r.squared

  return(list(pw = pweights,
              pwdelta_j = pwdelta_j,
              pwdelta = pwdelta,
              prog_Rsq = prog_Rsq,
              bal_Rsq = bal_Rsq))

}

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
  pw_full <- pw(data = data, covariates = covariates, treatment = treatment,
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

  pw_bw <- pw(data = data_bw, covariates = covariates, treatment = treatment,
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
