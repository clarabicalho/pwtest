#' Produces unweighted and prognosis-weighted test statistics with standard errors and p-values
#' @param data data.frame containing covariates, treatment assignment, and outcome variable
#' @param covariates character vector of names of placebo variables. This is not to be confused with `covs` argument passed onto `rdrobust`, which can be specified separately (`...`).
#' @param treatment name of variable indicating (binary) treatment assigned
#' @param running_var character string with running variable name.
#' @param outcome name of outcome variable
#' @param n_bootstraps numeric scalar indicating number of bootstraps to use for the resampling-based p-values
#' @param se_type character string. Determines which type of standard error should be returned. Takes the values of "bootstrap" if estimated from the bootstrap samples, or "conventional", "bias-adjusted", and "robust" standard errors returned by `rdrobust`.
#' @param rd_estimator character. Whether to use the conventional ("h") or the bias-corrected local-polynomial point estimator ("b"). See `rdrobust()` for more details. Defaults to conventional estimate ("h").
#' @param return_uwtest logical. Whether to return results from unweighted test as well.
#' @param ... arguments passed onto `rdrobust()` function
#' @import rdrobust rdrobust
#' @importFrom utils modifyList
#' @importFrom stats pnorm
#' @export

pwtest_rdd <- function(data,
                       covariates = c("X1", "X2", "X3"),
                       treatment = "Z",
                       running_var = NULL,
                       outcome = "Y",
                       n_bootstraps = 500,
                       se_type = "conventional",
                       rd_estimator = "h",
                       return_uwtest = FALSE,
                       ...) {

  # REVIEW write warning if specifying differing values for y and outcome, x and treatment
  # write message if set se_type = 'analytic' and rdd = TRUE

  if (is.null(running_var)) stop("Need to specify 'running_var'.")
  if ((length(se_type) != 1L | !se_type %in% c("bootstrap", "conventional", "bias-corrected", "robust"))) {
    stop("`se_type` must take be set to either 'bootstrap', 'conventional', 'bias-corrected', 'robust'.")
  }

  # restrict data to variables of interest
  data <- data[,c(treatment, covariates, outcome, running_var)] %>%
    tidyr::drop_na()

  # standardize data relative to entire study group (finite population)
  data <- data %>%
    mutate_at(.vars = c(outcome, covariates),
              .funs = scale) %>% as.data.frame()
  # observed statistics-------------------------------------

  # any arguments not supplied take `pw_delta_rdd` default values
  argg <- as.list(match.call())[-1]
  c <- NULL

  arg_formals <- formals(pwtest_rdd) ## formals with default arguments
  for (v in names(arg_formals)){
    if (!(v %in% names(argg)))
      ## if arg value is missing, add its default
      argg <- append(argg, arg_formals[v])
  }
  pwdelta_obs <- do.call("pw_delta_rdd", args = argg)

  # check that no NA coefficients
  # coeflabs <- grep("coef_", names(pwdelta_obs), value = TRUE)
  coef_na <- is.na(pwdelta_obs$pw)
  if (any(coef_na)) {
    vars <- names(which(coef_na))
    warning(paste0("The following variables return missing coefficients in regression Y[Z=0] ~ X[Z=0]: ", vars, ". Consider diagnosing or omitting those variables from the test."))
  }

  # bootstrapping -------------------------------------

  # resample from control group with replacement
  treat_i <- which(data[[treatment]] == 1)
  control_i <- which(data[[treatment]] == 0)
  samples <- get_bootstrap_samples(control_i, length(treat_i), n_bootstraps)

  # standardize bootstrap population relative to control group SD
  data_stdc <- std_data(data, c(outcome, covariates), treatment)

  # obtain delta distribution from bootstrap samples
  pwdelta_from_draw <- function(z) {
    bsample <- data_stdc[z, ]
    # change treatment condition from control to treatment for bootstrap treatment group
    bsample[[treatment]][1:length(treat_i)] <- 1

    bs_t <- bsample[[treatment]] == 1
    # invert the running variable for the bootstrap sample of treatment observations
    bsample[bs_t, running_var] <- -(bsample[bs_t, running_var])
    # cutoff takes default value as in rdrobust() if not specified by user
    cutoff <- ifelse("c" %in% names(argg), argg$c, 0)
    pwbstats <- tryCatch(
      expr = {
        do.call("pw_delta_rdd", args = modifyList(argg, list(data = bsample)))
      },
      error = function(e) {
        out <- vector(mode = "list", length = length(pwdelta_obs))
        names(out) <- names(pwdelta_obs)
        return(out)
      }
    )
    return(c(pwbstats))
  }

  bstats <- capture_warnings_apply(samples, 2, pwdelta_from_draw)

  # bootstrap sample checks ----------------------------------

  nc_bootstraps <- sum(!is.na(unlist(sapply(bstats, function(x) x[["pwdelta"]]))))
  nc_bootstraps <- min(nc_bootstraps, sum(!is.na(unlist(sapply(bstats, function(x) x[["uwdelta"]])))))

  # re-sample until effective is equal or greater to bootstrap sample size argument
  while (nc_bootstraps < as.integer(n_bootstraps)) {
    add_samples <- get_bootstrap_samples(control_i, length(treat_i), n_bootstraps-nc_bootstraps)

    # obtain delta distribution based on resamples
    add_bstats <- capture_warnings_apply(add_samples, 2, pwdelta_from_draw)

    bstats <- c(bstats, add_bstats)
    na_check <- sapply(bstats, function(e) !is.null(e[["pwdelta"]]))
    na_check2 <- sapply(bstats, function(e) !is.null(e[["uwdelta"]]))
    nc_bootstraps <- sum(!is.na(unlist(sapply(bstats, function(x) x[["pwdelta"]]))))

    bstats <- bstats[na_check & na_check2]
    nc_bootstraps <- min(sum(na_check), sum(na_check2))
  }

  if(length(bstats) > n_bootstraps) bstats <- bstats[1:n_bootstraps]

  # p-value (from bootstrap distribution)
  #REVIEW: maybe standardization different between observed and bootstrap???? p-value is 1 for Bohlke
  if(se_type == "bootstrap"){
    pwdelta_se <- sd(unlist(sapply(bstats, function(x) x$pwdelta)), na.rm = TRUE) * (n_bootstraps - 1) / n_bootstraps
    pw_delta_p <- sum(abs(unlist(sapply(bstats, function(x) x$pwdelta))) >= abs(pwdelta_obs$pwdelta), na.rm = TRUE) / n_bootstraps
  }

  if (se_type %in% c("conventional", "bias-corrected", "robust")) {
    pwdelta_se <- switch(se_type,
                         conventional = pwdelta_obs$rdrobust_output$se["Conventional", ],
                         `bias-corrected` = pwdelta_obs$rdrobust_output$se["Bias-Corrected", ],
                         robust = pwdelta_obs$rdrobust_output$se["Robust", ],
    )
    # p-value (from analytic SE)
    pw_delta_t <- pwdelta_obs$pwdelta / pwdelta_se
    pw_delta_p <- 2 * pnorm(-abs(pw_delta_t))
  }

  # output ------------------------------------------------------------------

  dii_dist <- do.call(rbind, lapply(bstats, function(l) l$dii))

  estimates <- data.frame(
    test_type = "continuity",
    n_bootstraps = n_bootstraps,
    uwdelta = pwdelta_obs$uwdelta,
    uwdelta_se = sd(apply(dii_dist, 2, sum)) * (n_bootstraps - 1) / n_bootstraps,
    uwdelta_p = sum(abs(rowSums(dii_dist)) >= abs(pwdelta_obs$uwdelta)) / n_bootstraps,
    pwdelta = pwdelta_obs$pwdelta,
    pwdelta_se = pwdelta_se,
    pwdelta_p = pw_delta_p
  )

  if(!return_uwtest){
    estimates <- estimates %>% dplyr::select(-uwdelta, -uwdelta_se, -uwdelta_p)
  }

  cov_table <- cbind(
    prognosis = pwdelta_obs$pw,
    dii = pwdelta_obs$dii,
    `p-value` = sapply(1:length(covariates), function(v){
      sum(abs(dii_dist[,v]) >= abs(pwdelta_obs$dii[v])) / nrow(dii_dist)
    }))

  estimates <- list(
    estimates = estimates,
    cov_table = cov_table,
    rdrobust_output = pwdelta_obs$rdrobust_output,
    prog_Rsq = pwdelta_obs$bal_Rsq,
    bal_Rsq = pwdelta_obs$bal_Rsq
  )

  return(estimates)
}
