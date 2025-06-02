#' Produces unweighted and prognosis-weighted test statistics with standard errors and p-values
#' @param data data.frame containing covariates, treatment assignment, and outcome variable
#' @param covariates character vector of names of placebo variables. This is not to be confused with `covs` argument passed onto `rdrobust` in the case of RD designs, which can be specified separately (`...`).
#' @param treatment name of variable indicating (binary) treatment assigned
#' @param running_var character string with running variable name. Ignored unless `rdd = TRUE`.
#' @param outcome name of outcome variable
#' @param nsims numeric scalar indicating number of bootstraps to use for the resampling-based p-values
#' @param oversample logical value. If `FALSE`, function gives resampling-based p-values using number of working bootstrap samples (which did not return an error) and prints a warning when that sample is smaller than the value in `nsims`. If `oversample = TRUE`, function continues to re-sample until number of working bootstrap samples reaches `nsims`.
#' @param se_type character string. Determines which type of standard error should be returned with the unweighted delta estimate. Takes the values of either `"analytic"` or `"bootstrap"` for non RDD cases, and can also take values "conventional", "bias-adjusted", and "robust" standard errors returned by `rdrobust` (when `rdd = TRUE`). If set to "analytic" in the latter, default results will use conventional standard errors.
#' @param simulation logical value.
#' @param rdd logical value. Whether test statistics are calculated using continuous RDD approach.
#' @param rd_estimator character. Whether to use the conventional ("h") or the bias-corrected local-polynomial point estimator ("b"). See `rdrobust()` for more details. Defaults to conventional estimate ("h").
#' @param equivalence logical value. If TRUE, carry out equivalence test as described in Hartman and Hidalgo (2018). For a full description of the bootstrapping procedure, see paper.
#' @param equiv_lower numerical value. If `equivalence = TRUE`, lower bound of the equivalence range as a multiplier for the standard deviation of potential outcomes under control. Takes 0.36 value by default (see Hartman and Hidalgo (2018)).
#' @param equiv_upper numerical value. If `equivalence = TRUE`, upper bound of the equivalence range as a multiplier for the standard deviation of potential outcomes under control. Takes 0.36 value by default (see Hartman and Hidalgo (2018)).
#' @param method character name of function to model (or train a model on) $Y^C(0)$, e.g., "glm", "ranger", "gxboost". Must be a function that returns an object of class compatible with `predict()`. Include additional function arguments in the call (`...`), otherwise default argument values for each method are supplied where default is specified. Requires package containing function to be installed and loaded separately.
#' @param predict_args named list. If using a non-linear model (e.g., "glm", "ranger", "xgb.train"), specify arguments passed onto `predict()` function for object of class returned by the function specified in `model` argument. For example, if `method = "ranger"`, this argument will defined the values passed onto to `predict.ranger()` to fit values of $Y^T(C)$.
#' @param ... arguments passed onto `rdrobust()` function
#' @import rdrobust
#' @importFrom stats pnorm
#' @export

pwtest <- function(data,
                   covariates = c("X1", "X2", "X3"),
                   treatment = "Z", running_var = NULL,
                   outcome = "Y",
                   nsims = 500, oversample = FALSE,
                   se_type = "analytic",
                   simulation = FALSE, rdd = FALSE,
                   rd_estimator = "h",
                   equivalence = FALSE,
                   equiv_lower = NULL,
                   equiv_upper = NULL,
                   method = "lm",
                   predict_args = NULL,
                   ...) {
  require(rdrobust)
  argg <- as.list(match.call())
  c <- NULL

  # REVIEW write warning if specifying differing values for y and outcome, x and treatment
  # write message if set se_type = 'analytic' and rdd = TRUE
  if (rdd & is.null(running_var)) stop("When using 'rdd = TRUE', need to specify 'running_var'.")

  if (!rdd & (length(se_type) != 1L |
              !se_type %in% c("analytic", "bootstrap"))) {
    stop("`se_type` must take be set to either 'analytic' or 'bootstrap'.")
  }

  if (rdd & (length(se_type) != 1L |
             !se_type %in% c("analytic", "bootstrap", "conventional", "bias-corrected", "robust"))) {
    stop("`se_type` must take be set to either 'analytic', 'bootstrap', 'conventional', 'bias-corrected', 'robust'.")
  }

  if(rdd & equivalence) stop("The equivalence test approach is not yet compatible with RDD designs.")

  # restrict data to variables of interest
  data <- data %>% dplyr::select(tidyselect::all_of(c(treatment, covariates, running_var, outcome)))

  # set default lower and upper bounds if null
  # reference value .36*SD of potential outcomes under treatment (see Hartman and Hidalgo 2018)
  if(equivalence & is.null(equiv_lower)) equiv_lower <- .36;
  if(equivalence & is.null(equiv_upper)) equiv_upper <- .36;

  # observed statistics-------------------------------------

  # any arguments not supplied take default values
  argg_def <- formals()
  missing_args <- setdiff(names(argg_def)[!names(argg_def) %in% names(argg)], "...")
  argg_add <- setNames(argg_def[missing_args], missing_args)
  argg <- c(argg, argg_add)
  argg$standardize <- TRUE
  argg$DIM <- uwdelta_obs$dim

  if (!rdd) {
    uwdelta_obs <- do.call("uw_delta", args = argg)
    # linear regression, default pw method
    if(method == "lm"){
      pwdelta_obs <- do.call("pw_delta", args = argg)
    } else {
      # non-linear methods
      # REVIEW NOT SURE WE SHOULD STANDARDIZE IN THESE CASES
      # standardize data relative to entire study group (the finite population)
      if(standardize){
        data <- data %>%
          dplyr::mutate_at(.vars = c(outcome, covariates),
                           .funs = stdr) %>% as.data.frame()
      }
      argg$data <- data[data[[treatment]]==0,]

      # SPECIFY/TRANSLATE ARGUMENT NAMES TO BE PASSED ONTO SUBSET OF
      # "COMPATIBLE" MODELS
      relabel_args <- function(argg, method){
        if(method == "glmnet"){
          argg$x <- as.matrix(argg$data[,covariates])
          argg$y <- argg$data[[outcome]]
        }
        if(method %in% c("xgboost", "xgb.train")){

        }
        if(method %in% c("gbart", "pbart", "lbart", "wbart", "mc.gbart")){
          argg$x.train <- argg$data[,covariates]
          argg$y.train <- argg$data[[outcome]]
        }
        return(argg)
      }

      argg <- relabel_args(argg, method)
      m_output <- do.call(method, args = argg)
      predict_args$object <- m_output
      predict_args$data <- data[data[[treatment]] == 1,]
      predict_args$newdata <- data[data[[treatment]] == 1,]
      m_fit_YT0 <- do.call("predict", args = predict_args)
      pwdelta_obs <- mean(m_fit_YT0, na.rm = TRUE) - mean(argg$data[[outcome]])
    }
  } else {
    # test of continuity (RD designs)
    pwdelta_obs <- do.call("pw_delta_rdd", args = argg)
  }

  # check that no NA coefficients
  # coeflabs <- grep("coef_", names(pwdelta_obs), value = TRUE)
  coef_na <- is.na(pwdelta_obs$pw)
  if (any(coef_na)) {
    vars <- names(which(coef_na))
    warning(paste0("The following variables return missing coefficients in regression Y[Z=0] ~ X[Z=0]: ", vars, ". Consider omitting those variables for pwdelta."))
  }

  # bootstrapping -------------------------------------

  # resample from control group with replacement
  treat_i <- which(data[[treatment]] == 1)
  control_i <- which(data[[treatment]] == 0)
  samples <- get_bootstrap_samples(control_i, length(treat_i), nsims)

  if(equivalence){
    samples2 <- get_bootstrap_samples(control_i, length(treat_i), nsims)
  }

  # standardize bootstrap population relative to control group SD
  data_stdc <- std_data(data, c(outcome, covariates), treatment)


  # obtain delta distribution from bootstrap samples
  delta_sim <- apply(samples, 2, function(z) {
    dat_uw <- data_stdc[z, ]
    # change treatment condition from control to treatment for bootstrap treatment group
    dat_uw[[treatment]][1:length(treat_i)] <- 1

    dat_pw <- data[z, ]
    # change treatment condition from control to treatment for bootstrap treatment group
    dat_pw[[treatment]][1:length(treat_i)] <- 1

    # permutation under equivalence test (first one-sided test using lower bound)
    if(equivalence){
      dat_uw[[outcome]][1:length(treat_i)] <- dat_uw[[outcome]][1:length(treat_i)] - equiv_lower
    }

    if (!rdd) {
      uwdelta_sim <- safe_delta(uw_delta, dat_uw, uwdelta_obs, covariates, treatment, outcome, standardize = FALSE)

      if(method == "lm"){
        pwdelta_sim <- safe_delta(pw_delta, dat_pw, pwdelta_obs, covariates, treatment, outcome, standardize = TRUE, DIM = uwdelta_sim$dim, simulation = simulation)

      } else {
        # non-linear methods
        # standardize data relative to entire study group (the finite population)
        if(standardize){
          dat_pw <- dat_pw %>%
            dplyr::mutate_at(.vars = c(outcome, covariates),
                             .funs = stdr) %>% as.data.frame()
        }

        predict_args$object <- m_output
        predict_args$data <- dat_pw[dat_pw[[treatment]] == 1,]
        predict_args$newdata <- dat_pw[dat_pw[[treatment]] == 1,]
        m_fit_YT0 <- do.call("predict", args = predict_args)

        predict_args$data <- dat_pw[dat_pw[[treatment]] == 0,]
        predict_args$newdata <- dat_pw[dat_pw[[treatment]] == 0,]
        m_fit_YC0 <- do.call("predict", args = predict_args)

        pwdelta_obs <- mean(m_fit_YT0, na.rm = TRUE) - mean(m_fit_YC0, na.rm = TRUE)
      }

      return(c(uwdelta_sim, pwdelta_sim))

    } else {
      bs_t <- dat_pw[[treatment]] == 1
      # invert the running variable for the bootstrap sample of treatment observations
      dat_pw[bs_t, running_var] <- -(dat_pw[bs_t, running_var])
      # cutoff takes default value as in rdrobust() if not specified by user
      cutoff <- ifelse("c" %in% names(argg), argg$c, 0)
      pwdelta_sim <- tryCatch(
        expr = {
          do.call("pw_delta_rdd", args = modifyList(argg, list(data = dat_pw,
                                                               standardize = TRUE,
                                                               simulation = FALSE)))
        },
        error = function(e) {
          out <- vector(mode = "list", length = length(pwdelta_obs))
          names(out) <- names(pwdelta_obs)
          return(out)
        }
      )

      return(c(pwdelta_sim))
    }
  })

  # repeat procedure for second sample (only applies in equivalence test)
  if(equivalence){

  delta_sim2 <- apply(samples2, 2, function(z) {
    dat_uw <- data_stdc[z, ]
    # change treatment condition from control to treatment for bootstrap treatment group
    dat_uw[[treatment]][1:length(treat_i)] <- 1

    dat_pw <- data[z, ]
    # change treatment condition from control to treatment for bootstrap treatment group
    dat_pw[[treatment]][1:length(treat_i)] <- 1

    # permutation under equivalence test (second one-sided test using upper bound)
    dat_uw[[outcome]][1:length(treat_i)] <- dat_uw[[outcome]][1:length(treat_i)] + equiv_upper

    uwdelta_sim <- safe_delta(uw_delta, dat_uw, uwdelta_obs, covariates, treatment, outcome, standardize = FALSE)
    pwdelta_sim <- safe_delta(pw_delta, dat_pw, pwdelta_obs, covariates, treatment, outcome, standardize = TRUE, DIM = uwdelta_sim$dim, simulation = simulation)
    return(c(uwdelta_sim, pwdelta_sim))

  })

  }

  # bootstrap sample checks ----------------------------------

  effective_sample <- min(
    sum(!is.na(sapply(delta_sim, function(x) x[["uwdelta"]]))),
    sum(!is.na(sapply(delta_sim, function(x) x[["pwdelta"]])))
  )

  # REVIEW implement resampling in cases
  # when effective sample of second distribution < nsims and oversample == TRUE
  if(equivalence){
    effective_sample2 <- min(
      sum(!is.na(sapply(delta_sim2, function(x) x[["pwdelta"]]))),
      sum(!is.na(sapply(delta_sim2, function(x) x[["pwdelta"]])))
    )
  }

  if (!identical(effective_sample, as.integer(nsims)) & !oversample) {
    warning(paste0("Effective bootstrap sample to calculate delta p-values is of size ", effective_sample, ". Consider changing argument `oversample` to `TRUE`"))
  }

  if (oversample &
      (!identical(effective_sample, as.integer(nsims)))) {
    # re-sample until effective is equal or greater to bootstrap sample size argument
    while (effective_sample < as.integer(nsims)) {
      samples <- get_bootstrap_samples(control_i, length(treat_i), nsims-effective_sample)

      # restrict bootstrap sample to complete observations

      # obtain delta distribution based on resamples
      # REVIEW: not considering equivalence test
      delta_sim_add <- apply(samples, 2, function(z) {
        dat_uw <- data_stdc[z, ]
        # change treatment condition from control to treatment for sample
        dat_uw[[treatment]][1:length(treat_i)] <- 1

        dat_pw <- data[z, ]
        # change treatment condition from control to treatment for sample
        dat_pw[[treatment]][1:length(treat_i)] <- 1

        if (!rdd) {
          uwdelta_sim <- safe_delta(uw_delta, dat_uw, uwdelta_obs, covariates, treatment, outcome, standardize = FALSE)
          pwdelta_sim <- safe_delta(pw_delta, dat_pw, pwdelta_obs, covariates, treatment, outcome, standardize = TRUE, DIM = uwdelta_sim$dim, simulation = simulation)
          return(c(uwdelta_sim, pwdelta_sim))
        } else {
          # RD DESIGN ############
          # cutoff takes default value as in rdrobust() if not specified by user
          cutoff <- ifelse("c" %in% names(argg), argg$c, 0)
          # invert the running variable for the bootstrap sample of treatment observations
          bs_t <- dat_pw[[treatment]] == 1
          dat_pw[bs_t, running_var] <- -(dat_pw[bs_t, running_var] - cutoff)

          pwdelta_sim <- tryCatch(
            expr = {
              pw_delta_rdd(dat_pw, covariates, running_var, treatment, outcome,
                           rd_estimator,
                           standardize = TRUE, simulation = simulation
              )
            },
            error = function(e) {
              out <- rep(NA, length(pwdelta_obs))
              names(out) <- names(pwdelta_obs)
              return(out)
            }
          )
          return(pwdelta_sim)
        }
      })

      delta_sim <- c(delta_sim, delta_sim_add)
      na_check <- sapply(delta_sim, function(e) !is.na(e[["uwdelta"]]) & !is.na(e[["pwdelta"]]))
      delta_sim <- delta_sim[na_check]
      effective_sample <- sum(na_check)
    }
    # if bootstrap sample greater than argument, randomly select samples to match sample size argument
    if (effective_sample > as.integer(nsims)) {
      delta_sim <- delta_sim[sample(1:length(delta_sim), nsims, replace = FALSE)]
    }
  }

  # oversample for upper bound distribution of equivalence test
  if (oversample &
      (!identical(effective_sample2, as.integer(nsims)))) {
    # re-sample until effective is equal or greater to bootstrap sample size argument
    while (effective_sample2 < as.integer(nsims)) {
      samples <- get_bootstrap_samples(control_i, length(treat_i), nsims-effective_sample2)

      # restrict bootstrap sample to complete observations

      # obtain delta distribution based on resamples
      delta_sim_add <- apply(samples, 2, function(z) {
        dat_uw <- data_stdc[z, ]
        # change treatment condition from control to treatment for sample
        dat_uw[[treatment]][1:length(treat_i)] <- 1

        dat_pw <- data[z, ]
        # change treatment condition from control to treatment for sample
        dat_pw[[treatment]][1:length(treat_i)] <- 1


        uwdelta_sim <- safe_delta(uw_delta, dat_uw, uwdelta_obs, covariates, treatment, outcome, standardize = FALSE)
        pwdelta_sim <- safe_delta(pw_delta, dat_pw, pwdelta_obs, covariates, treatment, outcome, standardize = TRUE, DIM = uwdelta_sim$dim, simulation = simulation)
        return(c(uwdelta_sim, pwdelta_sim))
      })

      delta_sim2 <- c(delta_sim2, delta_sim_add)
      na_check <- sapply(delta_sim2, function(e) !is.na(e[["uwdelta"]]) & !is.na(e[["pwdelta"]]))
      delta_sim2 <- delta_sim2[na_check]
      effective_sample2 <- sum(na_check)
    }
    # if bootstrap sample greater than argument, randomly select samples to match sample size argument
    if (effective_sample2 > as.integer(nsims)) {
      delta_sim2 <- delta_sim2[sample(1:length(delta_sim2), nsims, replace = FALSE)]
    }
  }

  # SEs as SD of sampling distribution
  uwdelta_se2 <- sd(unlist(sapply(delta_sim, function(x) x[["uwdelta"]])), na.rm = TRUE) * (nsims - 1) / nsims

  # p-value (from bootstrap distribution)
  #REVIEW: maybe standardization different between observed and bootstrap???? p-value is 1 for Bohlke
  if(!rdd | (rdd & se_type == "bootstrap")){
    pwdelta_se <- sd(unlist(sapply(delta_sim, function(x) x$pwdelta)), na.rm = TRUE) * (nsims - 1) / nsims
    pw_delta_p <- sum(abs(unlist(sapply(delta_sim, function(x) x$pwdelta))) >= abs(pwdelta_obs$pwdelta), na.rm = TRUE) / effective_sample

  }

  if (rdd & se_type %in% c("analytic", "conventional", "bias-corrected", "robust")) {
    pwdelta_se <- switch(se_type,
                         analytic = pwdelta_obs$rdrobust_output$se["Conventional", ],
                         conventional = pwdelta_obs$rdrobust_output$se["Conventional", ],
                         `bias-corrected` = pwdelta_obs$rdrobust_output$se["Bias-Corrected", ],
                         robust = pwdelta_obs$rdrobust_output$se["Robust", ],
    )
    # p-value (from analytic SE)
    pw_delta_t <- pwdelta_obs$pwdelta / pwdelta_se
    pw_delta_p <- 2 * pnorm(-abs(pw_delta_t))
  }

  # output ------------------------------------------------------------------

  if (!rdd) {
      if(method == "lm"){
        if(!equivalence){
        estimates <- c(
          uwdelta_obs[-length(uwdelta_obs)],
          uwdelta_se = switch(se_type,
                              analytic = uwdelta_obs$uwdelta_se,
                              bootstrap = uwdelta_se2
          ),
          uwdelta_p = get_pvalue_uw(delta_sim,  uwdelta_obs$uwdelta, effective_sample),
          # number of unique bootstrap samples (if different from nsims)
          n_bootstrap = effective_sample,
          pwdelta_obs,
          pwdelta_se = unname(pwdelta_se),
          pwdelta_p = get_pvalue_pw(delta_sim, pwdelta_obs$pwdelta, effective_sample),
          prog_Rsq = pwdelta_obs$prog_Rsq,
          prog_Rsq_mean = mean(delta_sim$prog_Rsq, na.rm = TRUE),
          bal_Rsq = pwdelta_obs$bal_Rsq,
          bal_Rsq_mean = mean(delta_sim$bal_Rsq, na.rm = TRUE)
        )
    } else {
      estimates <- c(
        uwdelta_obs[-length(uwdelta_obs)],
        uwdelta_se = switch(se_type,
                            analytic = uwdelta_obs$uwdelta_se,
                            bootstrap = uwdelta_se2),
        uwdelta_pl = sum(unlist(sapply(delta_sim, function(x) x[["uwdelta"]])) <= uwdelta_obs$uwdelta) / effective_sample,
        uwdelta_pu = sum(unlist(sapply(delta_sim2, function(x) x[["uwdelta"]])) >= uwdelta_obs$uwdelta) / effective_sample2,
        # number of unique bootstrap samples (if different from nsims)
        n_bootstrap_l = effective_sample,
        n_bootstrap_u = effective_sample2,
        pwdelta_obs,
        pwdelta_se = unname(pwdelta_se),
        pwdelta_pl = sum(unlist(sapply(delta_sim, function(x) x[["pwdelta"]])) <= pwdelta_obs$pwdelta) / effective_sample,
        pwdelta_pu = sum(unlist(sapply(delta_sim2, function(x) x[["uwdelta"]])) >= pwdelta_obs$pwdelta) / effective_sample2,
        prog_Rsq = pwdelta_obs$prog_Rsq,
        prog_Rsq_mean = mean(delta_sim$prog_Rsq, na.rm = TRUE),
        bal_Rsq = pwdelta_obs$bal_Rsq,
        bal_Rsq_mean = mean(delta_sim$bal_Rsq, na.rm = TRUE)
      )

      estimates <- c(estimates,
                     uwdelta_equivtest = ifelse(estimates[["uwdelta_pl"]] < .05 & estimates[["uwdelta_pu"]] < .05, "Reject", "Fail to reject"),
                     pwdelta_equivtest = ifelse(estimates[["pwdelta_pl"]] < .05 & estimates[["pwdelta_pu"]] < .05, "Reject", "Fail to reject"))
    }
      } else {
        #non-linear methods
        estimates <- list(
          uwdelta_obs[-length(uwdelta_obs)],
          uwdelta_se = switch(se_type,
                              analytic = uwdelta_obs$uwdelta_se,
                              bootstrap = uwdelta_se2
          ),
          uwdelta_p = get_pvalue_uw(delta_sim,  uwdelta_obs$uwdelta, effective_sample),
          # number of unique bootstrap samples (if different from nsims)
          n_bootstrap = effective_sample,
          pwdelta = pwdelta_obs,
          pwdelta_se = unname(pwdelta_se),
          pwdelta_p = get_pvalue_pw(delta_sim, pwdelta_obs$pwdelta, effective_sample),
          prog_Rsq = pwdelta_obs$prog_Rsq,
          prog_Rsq_mean = mean(delta_sim$prog_Rsq, na.rm = TRUE),
          bal_Rsq = pwdelta_obs$bal_Rsq,
          bal_Rsq_mean = mean(delta_sim$bal_Rsq, na.rm = TRUE)
        )
      }
  } else {
    estimates <- list(
      dii = pwdelta_obs$dii,
      uwdelta = pwdelta_obs$uwdelta,
      uwdelta_se = uwdelta_se2,
      uwdelta_p =
        sum(abs(unlist(sapply(delta_sim, function(x) x[["uwdelta"]])))>= abs(pwdelta_obs$uwdelta), na.rm = TRUE) / effective_sample,
      pw = pwdelta_obs$pw,
      pwdelta = pwdelta_obs$pwdelta,
      pwdelta_se = pwdelta_se,
      pwdelta_p = pw_delta_p,
      # number of unique bootstrap samples (if different from nsims)
      n_bootstrap = effective_sample,
      rdrobust_output = pwdelta_obs$rdrobust_output,
      prog_Rsq = pwdelta_obs$prog_Rsq,
      bal_Rsq = pwdelta_obs$bal_Rsq
    )
  }


  return(estimates)
}
