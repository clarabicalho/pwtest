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
                   rd_estimator = "h", ...) {
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

  # restrict data to variables of interest
  data <- data %>% dplyr::select(tidyselect::all_of(c(treatment, covariates, running_var, outcome)))

  # observed statistics-------------------------------------

  if (!rdd) {
    uwdelta_obs <- uw_delta(data, covariates, treatment, outcome,
                            standardize = TRUE,
                            simulation = simulation
    )
    pwdelta_obs <- pw_delta(data, covariates, treatment, outcome,
                            standardize = TRUE,
                            DIM = uwdelta_obs$dim,
                            simulation = simulation
    )
  } else {
    # any arguments not supplied take `pw_delta_rdd` default values
    argg$standardize <- TRUE
    argg_def <- formals()
    missing_args <- setdiff(names(argg_def)[!names(argg_def) %in% names(argg)], "...")
    argg_add <- setNames(argg_def[missing_args], missing_args)
    argg <- c(argg, argg_add)
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
  samples <- replicate(
    nsims,
    c(
      sample(control_i, length(treat_i), replace = TRUE),
      sample(control_i, length(control_i), replace = TRUE)
    )
  )

  # standardize bootstrap population relative to control group SD
  data_stdc <- data %>%
    mutate_at(
      .vars = c(outcome, covariates),
      .funs = function(x) {
        x_c <- x[data[[treatment]] == 0]
        sd_x <- sd(x_c, na.rm = TRUE) * (sum(!is.na(x_c)) - 1) / sum(!is.na(x_c))
        return((x - mean(x_c, na.rm = TRUE)) / sd_x)
      }
    ) %>%
    as.data.frame()

  # obtain delta distribution from bootstrap samples
  delta_sim <- apply(samples, 2, function(z) {
    dat_uw <- data_stdc[z, ]
    # change treatment condition from control to treatment for bootstrap treatment group
    dat_uw[[treatment]][1:length(treat_i)] <- 1

    dat_pw <- data[z, ]
    # change treatment condition from control to treatment for bootstrap treatment group
    dat_pw[[treatment]][1:length(treat_i)] <- 1

    if (!rdd) {
      uwdelta_sim <- tryCatch(
        expr = {
          uw_delta(dat_uw, covariates, treatment, outcome, standardize = FALSE)
        },
        error = function(e) {
          out <- vector(mode = "list", length = length(uwdelta_obs))
          names(out) <- names(uwdelta_obs)
          return(out)
          message("Error with UW delta estimation of bootstrap sample.")
        }
      )

      pwdelta_sim <- tryCatch(
        expr = {
          pw_delta(dat_pw, covariates, treatment, outcome,
                   standardize = TRUE,
                   DIM = uwdelta_sim$dim,
                   simulation = simulation
          )
        },
        error = function(e) {
          out <- vector(mode = "list", length = length(pwdelta_obs))
          names(out) <- names(pwdelta_obs)
          return(out)
          # REVIEW line below not printing with results
          message("Error with pw delta estimation of bootstrap sample.")
        }
      )

      return(c(uwdelta_sim, pwdelta_sim))
    } else {
      # cutoff takes default value as in rdrobust() if not specified by user
      cutoff <- ifelse("c" %in% names(argg), argg$c, 0)
      # invert the running variable for the bootstrap sample of treatment observations
      bs_t <- dat_pw[[treatment]] == 1
      dat_pw[bs_t, running_var] <- -(dat_pw[bs_t, running_var] - cutoff)
      arggn <- argg
      arggn$data <- dat_pw
      arggn$standardize <- TRUE # REVIEW: already standardized
      arggn$simulation <- FALSE
      # relative to control group

      # run simulated value of pw statistic
      pwdelta_sim <- tryCatch(
        expr = {
          do.call("pw_delta_rdd", args = arggn)
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

  # bootstrap sample checks ----------------------------------

  effective_sample <- min(
    sum(!is.na(sapply(delta_sim, function(x) x[["uwdelta"]]))),
    sum(!is.na(sapply(delta_sim, function(x) x[["pwdelta"]])))
  )

  if (!identical(effective_sample, as.integer(nsims)) & !oversample) {
    warning(paste0("Effective bootstrap sample to calculate delta p-values is of size ", effective_sample, ". Consider changing argument `oversample` to `TRUE`"))
  }

  if (oversample &
      (!identical(effective_sample, as.integer(nsims)))) {
    # re-sample until effective is equal or greater to bootstrap sample size argument
    while (effective_sample < as.integer(nsims)) {
      samples <- replicate(
        nsims-effective_sample,
        c(
          sample(control_i, length(treat_i), replace = TRUE),
          sample(control_i, length(control_i), replace = TRUE)
        )
      )

      # restrict bootstrap sample to complete observations

      # obtain delta distribution based on resamples
      delta_sim_add <- apply(samples, 2, function(z) {
        dat_uw <- data_stdc[z, ]
        # change treatment condition from control to treatment for sample
        dat_uw[[treatment]][1:length(treat_i)] <- 1

        dat_pw <- data[z, ]
        # change treatment condition from control to treatment for sample
        dat_pw[[treatment]][1:length(treat_i)] <- 1

        if (!rdd) {
          uwdelta_sim <- tryCatch(
            expr = {
              uw_delta(dat_uw, covariates, treatment, outcome, standardize = FALSE)
            },
            error = function(e) {
              out <- vector(mode = "list", length = length(uwdelta_obs))
              names(out) <- names(uwdelta_obs)
              return(out)
              message("Error with UW delta estimation of bootstrap sample.")
            }
          )

          pwdelta_sim <- tryCatch(
            expr = {
              pw_delta(dat_pw, covariates, treatment, outcome,
                       standardize = TRUE,
                       DIM = uwdelta_sim$dim,
                       simulation = simulation
              )
            },
            error = function(e) {
              out <- vector(mode = "list", length = length(pwdelta_obs))
              names(out) <- names(pwdelta_obs)
              return(out)
              # REVIEW below not printing with results
              message("Error with pw delta estimation of bootstrap sample.")
            }
          )

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
              # REVIEW line below not printing with results
              # message("Error with pw delta estimation of bootstrap sample.")
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
    estimates <- c(
      uwdelta_obs[-length(uwdelta_obs)],
      uwdelta_se = switch(se_type,
                          analytic = uwdelta_obs$uwdelta_se,
                          bootstrap = uwdelta_se2
      ),
      uwdelta_p =
        sum(abs(unlist(sapply(delta_sim, function(x) x[["uwdelta"]]))) >=
              abs(uwdelta_obs$uwdelta)) / effective_sample,
      # number of unique bootstrap samples (if different from nsims)
      n_bootstrap = effective_sample,
      pwdelta_obs,
      pwdelta_se = unname(pwdelta_se),
      pwdelta_p =
        sum(abs(unlist(sapply(delta_sim, function(x) x[["pwdelta"]]))) >= abs(pwdelta_obs$pwdelta), na.rm = TRUE) / effective_sample,
      prog_Rsq = pwdelta_obs$prog_Rsq,
      bal_Rsq = pwdelta_obs$bal_Rsq
    )
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
