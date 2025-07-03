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

#' Wrapper for pwtest with option for automatic model selection
#' @param data data.frame containing covariates, treatment assignment, and outcome variable
#' @param method character. Either "auto" for automatic model selection via contest() or "manual" for user-specified parameters
#' @param cv_folds integer. Number of cross-validation folds for contest() when method="auto"
#' @param debug logical. Whether to print detailed debugging information during contest()
#' @param ... Additional arguments passed to pwtest() 
#' @export

prognostic_balance <- function(data, 
                               method = c("auto", "manual"),
                               cv_folds = 3,
                               debug = FALSE,
                               ...) {
  
  method <- match.arg(method)
  
  # Extract pwtest parameters from ...
  pwtest_args <- list(...)
  
  # Set defaults for required pwtest parameters if not provided
  if (!"covariates" %in% names(pwtest_args)) {
    pwtest_args$covariates <- c("X1", "X2", "X3")
  }
  if (!"treatment" %in% names(pwtest_args)) {
    pwtest_args$treatment <- "Z"
  }
  if (!"outcome" %in% names(pwtest_args)) {
    pwtest_args$outcome <- "Y"
  }
  
  # Validate that required variables exist in data
  required_vars <- c(pwtest_args$treatment, pwtest_args$covariates, pwtest_args$outcome)
  missing_vars <- setdiff(required_vars, names(data))
  if (length(missing_vars) > 0) {
    stop("Missing variables in data: ", paste(missing_vars, collapse = ", "))
  }
  
  if (method == "manual") {
    # Manual mode: direct pwtest call with user parameters
    pwtest_args$data <- data
    result <- do.call(pwtest, pwtest_args)
    
    return(list(
      method = "manual",
      pwtest_output = result
    ))
    
  } else {
    # Auto mode: use contest() for model selection, then compare models
    
    if (debug) cat("=== Model selection diagnostics ===\n")
    
    # Extract control group data for contest()
    control_data <- data %>% 
      dplyr::filter(!!sym(pwtest_args$treatment) == 0) %>%
      dplyr::select(tidyselect::all_of(c(pwtest_args$outcome, pwtest_args$covariates)))
    
    if (nrow(control_data) == 0) {
      stop("No control group observations found. Check treatment variable coding.")
    }
    
    if (debug) {
      cat("Control group size:", nrow(control_data), "\n")
      cat("Covariates:", paste(pwtest_args$covariates, collapse = ", "), "\n")
      cat("Outcome:", pwtest_args$outcome, "\n")
    }
    
    # Run contest() to find best model
    if (debug) cat("\n--- Running model selection contest ---\n")
    contest_results <- contest(
      data_control = control_data,
      covariates = pwtest_args$covariates,
      outcome = pwtest_args$outcome,
      cv_folds = cv_folds,
      debug = debug
    )
    
    # Validate contest results
    validate_contest_output(contest_results)
    
    if (debug) {
      cat("\n--- Contest winner ---\n")
      cat("Best model:", contest_results$best_model_name, "\n")
      cat("Best R²:", round(contest_results$best_cv_rsq, 4), "\n")
      cat("Linear baseline R²:", round(contest_results$linear_cv_rsq, 4), "\n")
    } else {
      # Still show key info even without debug
      cat("Contest selected:", contest_results$best_model_name, 
          "(R² =", round(contest_results$best_cv_rsq, 4), ")\n")
      
      # Show model selection notes if relevant
      if (length(contest_results$model_selection_summary$overfitting_models) > 0) {
        cat("Note:", contest_results$model_selection_summary$explanation, "\n")
      }
    }
    
    # Prepare pwtest arguments for contest winner
    contest_pwtest_args <- pwtest_args
    contest_pwtest_args$data <- data
    contest_pwtest_args$spec_mod <- contest_results$best_spec
    contest_pwtest_args$engine <- contest_results$best_engine
    contest_pwtest_args$recipe <- contest_results$best_recipe
    
    # Run pwtest with contest winner
    if (debug) cat("\n--- Running pwtest with contest winner ---\n")
    pwtest_contest <- do.call(pwtest, contest_pwtest_args)
    
    # Prepare pwtest arguments for linear baseline
    linear_pwtest_args <- pwtest_args
    linear_pwtest_args$data <- data
    linear_pwtest_args$spec_mod <- linear_reg(mode = "regression", engine = "lm")
    linear_pwtest_args$engine <- "lm"
    linear_pwtest_args$recipe <- NULL  # No preprocessing for baseline
    
    # Run pwtest with linear baseline
    if (debug) cat("\n--- Running pwtest with linear baseline ---\n")
    pwtest_linear <- do.call(pwtest, linear_pwtest_args)
    
    # Compile results
    if (debug) cat("\n--- Compiling comparison results ---\n")
    
    results <- list(
      method = "auto",
      contest_results = list(
        model = contest_results$best_model_name,
        cv_rsq = contest_results$best_cv_rsq,
        p_value = pwtest_contest$pwdelta_p,
        pwdelta = pwtest_contest$pwdelta,
        pwdelta_se = pwtest_contest$pwdelta_se
      ),
      linear_results = list(
        model = "Linear Regression",
        cv_rsq = contest_results$linear_cv_rsq,
        p_value = pwtest_linear$pwdelta_p,
        pwdelta = pwtest_linear$pwdelta,
        pwdelta_se = pwtest_linear$pwdelta_se
      ),
      full_pwtest_outputs = list(
        contest_winner = pwtest_contest,
        linear_baseline = pwtest_linear
      ),
      contest_diagnostics = contest_results$model_selection_summary  # Always include for user info
    )
    
    # Add comparison summary
    results$summary <- list(
      improvement_in_rsq = contest_results$best_cv_rsq - contest_results$linear_cv_rsq,
      contest_more_significant = pwtest_contest$pwdelta_p < pwtest_linear$pwdelta_p,
      contest_stronger_effect = abs(pwtest_contest$pwdelta) > abs(pwtest_linear$pwdelta),
      recommendation = if (contest_results$best_cv_rsq > contest_results$linear_cv_rsq + 0.05) {
        paste("Use", contest_results$best_model_name, "- better performance: R^2 improvement is greater than .05")
      } else if (abs(pwtest_contest$pwdelta_p - pwtest_linear$pwdelta_p) > 0.1) {
        "Models give different conclusions - investigate further"
      } else {
        "Linear model sufficient - similar performance and interpretability"
      }
    )
    
    if (debug) {
      cat("\n--- Results Summary ---\n")
      cat("Contest winner p-value:", round(pwtest_contest$pwdelta_p, 4), "\n")
      cat("Linear baseline p-value:", round(pwtest_linear$pwdelta_p, 4), "\n")
      cat("R² improvement:", round(results$summary$improvement_in_rsq, 4), "\n")
      cat("Recommendation:", results$summary$recommendation, "\n")
    }
    
    # Add prognostic_balance class for potential method dispatch
    class(results) <- c("prognostic_balance", "list")
    
    return(results)
  }
}

#' Print method for prognostic_balance results
#' @param x prognostic_balance object
#' @param ... additional arguments
#' @export

print.prognostic_balance <- function(x, ...) {
  cat("=== Prognostic Balance Test Results ===\n\n")
  
  if (x$method == "manual") {
    cat("Mode: Manual (user-specified model)\n")
    cat("P-value:", round(x$pwtest_output$pwdelta_p, 4), "\n")
    cat("Effect size (δ):", round(x$pwtest_output$pwdelta, 4), "\n")
    cat("Standard error:", round(x$pwtest_output$pwdelta_se, 4), "\n")
    
  } else {
    cat("Mode: Automatic (contest-selected model)\n\n")
    
    cat("Contest Winner:\n")
    cat("  Model:", x$contest_results$model, "\n")
    cat("  CV R²:", round(x$contest_results$cv_rsq, 4), "\n")
    cat("  P-value:", round(x$contest_results$p_value, 4), "\n")
    cat("  Effect size (δ):", round(x$contest_results$pwdelta, 4), "\n\n")
    
    cat("Linear Baseline:\n")
    cat("  Model:", x$linear_results$model, "\n")
    cat("  CV R²:", round(x$linear_results$cv_rsq, 4), "\n")
    cat("  P-value:", round(x$linear_results$p_value, 4), "\n")
    cat("  Effect size (δ):", round(x$linear_results$pwdelta, 4), "\n\n")
    
    cat("Summary:\n")
    cat("  R² improvement:", round(x$summary$improvement_in_rsq, 4), "\n")
    cat("  Recommendation:", x$summary$recommendation, "\n")
    
    # Add model selection context if available  
    if ("contest_diagnostics" %in% names(x) && !is.null(x$contest_diagnostics)) {
      if (length(x$contest_diagnostics$overfitting_models) > 0) {
        cat("\nModel Selection:\n")
        cat(" ", x$contest_diagnostics$explanation, "\n")
      }
    }
  }
  
  cat("\n")
  invisible(x)
}
