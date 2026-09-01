# Clean and standardize data before calculating coefficients of prognosis
# regression on control units
# Returns vector of labelled covariate (prognosis) weights

#' @importFrom stats sd
stdr <- function(x){
  (x - mean(x, na.rm = TRUE))/(stats::sd(x, na.rm = TRUE)*(length(na.omit(x))-1)/length(na.omit(x)))
}

# Return difference in means from two-tailed t-test
diff_in_means <- function(data, covariates, control_i, treat_i, dim_only = FALSE){
  out <- sapply(covariates, function(x) {
  tryCatch({
    dim <- mean(data[treat_i, x], na.rm = TRUE) - mean(data[control_i, x], na.rm = TRUE)
    if(!dim_only){
      tres <- t.test(data[treat_i, x], data[control_i, x], na.rm = TRUE)
      return(c(dim = unname(dim), ttest_p = tres$p.value))
    } else {
      return(c(dim = unname(dim)))
    }
  }, error = function(e) {
    message(sprintf("Error in t.test for covariate '%s': %s", x, e$message))
    return(c(dim = NA, ttest_p = NA))
  })
  })
  if(!dim_only) out <- as.data.frame(t(out), row.names = covariates)
  else{
    out <- data.frame(dim = unname(out), row.names = covariates)
  }
  return(out)
}


capture_warnings_apply <- function(X, MARGIN, FUN, ...) {
  warnings <- character()  # store warning messages
  errors <- character()

  # wrapper around FUN to collect warnings
  wrapper <- function(...) {
    withCallingHandlers(
      tryCatch(
        FUN(...),
        error = function(e) {
          errors <<- c(errors, conditionMessage(e))
          return(NULL)  # return NA on error
        }
      ),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        # if (immediate_warning) {
        #   message("Warning: ", conditionMessage(w))
        # }
        invokeRestart("muffleWarning")
      }
    )
  }


  result <- apply(X, MARGIN, wrapper, ...)

  unique_warnings <- unique(warnings)
  unique_errors <- unique(errors)

  if (length(unique_warnings) > 0) {
    warning("Warnings:\n", paste(unique_warnings, collapse = "\n"))
  }
  if (length(unique_errors) > 0) {
    warning("Errors:\n", paste(unique_errors, collapse = "\n"))
  }

  result
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
pw <- function(data, covariates, treatment, outcome){

  # listwise deletion of observations with missing values
  z0 <- data[[treatment]] == 0

  # standardize control data for prognosis regression
  # data_c <- data_pw %>% dplyr::filter(z0) %>%
  #   dplyr::mutate_at(.vars = c(outcome, covariates),
  #                    .funs = stdr) %>% as.data.frame()

  # check if variance = 0 and return error
  temp_data <- data %>%
    dplyr::select(tidyselect::all_of(c(covariates, outcome)))
  check_var <- apply(temp_data, 2, stats::var, na.rm = TRUE)
  var_na <- names(check_var[is.na(check_var)])
  if(length(var_na)>=1L) stop(paste0("The following variables are constant in the control group among complete cases, so cannot be standardized: ",
                                     paste0(var_na, collapse = ", "),
                                     ". Consider an alternative, for example, excluding the covariate(s)."))
  # calculate prognosis weights
  X <- as.matrix(data[z0,covariates], ncol = length(covariates))
  Y <- as.matrix(data[z0,outcome], ncol = 1)
  weights <- stats::coef(stats::lm(Y ~ X))[-1] # remove intercept (0 in expectation)

  names(weights) <- covariates
  return(weights)
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
pw_rdd <- function(data, covariates, treatment, outcome, standardize = TRUE, simulation = FALSE){

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
  data_uwc <- data_uw[stats::complete.cases(data_uw[,covariates]),]
  Nc <- nrow(data_uwc) # changed from data_uw
  z0c <- data_uwc[[treatment]] == 0
  n1c <- sum(!z0c, na.rm = TRUE)
  n0c <- sum(z0c, na.rm = TRUE)

  # analytic standard error of unweighted delta
  if(length(covariates)>1L){

    varcov <- stats::cov(data_uwc[,covariates], use = "everything")
    sum_sigma2 <- sum(diag(varcov)*(nrow(data_uwc)-1)/nrow(data_uwc))
    sum_covs <- sum(varcov[upper.tri(varcov)]*(nrow(data_uwc)-1)/nrow(data_uwc))
    multiplier <- ((Nc^2)/(Nc-1))/(n1c*n0c)

    # analytic SE
    uwdelta_se <- sqrt(multiplier*sum(sum_sigma2, 2*sum_covs))

  } else { # when only one covariate
    uwdelta_se <- sqrt(stats::var(data_uwc[!z0c,covariates])/n1c +
                         stats::var(data_uwc[z0c,covariates])/n0c)
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

  pweights <- pw(data, covariates, treatment, outcome)
  pwdelta_j <- pweights*DIM
  # names(pwdelta_j) <- paste0("pw_", covariates)

  # prognosis weighted delta
  pwdelta <- sum(pweights*DIM, na.rm = TRUE)

  # R-squared from prognosis regression
  prog_mod_f <- paste(c(outcome, paste(covariates, collapse = " + ")), collapse = " ~ ")
  prog_mod <- stats::lm(formula = prog_mod_f, data = data[data[[treatment]] == 0, ])
  prog_Rsq <- summary(prog_mod)$r.squared

  # R-squared from balance regression
  bal_mod_f <- paste(c(treatment, paste(covariates, collapse = " + ")), collapse = " ~ ")
  bal_mod <- stats::lm(formula = bal_mod_f, data = data)
  bal_Rsq <- summary(bal_mod)$r.squared

  return(list(pw = pweights,
              pwdelta_j = pwdelta_j,
              pwdelta = pwdelta,
              prog_Rsq = prog_Rsq,
              bal_Rsq = bal_Rsq))

}

# pwtest helpers

#' Generate bootstrap samples
#' @param control_indices integer vector. Control-group row indices in the data.
#' @param treatment_n integer. Number of units assigned to treatment in data.
#' @param bootstrap_n integer. Number of bootstrap draws.
get_bootstrap_samples <- function(control_indices, treatment_n, bootstrap_n) {
  replicate(
    bootstrap_n,
    c(
      sample(control_indices, treatment_n, replace = TRUE),
      sample(control_indices, length(control_indices), replace = TRUE)
    )
  )
}

# Safely calculate delta metrics
safe_delta <- function(func, data, ref_obj, ...) {
  tryCatch(
    func(data, ...),
    error = function(e) {
      out <- vector(mode = "list", length = length(ref_obj))
      names(out) <- names(ref_obj)
      return(out)
      message("Error with delta estimation of bootstrap sample.")
    }
  )
}

# Standardize data relative to control group
std_data <- function(data, variables, treatment_col) {
  data %>%
    dplyr::mutate_at(
      .vars = variables,
      .funs = function(x) {
        x_c <- x[data[[treatment_col]] == 0]
        sd_x <- stats::sd(x_c, na.rm = TRUE) * (sum(!is.na(x_c)) - 1) / sum(!is.na(x_c))
        (x - mean(x_c, na.rm = TRUE)) / sd_x
      }
    ) %>%
    as.data.frame()
}

get_pvalue_uw <- function(sim_d, ref_d, bootstrap_n){
  sum(abs(sim_d) >= abs(ref_d)) / bootstrap_n
}

get_pvalue_pw <- function(sim_d, ref_d, bootstrap_n){
  sum(abs(unlist(sapply(sim_d, function(x) x[["pwdelta"]]))) >= abs(ref_d)) / bootstrap_n
}

# Calculate standard error adjustments
adjust_se <- function(values, nsims) {
  stats::sd(values, na.rm = TRUE) * (nsims - 1) / nsims
}

#' ML model selection and cross-validation wrapper
#' @param data data.frame containing control group observations
#' @param covariates character vector of covariate names
#' @param outcome character name of outcome variable
#' @param cv_folds integer. Number of cross-validation folds (default 5)
#' @param min_penalty_exp minimum penalty exponent (default -6)
#' @param max_penalty_exp maximum penalty exponent (default -1)
#' @param verbose logical. Whether to print detailed contest report
#' @export
contest <- function(data, covariates, outcome, cv_folds = 5,
                    min_penalty_exp = -6, max_penalty_exp = -1,
                    verbose = TRUE, subset_workflow = NULL) {

  if(!is.null(subset_workflow) &
     any(!subset_workflow %in% c("linear", "linear_poly", "linear_interact", "linear_poly_interact",
                                 "lasso_poly", "lasso_interact", "lasso_poly_interact",
                                 "random_forest", "gradient_boosting"))) stop("One or more values of `subset_workflow` invalid. See Details in ??contest.")
  
  if(verbose) cat("=== CONTEST: Model Selection Pipeline ===\n")

  # Step 1: Tune hyperparameters via CV
  tuned_models <- tune_hyperparams(data = data,
                                   outcome = outcome,
                                   covariates = covariates,
                                   cv_folds = cv_folds,
                                   min_penalty_exp = min_penalty_exp,
                                   max_penalty_exp = max_penalty_exp,
                                   verbose = verbose,
                                   subset_workflow = subset_workflow)

  # Step 2: Select best model based on full data performance
  best_model <- pick_winner(tuned_models = tuned_models,
                            data = data,
                            outcome = outcome,
                            covariates = covariates,
                            formula = NULL,
                            verbose = verbose)

  if(verbose) {
    cat("\n=== WINNER ===\n")
    cat("Model:", best_model$best_model_name, "\n")
    cat("CV R-squared:", round(best_model$cv_rsq, 4), "\n")
    cat("Full data R-squared:", round(best_model$full_rsq, 4), "\n")
  }

  # Return in format compatible with existing pwtest code
  contest_results <- list(
    # Core model components for pwtest integration
    best_spec = best_model$best_spec,
    best_engine = best_model$best_engine,
    best_recipe = best_model$best_recipe,

    # Model identification and performance
    best_model_name = best_model$best_model_name,
    best_cv_rsq = best_model$cv_rsq,
    best_full_data_rsq = best_model$full_rsq,

    # Hyperparameters
    best_hyperparameters = best_model$best_hyperparams,

    # Add contest class for compatibility
    class = c("contest_results", "list")
  )

  return(contest_results)
}

#' Tune hyperparameters via cross-validation
#' @param data data.frame containing control group observations
#' @param covariates character vector of covariate names
#' @param outcome character name of outcome variable
#' @param cv_folds integer. Number of cross-validation folds (default 5)
#' @param min_penalty_exp minimum penalty exponent
#' @param max_penalty_exp maximum penalty exponent
#' @param verbose logical. Whether to print detailed information
#' @importFrom rlang .data
#' @importFrom tibble tibble
#' @importFrom dplyr select mutate filter arrange desc semi_join pull
#' @importFrom tidyselect all_of
#' @importFrom tidyr drop_na
#' @importFrom purrr map map_dbl map_lgl map2 map_chr
#' @importFrom stats var as.formula predict
#' @importFrom rsample vfold_cv assessment analysis
#' @importFrom recipes recipe step_normalize step_poly step_corr step_interact all_predictors prep bake
#' @importFrom parsnip linear_reg set_engine set_mode rand_forest boost_tree fit
#' @importFrom workflows workflow add_recipe add_model extract_spec_parsnip extract_preprocessor
#' @importFrom tune tune fit_resamples tune_grid finalize_workflow select_best collect_metrics
#' @importFrom dials grid_regular penalty mtry trees learn_rate
#' @importFrom scales log10_trans
#' @importFrom yardstick metric_set rmse rsq
tune_hyperparams <- function(data, outcome, covariates, cv_folds = 5,
                             min_penalty_exp = -6, max_penalty_exp = -1,
                             verbose = TRUE, subset_workflow = NULL) {

  if(verbose) cat("\n-- Step 1: Hyperparameter Tuning via Cross-Validation --\n")

  # Data diagnostics
  analysis_data <- data %>% tidyr::drop_na()

  if(verbose) {
    cat("Data dimensions:", nrow(analysis_data), "x", ncol(analysis_data), "\n")
    cat("Outcome variance:", stats::var(analysis_data[[outcome]], na.rm = TRUE), "\n")
  }

  if (nrow(analysis_data) < cv_folds) {
    stop("Insufficient data for cross-validation. Need at least ", cv_folds, " observations, got ", nrow(analysis_data))
  }

  # Cross-validation setup
  cv_folds_obj <- rsample::vfold_cv(analysis_data, v = cv_folds, strata = NULL)

  # Recipe creation
  base_recipe <- recipes::recipe(stats::as.formula(paste(outcome, "~ .")), data = analysis_data)

  poly_recipe <- base_recipe %>%
    recipes::step_normalize(recipes::all_predictors()) %>%
    recipes::step_poly(recipes::all_predictors(), degree = 2) %>%
    recipes::step_corr(recipes::all_predictors(), threshold = 0.99)

  interact_recipe <- base_recipe %>%
    recipes::step_normalize(recipes::all_predictors()) %>%
    recipes::step_interact(terms = ~ recipes::all_predictors():recipes::all_predictors()) %>%
    recipes::step_corr(recipes::all_predictors(), threshold = 0.99)

  poly_interact_recipe <- base_recipe %>%
    recipes::step_normalize(recipes::all_predictors()) %>%
    recipes::step_poly(recipes::all_predictors(), degree = 2) %>%
    recipes::step_interact(terms = ~ recipes::all_predictors():recipes::all_predictors()) %>%
    recipes::step_corr(recipes::all_predictors(), threshold = 0.99)

  # Model specifications
  linear_spec <- parsnip::linear_reg() %>% parsnip::set_engine("lm")

  lasso_grid <- dials::grid_regular(
    dials::penalty(range = c(min_penalty_exp, max_penalty_exp), trans = scales::log10_trans()),
    levels = 15
  )

  lasso_spec <- parsnip::linear_reg(penalty = tune::tune(), mixture = 1) %>%
    parsnip::set_engine("glmnet")

  rf_grid <- dials::grid_regular(
    dials::mtry(range = c(1, min(length(covariates), 5))),
    dials::trees(range = c(50, 200)),
    levels = 3
  )

  gbm_grid <- dials::grid_regular(
    dials::mtry(range = c(1, min(length(covariates), 5))),
    dials::trees(range = c(50, 200)),
    dials::learn_rate(range = c(0.01, 0.3)),
    levels = 3
  )

  rf_spec <- parsnip::rand_forest(mtry = tune::tune(), trees = tune::tune()) %>%
    parsnip::set_engine("ranger") %>%
    parsnip::set_mode("regression")

  gbm_spec <- parsnip::boost_tree(mtry = tune::tune(), trees = tune::tune(), learn_rate = tune::tune()) %>%
    parsnip::set_engine("xgboost") %>%
    parsnip::set_mode("regression")

  # Create workflows
  targeted_workflows <- tibble::tibble(
    wflow_id = c("linear", "linear_poly", "linear_interact", "linear_poly_interact",
                 "lasso_poly", "lasso_interact", "lasso_poly_interact",
                 "random_forest", "gradient_boosting"),
    info = list(
      list(workflow = list(workflows::workflow() %>% workflows::add_recipe(base_recipe) %>% workflows::add_model(linear_spec))),
      list(workflow = list(workflows::workflow() %>% workflows::add_recipe(poly_recipe) %>% workflows::add_model(linear_spec))),
      list(workflow = list(workflows::workflow() %>% workflows::add_recipe(interact_recipe) %>% workflows::add_model(linear_spec))),
      list(workflow = list(workflows::workflow() %>% workflows::add_recipe(poly_interact_recipe) %>% workflows::add_model(linear_spec))),
      list(workflow = list(workflows::workflow() %>% workflows::add_recipe(poly_recipe) %>% workflows::add_model(lasso_spec))),
      list(workflow = list(workflows::workflow() %>% workflows::add_recipe(interact_recipe) %>% workflows::add_model(lasso_spec))),
      list(workflow = list(workflows::workflow() %>% workflows::add_recipe(poly_interact_recipe) %>% workflows::add_model(lasso_spec))),
      list(workflow = list(workflows::workflow() %>% workflows::add_recipe(base_recipe) %>% workflows::add_model(rf_spec))),
      list(workflow = list(workflows::workflow() %>% workflows::add_recipe(base_recipe) %>% workflows::add_model(gbm_spec)))
    )
  )

  robust_metrics <- yardstick::metric_set(yardstick::rmse, yardstick::rsq)

  # Fit all models with CV
  if(verbose) cat("\n-- Step 2: Fitting Models with Cross-Validation --\n")

  # Mid-step if user specifies subset of workflows
  if(!is.null(subset_workflow)){
    targeted_workflows <- targeted_workflows[targeted_workflows$wflow_id %in% subset_workflow,]
  }


  # Validate all workflows were created properly
  if(verbose) {
    cat("Validating workflow creation...\n")
    for(i in 1:nrow(targeted_workflows)) {
      wf_name <- targeted_workflows$wflow_id[i]
      wf_obj <- targeted_workflows$info[[i]]$workflow[[1]]
      tryCatch({
        model_spec <- workflows::extract_spec_parsnip(wf_obj)
        cat("  ", wf_name, "- Model:", class(model_spec)[1], "check\n")
      }, error = function(e) {
        cat("  ", wf_name, "- ERROR in workflow creation:", e$message, "\n")
      })
    }
  }
  
  safe_fit_model <- function(wf_info, wf_name, lasso_grid, rf_grid, gbm_grid, cv_folds_obj, robust_metrics, verbose) {
    wf <- wf_info$workflow[[1]]

    if(verbose) cat("Fitting model:", wf_name, "\n")

    result <- tryCatch({
      if (wf_name %in% c("linear", "linear_poly", "linear_interact", "linear_poly_interact")) {
        tune::fit_resamples(wf, resamples = cv_folds_obj, metrics = robust_metrics)
      } else if (wf_name %in% c("lasso_poly", "lasso_interact", "lasso_poly_interact")) {
        tune::tune_grid(wf, resamples = cv_folds_obj, grid = lasso_grid, metrics = robust_metrics)
      } else if (wf_name == "random_forest") {
        tune::tune_grid(wf, resamples = cv_folds_obj, grid = rf_grid, metrics = robust_metrics)
      } else if (wf_name == "gradient_boosting") {
        tune::tune_grid(wf, resamples = cv_folds_obj, grid = gbm_grid, metrics = robust_metrics)
      }
    }, error = function(e) {
      if(verbose) cat("  ERROR:", e$message, "\n")
      return(NULL)
    })

    return(result)
  }

  results <- targeted_workflows %>%
    dplyr::mutate(
      fit_results = purrr::map2(.data$info, .data$wflow_id, ~{
        safe_fit_model(.x, .y, lasso_grid, rf_grid, gbm_grid, cv_folds_obj, robust_metrics, verbose)
      })
    )

  # Extract hyperparameters and create tuned model specifications
  tuned_models <- results %>%
    dplyr::mutate(
      tuned_info = purrr::map2(.data$wflow_id, .data$fit_results, function(id, fit_result) {
        if (is.null(fit_result)) return(NULL)

        # Get workflow
        wf <- targeted_workflows %>%
          dplyr::filter(wflow_id == id) %>%
          dplyr::pull(info) %>%
          .[[1]] %>%
          .$workflow[[1]]

        # Extract metrics
        metrics <- tune::collect_metrics(fit_result)

        if (id %in% c("lasso_poly", "lasso_interact", "lasso_poly_interact",
                      "random_forest", "gradient_boosting")) {
          # Get best hyperparameters
          best_params <- tune::select_best(fit_result, metric = "rsq")

          # Finalize workflow with best hyperparameters
          final_wf <- tune::finalize_workflow(wf, best_params)

          # Extract components
          best_spec <- workflows::extract_spec_parsnip(final_wf)
          best_recipe <- workflows::extract_preprocessor(final_wf)

          cv_rsq <- metrics %>%
            dplyr::semi_join(best_params, by = intersect(names(metrics), names(best_params))) %>%
            dplyr::filter(.data$.metric == "rsq") %>%
            dplyr::pull(.data$mean)

          return(list(
            spec = best_spec,
            recipe = best_recipe,
            engine = best_spec$engine,
            hyperparams = best_params,
            cv_rsq = cv_rsq,
            wflow_id = id
          ))
        } else {
          # No hyperparameters to tune
          best_spec <- workflows::extract_spec_parsnip(wf)
          best_recipe <- workflows::extract_preprocessor(wf)

          cv_rsq <- metrics %>%
            dplyr::filter(.data$.metric == "rsq") %>%
            dplyr::pull(.data$mean)

          return(list(
            spec = best_spec,
            recipe = best_recipe,
            engine = best_spec$engine,
            hyperparams = NULL,
            cv_rsq = cv_rsq,
            wflow_id = id
          ))
        }
      })
    )

  # Return all tuned models
  return(list(
    tuned_models = tuned_models,
    grids = list(lasso = lasso_grid, rf = rf_grid, gbm = gbm_grid),
    cv_folds = cv_folds
  ))
}

#' Select best model based on control data (control data passed by pwtest)
#' @param tuned_models output from tune_hyperparams
#' @param data full control group data
#' @param outcome outcome variable name
#' @param covariates covariate names
#' @param formula model formula (optional)
#' @param verbose whether to print diagnostics
pick_winner <- function(tuned_models, data, outcome, covariates, formula = NULL, verbose = TRUE) {

  if(is.null(formula)) formula <- as.formula(paste0(outcome, "~."))

  if(verbose) cat("\n-- Step 3: Evaluating Models on Full Control Data --\n")

  # Prepare data
  analysis_data <- data %>% tidyr::drop_na()

  # Evaluate each model on full dataset
  performance <- tuned_models$tuned_models %>%
    dplyr::filter(!sapply(.data$tuned_info, is.null)) %>%
    dplyr::mutate(
      full_fit = purrr::map(.data$tuned_info, function(model_info) {
        if(verbose) cat("Evaluating", model_info$wflow_id, "on full data... ")

        tryCatch({
          # Create workflow with recipe and model
          wf <- workflows::workflow() %>%
            workflows::add_recipe(model_info$recipe) %>%
            workflows::add_model(model_info$spec)

          # Fit on full data
          fit <- parsnip::fit(wf, data = analysis_data)

          if(verbose) cat("Success\n")
          return(fit)
        }, error = function(e) {
          if(verbose) cat("ERROR:", e$message, "\n")
          return(NULL)
        })
      }),

      full_rsq = purrr::map_dbl(.data$full_fit, function(fit) {
        if(is.null(fit)) return(0)

        tryCatch({
          preds <- predict(fit, analysis_data)
          actual <- analysis_data[[outcome]]
          yardstick::rsq_vec(actual, preds$.pred)
        }, error = function(e) 0)
      })
    )

  # Select best based on full data R²
  best_idx <- which.max(performance$full_rsq)

  if(verbose) {
    cat("\nModel performance comparison:\n")
    performance %>%
      dplyr::mutate(cv_rsq = sapply(.data$tuned_info, function(x) x$cv_rsq)) %>%
      dplyr::select(.data$wflow_id, .data$cv_rsq, .data$full_rsq) %>%
      dplyr::arrange(dplyr::desc(.data$full_rsq)) %>%
      print()
  }

  best_model_info <- performance$tuned_info[[best_idx]]

  return(list(
    best_spec = best_model_info$spec,
    best_engine = best_model_info$engine,
    best_recipe = best_model_info$recipe,
    best_hyperparams = best_model_info$hyperparams,
    best_model_name = best_model_info$wflow_id,
    cv_rsq = best_model_info$cv_rsq,
    full_rsq = performance$full_rsq[best_idx]
  ))
}

#' Validate contest results for prognostic_balance integration
#' @param contest_output Output from contest() function
validate_contest_output <- function(contest_output) {
  required_components <- c("best_spec", "best_engine", "best_model_name",
                           "best_cv_rsq", "best_full_data_rsq")

  missing_components <- setdiff(required_components, names(contest_output))

  if (length(missing_components) > 0) {
    stop("Contest output missing required components: ",
         paste(missing_components, collapse = ", "))
  }

  # Additional validation
  if (is.null(contest_output$best_spec)) {
    stop("Contest failed to identify a valid model specification")
  }

  return(TRUE)
}
