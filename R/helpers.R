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

# pwtest helpers

# Generate bootstrap samples
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
    mutate_at(
      .vars = variables,
      .funs = function(x) {
        x_c <- x[data[[treatment_col]] == 0]
        sd_x <- sd(x_c, na.rm = TRUE) * (sum(!is.na(x_c)) - 1) / sum(!is.na(x_c))
        (x - mean(x_c, na.rm = TRUE)) / sd_x
      }
    ) %>%
    as.data.frame()
}

get_pvalue_uw <- function(sim_d, ref_d, bootstrap_n){
  sum(abs(unlist(sapply(sim_d, function(x) x[["uwdelta"]]))) >= abs(ref_d)) / bootstrap_n
}

get_pvalue_pw <- function(sim_d, ref_d, bootstrap_n){
  sum(abs(unlist(sapply(sim_d, function(x) x[["pwdelta"]]))) >= abs(ref_d)) / bootstrap_n
}

# Calculate standard error adjustments
adjust_se <- function(values, nsims) {
  sd(values, na.rm = TRUE) * (nsims - 1) / nsims
}

#' ML model selection and cross-validation for prognosis weighting with comprehensive debugging
#' @param data_control data.frame containing control group observations
#' @param covariates character vector of covariate names
#' @param outcome character name of outcome variable
#' @param cv_folds integer number of cross-validation folds (default 5, increased from 3)
#' @param debug logical whether to print detailed debugging information
#' @param min_penalty_exp minimum penalty exponent (default -6 for more conservative start)
#' @param max_penalty_exp maximum penalty exponent (default -1 for small datasets)
                        
contest <- function(data_control, covariates, outcome, cv_folds = 5, debug = TRUE, 
                    min_penalty_exp = -6, max_penalty_exp = -1) {
  
  # Load required libraries
  require(tidymodels)
  require(dplyr)
  require(tidyselect)
  require(purrr)
  require(tidyr)
  require(yardstick)
  require(dials)
  require(rsample)
  require(scales)
  
  if(debug) cat("=== DEBUGGING ===\n")
  
  # Input data diagnostics
  analysis_data <- data_control %>%
    dplyr::select(tidyselect::all_of(c(outcome, covariates))) %>%
    tidyr::drop_na()
  
  if(debug) {
    cat("\n--- Input data diagnostics ---\n")
    cat("Original data dimensions:", nrow(data_control), "x", ncol(data_control), "\n")
    cat("Analysis data dimensions:", nrow(analysis_data), "x", ncol(analysis_data), "\n")
    cat("Outcome variable:", outcome, "\n")
    cat("Outcome variance:", var(analysis_data[[outcome]], na.rm = TRUE), "\n")
    cat("Outcome range:", paste(range(analysis_data[[outcome]], na.rm = TRUE), collapse = " to "), "\n")
    
    cat("\nCovariate summary:\n")
    cov_summary <- analysis_data[, covariates] %>%
      summarise_all(list(
        n_unique = ~length(unique(.)),
        variance = ~var(., na.rm = TRUE),
        range_width = ~diff(range(., na.rm = TRUE))
      ))
    print(cov_summary)
    
    # Check for problematic variables
    constant_vars <- sapply(analysis_data[, covariates], function(x) length(unique(x)) == 1)
    if(any(constant_vars)) {
      cat("WARNING: Constant variables detected:", paste(names(constant_vars)[constant_vars], collapse = ", "), "\n")
    }
  }
  
  if (nrow(analysis_data) < cv_folds) {
    stop("Insufficient data for cross-validation. Need at least ", cv_folds, " observations, got ", nrow(analysis_data))
  }
  
  # CV Diagnostics
  cv_folds_obj <- vfold_cv(analysis_data, v = cv_folds, strata = NULL)
  
  if(debug) {
    cat("\n--- CV Diagnostics ---\n")
    fold_sizes <- purrr::map_dbl(cv_folds_obj$splits, function(split) nrow(rsample::assessment(split)))
    fold_train_sizes <- purrr::map_dbl(cv_folds_obj$splits, function(split) nrow(rsample::analysis(split)))
    
    cat("CV folds:", cv_folds, "\n")
    cat("Assessment (validation) fold sizes:", paste(fold_sizes, collapse = ", "), "\n")
    cat("Training fold sizes:", paste(fold_train_sizes, collapse = ", "), "\n")
    cat("Min assessment fold size:", min(fold_sizes), "\n")
    
    if (min(fold_sizes) < 10) {
      cat("WARNING: Very small CV folds detected. Consider reducing cv_folds.\n")
    }
    
    # Test outcome variance in folds
    fold_outcome_vars <- purrr::map_dbl(cv_folds_obj$splits, function(split) {
      fold_data <- rsample::analysis(split)
      var(fold_data[[outcome]], na.rm = TRUE)
    })
    cat("Outcome variance in training folds:", paste(round(fold_outcome_vars, 4), collapse = ", "), "\n")
  }
  
  # Recipe creation and diagnostics
  base_recipe <- recipe(as.formula(paste(outcome, "~ .")), data = analysis_data)
  
  # Create more conservative recipes with normalization
  poly_recipe <- base_recipe %>%
    step_normalize(all_predictors()) %>%
    step_poly(all_predictors(), degree = 2) %>%
    step_corr(all_predictors(), threshold = 0.95)
  
  interact_recipe <- base_recipe %>%
    step_normalize(all_predictors()) %>%
    step_interact(terms = ~ all_predictors():all_predictors()) %>%
    step_corr(all_predictors(), threshold = 0.95)
  
  poly_interact_recipe <- base_recipe %>%
    step_normalize(all_predictors()) %>%
    step_poly(all_predictors(), degree = 2) %>%
    step_interact(terms = ~ all_predictors():all_predictors()) %>%
    step_corr(all_predictors(), threshold = 0.95)
  
  if(debug) {
    cat("\n--- Recipe Diagnostics ---\n")
    
    # Test recipe preprocessing
    tryCatch({
      poly_prepped <- prep(poly_recipe, training = analysis_data)
      poly_baked <- bake(poly_prepped, new_data = NULL)
      cat("Base features:", length(covariates), "\n")
      cat("Features after poly recipe:", ncol(poly_baked) - 1, "\n")
      
      # Check for remaining issues
      feature_vars <- sapply(poly_baked[, -which(names(poly_baked) == outcome)], var, na.rm = TRUE)
      n_zero_var <- sum(feature_vars == 0, na.rm = TRUE)
      if(n_zero_var > 0) {
        cat("WARNING: Features with zero variance after preprocessing:", n_zero_var, "\n")
      }
      
      # Check correlation matrix
      cor_matrix <- cor(poly_baked[, -which(names(poly_baked) == outcome)], use = "complete.obs")
      max_cor <- max(abs(cor_matrix[upper.tri(cor_matrix)]), na.rm = TRUE)
      cat("Max correlation after step_corr:", round(max_cor, 4), "\n")
      
      # Test other recipes
      interact_prepped <- prep(interact_recipe, training = analysis_data)
      interact_baked <- bake(interact_prepped, new_data = NULL)
      cat("Features after interact recipe:", ncol(interact_baked) - 1, "\n")
      
      poly_interact_prepped <- prep(poly_interact_recipe, training = analysis_data)
      poly_interact_baked <- bake(poly_interact_prepped, new_data = NULL)
      cat("Features after poly_interact recipe:", ncol(poly_interact_baked) - 1, "\n")
      
    }, error = function(e) {
      cat("ERROR in recipe preprocessing:", e$message, "\n")
    })
  }
  
  # Model specifications
  linear_spec <- linear_reg() %>% set_engine("lm")
  
  # Create more conservative lasso penalty grid based on dataset size
  n_obs <- nrow(analysis_data)
  n_features_base <- length(covariates)
  
  if(debug) {
    cat("\n--- Lasso tuning ---\n")
    cat("Sample size:", n_obs, "\n")
    cat("Base features:", n_features_base, "\n")
  }
  
  # More conservative penalty grid
  lasso_grid <- grid_regular(
    dials::penalty(range = c(min_penalty_exp, max_penalty_exp), trans = scales::log10_trans()), 
    levels = 15
  )
  
  if(debug) {
    cat("Lasso penalty range: 10^", min_penalty_exp, " to 10^", max_penalty_exp, "\n")
    cat("Penalty grid points:", nrow(lasso_grid), "\n")
    cat("Sample penalties:", paste(signif(lasso_grid$penalty[1:5], 3), collapse = ", "), "...\n")
    cat("Penalty range:", paste(signif(range(lasso_grid$penalty), 3), collapse = " to "), "\n")
  }
  
  lasso_spec <- linear_reg(penalty = tune(), mixture = 1) %>% set_engine("glmnet")
  
  # Tree model grids
  rf_grid <- grid_regular(
    dials::mtry(range = c(1, min(length(covariates), 5))), 
    dials::trees(range = c(50, 200)), 
    levels = 3
  )
  
  gbm_grid <- grid_regular(
    dials::mtry(range = c(1, min(length(covariates), 5))),
    dials::trees(range = c(50, 200)),
    dials::learn_rate(range = c(0.01, 0.3)), 
    levels = 3
  )
  
  rf_spec <- rand_forest(mtry = tune(), trees = tune()) %>% 
    set_engine("ranger") %>% 
    set_mode("regression")
  
  gbm_spec <- boost_tree(mtry = tune(), trees = tune(), learn_rate = tune()) %>%
    set_engine("xgboost") %>%
    set_mode("regression")
  
  # Workflow creation
  targeted_workflows <- tibble(
    wflow_id = c("linear", "linear_poly", "linear_interact", "linear_poly_interact",
                 "lasso_poly", "lasso_interact", "lasso_poly_interact", 
                 "random_forest", "gradient_boosting"),
    info = list(
      list(workflow = list(workflow() %>% add_recipe(base_recipe) %>% add_model(linear_spec))),
      list(workflow = list(workflow() %>% add_recipe(poly_recipe) %>% add_model(linear_spec))),
      list(workflow = list(workflow() %>% add_recipe(interact_recipe) %>% add_model(linear_spec))),
      list(workflow = list(workflow() %>% add_recipe(poly_interact_recipe) %>% add_model(linear_spec))),
      list(workflow = list(workflow() %>% add_recipe(poly_recipe) %>% add_model(lasso_spec))),
      list(workflow = list(workflow() %>% add_recipe(interact_recipe) %>% add_model(lasso_spec))),
      list(workflow = list(workflow() %>% add_recipe(poly_interact_recipe) %>% add_model(lasso_spec))),
      list(workflow = list(workflow() %>% add_recipe(base_recipe) %>% add_model(rf_spec))),
      list(workflow = list(workflow() %>% add_recipe(base_recipe) %>% add_model(gbm_spec)))
    )
  )
  
  # Use robust metrics
  robust_metrics <- metric_set(rmse, rsq)
  
  # Model fitting
  if(debug) {
    cat("\n--- Model fitting ---\n")
    cat("Models to fit:", paste(targeted_workflows$wflow_id, collapse = ", "), "\n")
    cat("Number of workflows created:", nrow(targeted_workflows), "\n")
  }
  
  # Validate all workflows were created properly
  if(debug) {
    cat("Validating workflow creation...\n")
    for(i in 1:nrow(targeted_workflows)) {
      wf_name <- targeted_workflows$wflow_id[i]
      wf_obj <- targeted_workflows$info[[i]]$workflow[[1]]
      tryCatch({
        model_spec <- extract_spec_parsnip(wf_obj)
        cat("  ", wf_name, "- Model:", class(model_spec)[1], "✓\n")
      }, error = function(e) {
        cat("  ", wf_name, "- ERROR in workflow creation:", e$message, "\n")
      })
    }
  }
  
  # Enhanced fitting function with detailed error capture
  safe_fit_model <- function(wf_info, wf_name, lasso_grid, rf_grid, gbm_grid, cv_folds_obj, robust_metrics, debug) {
    wf <- wf_info$workflow[[1]]
    
    if(debug) cat("Fitting model:", wf_name, "\n")
    
    # Track warnings and errors
    warnings_captured <- character(0)
    errors_captured <- character(0)
    
    result <- tryCatch({
      withCallingHandlers({
        if (wf_name %in% c("linear", "linear_poly", "linear_interact", "linear_poly_interact")) {
          # No tuning needed
          fit_resamples(wf, resamples = cv_folds_obj, metrics = robust_metrics)
        } else if (wf_name %in% c("lasso_poly", "lasso_interact", "lasso_poly_interact")) {
          # Test single penalty first
          if(debug) {
            cat("  Testing single penalty value first...\n")
            test_wf <- wf %>% finalize_workflow(tibble(penalty = 0.01))
            test_fit <- tryCatch({
              fit(test_wf, data = analysis_data)
            }, error = function(e) {
              cat("  ERROR in single penalty test:", e$message, "\n")
              return(NULL)
            })
            
            if(!is.null(test_fit)) {
              test_preds <- predict(test_fit, analysis_data)
              pred_var <- var(test_preds$.pred, na.rm = TRUE)
              pred_range <- range(test_preds$.pred, na.rm = TRUE)
              cat("  Single penalty test - Prediction variance:", round(pred_var, 6), 
                  "Range:", paste(round(pred_range, 4), collapse = " to "), "\n")
              
              if(pred_var < 1e-10) {
                cat("  WARNING: Near-constant predictions detected\n")
              }
            }
          }
          
          # Tune with grid
          tune_grid(wf, resamples = cv_folds_obj, grid = lasso_grid, metrics = robust_metrics)
        } else if (wf_name == "random_forest") {
          tune_grid(wf, resamples = cv_folds_obj, grid = rf_grid, metrics = robust_metrics)
        } else if (wf_name == "gradient_boosting") {
          tune_grid(wf, resamples = cv_folds_obj, grid = gbm_grid, metrics = robust_metrics)
        }
      }, warning = function(w) {
        warning_msg <- w$message
        warnings_captured <<- c(warnings_captured, warning_msg)
        
        if(debug) {
          if (grepl("correlation computation is required", warning_msg)) {
            cat("  WARNING: Constant predictions detected (over-regularized)\n")
          } else if (grepl("prediction from rank-deficient fit", warning_msg)) {
            cat("  WARNING: Rank deficient fit (feature explosion)\n")
          } else {
            cat("  WARNING:", warning_msg, "\n")
          }
        }
        invokeRestart("muffleWarning")
      })
    }, error = function(e) {
      errors_captured <<- c(errors_captured, e$message)
      if(debug) cat("  ERROR:", e$message, "\n")
      return(NULL)
    })
    
    # Attach diagnostic info
    if(!is.null(result)) {
      attr(result, "warnings") <- warnings_captured
      attr(result, "errors") <- errors_captured
    }
    
    if(debug) {
      if(is.null(result)) {
        cat("  FAILED:", wf_name, "\n")
      } else {
        cat("  SUCCESS:", wf_name, "- Warnings:", length(warnings_captured), "Errors:", length(errors_captured), "\n")
      }
    }
    
    return(result)
  }
  
  # Fit all models
  results <- targeted_workflows %>%
    dplyr::mutate(
      fit_results = purrr::map2(info, wflow_id, ~{
        if(debug) cat("Processing workflow:", .y, "\n")
        safe_fit_model(.x, .y, lasso_grid, rf_grid, gbm_grid, cv_folds_obj, robust_metrics, debug)
      })
    )
  
  if(debug) {
    cat("Completed fitting. Results summary:\n")
    results_summary <- results %>%
      dplyr::mutate(
        fit_success = purrr::map_lgl(fit_results, ~!is.null(.x))
      ) %>%
      dplyr::select(wflow_id, fit_success)
    print(results_summary)
  }
  
  # Results extraction
  if(debug) cat("\n--- Results ---\n")
  
  model_performance <- results %>%
    dplyr::mutate(
      model_warnings = purrr::map_chr(fit_results, function(x) {
        if (is.null(x)) return("failed")
        warnings <- attr(x, "warnings")
        if (is.null(warnings) || length(warnings) == 0) return("none")
        return(paste(length(warnings), "warnings"))
      }),
      
      best_rmse = purrr::map_dbl(fit_results, function(x) {
        if (is.null(x)) return(Inf)
        
        tryCatch({
          metrics <- collect_metrics(x)
          if ("penalty" %in% names(metrics) || "mtry" %in% names(metrics)) {
            # Tuned model - select best based on rsq
            best_params <- select_best(x, metric = "rsq")
            metrics %>%
              dplyr::semi_join(best_params, by = intersect(names(metrics), names(best_params))) %>%
              dplyr::filter(.metric == "rmse") %>%
              dplyr::pull(mean)
          } else {
            # Non-tuned model
            metrics %>%
              dplyr::filter(.metric == "rmse") %>%
              dplyr::pull(mean)
          }
        }, error = function(e) {
          if(debug) cat("Error extracting RMSE for model:", e$message, "\n")
          return(Inf)
        })
      }),
      
      best_rsq = purrr::map_dbl(fit_results, function(x) {
        if (is.null(x)) return(0)
        
        tryCatch({
          metrics <- collect_metrics(x)
          if ("penalty" %in% names(metrics) || "mtry" %in% names(metrics)) {
            # Tuned model
            best_params <- select_best(x, metric = "rsq")
            metrics %>%
              dplyr::semi_join(best_params, by = intersect(names(metrics), names(best_params))) %>%
              dplyr::filter(.metric == "rsq") %>%
              dplyr::pull(mean)
          } else {
            # Non-tuned model
            metrics %>%
              dplyr::filter(.metric == "rsq") %>%
              dplyr::pull(mean)
          }
        }, error = function(e) {
          if(debug) cat("Error extracting R-squared for model:", e$message, "\n")
          return(0)
        })
      })
    )
  
  if(debug) {
    cat("\nModel Performance Summary:\n")
    perf_summary <- model_performance %>%
      dplyr::select(wflow_id, best_rmse, best_rsq, model_warnings) %>%
      dplyr::arrange(desc(best_rsq))
    print(perf_summary)
  }
  
  # Select best model
  best_model_idx <- which.max(model_performance$best_rsq)
  
  if (all(model_performance$best_rsq == 0)) {
    stop("All models failed to fit. Check data for constant variables, perfect collinearity, or insufficient variation.")
  }
  
  best_model_name <- model_performance$wflow_id[best_model_idx]
  best_workflow <- model_performance$info[[best_model_idx]]$workflow[[1]]
  best_fit_results <- model_performance$fit_results[[best_model_idx]]
  
  if(debug) {
    cat("\nBest model:", best_model_name, "\n")
    cat("Best R-squared:", round(model_performance$best_rsq[best_model_idx], 4), "\n")
    cat("Best RMSE:", round(model_performance$best_rmse[best_model_idx], 4), "\n")
  }
  
  # Get final model specification
  if (best_model_name %in% c("lasso_poly", "lasso_interact", "lasso_poly_interact", 
                             "random_forest", "gradient_boosting")) {
    best_params <- select_best(best_fit_results, metric = "rsq")
    final_workflow <- finalize_workflow(best_workflow, best_params)
    best_spec <- extract_spec_parsnip(final_workflow)
    
    if(debug && "penalty" %in% names(best_params)) {
      cat("Best lasso penalty:", best_params$penalty, "\n")
    }
    if(debug && "mtry" %in% names(best_params)) {
      cat("Best mtry:", best_params$mtry, "\n")
    }
    if(debug && "learn_rate" %in% names(best_params)) {
      cat("Best learning rate:", best_params$learn_rate, "\n")
    }
  } else {
    final_workflow <- best_workflow
    best_spec <- extract_spec_parsnip(final_workflow)
  }
  
  best_engine <- best_spec$engine
  
  # Extract recipe
  best_recipe <- NULL
  if (best_model_name %in% c("linear_poly", "linear_interact", "linear_poly_interact",
                             "lasso_poly", "lasso_interact", "lasso_poly_interact")) {
    best_recipe <- tryCatch({
      extract_preprocessor(model_performance$info[[best_model_idx]]$workflow[[1]])
    }, error = function(e) {
      if(debug) cat("Could not extract recipe, recreating...\n")
      # Recreate based on model name
      if (best_model_name %in% c("linear_poly", "lasso_poly")) {
        poly_recipe
      } else if (best_model_name %in% c("linear_interact", "lasso_interact")) {
        interact_recipe
      } else if (best_model_name %in% c("linear_poly_interact", "lasso_poly_interact")) {
        poly_interact_recipe
      }
    })
  }
  
  # Get linear baseline performance for comparison
  linear_idx <- which(model_performance$wflow_id == "linear")
  if (length(linear_idx) > 0 && model_performance$best_rsq[linear_idx] != 0) {
    linear_cv_rmse <- model_performance$best_rmse[linear_idx]
    linear_cv_rsq <- model_performance$best_rsq[linear_idx]
  } else {
    # Fallback if linear model failed
    linear_cv_rmse <- model_performance$best_rmse[best_model_idx]
    linear_cv_rsq <- model_performance$best_rsq[best_model_idx]
  }
  
  if(debug) cat("\n=== Results ===\n")
  
  # Return comprehensive results object optimized for prognostic_balance() integration
  contest_results <- list(
    # Core model components for pwtest() integration
    best_spec = best_spec,
    best_engine = best_engine,
    best_recipe = best_recipe,
    
    # Model identification and performance
    best_model_name = best_model_name,
    best_cv_rmse = model_performance$best_rmse[best_model_idx],
    best_cv_rsq = model_performance$best_rsq[best_model_idx],
    
    # Linear baseline for comparison
    linear_cv_rmse = linear_cv_rmse,
    linear_cv_rsq = linear_cv_rsq,
    
    # Diagnostics and debugging
    debug_info = if(debug) {
      list(
        data_summary = list(
          n_obs = nrow(analysis_data),
          n_features = length(covariates),
          outcome_var = var(analysis_data[[outcome]], na.rm = TRUE)
        ),
        cv_info = list(
          n_folds = cv_folds,
          min_fold_size = min(purrr::map_dbl(cv_folds_obj$splits, ~nrow(rsample::assessment(.x))))
        ),
        penalty_grid = lasso_grid,
        model_fitting_summary = results_summary
      )
    } else NULL,
    
    # Complete model performance for advanced users
    all_model_performance = model_performance %>%
      dplyr::select(wflow_id, best_rmse, best_rsq, model_warnings) %>%
      dplyr::arrange(desc(best_rsq))
  )
  
  # Add contest class for method dispatch if needed
  class(contest_results) <- c("contest_results", "list")
  
  return(contest_results)
}

#' Validate contest results for prognostic_balance integration
#' @param contest_output Output from contest() function
validate_contest_output <- function(contest_output) {
  required_components <- c("best_spec", "best_engine", "best_model_name", 
                           "best_cv_rsq", "linear_cv_rsq")
  
  missing_components <- setdiff(required_components, names(contest_output))
  
  if (length(missing_components) > 0) {
    stop("Contest output missing required components: ", 
         paste(missing_components, collapse = ", "))
  }
  
  # Additional validation
  if (is.null(contest_output$best_spec)) {
    stop("Contest failed to identify a valid model specification")
  }
  
  if (contest_output$best_cv_rsq <= 0) {
    warning("Contest winner has very low predictive performance (R² ≤ 0)")
  }
  
  return(TRUE)
}
