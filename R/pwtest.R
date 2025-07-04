#' Produces unweighted and prognosis-weighted test statistics with standard errors and p-values
#' @param data data.frame containing covariates, treatment assignment, and outcome variable
#' @param covariates character vector of names of placebo variables. This is not to be confused with `covs` argument passed onto `rdrobust` in the case of RD designs, which can be specified separately (`...`).
#' @param treatment name of variable indicating (binary) treatment assigned
#' @param outcome name of outcome variable
#' @param n_bootstraps numeric scalar indicating number of bootstraps to use for the resampling-based p-values
#' @param model_spec an object of class `model_spec` defining the functional form of the prediction model to fit $Y^C(0)$ on covariate matrix. See https://www.tidymodels.org/find/parsnip/#models for all available models.
#' @param engine character. Engine type used for the model specified in `model_spec`. See https://www.tidymodels.org/find/parsnip/#models for all available engines for each model. If not defined, uses default engine defined by `model_spec` type in package `parsnip`, if available.
#' @param formula object of class formula or character describing the model to fit `model_spec` on control sample. Defaults to regressing outcome on full set of covariates defined by `covariates`.
#' @param cv_auto logical Whether to perform automatic cross-validation and run test on contest winner. If `TRUE`, `model_spec` and `engine` will be overwritten by contest winner.
#' @param ... additional cross-validation arguments when `cv_auto = TRUE`. See `?pick_winner`.
#' @importFrom parsnip linear_reg fit set_engine
#' @importFrom yardstick metrics
#' @importFrom stats predict sd setNames as.formula t.test
#' @importFrom dplyr bind_cols select mutate_at
#' @importFrom rlang sym
#' @export
pwtest <- function(data,
                   covariates = c("X1", "X2", "X3"),
                   treatment = "Z",
                   outcome = "Y",
                   n_bootstraps = 500, # in models other than lm, the bootstrap sample also defines the training sample
                   model_spec = linear_reg(mode = "regression", engine = "lm"),
                   engine = NULL,
                   formula = NULL,
                   cv_auto = FALSE,
                   ...) {

  # extract arguments (default and user-supplied) from parent function
  argg <- as.list(match.call())[-1]
  arg_formals <- formals(pwtest) ## formals with default arguments
  for (v in names(arg_formals)){
    if (!(v %in% names(argg)))
      ## if arg value is missing, add its default
      argg <- append(argg, arg_formals[v])
  }

  # default regression model
  if(is.null(formula)) formula <- as.formula(paste0(outcome, "~."))

  # define model engine
  if(!is.null(engine)) model_spec <- model_spec %>% set_engine(engine)

  # restrict data to variables of interest
  data <- data[,c(treatment, covariates, outcome)]

  # standardize data relative to entire study group (finite population)
  data <- data %>%
    mutate_at(.vars = c(outcome, covariates),
              .funs = scale) %>% as.data.frame()

  treat_i <- which(data[[treatment]] == 1)
  control_i <- which(data[[treatment]] == 0)

  # Validate required variables exist -----------------------
  required_vars <- c(treatment, covariates, outcome)
  missing_vars <- setdiff(required_vars, names(data))
  if (length(missing_vars) > 0) {
    stop("Missing variables in data: ", paste(missing_vars, collapse = ", "))
  }

  control_data <- data[control_i, c(outcome, covariates)]
  if (nrow(control_data) == 0) {
    stop("No control group observations found. Check treatment variable coding.")
  }

  # run cross validation (cv_auto = TRUE) -------------------
  if(cv_auto){
    # pass on arguments from parent function
    argg_contest <- c(argg[names(argg) %in% names(formals(pick_winner))])
    argg_contest$data <- control_data

    # pick winner and output model specs
    winner <- do.call(pick_winner, argg_contest)
    argg$model_spec <- winner$model_spec
    argg$engine <- winner$engine
    # argg$recipe <- winner$recipe
  }

  # observed statistics -------------------------------------

  # data to fit Yc(0)
  datc <- data[control_i, c(outcome, covariates)]
  # data to predict Yt(0)
  datt <- data[treat_i, c(outcome, covariates)]
  fit_Yc <- model_spec %>% fit(formula = formula, data = datc)
  # extract standard metrics
  fit_metrics <- predict(fit_Yc, datc) %>%
    bind_cols(Y = as.vector(datc$Y)) %>% yardstick::metrics(Y, .pred)
  predict_Yt <-  predict(fit_Yc, datt)
  # observed \overline{\widehat{Y^T(0)}} - \overline{Y^C(0)}
  pwdelta_obs <-  mean(predict_Yt$.pred) - mean(datc$Y)

  # bootstrapping -------------------------------------

  # resample from control group with replacement
  samples <- get_bootstrap_samples(control_i, length(treat_i), n_bootstraps)

  # standardize bootstrap population relative to control group SD
  data_stdc <- std_data(data, c(outcome, covariates), treatment)

  # obtain delta distribution from bootstrap samples
  bstats <- apply(samples, 2, function(z) {
    bsample <- data_stdc[z, ]
    bsample[[treatment]][1:length(treat_i)] <- 1

    predict_Yc <- predict(fit_Yc, bsample[bsample[[treatment]] == 0, ])
    predict_Yt <-  predict(fit_Yc, bsample[bsample[[treatment]] == 1, ])
    dY <- mean(predict_Yt$.pred) - mean(predict_Yc$.pred)

    return(list(pwdelta=dY))
  })

  # check for completeness and draw additional bootstrap samples if necessary
  nc_bootstraps <- sum(!is.na(sapply(bstats, function(x) x$pwdelta)))

  while (nc_bootstraps < as.integer(n_bootstraps)) {
    add_samples <- get_bootstrap_samples(control_i, length(treat_i), n_bootstraps-nc_bootstraps)

    # obtain delta distribution from bootstrap samples
    add_bstats <- apply(add_samples, 2, function(z) {
      bsample <- data_stdc[z, ]
      bsample[[treatment]][1:length(treat_i)] <- 1

      predict_Yc <- predict(fit_Yc, bsample[bsample[[treatment]] == 0, ])
      predict_Yt <-  predict(fit_Yc, bsample[bsample[[treatment]] == 1, ])
      dY <- mean(predict_Yt$.pred) - mean(predict_Yc$.pred)

      return(list(pwdelta=dY))
    })

    bstats <- c(bstats, add_bstats)
    nc_bootstraps <- sum(!is.na(sapply(bstats, function(x) x$pwdelta)))
  }

  # SE bootstrap distribution
  pwdelta_se <- sd(unlist(sapply(bstats, function(x) x$pwdelta)), na.rm = TRUE) * (n_bootstraps - 1) / n_bootstraps %>% unname

  # always estimate pwdelta from linear model unless model_spec is lm -------

  # modify relevant arguments to fit linear regression
  if(!"linear_reg" %in% class(model_spec)){
    argg$model_spec <- linear_reg(mode = "regression", engine = "lm")
    argg$engine <- NULL; argg$formula <- NULL
    pwtest_lm <- do.call("pwtest", args = argg)
  } else {
    pwtest_lm <- NULL
  }


  # difference in means -----------------------------------------------------
  dim_ests <- sapply(covariates, function(x){
    tres <- t.test(data[treat_i, x], data[control_i, x], na.rm = TRUE)
    dim <- tres$estimate[1] - tres$estimate[2]
    return(c(dim = unname(dim), ttest_p = tres$p.value))
  })

  # prognosis ---------------------------------------------------------------
  if(is.null(pwtest_lm)) prognosis <- coef(fit_Yc$fit)
  else prognosis <- coef(pwtest_lm$fit_obj[[1]]$fit)

  cov_table <- t(rbind(prognosis = prognosis[covariates], dim_ests))

  # balance R-squared -------------------------------------------------------

  bal_formula <- as.formula(paste0(treatment, "~", paste(covariates, collapse = "+")))
  bal_rsqr <- summary(lm(bal_formula, data = data))$r.squared

  # output ------------------------------------------------------------------

  #pwdelta estimates
  est_data <- data.frame(
    method = class(model_spec)[1],
    engine = ifelse(is.null(engine), "default", engine),
    formula = deparse(formula),
    n_bootstraps = n_bootstraps,
    pwdelta = pwdelta_obs,
    pwdelta_se = pwdelta_se,
    pwdelta_p = get_pvalue_pw(bstats, pwdelta_obs, n_bootstraps)
  )

  # append estimates from linear regression if missing
  if(!is.null(pwtest_lm)) {
    est_data <- rbind(est_data, pwtest_lm$estimates)

    estimates <- list(
      estimates = est_data,
      cov_table = cov_table,
      fit_obj = list(pwtest_lm$fit_obj[[1]], fit_Yc),
      fit_metrics = rbind(
        pwtest_lm$fit_metrics,
        cbind(model = class(model_spec)[1], fit_metrics)
      ),
      prog_rsqr = pwtest_lm$fit_metrics %>% filter(.metric == "rsq") %>% pull(.estimate),
      bal_rsqr = bal_rsqr
    )

  } else {
    estimates <- list(
      estimates = est_data,
      cov_table = cov_table,
      fit_obj = list(fit_Yc),
      fit_metrics = cbind(model = class(model_spec)[1], fit_metrics),
      prog_Rsq = fit_metrics %>% filter(.metric == "rsq") %>% pull(.estimate),
      bal_Rsq = bal_rsqr
    )
  }

  return(estimates)
}

#' High-level wrapper for prognostic balance testing with automatic model selection
#' @param data data.frame containing covariates and outcome variable and subset to control group units
#' @param covariates character vector of names of placebo variables
#' @param outcome name of outcome variable
#' @param cv_folds integer. Number of cross-validation folds (see `?contest()`)
#' @param verbose logical. Whether to print detailed information during contests (see `?contest()`)
#' @export
pick_winner <- function(data,
                        outcome,
                        covariates,
                        cv_folds = 3,
                        verbose = FALSE) {

  ###########################################################################
  # AUTO MODE: Contest selection
  ###########################################################################

  argg <- as.list(match.call())[-1]

  # Run contest
  if (verbose) cat("=== Model selection diagnostics ===\n")
  contest_results <- do.call(contest, argg)

  if (!verbose) {
    cat("Contest selected:", contest_results$best_model_name,
        "(R-squared = ", round(contest_results$best_cv_rsq, 4), ")\n")
  }

  # Prepare pwtest with contest winner
  model_spec <- contest_results$best_spec
  engine <- contest_results$best_engine
  if (!is.null(contest_results$best_recipe)) {
    recipe <- contest_results$best_recipe
  }

  if (verbose) cat("\n--- Returning contest winner ---\n")

  # Return pwtest result directly
  return(list(model_spec = model_spec, engine = engine, recipe = recipe))

}
