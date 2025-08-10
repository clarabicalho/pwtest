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
#' @param return_uwtest logical. Whether to return results from unweighted test as well.
#' @param ... additional cross-validation arguments when `cv_auto = TRUE`. See `?pick_winner`.
#' @importFrom parsnip linear_reg fit set_engine
#' @importFrom yardstick metrics
#' @importFrom stats predict sd setNames as.formula t.test
#' @importFrom dplyr bind_cols select mutate_at
#' @importFrom rlang sym :=
#' @export
pwtest <- function(data,
                   covariates = c("X1", "X2", "X3"),
                   treatment = "Z",
                   outcome = "Y",
                   n_bootstraps = 500, # in models other than lm, the bootstrap sample also defines the training sample
                   model_spec = NULL,
                   engine = NULL,
                   formula = NULL,
                   cv_auto = FALSE,
                   return_uwtest = FALSE,
                   ...) {

  if(cv_auto & !is.null(model_spec)) {
    warning("Argument `model_spec` will be overwritten by contest winner since `cv_auto = TRUE`. To use a specific model, set `cv_auto = FALSE`.")
  }
  if(!cv_auto & (!is.null(model_spec) & all(class(model_spec) != "model_spec"))) {
    stop("Argument `model_spec` must be an object of class `model_spec` from package `parsnip`. See https://www.tidymodels.org/find/parsnip/#models for more details.")
  }
  if(!is.null(formula)){
    if(is.character(formula)) formula <- as.formula(formula)
    else{
      stop("Argument `formula` must be either of class character or formula.")
    }
  }
  if(any(!is.character(outcome), !is.character(treatment), !is.character(covariates))){
    stop("Missing outcome, treatment, or covariates. All arguments must be character vectors of variable names.")
  }

  # Validate required variables exist -----------------------
  required_vars <- c(treatment, covariates, outcome)
  missing_vars <- setdiff(required_vars, names(data))
  if (length(missing_vars) > 0) {
    stop("Missing variables in data: ", paste(missing_vars, collapse = ", "))
  }

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
  # index control and treatment rows
  treat_i <- which(data[[treatment]] == 1)
  control_i <- which(data[[treatment]] == 0)

  control_data <- data[control_i, c(outcome, covariates)]
  if (nrow(control_data) == 0) {
    stop("No control group observations found. Check treatment variable coding.")
  }

  # run cross validation (cv_auto = TRUE) -------------------
  if(cv_auto){
    # pass on arguments from parent function
    argg_contest <- c(argg[names(argg) %in% names(formals(contest))])
    argg_contest$data <- control_data

    # run contest to get best model based on full control data performance
    winner <- do.call(contest, argg_contest)
    model_spec <- winner$best_spec
    engine <- winner$best_engine
    argg$recipe <- winner$best_recipe  # May need this later for recipe-based models
  }


  # observed statistics -------------------------------------


  # standardize data relative to control group SD
  data <- std_data(data, c(outcome, covariates), treatment)

  # difference in means -----------------------------------------------------
  cov_ttest <- diff_in_means(data, covariates, control_i, treat_i)
  DIM <- cov_ttest$dim
  # if any missing/errors
  if(any(is.na(DIM))) stop(paste0("Difference in means could not be calculated for the following covariate(s): ",
                                  paste0(names(DIM)[is.na(DIM)], collapse = ", ")))

  # prognosis ---------------------------------------------------------------
  pweights <- pw(data, covariates, treatment, outcome)

  if(any(is.na(pweights))) warning(paste0("Covariate(s) ",
                                           paste(names(pweights)[is.na(pweights)], collapse = ", "),
                                           " have missing prognosis coefficients and receive weight of 0 in the prognosis-weighted test statistic."))
  cov_table <- cbind(prognosis = pweights[covariates], cov_ttest)
  pwdelta_obs <- sum(cov_table$dim[!is.na(cov_table$prognosis)] * cov_table$prognosis[!is.na(cov_table$prognosis)])

  if(!is.null(model_spec)){
    # data to fit Yc(0)
    datc <- data[control_i, c(outcome, covariates)]
    # data to predict Yt(0)
    datt <- data[treat_i, c(outcome, covariates)]
    fit_Yc <- model_spec %>% fit(formula = formula, data = datc)
    # extract standard metrics
    fit_metrics <- predict(fit_Yc, datc) %>%
      bind_cols(!!outcome := as.vector(datc[[outcome]])) %>%
      yardstick::metrics(outcome, .pred)
    predict_Yt <-  predict(fit_Yc, datt)
    predict_Yc <-  predict(fit_Yc, datc)
    # observed \overline{\widehat{Y^T(0)}} - \overline{Y^C(0)}
    pwdelta_obs <-  mean(predict_Yt$.pred, na.rm = TRUE) -
      mean(predict_Yc$.pred, na.rm = TRUE)
  }

  # bootstrapping -------------------------------------

  # resample from control group with replacement
  samples <- get_bootstrap_samples(control_i, length(treat_i), n_bootstraps)

  # obtain delta distribution from bootstrap samples
  pwdelta_from_draw <- function(z) {

    bsample <- data[z, ]
    bsample[[treatment]][1:length(treat_i)] <- 1

    if(is.null(model_spec)){
      # difference in means
      DIM <- diff_in_means(bsample, covariates, control_i, treat_i)$dim
      # if any missing/errors
      if(any(is.na(DIM))) stop(paste0("Difference in means could not be calculated for the following covariate(s): ",
                                      paste0(names(DIM)[is.na(DIM)], collapse = ", ")))

      # prognosis coefficients on standardized data
      bpweights <- pw(bsample, covariates, treatment, outcome)

      # replicate the delta statistic using the same covariate set
      # as the observed data set
      dY <- sum(DIM[!is.na(pweights)]*bpweights[!is.na(pweights)])

    } else {
      fit_Yc <- model_spec %>%
        fit(formula = formula,
            data = bsample[bsample[[treatment]]==0,])

      predict_Yc <- predict(fit_Yc, bsample[bsample[[treatment]] == 0, ])
      predict_Yt <-  predict(fit_Yc, bsample[bsample[[treatment]] == 1, ])
      dY <- mean(predict_Yt$.pred, na.rm = TRUE) -
        mean(predict_Yc$.pred, na.rm = TRUE)
    }

    return(list(pwdelta=dY))
  }

  bstats <- capture_warnings_apply(samples, 2, pwdelta_from_draw)

  # Identify which elements have NA values
  has_na <- lapply(bstats, function(x) any(is.na(x$pwdelta)))

  # Filter the bootstraps to remove when returns NA statistics
  bstats <- bstats[!unlist(has_na)]

  # check for completeness and draw additional bootstrap samples if necessary
  nc_bootstraps <- length(bstats)

  while (nc_bootstraps < as.integer(n_bootstraps)) {
    message("NAs in bootstrap sample statistics. Resampling...")
    add_samples <- get_bootstrap_samples(control_i, length(treat_i), n_bootstraps-nc_bootstraps)


    # obtain delta distribution from bootstrap samples
    add_bstats <- capture_warnings_apply(add_samples, 2, pwdelta_from_draw)

    bstats <- c(bstats, add_bstats)

    # Identify which elements have NA values
    has_na <- lapply(bstats, function(x) is.na(x$pwdelta))

    # Filter the bootstraps to remove when returns NA statistics
    bstats <- bstats[!unlist(has_na)]
    if(length(bstats) > n_bootstraps) bstats <- bstats[1:n_bootstraps]
    # check for completeness and draw additional bootstrap samples if necessary
    nc_bootstraps <- length(bstats)

  }

  # SE bootstrap distribution
  pwdelta_se <- sd(unlist(sapply(bstats, function(x) x$pwdelta)), na.rm = TRUE) * (n_bootstraps - 1) / n_bootstraps %>% unname

  # always estimate pwdelta from default sum of weighted difference in means ----

  # modify relevant arguments to fit linear regression
  if(!is.null(model_spec)){
    argg$model_spec <- NULL
    argg$engine <- NULL; argg$formula <- NULL
    argg$recipe <- NULL; argg$cv_auto <- FALSE
    argg$verbose <- FALSE
    argg$subset_workflow <- NULL
    pwtest_lr <- do.call("pwtest", args = argg)
  } else {
    pwtest_lr <- NULL
  }

  # prognosis R-squared -------------------------------------------------------

  prog_formula <- as.formula(paste0(outcome, "~", paste(covariates, collapse = "+")))
  prog_lr <- lm(prog_formula, data = data[control_i,])
  prog_rsqr <- summary(prog_lr)$r.squared

  # balance R-squared -------------------------------------------------------

  bal_formula <- as.formula(paste0(treatment, "~", paste(covariates, collapse = "+")))
  bal_rsqr <- summary(lm(bal_formula, data = data))$r.squared

  # output ------------------------------------------------------------------

  #pwdelta estimates
  est_data <- data.frame(
    method = ifelse(is.null(model_spec), "default", class(model_spec)[1]),
    engine = ifelse(is.null(engine), "default", engine),
    formula = deparse(formula),
    n_bootstraps = n_bootstraps,
    pwdelta = pwdelta_obs,
    pwdelta_se = pwdelta_se,
    pwdelta_p = get_pvalue_pw(bstats, pwdelta_obs, n_bootstraps)
  )

  # unweighted test
  if(return_uwtest){
    uwdelta <- uwtest(data = data, covariates = covariates, treatment = treatment,
                      outcome = outcome, n_bootstraps = n_bootstraps)
    est_data <- cbind(est_data, uwdelta[,-1])
  }

  # append estimates from linear regression if missing
  if(!is.null(pwtest_lr)) {
    est_data <- rbind(est_data, pwtest_lr$estimates)

    estimates <- list(
      estimates = est_data,
      cov_table = cov_table,
      fit_obj = list(prog_lr, fit_Yc),
      fit_metrics = rbind(pwtest_lr$fit_metrics,
        cbind(model = class(model_spec)[1], fit_metrics)
      ),
      prog_rsqr = prog_rsqr,
      bal_rsqr = bal_rsqr
    )

  } else {
    estimates <- list(
      estimates = est_data,
      cov_table = cov_table,
      fit_obj = list(prog_lr),
      fit_metrics = data.frame(model = "default",
                               .metric = "rsq",
                               .estimator = "standard",
                               .estimate = prog_rsqr),
      prog_Rsq = prog_rsqr,
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
#' @param min_penalty_exp minimum penalty exponent (default -6 for more conservative start)
#' @param max_penalty_exp maximum penalty exponent (default -1 for small data sets)
#' @param subset_workflow character vector containing the labels of workflows to subset the contest on. See Details for more information.
#' @param verbose logical. Whether to print detailed information during contests (see `?contest()`)
#' @export
pick_winner <- function(data,
                        outcome,
                        covariates,
                        cv_folds = 5,
                        min_penalty_exp = -1,
                        max_penalty_exp = 6,
                        subset_workflow = NULL,
                        verbose = FALSE
) {

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
