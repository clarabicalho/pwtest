#' Produces unweighted and prognosis-weighted test statistics with standard errors and p-values
#' @param data data.frame containing covariates, treatment assignment, and outcome variable
#' @param covariates character vector of names of placebo variables. This is not to be confused with `covs` argument passed onto `rdrobust` in the case of RD designs, which can be specified separately (`...`).
#' @param treatment name of variable indicating (binary) treatment assigned
#' @param outcome name of outcome variable
#' @param n_bootstraps numeric scalar indicating number of bootstraps to use for the resampling-based p-values
#' @param equivalence logical value. If TRUE, carry out equivalence test as described in Hartman and Hidalgo (2018). For a full description of the bootstrapping procedure, see paper.
#' @param equiv_bounds numerical vector. If `equivalence = TRUE`, lower and upper bounds respectively of the equivalence range as a multiplier for the standard deviation of potential outcomes under control. Both values are set to 0.36 if set to `NULL` (see Hartman and Hidalgo (2018)).
#' @param model_spec an object of class `model_spec` defining the functional form of the prediction model to fit $Y^C(0)$ on covariate matrix. See https://www.tidymodels.org/find/parsnip/#models for all available models.
#' @param engine character. Engine type used for the model specified in `model_spec`. See https://www.tidymodels.org/find/parsnip/#models for all available engines for each model. If not defined, uses default engine defined by `model_spec` type in package `parsnip`, if available.
#' @param formula object of class formula or character describing the model to fit `model_spec` on control sample. Defaults to regressing outcome on full set of covariates defined by `covariates`.
#' @import tidymodels
#' @export

pwtest <- function(data,
                   covariates = c("X1", "X2", "X3"),
                   treatment = "Z",
                   outcome = "Y",
                   n_bootstraps = 500, # in models other than lm, the bootstrap sample also defines the training sample
                   spec_mod = linear_reg(mode = "regression", engine = "lm"),
                   engine = NULL,
                   formula = NULL) {

  if(is.null(formula)) formula <- as.formula(paste0(outcome, "~."))

  # restrict data to variables of interest
  data <- data %>% dplyr::select(tidyselect::all_of(c(treatment, covariates, outcome)))

  # observed statistics-------------------------------------

  # standardize data relative to entire study group (the finite population)
  data <- data %>%
    dplyr::mutate_at(.vars = c(outcome, covariates),
                     .funs = scale) %>% as.data.frame()

  # define model engine
  if(!is.null(engine)) spec_mod <- spec_mod %>% set_engine(engine)

  # train data, Yc(0)
  datc <- data %>% filter(!!sym(treatment) == 0) %>% dplyr::select(all_of(c(outcome, covariates)))

  # test data, Yt(0)
  datt <- data %>% filter(!!sym(treatment) == 1) %>% dplyr::select(all_of(c(outcome, covariates)))

  fit_Yc <- spec_mod %>% fit(formula = formula, data = datc)

  # extract standard metrics
  fit_metrics <- predict(fit_Yc, datc) %>%
    bind_cols(Y = as.vector(datc$Y)) %>% metrics(Y, .pred)

  predict_Yt <-  predict(fit_Yc, datt)

  # observed \overline{Y^T(0)} - \overline{Y^C(0)}
  pwdelta_obs <-  mean(predict_Yt$.pred) - mean(datc$Y)

  # bootstrapping -------------------------------------

  # resample from control group with replacement
  treat_i <- which(data[[treatment]] == 1)
  control_i <- which(data[[treatment]] == 0)
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

  # SE and p-value from bootstrap distribution
  pwdelta_se <- sd(unlist(sapply(bstats, function(x) x$pwdelta)), na.rm = TRUE) * (n_bootstraps - 1) / n_bootstraps
  pw_delta_p <- sum(abs(unlist(sapply(bstats, function(x) x$pwdelta))) >= abs(pwdelta_obs), na.rm = TRUE) / n_bootstraps

  # output ------------------------------------------------------------------
  estimates <- list(
    n_bootstraps = n_bootstraps,
    pwdelta = pwdelta_obs,
    pwdelta_se = unname(pwdelta_se),
    pwdelta_p = get_pvalue_pw(bstats, pwdelta_obs, n_bootstraps),
    fit_obj = fit_Yc,
    fit_metrics = fit_metrics
  )

  return(estimates)
}


# Example:
library(tidymodels)
source("~/Library/CloudStorage/Dropbox/Github-Cloud/prognostic_balance/code/simulations/R/simulate_experiment_poly_4.R")

data <- simulate_experiment(N = 1000, imbalance = rep(.2, 3), prognosis = c(.3, .2, .1), poly = c(2, 1))

pwstat <- pwtest(data, covariates = c("X1", "X2", "X1_2"), spec_mod = rand_forest(mode = "regression"), n_bootstraps = 200)
pwstat2 <- pwtest(data, covariates = c("X1", "X2", "X1_2"), spec_mod = linear_reg(mode = "regression"), n_bootstraps = 200)

pwstat2 <- pwtest(data, covariates = c("X1", "X2", "X1_2"),
                  spec_mod = linear_reg(mode = "regression", penalty = .5),
                  engine = "glmnet", n_bootstraps = 200)
