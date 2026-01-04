# data <- simulate_experiment_complex(N = 500, target_r2 = .40,
#                                     imbalance = rep(.1, 5))
#
# data$W <- rnorm(500)
#
# test <- pwtest(data = data, covariates = c("X1", "X2", "X1_sq",
#                                            "X2_sq", "X1_X2"),
#                outcome = "Y", treatment = "Z", cv_folds = 3,
#                model_spec = NULL, cv_auto = TRUE,
#                subset_workflow = "random_forest",
#                verbose = FALSE)
#
# test <- pwtest_rdd(data = data, covariates = c("X1", "X2", "X1_sq",
#                                                "X2_sq", "X1_X2"),
#                    outcome = "Y", treatment = "Z", running_var = "W")
#
# output <- contest(data = data,
#                   covariates = c("X1", "X2", "X1_sq",
#                                  "X2_sq", "X1_X2"),
#                   outcome = "Y0", cv_folds = 3,
#                   subset_workflow = "random_forest")
