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

