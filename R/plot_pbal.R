#' Generates plot showing standardized covariate difference-in-means and prognosis from prognosis-weighted test results
#' @param pwtest_output Output from `pwtest()` function.
#' @param label_option character. Label covariates with most extreme values ("minmax" is default), all covariates ("all"), or none ("none").
#' @details The plot offers a visual diagnostic for covariate-by-covariate standardized difference in means and prognosis weights from standardized coefficients of prognostic regression of $Y^C(0)$ on set of covariates. The red dot indicates the value of the R-squared from the prognosis regression on the x-axis and the balance regression on the y-axis.
#' @import ggplot2
#' @export

plot_pbal <- function(pwtest_output, label_option = "minmax"){
  pdat <- as.data.frame(pwtest_output$cov_table)

    p <- ggplot(data = pdat) +
      geom_vline(xintercept = 0, linetype = "dashed", alpha = .6) +
      geom_hline(yintercept = 0, linetype = "dashed", alpha = .6) +
      geom_point(aes(x = prognosis, y = dim)) +
      geom_point(aes(x = pwtest_output$prog_Rsq, y = pwtest_output$bal_Rsq),
                 color = "red") +
      labs(y = "Covariate difference in means",
           x = "Covariate prognosis") +
      theme_classic()

    if(label_option != "none"){

    if(label_option == "all") sel_label <- 1:nrow(pdat)

    if(label_option == "minmax"){
      sel_label <- which(
        pdat$dim %in% c(min(pdat$dim, na.rm = TRUE), max(pdat$dim, na.rm = TRUE)) |
          pdat$prognosis %in% c(min(pdat$prognosis, na.rm = TRUE),
                                max(pdat$prognosis, na.rm = TRUE)))
      }

    p <- p + geom_text(data = pdat[sel_label,],
                       aes(x = prognosis, y = dim,
                           label = rownames(pdat), vjust = "inward",
                           hjust = "inward"), size = 3)
  }

  return(p)
}

