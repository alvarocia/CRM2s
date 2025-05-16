#' @title 3+3 Design Simulation for Logistic model with 2 parameters
#' @description Simulates one trial using a 3+3 dose-escalation method for estimating the Maximum Tolerated Dose (MTD).
#'
#' @param theta True value of the vector of parameters for the dose-toxicity curve. Default is c(-3,2).
#' @param theta_0 Nominal value for the vector of parameters. Default is c(-3.1,1.8).
#' @param n_initial Number of patients per dose level. Default is 3.
#' @param p_tox_init Initial toxicity probability for computing starting dose. Default is 0.02.
#' @param delta_dosis Step size for dose escalation. Default is 0.05.
#' @param seed Random seed for reproducibility. Default is 1234.
#' @param show_plot Logical. If TRUE, plot dose level vs. patient index. Default is FALSE.
#'
#' @return A list with:
#'   \item{n_toxicities}{Total number of toxicities observed}
#'   \item{mtd_estimated}{Estimated MTD (last safe dose)}
#'   \item{n_patients}{Total number of patients enrolled}
#'   \item{x}{Dose levels assigned}
#'   \item{y}{Observed toxicity outcomes (1 = toxic, 0 = non-toxic)}
#'
#' @examples
#' res <- logistic_3_3()
#' print(res$mtd_estimated)
#'
#' @export
logistic_3_3 <- function(
  theta_0 = c(-3,2),
  theta = c(-3.1,1.8),
  n_initial = 3,
  p_tox_init = 0.05,
  delta_dosis = 0.1,
  seed = 1234,
  show_plot = FALSE
) {
  set.seed(seed)

  # Compute initial dose and sequence
  x0 <- (log(p_tox_init/(1-p_tox_init))-theta_0[1])/theta_0[2]

  d_initial <- c()
  for (k in 0:floor((1-x0)/delta_dosis + 5)) {
    new_dose <- x0 + k * delta_dosis
    d_initial <- c(d_initial, new_dose)
  }

  y <- NULL
  y_old <- NULL
  y_new <- NULL
  d_new <- d_initial[1]
  d_old <- NULL
  x <- NULL
  j <- 1

  # Main 3+3 escalation loop
  while (((sum(y_old) + sum(y_new)) < 2) & j < 100) {
    y_old <- y_new

    if (sum(y_new) == 1) {
      d_old <- d_new
      d_new <- d_new
    }
    if (sum(y_new) == 0) {
      d_old <- d_new
      d_new <- d_new + delta_dosis
    }

    prob_tox <- 1/(1+exp(-theta[1]-theta[2]*d_new))
    x <- c(x, rep(d_new, n_initial))
    y_new <- rbinom(n_initial, size = 1, prob = prob_tox)
    y <- c(y, y_new)
    j <- j + 1
  }

  mtd_est <- d_new - delta_dosis
  n_tox <- sum(y)
  n_patients <- (j - 1) * n_initial

  if (show_plot) {
    plot(1:length(x), x,
         pch = ifelse(y == 1, 4, 1),
         xlab = "Patient number",
         ylab = "Dose level",
         main = "3+3 Design",
         ylim = c(min(x) * 0.9, max(x) * 1.1))

    legend("bottomright",
           legend = c("Toxicity", "No toxicity"),
           pch = c(4, 1),
           col = c("black", "black"),
           pt.cex = 1.2,
           bty = "n")
  }

  return(list(
    n_toxicities = n_tox,
    mtd_estimated = mtd_est,
    n_patients = n_patients,
    x = x,
    y = y
  ))
}
