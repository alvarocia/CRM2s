#' @title 3+3 Design Simulation for power model
#' @description Simulates one trial using a 3+3 dose-escalation method for estimating the Maximum Tolerated Dose (MTD).
#'
#' @param theta True value for the dose-toxicity curve. Default is 3.
#' @param theta_0 Nominal value for the dose-toxicity curve. Default is 2.7.
#' @param n_initial Number of patients per dose level. Default is 3.
#' @param p_tox_init Initial toxicity probability for computing starting dose. Default is 0.02.
#' @param delta_dosis Step size for dose escalation. Default is 0.05.
#' @param seed Random seed for reproducibility. Default is 1234.
#' @param show_plot Logical. If TRUE, plot dose level vs. patient index. Default is FALSE.
#'
#' @return A list with:
#' \describe{
#'   \item{n_toxicities}{Total number of toxicities observed.}
#'   \item{mtd_estimated}{Estimated MTD (last safe dose).}
#'   \item{n_patients}{Total number of patients enrolled.}
#'   \item{x}{Dose levels assigned.}
#'   \item{y}{Observed toxicity outcomes (1 = toxic, 0 = non-toxic).}
#' }
#'
#' @examples
#' res <- power_3_3()
#' print(res$mtd_estimated)
#'
#' @export
power_3_3 <- function(
  theta = 3,
  theta_0 = 2.7,
  n_initial = 3,
  p_tox_init = 0.02,
  delta_dosis = 0.05,
  seed = 1234,
  show_plot = FALSE
) {
  set.seed(seed)

  # Compute initial dose and sequence
  x0 <- p_tox_init^(1/theta_0)

  d_initial <- c()
  new_dose <- NULL
  for (k in 0:floor((1-x0)/delta_dosis + 5)) {
    if ((x0 + k * delta_dosis) < 1) {
      new_dose <- x0 + k * delta_dosis
    }
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

    prob_tox <- d_new^theta
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
         main = "",
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
