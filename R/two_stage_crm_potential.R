#' @title C-optimal based two-stage Continual Reassessment Method (CRMC) Simulation for potential model
#' @description Performs one simulation run of a two-stage CRM design for estimating the Maximum Tolerated Dose (MTD) in phase I trials.
#'
#' @param p0 Target toxicity probability. Default is 0.4.
#' @param theta True value for the dose-toxicity curve. Default is 3.
#' @param theta_0 Nominal value of the dose-toxicity curve parameter (initial guess). Default is 2.7.
#' @param N Total number of patients (including both stages). Default is 24.
#' @param n_initial Number of patients per dose in stage 1. Default is 3.
#' @param q_0 Toxicity probability at first dose. Default is 0.02.
#' @param q_2 Fraction of patients in stage 1. Default is 0.5.
#' @param q_1 Target probability of observing at least one toxicity during stage 1. Default is 0.9.
#' @param fixed_optimal_dose Reference dose for optimal design estimation. Default is 0.2032.
#' @param show_plot Logical. If TRUE, plots the dose levels for patients in the trial. Default is FALSE.
#' @param seed Random seed for reproducibility. Default is 1234.
#'
#' @return A list with:
#' \describe{
#'   \item{n_toxicities}{Total number of toxicities observed.}
#'   \item{mtd_estimated}{Estimated Maximum Tolerated Dose (MTD). If no toxicity is observed in the first stage, the MTD is set to the largest dose level used.}
#'   \item{mle_theta}{Estimated value of the dose-toxicity parameter \eqn{\theta}. Set to \code{NA} if no toxicity is detected in stage 1.}
#'   \item{x}{Vector of dose levels administered.}
#'   \item{y}{Vector of toxicity outcomes (1 = toxic, 0 = non-toxic).}
#' }
#' @note
#' If no toxicity is observed during the first stage of the trial (i.e., \code{sum(y) == 0}), the simulation is terminated.
#' A warning is issued, and the MTD is conservatively estimated as the highest dose level reached.
#' The value of \code{mle_theta} is set to \code{NA} in this case.
#'
#' @examples
#' result <- two_stage_crm_potential(show_plot = TRUE)
#' print(result$mtd_estimated)
#'
#' @export
two_stage_crm_potential <- function(
  p0 = 0.4,
  theta = 3,
  theta_0 = 2.7,
  N = 24,
  n_initial = 3,
  q_0 = 0.02,
  q_2 = 0.5,
  q_1 = 0.9,
  fixed_optimal_dose = 0.2032,
  show_plot = FALSE,
  seed = 1234
) {
  # Step 1: Compute initial doses
  x0 <- q_0^(1 / theta_0)
  k_max <- floor(q_2 * N / n_initial)

  prob_h <- function(h) {
    prob <- 1
    for (k in 0:(k_max - 1)) {
      xk <- x0 + k * h
      prob <- prob * (1 - xk^theta_0)^n_initial
    }
    return(1 - prob - q_1)
  }

  h_opt <- uniroot(prob_h, lower = 0, upper = (1 - x0) / (k_max - 1))$root

  d_initial <- c()
  new_dose <- NULL
  for (k in 0:floor(N / n_initial)) {
    if ((x0 + k * h_opt) < 1) {
      new_dose <- x0 + k * h_opt
    }
    d_initial <- c(d_initial, new_dose)
  }
  # Step 2: Begin simulation for one run
  set.seed(seed)
  y <- NULL
  x <- NULL
  j <- 1

  while ((sum(y) < 1 | j < 4) &  (length(y)+n_initial)<(N+1)) {
    response_probabilities <- d_initial[j]^theta
    x <- c(x, rep(d_initial[j], n_initial))
    y <- c(y, rbinom(n_initial, size = 1, prob = response_probabilities))
    j <- j + 1
  }

  # After initial patient stage:
  if (sum(y) == 0) {
    warning("None toxicity detected in 1st stage. MTD estimated as the largest dose level")
    return(list(
      n_toxicities = 0,
      mtd_estimated = max(x),
      mle_theta = NA,
      x = x,
      y = y
    ))
  }

  result <- optim(par = theta_0, fn = log_likelihood, x = x, y = y,
                  method = "Brent", lower = 0.1, upper = theta_0+5)
  mle_theta <- result$par
  mtd_estimated <- p0^(1 / mle_theta)
  mtd_est_vector <- rep(mtd_estimated, length(x))

  for (i in 1:(N - (j - 1) * n_initial)) {
    opt_dose_estimated <- fixed_optimal_dose^(1 / mle_theta)
    response_probability <- opt_dose_estimated^theta

    x <- c(x, opt_dose_estimated)
    y <- c(y, rbinom(1, size = 1, prob = response_probability))

    result <- optim(par = mle_theta, fn = log_likelihood, x = x, y = y,
                    method = "Brent", lower = 0.1, upper = theta_0+5)
    mle_theta <- result$par
    mtd_estimated <- p0^(1 / mle_theta)
    mtd_est_vector <- c(mtd_est_vector, mtd_estimated)
  }

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

  return(list(n_toxicities = sum(y), mtd_estimated = mtd_estimated, mle_theta = mle_theta, x = x, y = y))
}

# Internal function (not exported)
log_likelihood <- function(theta, x, y) {
  log_likelihood_value <- sum(y * log(x^theta) + (1 - y) * log(1 - x^theta))
  return(-log_likelihood_value)
}
