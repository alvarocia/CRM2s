#' @title C-optimal based two-stage Continual Reassessment Method (CRM2s) Simulation for logistic model with 2 parameters
#' @description Performs one simulation run of a two-stage CRM design for estimating the Maximum Tolerated Dose (MTD) in phase I trials.
#'
#' @param p0 Target toxicity probability. Default is 0.4.
#' @param theta True value for the vector of parameters for the dose-toxicity curve. Default is c(-3,2).
#' @param theta_0 Nominal value of the vector of parameters. Default is c(-3.1, 1.8).
#' @param N Total number of patients (including both stages). Default is 24.
#' @param n_initial Number of patients per dose in stage 1. Default is 3.
#' @param q_0 Toxicity probability at first dose. Default is 0.02.
#' @param q_2 Fraction of patients in stage 1. Default is 0.5.
#' @param q_1 Target probability of observing at least one toxicity during stage 1. Default is 0.9.
#' @param lim_sup_prob Maximum acceptable probability of toxicity. Default is 0.7.
#' @param show_plot Logical. If TRUE, plots the dose levels for patients in the trial. Default is FALSE.
#' @param seed Random seed for reproducibility. Default is 1234.
#'
#' @return A list with:
#' \describe{
#'   \item{n_toxicities}{Total number of toxicities observed.}
#'   \item{mtd_estimated}{Estimated Maximum Tolerated Dose (MTD). If no toxicity is observed in the first stage, the MTD is set to the largest dose level used.}
#'   \item{mle_theta}{Vector of estimated parameters \eqn{\theta = (\alpha, \beta)} for the logistic dose-toxicity model. Set to \code{NA} if no toxicity is detected in stage 1.}
#'   \item{x}{Vector of dose levels administered.}
#'   \item{y}{Vector of toxicity outcomes (1 = toxic, 0 = non-toxic).}
#' }
#' @note
#' If no toxicity is observed during the first stage of the trial (i.e., \code{sum(y) == 0}), the simulation is terminated.
#' A warning is issued, and the MTD is conservatively estimated as the highest dose level reached.
#' The value of \code{mle_theta} is set to \code{NA} in this case.
#'
#' @examples
#' result <- two_stage_crm_logistic(show_plot = TRUE)
#' print(result$mtd_estimated)
#'
#' @export
two_stage_crm_logistic <- function(
  p0 = 0.4,
  theta = c(-3,2),
  theta_0 = c(-3.1,1.8),
  N = 24,
  n_initial = 3,
  q_0 = 0.05,
  q_2 = 0.4,
  q_1 = 0.9,
  lim_sup_prob = 0.7,
  show_plot = FALSE,
  seed = 1234
) {
  if (theta[2] <= 0) {
    stop(sprintf("Invalid value: beta = %.2f must be > 0", theta[2]))
  }
  if (theta_0[2] <= 0) {
    stop(sprintf("Invalid value: beta_0 = %.2f must be > 0", theta_0[2]))
  }
  if (log(p0/(1-p0))<= theta[1]) {
    stop(sprintf("Invalid value: alpha = %.2f must be < log(p0/(1-p0))", theta[1]))
  }
  if (log(p0/(1-p0))<= theta_0[1]) {
    stop(sprintf("Invalid value: alpha_0 = %.2f must be < log(p0/(1-p0))", theta_0[1]))
  }

  # Step 1: Compute initial doses
  x0 <- (log(q_0/(1-q_0))-theta_0[1])/theta_0[2]
  x_max <- (log(lim_sup_prob/(1-lim_sup_prob))-theta_0[1])/theta_0[2]
  k_max <- floor(q_2 * N / n_initial)

  prob_h <- function(h) {
    prob <- 1
    for (k in 0:(k_max - 1)) {
      xk <- x0 + k * h
      prob <- prob * (1 - 1/(1+exp(-theta_0[1]-theta_0[2]*xk)))^n_initial
    }
    return(1 - prob - q_1)
  }

  h_opt <- uniroot(prob_h, lower = 0, upper = 20)$root

  d_initial <- c()
  new_dose <- NULL
  for (k in 0:floor(N / n_initial)) {
    if (x0 + k * h_opt < x_max){
      new_dose <- x0 + k * h_opt
    }
    else {
      new_dose <- x_max
    }
    d_initial <- c(d_initial, new_dose)
  }
  # Step 2: Begin simulation for one run
  set.seed(seed)
  y <- NULL
  x <- NULL
  j <- 1

  while ((sum(y) < 1 | j < 4) &  (length(y)+n_initial)<(N+1)) {
    response_probabilities <- 1/(1+exp(-theta[1]-theta[2]*d_initial[j]))
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

  # model <- glm(y ~ x, family = binomial(link = "logit"))
  #
  # mle_theta <- as.numeric(coef(model))

  neg_loglik <- function(par, xx, yy) {
    b0 <- par[1]
    b1 <- par[2]
    eta <- b0 + b1 * xx
    p <- 1 / (1 + exp(-eta))
    # Bound p away from 0 and 1 to avoid log(0)
    eps <- 1e-8
    p <- pmax(pmin(p, 1 - eps), eps)
    -sum(yy * log(p) + (1 - yy) * log(1 - p))
  }

  # Fit with constraints:
  result <- optim(
    par = theta_0,                     # initial values
    fn = neg_loglik,
    xx = x,
    yy = y,
    method = "L-BFGS-B",
    lower = c(-Inf, 1e-08),
    upper = c(log(p0/(1-p0))-1e-08, Inf)
  )

  # Extract parameter estimates
  mle_theta <- result$par
  mtd_estimated <- (log(p0/(1-p0))-mle_theta[1])/mle_theta[2]
  mtd_est_vector <- rep(mtd_estimated, length(x))

  if ((x0<=mtd_estimated)& (mtd_estimated<=x_max)){
    opt_dose_estimated <- mtd_estimated
  }
  else if (mtd_estimated>x_max){
    opt_dose_estimated <- x_max
  }
  else {
    opt_dose_estimated <- x0
  }


  for (i in 1:(N - (j - 1) * n_initial)) {

    response_probability <- 1/(1+exp(-theta[1]-theta[2]*opt_dose_estimated))

    x <- c(x, opt_dose_estimated)
    y <- c(y, rbinom(1, size = 1, prob = response_probability))

    model <- suppressWarnings(glm(y ~ x, family = binomial(link = "logit")))

    mle_theta <- as.numeric(coef(model))
    mtd_estimated <- (log(p0/(1-p0))-mle_theta[1])/mle_theta[2]


    if ((x0<=mtd_estimated)& (mtd_estimated<=x_max)){
      opt_dose_estimated <- mtd_estimated
    }
    else if (mtd_estimated>x_max){
      opt_dose_estimated <- x_max
    }
    else {
      opt_dose_estimated <- x0
    }


    mtd_est_vector <- c(mtd_est_vector, opt_dose_estimated)
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

  return(list(n_toxicities = sum(y), mtd_estimated = opt_dose_estimated, mle_theta = mle_theta, x = x, y = y))
}
