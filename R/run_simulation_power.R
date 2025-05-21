#' @title Run Simulation Comparison Between 3+3 and CRM2s for the power model
#' @description Compares the MTD estimation and toxicity count between the power 3+3 method and the two-stage CRM using multiple replications.
#'
#' @param num_rep Number of replications to run. Default is 1000.
#' @param seed Base random seed for reproducibility. Default is 1234.
#' @param save_plot Logical. If TRUE, saves comparison plots as PDF. Default is FALSE.
#' @param p0 Target toxicity probability. Default is 0.4.
#' @param theta_0 Nominal value of theta used in CRM and 3+3 escalation models. Default is 2.7.
#' @param theta True theta used for generating the MTD reference. Default is 3.
#' @param N Total number of patients in the CRM design. Default is 24.
#' @param n_initial Number of patients per dose level in the CRM design (not used in 3+3). Default is 1.
#' @param q_0 Initial toxicity probability for CRM design. Default is 0.05.
#' @param q_2 Fraction of patients in CRM stage 1. Default is 0.5.
#' @param q_1 Target probability of observing at least one toxicity in CRM stage 1. Default is 0.9.
#' @param p_tox_init_3_3 Initial toxicity probability for the 3+3 model. Default is 0.02.
#' @param delta_dosis_3_3 Step size for dose escalation in the 3+3 model. Default is 0.055.
#' @param fixed_optimal_dose Reference dose for optimal CRM design. Default is 0.2032.
#'
#' @return A \code{data.frame} with one row per method ("3+3" and "CRM2s") and the following columns:
#' \describe{
#'   \item{method}{Design used ("3+3" or "CRM2s")}
#'   \item{mean_pat}{Mean of patients}
#'   \item{median_pat}{Median of patients}
#'   \item{mean_mtd}{Mean of the estimated MTDs}
#'   \item{sd_mtd}{Standard deviation of the estimated MTDs}
#'   \item{median_mtd}{Median of the estimated MTDs}
#'   \item{min_mtd}{Minimum of the estimated MTDs}
#'   \item{q1_mtd}{First quartile (Q1) of the estimated MTDs}
#'   \item{q3_mtd}{Third quartile (Q3) of the estimated MTDs}
#'   \item{max_mtd}{Maximum of the estimated MTDs}
#'   \item{siqr_mtd}{Semi-interquartile range of the estimated MTDs}
#'   \item{mean_tox}{Mean number of toxicities}
#'   \item{sd_tox}{Standard deviation of the number of toxicities}
#'   \item{median_tox}{Median number of toxicities}
#'   \item{min_tox}{Minimum number of toxicities}
#'   \item{q1_tox}{First quartile (Q1) of the number of toxicities}
#'   \item{q3_tox}{Third quartile (Q3) of the number of toxicities}
#'   \item{max_tox}{Maximum number of toxicities}
#'   \item{siqr_tox}{Semi-interquartile range of the number of toxicities}
#' }
#' @import lattice
#' @examples
#' df <- run_simulation_power(num_rep = 100)
#' head(df)
#' @export
run_simulation_power <- function(num_rep = 1000,
                                     seed = 1234,
                                     save_plot = FALSE,
                                     p0 = 0.4,
                                     theta_0 = 2.7,
                                     theta = 3,
                                     N = 24,
                                     n_initial = 1,
                                     q_0 = 0.05,
                                     q_2 = 0.5,
                                     q_1 = 0.9,
                                     p_tox_init_3_3 = 0.02,
                                     delta_dosis_3_3 = 0.055,
                                     fixed_optimal_dose = 0.2032) {
  mtd_true <- p0^(1 / theta)

  # Run 3+3 and CRM simulations
  n_tox_3_3 <- numeric(num_rep)
  mtd_3_3_estimated <- numeric(num_rep)
  mtd_proposal_estimated <- numeric(num_rep)
  n_tox_proposal <- numeric(num_rep)
  n_patients_3_3 <- numeric(num_rep)

  for (i in 1:num_rep) {
    results_3_3 <- power_3_3(seed = seed + i,
                                 theta = theta,
                                 theta_0 = theta_0,
                                 n_initial = 3,
                                 p_tox_init = p_tox_init_3_3,
                                 delta_dosis = delta_dosis_3_3,
                                 show_plot = FALSE)
    n_tox_3_3[i] <- results_3_3$n_toxicities
    mtd_3_3_estimated[i] <- results_3_3$mtd_estimated
    n_patients_3_3[i] <- results_3_3$n_patients

    results_crm <- two_stage_crm_power(seed = seed + i,
                                           p0 = p0,
                                           n_initial = n_initial,
                                           N = N,
                                           theta = theta,
                                           theta_0 = theta_0,
                                           q_0 = q_0,
                                           q_2 = q_2,
                                           q_1 = q_1,
                                           fixed_optimal_dose = fixed_optimal_dose,
                                           show_plot = FALSE)
    mtd_proposal_estimated[i] <- results_crm$mtd_estimated
    n_tox_proposal[i] <- results_crm$n_toxicities
  }

  # Run two-stage CRM simulations


  # Visualization using lattice
  df <- data.frame(
    value = c(mtd_3_3_estimated, mtd_proposal_estimated),
    group = factor(rep(c("3+3", "CRM2s"), each = num_rep))
  )

  my_breaks <- seq(min(df$value, na.rm = TRUE)-delta_dosis_3_3/2,
                   max(df$value, na.rm = TRUE)+3*delta_dosis_3_3/2,
                   by = delta_dosis_3_3)

  if (save_plot) {
    if (!dir.exists("PLOTS")) {
      dir.create("PLOTS")
    }
    pdf("PLOTS/compare_MTD_power.pdf", width = 5.6, height = 4.2)
  }

  # --- Plot 1: Estimated MTD ---
  print(lattice::bwplot(value ~ group, data = df, coef = 0, pch = "|",
                  main = "",
                  xlab = "Method",
                  ylab = "Estimated MTD",
                  ylim = range(0.7 * min(df$value, na.rm = TRUE),
                               1.1 * max(df$value, na.rm = TRUE)),
                  par.settings = list(
                    box.rectangle = list(col = "black"),
                    box.umbrella = list(col = "transparent"),
                    plot.symbol = list(col = "black"),
                    box.dot = list(col = "black"),
                    superpose.symbol = list(col = "black")
                  ),
                  panel = function(...) {
                    panel.hanoi(breaks = my_breaks, col = "grey", ...)
                    lattice::panel.bwplot(...)
                    lattice::panel.abline(h = mtd_true, lty = 3, col = "red", lwd = 2)
                  }))
  if (save_plot) {
    dev.off()
  }

  df_2 <- data.frame(
    value = c(n_tox_3_3, n_tox_proposal),
    group = factor(rep(c("3+3", "CRM2s"), each = num_rep))
  )

  my_breaks <- seq(0.5, 20.5, by = 1)

  if (save_plot) {
    pdf("PLOTS/compare_tox_power.pdf", width = 5.6, height = 4.2)
  }

  print(lattice::bwplot(value ~ group, data = df_2, coef = 0, pch = "|",
                  main = "",
                  xlab = "Method",
                  ylab = "Number of toxicities",
                  ylim = range( min(df_2$value, na.rm = TRUE)-1,
                               max(df_2$value, na.rm = TRUE)+1),
                  par.settings = list(
                    box.rectangle = list(col = "black"),
                    box.umbrella = list(col = "transparent"),
                    plot.symbol = list(col = "black"),
                    box.dot = list(col = "black"),
                    superpose.symbol = list(col = "black")
                  ),
                  panel = function(...) {
                    panel.hanoi(breaks = my_breaks, col = "grey", ...)
                    lattice::panel.bwplot(...)
                  }))
  if (save_plot) {
    dev.off()
  }

  result_df <- data.frame(
    method = c("3+3", "CRM2s"),
    mean_pat   = c(mean(n_patients_3_3), N),
    median_pat = c(median(n_patients_3_3), N),
    mean_mtd   = c(mean(mtd_3_3_estimated), mean(mtd_proposal_estimated)),
    sd_mtd    = c(sd(mtd_3_3_estimated), sd(mtd_proposal_estimated)),
    median_mtd = c(median(mtd_3_3_estimated), median(mtd_proposal_estimated)),
    min_mtd    = c(min(mtd_3_3_estimated), min(mtd_proposal_estimated)),
    q1_mtd     = c(as.numeric(quantile(mtd_3_3_estimated, 0.25)),
                   as.numeric(quantile(mtd_proposal_estimated, 0.25))),
    q3_mtd     = c(as.numeric(quantile(mtd_3_3_estimated, 0.75)),
                   as.numeric(quantile(mtd_proposal_estimated, 0.75))),
    max_mtd    = c(max(mtd_3_3_estimated), max(mtd_proposal_estimated)),
    siqr_mtd    = c(IQR(mtd_3_3_estimated)/2, IQR(mtd_proposal_estimated)/2),

    mean_tox   = c(mean(n_tox_3_3), mean(n_tox_proposal)),
    sd_tox   = c(sd(n_tox_3_3), sd(n_tox_proposal)),
    median_tox = c(median(n_tox_3_3), median(n_tox_proposal)),
    min_tox    = c(min(n_tox_3_3), min(n_tox_proposal)),
    q1_tox     = c(as.numeric(quantile(n_tox_3_3, 0.25)),
                   as.numeric(quantile(n_tox_proposal, 0.25))),
    q3_tox     = c(as.numeric(quantile(n_tox_3_3, 0.75)),
                   as.numeric(quantile(n_tox_proposal, 0.75))),
    max_tox    = c(max(n_tox_3_3), max(n_tox_proposal)),
    siqr_tox    = c(IQR(n_tox_3_3)/2, IQR(n_tox_proposal)/2)
  )

  return(result_df)
}
