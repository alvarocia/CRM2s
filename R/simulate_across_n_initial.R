#' @title Simulation Study for Different n_initial Values (Potential and Logistic Models)
#' @description Runs simulations using run_simulation_potential() and run_simulation_logistic() with varying n_initial values and summarizes key results.
#'
#' @param num_rep Number of repetitions for each simulation. Default is 500.
#' @param seed Random seed. Default is 1234.
#'
#' @return A list of two data frames: \code{$potential} and \code{$logistic}, each containing results for CRMC at n_initial = 1:4 and 3+3.
#' @export
#' @examples
#' result_list <- simulate_across_n_initial()
#' head(result_list$potential)
simulate_across_n_initial <- function(num_rep = 500, seed = 1234) {
  # ----------------------------
  # POTENTIAL MODEL
  # ----------------------------
  crm_potential_results <- lapply(1:4, function(n_init) {
    df <- run_simulation_potential(
      num_rep = num_rep,
      seed = seed,
      n_initial = n_init,
      save_plot = FALSE
    )
    crm_row <- df[df$method == "CRMC", ]
    crm_row$n_initial <- n_init
    return(crm_row)
  })
  crm_potential_df <- do.call(rbind, crm_potential_results)

  df_3_3_pot <- run_simulation_potential(
    num_rep = num_rep,
    seed = seed,
    n_initial = 3,
    save_plot = FALSE
  )
  row_3_3_pot <- df_3_3_pot[df_3_3_pot$method == "3+3", ]
  row_3_3_pot$n_initial <- NA

  df_potential <- rbind(crm_potential_df, row_3_3_pot)
  df_potential <- df_potential[, c("method", "n_initial", "mean_pat", "median_pat",
                                   "mean_mtd", "var_mtd", "median_mtd", "iqr_mtd",
                                   "mean_tox", "var_tox", "median_tox", "iqr_tox")]

  # ----------------------------
  # LOGISTIC MODEL
  # ----------------------------
  crm_logistic_results <- lapply(1:4, function(n_init) {
    df <- run_simulation_logistic(
      num_rep = num_rep,
      seed = seed,
      n_initial = n_init,
      save_plot = FALSE
    )
    crm_row <- df[df$method == "CRMC", ]
    crm_row$n_initial <- n_init
    return(crm_row)
  })
  crm_logistic_df <- do.call(rbind, crm_logistic_results)

  df_3_3_log <- run_simulation_logistic(
    num_rep = num_rep,
    seed = seed,
    n_initial = 3,
    save_plot = FALSE
  )
  row_3_3_log <- df_3_3_log[df_3_3_log$method == "3+3", ]
  row_3_3_log$n_initial <- NA

  df_logistic <- rbind(crm_logistic_df, row_3_3_log)
  df_logistic <- df_logistic[, c("method", "n_initial", "mean_pat", "median_pat",
                                 "mean_mtd", "var_mtd", "median_mtd", "iqr_mtd",
                                 "mean_tox", "var_tox", "median_tox", "iqr_tox")]

  # Return both results
  return(list(
    potential = df_potential,
    logistic = df_logistic
  ))
}
