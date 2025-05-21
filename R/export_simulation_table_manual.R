#' @title Export Simulation Summary to LaTeX
#' @description Writes a LaTeX table with MTD and toxicity summaries using cat(), without external dependencies.
#'
#' @param df A data frame with simulation summaries for one model (e.g., \code{result_list$power}) as returned by \code{simulate_across_n_initial()}.
#' @param file_path Path to the output .tex file.
#' @export
#' @examples
#' result_list <- simulate_across_n_initial()
#' export_simulation_table_manual(result_list$power, "table_power.tex")
#' export_simulation_table_manual(result_list$logistic, "table_logistic.tex")
export_simulation_table_manual <- function(df, file_path = "RESULTS/simulation_results.tex") {
  con <- file(file_path, open = "wt")

  cat("\\begin{tabular}{llcccc}\n", file = con)
  cat("\\toprule\n", file = con)
  cat("Method & $n_{\\text{initial}}$ & Mean MTD (SD)  & Median MTD (SIQR)  & Mean Tox (SD)  & Median Tox (SIQR)  \\\\\n", file = con)
  cat("\\midrule\n", file = con)

  df_crm <- df[df$method == "CRM2s", ]
  df_3_3 <- df[df$method == "3+3", ]

  cat(sprintf("\\multirow{%d}{*}{CRM2s}", nrow(df_crm)), file = con)

  for (i in seq_len(nrow(df_crm))) {
    row <- df_crm[i, ]
    prefix <- if (i > 1) " &" else " &"
    cat(sprintf(
      "%s %d & %s & %s & %s & %s \\\\\n",
      prefix,
      row$n_initial,
      sprintf("%.3f (%.3f)", row$mean_mtd, row$sd_mtd),
      sprintf("%.3f (%.3f)", row$median_mtd, row$siqr_mtd),
      sprintf("%.3f (%.3f)", row$mean_tox, row$sd_tox),
      sprintf("%.1f (%.2f)", row$median_tox, row$siqr_tox)
    ), file = con)
  }

  cat("\\midrule\n", file = con)
  cat(sprintf(
    "3+3 & -- & %s & %s & %s & %s \\\\\n",
    sprintf("%.3f (%.3f)", df_3_3$mean_mtd, df_3_3$sd_mtd),
    sprintf("%.3f (%.3f)", df_3_3$median_mtd, df_3_3$siqr_mtd),
    sprintf("%.3f (%.3f)", df_3_3$mean_tox, df_3_3$sd_tox),
    sprintf("%.1f (%.2f)", df_3_3$median_tox, df_3_3$siqr_tox)
  ), file = con)

  cat("\\bottomrule\n", file = con)
  cat("\\end{tabular}\n", file = con)

  close(con)
  message("✅ LaTeX table written to ", normalizePath(file_path))
}
