#' @title Export Simulation Summary to LaTeX
#' @description Writes a LaTeX table with MTD and toxicity summaries using cat(), without external dependencies.
#'
#' @param df A data frame from simulate_across_n_initial()[[model]].
#' @param file_path Path to the output .tex file.
#' @export
#' @examples
#' result_list <- simulate_across_n_initial()
#' export_simulation_table_manual(result_list$potential, "table_potential.tex")
export_simulation_table_manual <- function(df, file_path = "simulation_results.tex") {
  con <- file(file_path, open = "wt")

  cat("\\begin{tabular}{llcccc}\n", file = con)
  cat("\\toprule\n", file = con)
  cat("Method & $n_{\\text{initial}}$ & Mean (Var) MTD & Median (IQR) MTD & Mean (Var) Tox & Median (IQR) Tox \\\\\n", file = con)
  cat("\\midrule\n", file = con)

  df_crm <- df[df$method == "2stage", ]
  df_3_3 <- df[df$method == "3+3", ]

  cat(sprintf("\\multirow{%d}{*}{CRM}", nrow(df_crm)), file = con)

  for (i in seq_len(nrow(df_crm))) {
    row <- df_crm[i, ]
    prefix <- if (i > 1) " &" else " &"
    cat(sprintf(
      "%s %d & %s & %s & %s & %s \\\\\n",
      prefix,
      row$n_initial,
      sprintf("%.3f (%.3f)", row$mean_mtd, row$var_mtd),
      sprintf("%.3f (%.3f)", row$median_mtd, row$iqr_mtd),
      sprintf("%.3f (%.3f)", row$mean_tox, row$var_tox),
      sprintf("%.3f (%.3f)", row$median_tox, row$iqr_tox)
    ), file = con)
  }

  cat("\\midrule\n", file = con)
  cat(sprintf(
    "3+3 & -- & %s & %s & %s & %s \\\\\n",
    sprintf("%.3f (%.3f)", df_3_3$mean_mtd, df_3_3$var_mtd),
    sprintf("%.3f (%.3f)", df_3_3$median_mtd, df_3_3$iqr_mtd),
    sprintf("%.3f (%.3f)", df_3_3$mean_tox, df_3_3$var_tox),
    sprintf("%.3f (%.3f)", df_3_3$median_tox, df_3_3$iqr_tox)
  ), file = con)

  cat("\\bottomrule\n", file = con)
  cat("\\end{tabular}\n", file = con)

  close(con)
  message("✅ LaTeX table written to ", normalizePath(file_path))
}
