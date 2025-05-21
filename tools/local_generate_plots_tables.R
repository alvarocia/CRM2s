# Local document. Generate plots and tables
result_list <- simulate_across_n_initial(num_rep=1000)
export_simulation_table_manual(result_list$power, "RESULTS/table_power.tex")
export_simulation_table_manual(result_list$logistic, "RESULTS/table_logistic.tex")


run_simulation_logistic(save_plot = TRUE, n_initial = 2)
run_simulation_power(save_plot = TRUE, n_initial = 2)
