# CRM2s: A two-stage algorithm for the Continuous Reassessment Method

**CRM2s** is an R package designed to simulate and compare dose-escalation methods in phase I clinical trials. It includes both the traditional **3+3 design** and a **two-stage Continual Reassessment Method (CRM)**, implemented for two types of dose-toxicity models: **logistic** and **power**.

------------------------------------------------------------------------

## 📦 Features

-   Simulate individual clinical trials using:
    -   3+3 design (logistic or power model)
    -   CRM2s: Two-stage CRM based on C-optimality
-   Compare methods based on MTD estimation and toxicity outcomes
-   Visualize patient-level dose allocations and outcomes

------------------------------------------------------------------------

## 🛠 Installation

### From GitHub

``` r
devtools::install_github("alvarocia/CRM2s")
```

------------------------------------------------------------------------

## 🚀 Getting Started

### Run a two-stage CRM simulation (power model)

``` r
library(CRM2s)

res <- two_stage_crm_power(show_plot = TRUE)
print(res$mtd_estimated)
```

### Compare 3+3 and CRM2s designs under the logistic model

``` r
result <- run_simulation_logistic(num_rep = 100)
print(result)
```

### Export simulation results to a LaTeX file

``` r
res <- simulate_across_n_initial()
export_simulation_table_manual(res$logistic, "RESULTS/table_logistic.tex")
```

------------------------------------------------------------------------

## 📊 Function Overview

| Function                     | Purpose                                                           |
|-----------------------------|-------------------------------------------------------------------|
| `logistic_3_3()`            | Run one 3+3 trial under logistic dose-toxicity model              |
| `power_3_3()`               | Run one 3+3 trial under power model                               |
| `two_stage_crm_logistic()`  | Run one two-stage CRM simulation using logistic model             |
| `two_stage_crm_power()`     | Run one two-stage CRM simulation using power model                |
| `run_simulation_logistic()` | Run multiple simulations (logistic) and compare designs           |
| `run_simulation_power()`    | Run multiple simulations (power) and compare designs              |
| `simulate_across_n_initial()` | Vary `n_initial` and summarize results across configurations    |
| `export_simulation_table_manual()` | Export formatted LaTeX summary table                      |

---

## 📁 Package Structure

```         
CRM2s/
├── R/
│   ├── power_3_3.R
│   ├── logistic_3_3.R
│   ├── two_stage_crm_*.R
│   ├── run_simulation_*.R
│   ├── simulate_across_n_initial.R
│   └── plot_functions.R
├── tools/
├── PLOTS/
├── RESULTS/
├── tests/
├── CRM2s.pdf
├── DESCRIPTION
├── NAMESPACE
├── LICENSE.md
└── README.md
```

---

## 📄 License

This package is released under the terms of the [GNU General Public License v3.0](LICENSE.md).

------------------------------------------------------------------------

## 👤 Author

Developed by Alvaro Cia-Mina For contributions or issues, visit: <https://github.com/alvarocia/CRM2s>
