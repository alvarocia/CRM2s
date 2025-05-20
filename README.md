# CRMC: C-optimality based Continual Reassessment Method

**CRMC** is an R package designed to simulate and compare dose-escalation methods in phase I clinical trials. It includes both the traditional **3+3 design** and a **two-stage Continual Reassessment Method (CRM)**, implemented for two types of dose-toxicity models: **logistic** and **potential**.

------------------------------------------------------------------------

## 📦 Features

-   Simulate individual clinical trials using:
    -   3+3 design (logistic or potential model)
    -   CRMC: Two-stage CRM based on C-optimality
-   Compare methods based on MTD estimation and toxicity outcomes
-   Visualize patient-level dose allocations and outcomes

------------------------------------------------------------------------

## 🛠 Installation

### From GitHub

``` r
devtools::install_github("alvarocia/CRMC")
```

------------------------------------------------------------------------

## 🚀 Getting Started

### Run a two-stage CRM simulation (potential model)

``` r
library(CRMC)

res <- two_stage_crm_potential(show_plot = TRUE)
print(res$mtd_estimated)
```

### Compare 3+3 and CRM designs under the logistic model

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

| Function | Purpose |
|--------------------------|----------------------------------------------|
| `logistic_3_3()` | Run one 3+3 trial under logistic dose-toxicity model |
| `potential_3_3()` | Run one 3+3 trial under potential model |
| `two_stage_crm_logistic()` | Run one two-stage CRM simulation using logistic model |
| `two_stage_crm_potential()` | Run one two-stage CRM simulation using potential model |
| `run_simulation_logistic()` | Run multiple simulations (logistic) and compare designs |
| `run_simulation_potential()` | Run multiple simulations (potential) and compare designs |
| `simulate_across_n_initial()` | Vary `n_initial` and summarize results across configurations |
| `export_simulation_table_manual()` | Export formatted LaTeX summary table |

------------------------------------------------------------------------

## 📁 Package Structure

```         
CRMC/
├── R/
│   ├── potential_3_3.R
│   ├── logistic_3_3.R
│   ├── two_stage_crm_*.R
│   ├── run_simulation_*.R
│   ├── simulate_across_n_initial.R
│   └── plot_functions.R
├── PLOTS/
├── RESULTS/
├── tests/
├── CRMC.pdf
├── DESCRIPTION
├── NAMESPACE
└── README.md
```

------------------------------------------------------------------------

## 📄 License

This package is released under the GLP License.

------------------------------------------------------------------------

## 👤 Author

Developed by Alvaro Cia-Mina For contributions or issues, visit: <https://github.com/alvarocia/CRMC>
