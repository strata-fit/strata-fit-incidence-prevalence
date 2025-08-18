# STRATA-FIT | Incident and Prevalence Analysis (R)

This repository contains the R scripts used to calculate the **incidence** and **prevalence** of Difficult-to-Treat Rheumatoid Arthritis (D2T RA) as part of the **STRATA-FIT** federated analytics framework.

## 📦 Purpose

The code in this repository enables each participating center in the STRATA-FIT project to:

- ✅ Load and validate their local data structure
- 📊 Calculate **prevalence** of D2T RA over time
- 📈 Estimate **incidence** using survival analysis techniques (e.g., Kaplan-Meier)
- 🔍 Perform **sensitivity analyses** on multiple D2T definitions
- 🛠 Prepare harmonized, locally validated output that can be shared in a **privacy-preserving federated setting**

These steps ensure consistency in local preprocessing before contributions are securely transferred to the **federated analytics environment**.

## 🧪 Key Features

- Prevalence analysis stratified by calendar year and D2T definition
- Kaplan-Meier survival analysis using imputed datasets
- Support for multiple definitions of D2T RA (base and sensitivity variants)
- Clean visualizations for local review
- Exportable data summaries for federation

## 🗂 Repository Structure

```
/scripts
  ├── 01_prevalence_analysis.R
  ├── 02_incidence_analysis.R
  ├── utils/
  │     └── helper_functions.R
/data
  └── (User-provided input files - not committed)
/outputs
  └── (Local output summaries, plots, tables)
/README.md
```

> 📌 All patient-level data **remains local** and is never committed or shared. Only aggregated summaries are exported for federation.

## 🚀 Getting Started

### 1. Clone the repository

```bash
git clone https://github.com/YOUR-USERNAME/strata-fit-prevalence.git
cd strata-fit-prevalence
```

### 2. Install required R packages

Open R and run:

```r
install.packages(c("dplyr", "ggplot2", "survival", "survminer", "readr", "purrr"))
```

### 3. Run local checks and analysis

Update your `.R` scripts to load your center’s data, then execute:

```r
source("scripts/01_prevalence_analysis.R")
source("scripts/02_incidence_analysis.R")
```

The output plots and tables will appear in the `/outputs` folder.

## 🧠 Background

This code supports the **federated survival modeling** goals of the STRATA-FIT initiative, which aims to analyze D2T RA data across multiple centers without ever sharing patient-level data.

All logic is aligned with the federated model — only minimal, pre-aggregated outputs (such as summary statistics or survival probabilities) are generated for sharing.

## 🤝 Contributing

Each center is encouraged to:

- Adapt input paths and formats as needed
- Validate outputs visually and statistically
- Submit any improvements or extensions via pull request

For substantial changes, please open an issue first to discuss proposed modifications.

## 📄 License

This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.

---

🧬 **STRATA-FIT** — Stratification of Patients with Difficult-to-Treat Rheumatoid Arthritis Using Federated Learning  
