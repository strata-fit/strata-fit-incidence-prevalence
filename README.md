# STRATA-FIT | Incident and Prevalence Analysis (R)

This repository contains the R scripts used to calculate the **incidence** and **prevalence** of Difficult-to-Treat Rheumatoid Arthritis (D2T RA) as part of the **STRATA-FIT** federated analytics framework.

## Purpose

The code in this repository enables each participating center in the STRATA-FIT project to:

-  Load and validate their local data structure
-  Calculate **prevalence** of D2T RA over time
-  Estimate **incidence** using survival analysis techniques (e.g., Kaplan-Meier)
-  Perform **sensitivity analyses** on multiple D2T definitions
-  Prepare harmonized, locally validated output that can be shared in a **privacy-preserving federated setting**

These steps ensure consistency in local preprocessing before contributions are securely transferred to the **federated analytics environment**.




> All patient-level data **remains local** and is never committed or shared. Only aggregated summaries are exported and shared with WP3. 

##  Getting Started

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
source("V5_Incidence&Prevalence_script.R")
```

It is important to run the code in the followin sequence if you wish you run it separatly and no in the full version. 

```r
source("V5_Preprocessing&Imputation.R")
source("V5_D2TDefinitions.R")
source("V5_Incidence&Sensitivity.R")
source("V5_Persistance.R")
source("V5_Prevalence&CoxRegression.R")
```


The output plots and tables will appear in the `/output' on the side. Please upload these files into an excel or within the Synopsis Paper shared. 

## Background

This code supports the **federated survival modeling** goals of the STRATA-FIT initiative, which aims to analyze D2T RA data across multiple centers without ever sharing patient-level data.

All logic is aligned with the federated model — only minimal, pre-aggregated outputs (such as summary statistics or survival probabilities) are generated for sharing.

## Contributing

Each center is encouraged to:

- Adapt input paths and formats as needed
- Validate outputs visually and statistically
- Submit any improvements or extensions via pull request

For substantial changes, please open an issue first to discuss proposed modifications.


---

 **STRATA-FIT** — Stratification of Patients with Difficult-to-Treat Rheumatoid Arthritis Using Federated Learning  
