###############################################################################
# Incidence and Prevalence of D2T RA 
# VERSION 6.0
# P.M.J. Welsing & C. Ripepi Jul 2025
###############################################################################
rm(list = ls())

library(dplyr)
library(mice)
library(ggplot2)
library(ggsurvfit)
library(survival)
library(haven)
library(survminer)
library(readr)
library(tidyverse)
library(zoo)
library(lattice)
library(mitml)
library(pan)
library(data.table)
library(broom)
library(patchwork)


###############################################################################
#  Data Selection 
###############################################################################
data_node<- read_csv("STRATA-FIT20250925.csv") # Upload your own data 

data_node <- data_node %>%
  mutate(tsDMARD = case_when(
    is.na(tsDMARD) ~ NA_real_,
    tsDMARD >= 1   ~ 1,
    TRUE           ~ 0
  ))

#Convert CRP from mg/L to mg/dL if needed 
#data_node <- data_node %>%
#  mutate(CRP = CRP / 10)  # 1 mg/dL = 10 mg/L

# Summary per patient
patient_summary <- data_node %>%
  filter(!is.na(Visit_months_from_diagnosis)) %>%
  group_by(pat_ID) %>%
  summarise(
    Female = first(Sex == 1),
    RF_positive = first(RF_positivity == 1),
    aCCP_positive = first(anti_CCP == 1),
    Age_diagnosis = first(Age_diagnosis),
    Year_diagnosis = first(Year_diagnosis),
    FU_years = max(Visit_months_from_diagnosis, na.rm = TRUE) / 12,
    n_visits = n(),
    DAS28_mean = mean(DAS28, na.rm = TRUE),
    SJC28_mean = mean(SJC28, na.rm = TRUE),
    TJC28_mean = mean(TJC28, na.rm = TRUE),
    ESR_mean = mean(ESR, na.rm = TRUE),
    CRP_mean = mean(CRP, na.rm = TRUE),
    VAS_patient_mean = mean(Pat_global, na.rm = TRUE),
    VAS_physician_mean = mean(Ph_global, na.rm = TRUE),
    # GC variables
    GC_pct = mean(GC == 1, na.rm = TRUE) * 100,
    GC_highdose_pct = mean(GC_dose >= 7.5, na.rm = TRUE) * 100,
    # Previous DMARDs
    N_prev_csDMARD = first(N_prev_csDMARD),
    N_prev_bDMARD = first(N_prev_bDMARD),
    N_prev_tsDMARD = first(N_prev_tsDMARD)
  ) %>%
  mutate(
    visits_per_year = ifelse(is.na(FU_years) | FU_years == 0, NA, n_visits / FU_years)
  )

# --- Baseline features & "followed from diagnosis (0–6 months)" ---
dx_window_months <- 6

baseline <- data_node %>%
  filter(!is.na(Visit_months_from_diagnosis)) %>%
  arrange(pat_ID, Visit_months_from_diagnosis) %>%
  group_by(pat_ID) %>%
  summarise(
    first_visit_months = first(Visit_months_from_diagnosis),
    Disease_duration_first_mo = first(Visit_months_from_diagnosis),   # months from dx to 1st visit
    DAS28_first = if (all(is.na(DAS28))) NA_real_ else first(na.omit(DAS28)),
    Symptom_duration_first_mo = {  # Symptom duration column best-effort (only if you have one)
      cols <- names(pick(everything()))
      val <- if ("Symptom_duration_first_mo" %in% cols) {
        Symptom_duration_first_mo} else if ("Symptom_duration_months" %in% cols) {
        Symptom_duration_months} else if ("Symptom_duration_years" %in% cols) {
        Symptom_duration_years * 12} else {
        NA_real_}
      if (all(is.na(val))) NA_real_ else first(na.omit(val))
    },
    N_prev_csDMARD = first(N_prev_csDMARD),
    N_prev_bDMARD  = first(N_prev_bDMARD),
    N_prev_tsDMARD = first(N_prev_tsDMARD),
    .groups = "drop"
  ) %>%
  mutate(
    Followed_from_dx = !is.na(first_visit_months) &
      first_visit_months >= 0 &
      first_visit_months <= dx_window_months
  )

# Percentages Yes/No
follow_yes <- sum(baseline$Followed_from_dx %in% TRUE)
follow_no  <- sum(baseline$Followed_from_dx %in% FALSE)
follow_den <- follow_yes + follow_no
pct_yes <- if (follow_den > 0) round(follow_yes / follow_den * 100, 1) else NA_real_
pct_no  <- if (follow_den > 0) round(follow_no  / follow_den * 100, 1) else NA_real_

# Subgroups
yes_grp <- baseline %>% dplyr::filter(Followed_from_dx)
no_grp  <- baseline %>% dplyr::filter(!Followed_from_dx)

# --- Build summary table --- Addition of Medication and distiction of follow-up groups
Table1 <- list(
  "Female, n (%)" = {
    tab <- table(patient_summary$Female)
    paste0(tab[2], " (", round(tab[2] / sum(tab) * 100, 1), "%)") },
  "RF-positive, n (%)" = {
    tab <- table(patient_summary$RF_positive)
    paste0(tab[2], " (", round(tab[2] / sum(tab) * 100, 1), "%)") },
  "aCCP-positive, n (%)" = {
    tab <- table(patient_summary$aCCP_positive)
    paste0(tab[2], " (", round(tab[2] / sum(tab) * 100, 1), "%)")},
  "Age at diagnosis" = {
    m <- mean(patient_summary$Age_diagnosis, na.rm = TRUE)
    s <- sd(patient_summary$Age_diagnosis, na.rm = TRUE)
    paste0(round(m, 1), " (", round(s, 1), ")")},
  "Initial calendar year of follow-up, median (Q1, Q3)" = {
    q <- quantile(patient_summary$Year_diagnosis, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
    paste0(round(q[2], 0), " (", round(q[1], 0), "–", round(q[3], 0), ")")},
  "Follow-up (years)" = {
    q <- quantile(patient_summary$FU_years, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
    paste0(round(q[2], 1), " (", round(q[1], 1), "–", round(q[3], 1), ")")},
  "# visits per patient per year" = {
    m <- mean(patient_summary$visits_per_year, na.rm = TRUE)
    s <- sd(patient_summary$visits_per_year, na.rm = TRUE)
    paste0(round(m, 1), " (", round(s, 1), ")")},
  "DAS28" = {
    m <- mean(patient_summary$DAS28_mean, na.rm = TRUE)
    s <- sd(patient_summary$DAS28_mean, na.rm = TRUE)
    paste0(round(m, 2), " (", round(s, 2), ")")},
  "28SJC" = {
    m <- mean(patient_summary$SJC28_mean, na.rm = TRUE)
    s <- sd(patient_summary$SJC28_mean, na.rm = TRUE)
    paste0(round(m, 2), " (", round(s, 2), ")") },
  "28TJC" = {
    m <- mean(patient_summary$TJC28_mean, na.rm = TRUE)
    s <- sd(patient_summary$TJC28_mean, na.rm = TRUE)
    paste0(round(m, 2), " (", round(s, 2), ")")},
  "ESR, median (Q1, Q3)" = {
    q <- quantile(patient_summary$ESR_mean, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
    paste0(round(q[2], 1), " (", round(q[1], 1), "–", round(q[3], 1), ")")},
  "CRP, median (Q1, Q3)" = {
    q <- quantile(patient_summary$CRP_mean, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
    paste0(round(q[2], 1), " (", round(q[1], 1), "–", round(q[3], 1), ")")},
  "VAS patient" = {
    m <- mean(patient_summary$VAS_patient_mean, na.rm = TRUE)
    s <- sd(patient_summary$VAS_patient_mean, na.rm = TRUE)
    paste0(round(m, 1), " (", round(s, 1), ")")},
  "VAS physician" = {
    m <- mean(patient_summary$VAS_physician_mean, na.rm = TRUE)
    s <- sd(patient_summary$VAS_physician_mean, na.rm = TRUE)
    paste0(round(m, 1), " (", round(s, 1), ")")  },
  "GC use, % of visits" = {
    m <- mean(patient_summary$GC_pct, na.rm = TRUE)
    s <- sd(patient_summary$GC_pct, na.rm = TRUE)
    paste0(round(m, 1), "% (", round(s, 1), "%)")},
  "GC ≥7.5mg, % of visits" = {
    m <- mean(patient_summary$GC_highdose_pct, na.rm = TRUE)
    s <- sd(patient_summary$GC_highdose_pct, na.rm = TRUE)
    paste0(round(m, 1), "% (", round(s, 1), "%)")},
  "Previous csDMARDs, mean (SD)" = {
    m <- mean(patient_summary$N_prev_csDMARD, na.rm = TRUE)
    s <- sd(patient_summary$N_prev_csDMARD, na.rm = TRUE)
    paste0(round(m,1), " (", round(s,1), ")")},
  "Previous bDMARDs, mean (SD)" = {
    m <- mean(patient_summary$N_prev_bDMARD, na.rm = TRUE)
    s <- sd(patient_summary$N_prev_bDMARD, na.rm = TRUE)
    paste0(round(m,1), " (", round(s,1), ")")},
  "Previous tsDMARDs, mean (SD)" = {
    m <- mean(patient_summary$N_prev_tsDMARD, na.rm = TRUE)
    s <- sd(patient_summary$N_prev_tsDMARD, na.rm = TRUE)
    paste0(round(m,1), " (", round(s,1), ")")}
)
Table1 <- c(
  "Total number of patients (n)" = n_distinct(data_node$pat_ID),
  Table1
)  # Adding total cohort number 

# --- Append dynamic-name rows safely ---
nm_yes <- sprintf("Followed from diagnosis (0-%d months), n (%%) - Yes", dx_window_months)
nm_no  <- sprintf("Followed from diagnosis (0-%d months), n (%%) - No",  dx_window_months)

Table1[[nm_yes]] <- if (!is.na(pct_yes)) sprintf("%d (%.1f%%)", follow_yes, pct_yes) else "NA"
Table1[[nm_no ]] <- if (!is.na(pct_no )) sprintf("%d (%.1f%%)", follow_no,  pct_no ) else "NA"

Table1[["Symptom duration at diagnosis (months), median (Q1, Q3) [Since Diagnosis]"]] <- {
  if (nrow(yes_grp) == 0) "NA" else {
    q <- quantile(yes_grp$Symptom_duration_first_mo, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
    sprintf("%.1f (%.1f-%.1f)", q[2], q[1], q[3])
  }
}
Table1[["DAS28 at first visit, mean (SD) [Since Diagnosis]"]] <- {
  if (nrow(yes_grp) == 0) "NA" else {
    m <- mean(yes_grp$DAS28_first, na.rm = TRUE); s <- sd(yes_grp$DAS28_first, na.rm = TRUE)
    sprintf("%.2f (%.2f)", m, s)
  }
}
Table1[["Disease duration at 1st visit (months), median (Q1, Q3) [Intake]"]] <- {
  if (nrow(no_grp) == 0) "NA" else {
    q <- quantile(no_grp$Disease_duration_first_mo, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
    sprintf("%.1f (%.1f-%.1f)", q[2], q[1], q[3])
  }
}
Table1[["Previous csDMARDs, mean (SD) [Intake]"]] <- {
  if (nrow(no_grp) == 0) "NA" else {
    m <- mean(no_grp$N_prev_csDMARD, na.rm = TRUE); s <- sd(no_grp$N_prev_csDMARD, na.rm = TRUE)
    sprintf("%.1f (%.1f)", m, s)
  }
}
Table1[["Previous bDMARDs, mean (SD) [Intake]"]] <- {
  if (nrow(no_grp) == 0) "NA" else {
    m <- mean(no_grp$N_prev_bDMARD, na.rm = TRUE); s <- sd(no_grp$N_prev_bDMARD, na.rm = TRUE)
    sprintf("%.1f (%.1f)", m, s)
  }
}
Table1[["Previous tsDMARDs, mean (SD) [Intake]"]] <- {
  if (nrow(no_grp) == 0) "NA" else {
    m <- mean(no_grp$N_prev_tsDMARD, na.rm = TRUE); s <- sd(no_grp$N_prev_tsDMARD, na.rm = TRUE)
    sprintf("%.1f (%.1f)", m, s)
  }
}
#### 
# --- Missingness of disease activity markers (per-visit) ---

# Variables you want missingness for
disease_markers <- c("DAS28", "SJC28", "TJC28", "ESR", "CRP", "Pat_global", "Ph_global")

# Keep only visits with at least one disease activity marker present
vis_with_marker <- data_node %>%
  dplyr::filter(
    dplyr::if_any(dplyr::all_of(disease_markers), ~ !is.na(.x))
  )

# Percent missing within those visits
missingness <- vis_with_marker %>%
  summarise(
    DAS28       = mean(is.na(DAS28))       * 100,
    `28SJC`     = mean(is.na(SJC28))       * 100,
    `28TJC`     = mean(is.na(TJC28))       * 100,
    ESR         = mean(is.na(ESR))         * 100,
    CRP         = mean(is.na(CRP))         * 100,
    `VAS patient`   = mean(is.na(Pat_global)) * 100,
    `VAS physician` = mean(is.na(Ph_global))  * 100
  )
# --- Add Missingness* (%) block to Table 1 ---

# Header row (just a label, optional)
Table1[["Missingness* (%)"]] <- ""

Table1[["  DAS28"]] <- sprintf("%.1f", missingness$DAS28)
Table1[["  28SJC"]] <- sprintf("%.1f", missingness$`28SJC`)
Table1[["  28TJC"]] <- sprintf("%.1f", missingness$`28TJC`)
Table1[["  ESR"]]   <- sprintf("%.1f", missingness$ESR)
Table1[["  CRP"]]   <- sprintf("%.1f", missingness$CRP)
Table1[["  VAS patient"]]   <- sprintf("%.1f", missingness$`VAS patient`)
Table1[["  VAS physician"]] <- sprintf("%.1f", missingness$`VAS physician`)

# Convert to data frame
summary_table1 <- data.frame(Variable = names(Table1), Value = unlist(Table1), row.names = NULL)
print(summary_table1)

# Reorder data_node
data_node <- data_node[order(data_node$pat_ID, data_node$Visit_months_from_diagnosis),]

# Cleanup and patient count
rm(patient_summary, summary_table1)
n_distinct(data_node$pat_ID)

###############################################################################
#  Create derived variables needed for D2T RA definition and analyses 
###############################################################################
data_node <- data_node %>%
  group_by(pat_ID) %>%
  arrange(Visit_months_from_diagnosis, .by_group = TRUE) %>%
  mutate(
    cum_csDMARD1 = cumsum(!duplicated(csDMARD1) & !is.na(csDMARD1)),
    cum_csDMARD2 = cumsum(!duplicated(csDMARD2) & !is.na(csDMARD2)),
    cum_csDMARD3 = cumsum(!duplicated(csDMARD3) & !is.na(csDMARD3)),
    cum_bDMARD = cumsum(!duplicated(bDMARD) & !is.na(bDMARD)),
    cum_tsDMARD = cumsum(!duplicated(tsDMARD) & !is.na(tsDMARD)),
    cum_btsDMARD = cum_bDMARD + cum_tsDMARD +
      coalesce(as.integer(N_prev_bDMARD), 0L) +
      coalesce(as.integer(N_prev_tsDMARD), 0L),
    cum_GC = cumsum(!duplicated(GC) & !is.na(GC)),
    cum_csDMARD = max(cum_csDMARD1, cum_csDMARD2, cum_csDMARD3, N_prev_csDMARD, na.rm= TRUE),
    cum_btsDMARDmin = cummin(cum_btsDMARD)) %>% # Define cum_btsDMARDmin (cum_btsDMARDs at start follow-up) for left/period censoring 
  ungroup()

#define current calender year variable
data_node$current_year <- floor(data_node$Year_diagnosis + (data_node$month_diagnosis/12) +  (data_node$Visit_months_from_diagnosis)/12)

data_node <- data_node %>%
  group_by(pat_ID) %>%
  arrange(Visit_months_from_diagnosis) %>%
  mutate(
    # mark when cum_btsDMARD increases (a new counted b/tsDMARD start)
    new_DMARD_start = if_else( cum_btsDMARD > lag(cum_btsDMARD, default = 0),
                               Visit_months_from_diagnosis, as.numeric(NA))) %>%
  fill(new_DMARD_start, .direction = "down") %>%        # carry start time forward
  mutate(delta_DMARD_time = Visit_months_from_diagnosis - new_DMARD_start) %>%
  select(-new_DMARD_start) %>%                           # drop helper
  ungroup()

#define current age / age category
data_node$current_age <- round(data_node$Age_diagnosis + (data_node$Visit_months_from_diagnosis/12) + 0.5) 
#current age categories
data_node$current_age_cat <- factor(
  ifelse(data_node$current_age < 40, '< 40',
         ifelse(data_node$current_age < 50, '40-50',
                ifelse(data_node$current_age < 60, '50-60',
                       ifelse(data_node$current_age < 70, '60-70',
                              ifelse(data_node$current_age < 80, '70-80',
                                     ifelse(data_node$current_age >= 80, '> 80', NA)))))),
  levels = c('< 40', '40-50', '50-60', '60-70', '70-80', '> 80'))

# define MOA change & Identify increases in cumulative MOA
data_node <- data_node %>%
  group_by(pat_ID) %>%
  arrange(Visit_months_from_diagnosis) %>%
  mutate(changed_MOA = ifelse(cum_btsDMARD > lag(cum_btsDMARD, default = 0), 1, 0)) %>%
  ungroup()

###############################################################################
#  Create imputed DAS28, components, and VAS_phys
###############################################################################
imputation_vars <- c("DAS28", "ESR", "CRP", "TJC28", "SJC28", "Pat_global", "Ph_global")  # Variables to be imputed
auxiliary_vars <- c("Visit_months_from_diagnosis", "Age_diagnosis", "Sex", "changed_MOA")  # Auxiliary variables (used for imputation, not imputed)
data_for_imputation <- data_node %>%
  select(all_of(c(imputation_vars, auxiliary_vars)))

#Create Predictor Matrix (Defines relationships between variables)
predictor_matrix <- make.predictorMatrix(data_for_imputation)

# Set auxiliary variables as predictors (but not to be imputed)
predictor_matrix[auxiliary_vars, ] <- 1  # Auxiliary variables predict other variables
predictor_matrix[, auxiliary_vars] <- 1  # Auxiliary variables can be used for imputation
predictor_matrix[auxiliary_vars, auxiliary_vars] <- 0  # Auxiliary variables themselves are not imputed

# Define Imputation Methods
method_vector <- rep("", ncol(data_for_imputation))  # Default to no imputation
names(method_vector) <- colnames(data_for_imputation)  # Ensure names match dataset columns

# Apply "pmm" method only to variables that need imputation
method_vector[imputation_vars] <- "pmm"  # Predictive Mean Matching for imputation

# Perform Multiple Imputation Using MICE
imputed_data <- mice(
  data_for_imputation,
  m = 10,                   # Number of imputed datasets
  method = method_vector,   # Apply "pmm" to imputed variables only
  predictorMatrix = predictor_matrix, 
  maxit = 10,               # More iterations for better convergence
  printFlag = TRUE  )        # Show progress

#Extract One of the Imputed Datasets 
df_imputed <- complete(imputed_data, 1)
df_imputed2 <- complete(imputed_data, 2)
df_imputed3 <- complete(imputed_data, 3)
df_imputed4 <- complete(imputed_data, 4)
df_imputed5 <- complete(imputed_data, 5)
df_imputed6 <- complete(imputed_data, 6)
df_imputed7 <- complete(imputed_data, 7)
df_imputed8 <- complete(imputed_data, 8)
df_imputed9 <- complete(imputed_data, 9)
df_imputed10 <- complete(imputed_data, 10)

# Apply the operations to each imputed dataset
imputed_list <- lapply(1:10, function(i) {
  df <- complete(imputed_data, i) %>%
    select(c("ESR", "CRP", "TJC28", "SJC28", "Pat_global", "Ph_global")) %>%
    rename(ESR_imp = "ESR",
           CRP_imp = "CRP",
           TJC28_imp = "TJC28",
           SJC28_imp = "SJC28",
           Pat_global_imp = "Pat_global",
           Ph_global_imp = "Ph_global" )
  df$DAS28_imp <- 0.56 * sqrt(df$TJC28_imp) + 0.28 * sqrt(df$SJC28_imp) +
    0.70 * log(df$ESR_imp) + 0.014 * df$Pat_global_imp 
  df$CDAI <- df$TJC28_imp + df$SJC28_imp + (df$Pat_global_imp / 10) + (df$Ph_global_imp / 10)   # Calcualtion of the CDAI
  return(df)})

names(imputed_list) <- paste0("df_imputed", 1:10)

# Wrap post-imputation processing into a loop
processed_list <- lapply(1:10, function(i) {
  # Combine with original variables (e.g., IDs, visit time, diagnosis year)
  data_node <- cbind(data_node, imputed_list[[i]])
  
  # Sort by patient ID and visit time
  data_node <- data_node[order(data_node$pat_ID, data_node$Visit_months_from_diagnosis), ]
  
  # Add rolling averages
  data_node <- data_node %>%
    group_by(pat_ID) %>%
    mutate(
      rol_av_DAS28 = rollmean(DAS28_imp, k = 2, fill = NA, align = "right"),
      rol_av_DAS28_ni =rollmean(DAS28, k = 2, fill = NA, align = "right"),
      rol_av_Pat_global = rollmean(Pat_global_imp, k = 2, fill = NA, align = "right"),
      rol_av_Ph_global = rollmean(Ph_global_imp, k = 2, fill = NA, align = "right"),
      rol_av_CRP = rollmean(CRP_imp, k = 2, fill = NA, align = "right")
    ) %>%
    ungroup()
  
  # Compute follow-up from year 2000
  data_node$FU_2000 <- ifelse(
    data_node$Year_diagnosis >= 2000,
    data_node$Visit_months_from_diagnosis,
    data_node$Visit_months_from_diagnosis - ((2000 - data_node$Year_diagnosis) * 12)
  )
  
  return(data_node)
})

names(processed_list) <- paste0("data_node", 1:10)

###############################################################################
#  D2T Criteria 
###############################################################################

#Criterion 1 (treatment failure definition): 
#     - Requires patients to have >2 btsDMARDs, or exactly 2 with ≥6 months follow-up. The purpose is to demonstrate that the second treatment has failed.
#     - This ensures a minimum exposure/number of drugs before further classification.

#Criterion 2 (change of mechanism of action / MOA definition):
#     - Applied only after Criterion 1 is satisfied.
#     - Checks whether there is a change in MOA (e.g., from 3rd to 4th drug). 'changed_MOA' is a binary indicator: it records whether a change occurs at that visit, not how many total changes happened.
#     - No need to re-impose the “≥2 btsDMARD” rule here, since it is already guaranteed by Criterion 1.
#     - If disease activity data are missing, a change in MOA serves as a proxy for increased disease activity.

#Sensitivity with CDAI replacing the DAS28 disease activity index (sens5)

# Apply D2T logic to all 10 datasets
final_list <- lapply(processed_list, function(data_node) {
  data_node <- data_node %>%
    mutate(
      # Criterion 1.  Added time criterium to show failure of the secondary btsDMARD
      D2T_crit1 = ifelse((cum_btsDMARD > 2) | (cum_btsDMARD == 2 & delta_DMARD_time >= 6) , 1, 0),
      D2T_crit1a = ifelse(( cum_btsDMARD > 3) | ( cum_btsDMARD == 3 & delta_DMARD_time >= 6), 1, 0),
      D2T_crit1b = ifelse((( cum_btsDMARD > 2) | (cum_btsDMARD == 2 & delta_DMARD_time >= 6)) & FU_2000 < 120, 1, 0),
      
      
      # Criterion 2.  # Added Elevated CRP as a description of high disease activity     ## ADDITION OF CALCULATED CDAI 
      D2T_crit2 = ifelse((DAS28_imp > 3.2 | CRP_imp > 1.0 |  (changed_MOA == 1 & cum_btsDMARD > 1)) , 1, 0),
      D2T_crit2a = ifelse((rol_av_DAS28 > 3.2 | rol_av_CRP > 1.0 | (changed_MOA == 1 & cum_btsDMARD > 1) ) , 1, 0),
      D2T_crit2b = ifelse((DAS28 > 3.2 | CRP_imp > 1.0 | (changed_MOA == 1 & cum_btsDMARD > 1)) , 1, 0),
      D2T_crit2sens1 = ifelse(DAS28_imp > 3.2 | CRP_imp > 1.0 , 1, 0),
      D2T_crit2sens2 = ifelse(DAS28 > 3.2 | CRP > 1.0 | (changed_MOA == 1 & cum_btsDMARD > 1), 1, 0),
      D2T_crit2sens5 = ifelse((CDAI > 10 | rol_av_CRP > 1.0 | (changed_MOA == 1 & cum_btsDMARD > 1)) , 1, 0),
      D2T_crit2sens6 = ifelse(rol_av_DAS28_ni > 3.2 | CRP > 1.0 | (changed_MOA == 1 & cum_btsDMARD > 1), 1, 0),
      
      
      # Criterion 3
      D2T_crit3 = ifelse(Pat_global_imp > 50 | Ph_global_imp > 50, 1, 0),
      D2T_crit3a = ifelse(Pat_global_imp > 75 | Ph_global_imp > 75, 1, 0),
      D2T_crit3b = ifelse(Pat_global > 50 | Ph_global > 50, 1, 0),
      
      # D2T steps. (Seven steps with more restrictive criteria)
      D2T_step0 = ifelse(D2T_crit1 == 1, 1, 0),
      D2T_step1 = ifelse(D2T_crit1 == 1 & D2T_crit2 == 1, 1, 0),
      D2T_step2 = ifelse(D2T_crit1 == 1 & D2T_crit2 == 1 & D2T_crit3 == 1, 1, 0),
      D2T_step3 = ifelse(D2T_crit1 == 1 & D2T_crit2a & D2T_crit3 == 1, 1, 0),
      D2T_step4 = ifelse(D2T_crit1 == 1 & D2T_crit2a & D2T_crit3a == 1, 1, 0),
      D2T_step5 = ifelse(D2T_crit1 == 1 & D2T_crit2a & D2T_crit3a == 1 & FU_2000 <= 120, 1, 0),
      D2T_step6 = ifelse(D2T_crit1 == 1 & D2T_crit2a & D2T_crit3a == 1 & FU_2000 <= 60, 1, 0),
      
      # Sensitivity analyses
      D2T_RA_sens1 = ifelse(D2T_crit1 == 1 & D2T_crit2sens1 == 1 & D2T_crit3 == 1, 1, 0),
      D2T_RA_sens2 = ifelse(D2T_crit1 == 1 & D2T_crit2sens2 == 1 & D2T_crit3b == 1, 1, 0),
      D2T_RA_sens3 = ifelse(D2T_crit1 == 1 & D2T_crit2a  & D2T_crit3 == 1, 1, 0),
      D2T_RA_sens4 = ifelse(D2T_crit1a == 1 & D2T_crit2a & D2T_crit3a == 1, 1, 0),   # 3rd MOA as a sensitivity analysis -- Removed later
      D2T_RA_sens5 = ifelse(D2T_crit1 == 1 & D2T_crit2sens5 == 1 & D2T_crit3 == 1, 1, 0),   # CDAI Instead of DAS
      D2T_RA_sens6 = ifelse(D2T_crit1 == 1 & D2T_crit2sens6 == 1 & D2T_crit3b == 1, 1, 0),   # No imputation or w/rolling das28
    ) %>%
    
    # Ever-D2T variables (per patient)
    group_by(pat_ID) %>%
    mutate(across(
      starts_with("D2T_step"),
      ~ ifelse(all(is.na(.x)), NA, max(.x, na.rm = TRUE)),
      .names = "{.col}_Ever"
    )) %>%
    ungroup()
  return(data_node)
})

###############################################################################
#  Incidence D2T RA time to first D2T RA
################################################################################
calculate_incidence_per_dataset <- function(data_node, defs = paste0("D2T_step", 0:6), time_point = 60) {
  
  get_incidence_data <- function(def, data_node) {
    tmp <- data_node %>%
      group_by(pat_ID) %>%
      summarise(
        definition = def,
        event = ifelse(all(is.na(.data[[def]])), NA, max(.data[[def]], na.rm = TRUE)),
        .groups = "drop"
      )
    
    tmp_tte <- data_node %>%
      group_by(pat_ID) %>%
      summarise(
        TTE = if (any(.data[[def]] == 1, na.rm = TRUE)) {
          min(Visit_months_from_diagnosis[.data[[def]] == 1], na.rm = TRUE)
        } else {
          NA_real_
        },
        .groups = "drop"
      )
    
    tmp_fu <- data_node %>%
      group_by(pat_ID) %>%
      summarise(
        minFU = min(Visit_months_from_diagnosis, na.rm = TRUE),
        maxFU = max(Visit_months_from_diagnosis, na.rm = TRUE),
        Year_diagnosis = first(Year_diagnosis),
        cum_btsDMARDmin = first(cum_btsDMARDmin, na.rm = TRUE),
        Age_diagnosis = first(Age_diagnosis),
        anti_CCP = first(anti_CCP),
        RF_positivity = first(RF_positivity),
        .groups = "drop"
      )
    
    tmp <- tmp %>%
      left_join(tmp_tte, by = "pat_ID") %>%
      left_join(tmp_fu, by = "pat_ID") %>%
      mutate(
        TTE = ifelse(is.na(TTE) | is.infinite(TTE), maxFU, TTE),
        start_date = pmax(2006, Year_diagnosis),
        TTE_2006 = pmax(0, (Year_diagnosis + TTE / 12 - start_date) * 12),
        TTE_start_fu = pmax(0, (Year_diagnosis + minFU / 12 - start_date) * 12),
        cens = case_when(
          cum_btsDMARDmin > 2 | (cum_btsDMARDmin > 1 & event == 1) ~ "interval",
          event == 0 ~ "right",
          TRUE ~ "no"
        ),
        time = case_when(
          cens == "interval" ~ 0,
          TRUE ~ TTE_2006
        ),
        time2 = case_when(
          cens == "interval" ~ TTE_start_fu,
          cens == "right" ~ Inf,
          TRUE ~ TTE_2006
        ),
        status = case_when(
          cens == "interval" ~ 2,
          cens == "no" ~ 1,
          TRUE ~ 0
        )
      )
    
    return(tmp)
  }
  # Collect all definitions into one dataset
  data_inc <- map_dfr(defs, get_incidence_data, data_node = data_node)
  
  # Filter for right-censored and exact events
  data_km <- data_inc %>%
    filter(status %in% c(0, 1)) %>%
    mutate(definition = factor(definition, levels = defs))
  
  # Fit Kaplan-Meier model
  fit_km <- survfit(Surv(time, status) ~ definition, data = data_km)
  
  # Extract 5-year summary
  summary_5yr <- summary(fit_km, times = time_point)
  
  df <- data.frame(
    Definition = gsub("^definition=", "", summary_5yr$strata),
    Time = summary_5yr$time,
    Survival = summary_5yr$surv,
    Lower_CI = summary_5yr$lower,
    Upper_CI = summary_5yr$upper
  ) %>%
    mutate(
      Incidence = 1 - Survival,
      Incidence_CI = sprintf("%.1f%% (%.1f–%.1f%%)",
                             (1 - Survival) * 100,
                             (1 - Upper_CI) * 100,
                             (1 - Lower_CI) * 100)
    )
  
  return(df)
}

# Run on each of the final imputed datasets
incidence_results_list <- lapply(final_list, calculate_incidence_per_dataset)

# Add identifier
for (i in seq_along(incidence_results_list)) {
  incidence_results_list[[i]]$dataset <- paste0("Imputation_", i)
}
# Combine all results
incidence_all <- bind_rows(incidence_results_list)


# Function to apply Rubin's Rules to one definition
rubins_rule_pool <- function(data_subset) {
  m <- nrow(data_subset)
  q_bar <- mean(data_subset$Incidence)  # Pooled point estimate
  u_bar <- mean((data_subset$Upper_CI - data_subset$Lower_CI)^2 / (4 * 1.96^2))  # Average within-imputation variance
  b <- var(data_subset$Incidence)  # Between-imputation variance
  t_var <- u_bar + (1 + 1/m) * b    # Total variance
  se <- sqrt(t_var)
  lower <- q_bar - 1.96 * se
  upper <- q_bar + 1.96 * se
  
  tibble(
    Pooled_Incidence = q_bar,
    Lower_CI = pmax(0, lower),
    Upper_CI = pmin(1, upper),
    Incidence_CI = sprintf("%.1f%% (%.1f–%.1f%%)", q_bar*100, lower*100, upper*100)
  )
}

# Apply Rubin’s Rules to each D2T definition
Table2_5yrIncidence  <- incidence_all %>%
  group_by(Definition) %>%
  group_modify(~ rubins_rule_pool(.x)) %>%
  ungroup()



### Plot the Incidence Curve 

# ---- Step 1: Define incidence extraction function ----
get_incidence_data <- function(def, data_node) {
  tmp <- data_node %>%
    group_by(pat_ID) %>%
    summarise(
      definition = def,
      event = ifelse(all(is.na(.data[[def]])), NA, max(.data[[def]], na.rm = TRUE)),
      .groups = "drop"
    )
  
  tmp_tte <- data_node %>%
    group_by(pat_ID) %>%
    summarise(
      TTE = if (any(.data[[def]] == 1, na.rm = TRUE)) {
        min(Visit_months_from_diagnosis[.data[[def]] == 1], na.rm = TRUE)
      } else {
        NA_real_
      },
      .groups = "drop"
    )
  
  tmp_fu <- data_node %>%
    group_by(pat_ID) %>%
    summarise(
      minFU = min(Visit_months_from_diagnosis, na.rm = TRUE),
      maxFU = max(Visit_months_from_diagnosis, na.rm = TRUE),
      Year_diagnosis = first(Year_diagnosis),
      Age_diagnosis = first(Age_diagnosis),
      RF_positivity = first(RF_positivity),
      anti_CCP = first(anti_CCP),
      cum_btsDMARDmin = first(cum_btsDMARDmin, na.rm = TRUE),
      .groups = "drop"
    )
  
  tmp <- tmp %>%
    left_join(tmp_tte, by = "pat_ID") %>%
    left_join(tmp_fu, by = "pat_ID") %>%
    mutate(
      TTE = ifelse(is.na(TTE) | is.infinite(TTE), maxFU, TTE),
      start_date = pmax(2006, Year_diagnosis),
      TTE_2006 = pmax(0, (Year_diagnosis + TTE / 12 - start_date) * 12),
      TTE_start_fu = pmax(0, (Year_diagnosis + minFU / 12 - start_date) * 12),
      cens = case_when(
        cum_btsDMARDmin > 2 | (cum_btsDMARDmin > 1 & event == 1) ~ "interval",
        event == 0 ~ "right",
        TRUE ~ "no"
      ),
      time = case_when(
        cens == "interval" ~ 0,
        TRUE ~ TTE_2006
      ),
      time2 = case_when(
        cens == "interval" ~ TTE_start_fu,
        cens == "right" ~ Inf,
        TRUE ~ TTE_2006
      ),
      status = case_when(
        cens == "interval" ~ 2,
        cens == "no" ~ 1,
        TRUE ~ 0
      )
    )
  
  # --- NEW: truncate follow-up at 15 years (180 months) ---
  max_months <- 15 * 12
  tmp <- tmp %>%
    mutate(
      status = ifelse(time > max_months, 0L, status),  # events after 15y become censored
      time   = pmin(time, max_months)                  # cap time at 15y
    )

  return(tmp)
}

# ---- Step 2: Run survival for each dataset and definition ----
defs <- paste0("D2T_step", 0:6)

get_surv_df <- function(data_node, def) {
  data_def <- get_incidence_data(def, data_node) %>%
    filter(status %in% c(0, 1))
  
  fit <- survfit(Surv(time, status) ~ 1, data = data_def)
  
  tibble(
    time = fit$time,
    surv = 1 - fit$surv,  # cumulative incidence
    dataset = NA,
    definition = def
  )
}


all_surv_data <- map2_dfr(
  final_list,
  paste0("Imputation_", 1:10),
  function(data_node, ds_name) {
    map_dfr(defs, function(def) {
      df <- get_surv_df(data_node, def)
      df$dataset <- ds_name
      return(df)
    })
  }
)

# ---- Step 3: Pool over time using Rubin-style averaging ----
all_surv_data <- all_surv_data %>%
  mutate(time_rounded = round(time))

pooled_curve <- all_surv_data %>%
  group_by(definition, time_rounded) %>%
  summarise(
    mean_inc = mean(surv, na.rm = TRUE),
    sd_inc = sd(surv, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    lower = pmax(0, mean_inc - 1.96 * sd_inc),
    upper = pmin(1, mean_inc + 1.96 * sd_inc)
  )

# ---- Step 4: Plot pooled cumulative incidence ----
def_labels <- c(
  "D2T_step0" = "1. DMARDs ≥2",
  "D2T_step1" = "2. + DAS28 ≥ 3.2",
  "D2T_step2" = "3. + VAS >50",
  "D2T_step3" = "4. Rolling DAS28",
  "D2T_step4" = "5. VAS >75",
  "D2T_step5" = "6. ≤10y after 2000",
  "D2T_step6" = "7. ≤5y after 2000"
)

Inc_plot <- ggplot(pooled_curve, aes(x = time_rounded / 12, y = mean_inc, color = definition)) +
  geom_line(linewidth = 1.2) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = definition), alpha = 0.2, color = NA) +
  scale_color_manual(values = scales::hue_pal()(length(defs)), labels = def_labels) +
  scale_fill_manual(values = scales::hue_pal()(length(defs)), labels = def_labels) +
  labs(
    x = "Years from Risk Start",
    y = "Cumulative Incidence",
    color = "Definition",
    fill = "Definition",
    title = "Cumulative Incidence of D2T RA Over Time per definition in the UMCU"
  ) +
  scale_x_continuous(breaks = seq(0, 20, by = 2), expand = c(0, 0)) + 
  coord_cartesian(xlim = c(0, 15)) +  # Limit at 15 years of follow-up
  scale_y_continuous(limits = c(0, 0.30), expand = c(0, 0)) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "right")

Inc_plot

################################################################
# Baseline Characteristics at different D2T Definitions 
################################################################
defs <- paste0("D2T_step", 0:6)

baseline_list <- lapply(final_list, function(data_node) {
  purrr::map_dfr(defs, function(def) {
    data_sub <- data_node %>%
      group_by(pat_ID) %>%
      summarise(
        definition    = def,
        D2T_status    = if (all(is.na(.data[[def]]))) NA else max(.data[[def]], na.rm = TRUE),
        Sex           = first(Sex),
        RF            = if (all(is.na(RF_positivity))) NA else max(RF_positivity, na.rm = TRUE),
        aCCP          = if (all(is.na(anti_CCP)))     NA else max(anti_CCP,     na.rm = TRUE),
        Age_diag      = first(Age_diagnosis),
        Year_diagnosis = first(Year_diagnosis),
        .groups = "drop"
      ) %>%
      filter(D2T_status == 1)
    
    # handle empty subsets so downstream means/SDs don't become NaN/NA
    if (nrow(data_sub) == 0) {
      return(tibble(
        Definition          = def,
        N                   = 0,
        Female_pct          = NA_real_,
        RF_pos_pct          = NA_real_,
        aCCP_pos_pct        = NA_real_,
        Age_mean            = NA_real_,
        Age_sd              = NA_real_,
        diag_lt2006_n       = 0,
        diag_lt2006_pct     = NA_real_,
        diag_2006_2010_n    = 0,
        diag_2006_2010_pct  = NA_real_,
        diag_2010_2015_n    = 0,
        diag_2010_2015_pct  = NA_real_,
        diag_gt2015_n       = 0,
        diag_gt2015_pct     = NA_real_
      ))
    }
    
    total_n <- nrow(data_sub)
    
    # Diagnosis year categories
    diag_lt2006_n      <- sum(data_sub$Year_diagnosis < 2006, na.rm = TRUE)
    diag_2006_2010_n   <- sum(data_sub$Year_diagnosis >= 2006 & data_sub$Year_diagnosis < 2010, na.rm = TRUE)
    diag_2010_2015_n   <- sum(data_sub$Year_diagnosis >= 2010 & data_sub$Year_diagnosis <= 2015, na.rm = TRUE)
    diag_gt2015_n      <- sum(data_sub$Year_diagnosis > 2015, na.rm = TRUE)
    
    
tibble(
      Definition          = def,
      N                   = total_n,
      Female_pct          = mean(data_sub$Sex == 1,   na.rm = TRUE) * 100,
      RF_pos_pct          = mean(data_sub$RF  == 1,   na.rm = TRUE) * 100,
      aCCP_pos_pct        = mean(data_sub$aCCP == 1,  na.rm = TRUE) * 100,
      Age_mean            = mean(data_sub$Age_diag,   na.rm = TRUE),
      Age_sd              = sd(  data_sub$Age_diag,   na.rm = TRUE), 
      diag_lt2006_n       = diag_lt2006_n,
      diag_lt2006_pct     = 100 * diag_lt2006_n    / total_n,
      diag_2006_2010_n    = diag_2006_2010_n,
      diag_2006_2010_pct  = 100 * diag_2006_2010_n / total_n,
      diag_2010_2015_n    = diag_2010_2015_n,
      diag_2010_2015_pct  = 100 * diag_2010_2015_n / total_n,
      diag_gt2015_n       = diag_gt2015_n,
      diag_gt2015_pct     = 100 * diag_gt2015_n    / total_n
    )
  })
})


# ---- Combine into one dataframe ----
baseline_all <- bind_rows(baseline_list, .id = "Imputation")

# ---- Pooled summary by Rubin’s Rules (simple averaging) ----
Table3 <- baseline_all %>%
  group_by(Definition) %>%
  summarise(
    N_mean               = mean(N, na.rm = TRUE),
    Female_pct_mean      = mean(Female_pct,   na.rm = TRUE),
    RF_pos_pct_mean      = mean(RF_pos_pct,   na.rm = TRUE),
    aCCP_pos_pct_mean    = mean(aCCP_pos_pct, na.rm = TRUE),
    Age_mean             = mean(Age_mean,     na.rm = TRUE),
    Age_sd_pooled        = sqrt(mean(Age_sd^2, na.rm = TRUE)),
    diag_lt2006_n_mean      = mean(diag_lt2006_n,      na.rm = TRUE),
    diag_lt2006_pct_mean    = mean(diag_lt2006_pct,    na.rm = TRUE),
    diag_2006_2010_n_mean   = mean(diag_2006_2010_n,   na.rm = TRUE),
    diag_2006_2010_pct_mean = mean(diag_2006_2010_pct, na.rm = TRUE),
    diag_2010_2015_n_mean   = mean(diag_2010_2015_n,   na.rm = TRUE),
    diag_2010_2015_pct_mean = mean(diag_2010_2015_pct, na.rm = TRUE),
    diag_gt2015_n_mean      = mean(diag_gt2015_n,      na.rm = TRUE),
    diag_gt2015_pct_mean    = mean(diag_gt2015_pct,    na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Diagnosis_lt2006      = sprintf("%.0f (%.1f%%)", diag_lt2006_n_mean,      diag_lt2006_pct_mean),
    Diagnosis_2006_2010   = sprintf("%.0f (%.1f%%)", diag_2006_2010_n_mean,   diag_2006_2010_pct_mean),
    Diagnosis_2010_2015   = sprintf("%.0f (%.1f%%)", diag_2010_2015_n_mean,   diag_2010_2015_pct_mean),
    Diagnosis_gt2015      = sprintf("%.0f (%.1f%%)", diag_gt2015_n_mean,      diag_gt2015_pct_mean)
  )

print(Table3)


################################################################
# Average Disease AFTER D2T RA 
################################################################
# Step 1: Combine the 10 imputed + processed datasets
pooled_data <- bind_rows(final_list, .id = "imputation_id")

# Definitions (D2T steps)
defs <- paste0("D2T_step", 0:6)

# Apply to pooled dataset
Table3b <- purrr::map_dfr(defs, function(def) {
  
  # Step 1: Get first D2T event per patient for this definition
  D2T_times <- pooled_data %>%
    filter(.data[[paste0(def, "_Ever")]] == 1) %>%
    group_by(pat_ID) %>%
    summarise(
      definition = def,
      D2T_time = min(Visit_months_from_diagnosis[.data[[def]] == 1], na.rm = TRUE),
      .groups = "drop"
    )
  
  # Step 2: Join and keep post-D2T visits
  Table3b <- pooled_data %>%
    semi_join(D2T_times, by = "pat_ID") %>%
    mutate(definition = def) %>%
    left_join(D2T_times, by = c("pat_ID", "definition")) %>%
    filter(Visit_months_from_diagnosis >= D2T_time)
  
  # Step 3: Summarise post-D2T disease activity
  tibble(
    Definition = def,
    DAS28_mean = mean(Table3b$DAS28, na.rm = TRUE),
    DAS28_sd = sd(Table3b$DAS28, na.rm = TRUE),
    TJC28_mean = mean(Table3b$TJC28, na.rm = TRUE),
    TJC28_sd = sd(Table3b$TJC28, na.rm = TRUE),
    SJC28_mean = mean(Table3b$SJC28, na.rm = TRUE),
    SJC28_sd = sd(Table3b$SJC28, na.rm = TRUE),
    ESR_mean = mean(Table3b$ESR, na.rm = TRUE),
    ESR_sd = sd(Table3b$ESR, na.rm = TRUE),
    CRP_mean = mean(Table3b$CRP, na.rm = TRUE),
    CRP_sd = sd(Table3b$CRP, na.rm = TRUE)
  )
})

Table3b <- Table3b %>% mutate(
  DAS28 =  paste0(round(DAS28_mean,2) , " (", round(DAS28_sd, 2), ") "),
  TJC28 =  paste0(round(TJC28_mean,2) , " (", round(TJC28_sd, 2), ") "),
  SJC28 =  paste0(round(SJC28_mean,2) , " (", round(SJC28_sd, 2), ") "),
  ESR =  paste0(round(ESR_mean,2) , " (", round(ESR_sd, 2), ") "),
  CRP =  paste0(round(CRP_mean,2) , " (", round(CRP_sd, 2), ") ")) %>%
  select(Definition, DAS28, TJC28, SJC28, ESR, CRP)

print(Table3b)



################################################################
# Persistence
################################################################

# ---- Step 0: Filter cohort (Year_diagnosis >= 2006) ----
eligible_ids <- pooled_data %>%
  group_by(pat_ID) %>%
  summarise(Year_diagnosis = first(Year_diagnosis), .groups = "drop") %>%
  filter(Year_diagnosis >= 2006) %>%
  pull(pat_ID)

pooled_data2 <- pooled_data %>%
  filter(pat_ID %in% eligible_ids)

# ---- Step 1: Calculate Persistence ----
persistence_df <- purrr::map_dfr(defs, function(def) {
  
  d2t_first <- pooled_data2 %>%
    filter(.data[[def]] == 1) %>%
    group_by(pat_ID) %>%
    summarise(
      first_D2T_time = min(Visit_months_from_diagnosis, na.rm = TRUE),
      definition = def,
      .groups = "drop"
    )
  
  post_d2t_data <- pooled_data2 %>%
    semi_join(d2t_first, by = "pat_ID") %>%
    inner_join(d2t_first, by = "pat_ID") %>%
    filter(Visit_months_from_diagnosis >= first_D2T_time) %>%
    mutate(definition = def, is_D2T = .data[[def]] == 1)
  
  persistence <- post_d2t_data %>%
    group_by(pat_ID, definition) %>%
    summarise(
      total_visits = n(),
      d2t_visits = sum(is_D2T, na.rm = TRUE),
      persistence = d2t_visits / total_visits,
      .groups = "drop"
    ) %>%
    filter(total_visits > 1)  # Optional filter
  
  return(persistence)
})

# ---- Step 2: Summary of Persistence ----
persistence_summary <- persistence_df %>%
  group_by(definition) %>%
  summarise(
    mean_persistence   = mean(persistence, na.rm = TRUE),
    sd_persistence     = sd(persistence, na.rm = TRUE),
    median_persistence = median(persistence, na.rm = TRUE),
    Q1                 = quantile(persistence, 0.25, na.rm = TRUE),
    Q3                 = quantile(persistence, 0.75, na.rm = TRUE),
    min_persistence    = min(persistence, na.rm = TRUE),
    max_persistence    = max(persistence, na.rm = TRUE),
    N_patients         = n(),
    .groups = "drop"
  )

# ---- Step 3: Follow-up Summary per Definition ----
followup_df <- purrr::map_dfr(defs, function(def) {
  d2t_first <- pooled_data2 %>%
    filter(.data[[def]] == 1) %>%
    group_by(pat_ID) %>%
    summarise(
      first_D2T_time = min(Visit_months_from_diagnosis, na.rm = TRUE),
      definition = def,
      .groups = "drop"
    )
  
  post_d2t <- pooled_data2 %>%
    semi_join(d2t_first, by = "pat_ID") %>%
    inner_join(d2t_first, by = "pat_ID") %>%
    filter(Visit_months_from_diagnosis >= first_D2T_time) %>%
    mutate(definition = def)
  
  patient_summary <- post_d2t %>%
    group_by(pat_ID, definition) %>%
    summarise(
      followup_months = max(Visit_months_from_diagnosis) - first(first_D2T_time),
      n_visits = n_distinct(Visit_months_from_diagnosis),  # 
      .groups = "drop"
    )
  
  return(patient_summary)
})

followup_summary <- followup_df %>%
  group_by(definition) %>%
  summarise(
    mean_followup_months   = mean(followup_months, na.rm = TRUE),
    sd_followup_months     = sd(followup_months, na.rm = TRUE),
    median_followup_months = median(followup_months, na.rm = TRUE),
    
    mean_n_visits          = mean(n_visits, na.rm = TRUE),
    sd_n_visits            = sd(n_visits, na.rm = TRUE),
    median_n_visits        = median(n_visits, na.rm = TRUE),
    Q1_n_visits            = quantile(n_visits, 0.25, na.rm = TRUE),
    Q3_n_visits            = quantile(n_visits, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

# ---- Final Combined Summary Table (Optional) ----
Table3c <- persistence_summary %>%
  left_join(followup_summary, by = "definition")

Table3c <- Table3c %>% mutate(
  persistence =  paste0(round(median_persistence,2) , " (", round(Q1, 2), "-",round(Q3, 2),") "),
  FU =  paste0(round(mean_followup_months,2) , " (", round(sd_followup_months, 2), ") "),
  visits =  paste0(round(mean_n_visits,2) , " (", round(sd_n_visits, 2), ") "),
  median_visits =  paste0(round(median_n_visits,2) , " (", round(Q1_n_visits, 2), "-",round(Q3_n_visits, 2), ") ")) %>%
  select(definition, persistence, FU, visits, median_visits)

# View the result
print(Table3c)

df_persist <- Table3c %>%
  select(definition, persistence, FU, visits) %>%
  pivot_longer(cols = -definition, names_to = "metric", values_to = "value")

p_persist <- ggplot(df_persist, aes(definition, as.numeric(str_extract(value, "^[0-9\\.]+")), 
                                    group = metric, color = metric)) +
  geom_line() + geom_point() +
  facet_wrap(~metric, scales = "free_y") +
  labs(x = NULL, y = NULL, title = "Persistence & follow-up") +
  theme_minimal()

p_persist 

################################################################
# Sensitvity analysis 
################################################################
defs_all <- c("D2T_step2","D2T_step3", "D2T_RA_sens1", "D2T_RA_sens2", "D2T_RA_sens3", "D2T_RA_sens4", "D2T_RA_sens5", "D2T_RA_sens6")

#centralize labels/colors and reuse across both plots & risk table
def_labels <- c(
  #D2T_step2 = "Base Definition",
  D2T_step3    = "Base Definition (Rolling Average)",
  # D2T_RA_sens1 = "No 3rd MOA for DA",   # uncomment if you include sens1
  #D2T_RA_sens2 = "1. No Imputation",
  D2T_RA_sens3 = " 2. Diagnosed >2006",
  D2T_RA_sens4 = " 3. Failed 3 MOAs for Crit1",
  D2T_RA_sens5 = "4. CDAI for DA",
  D2T_RA_sens6 = "5. No Imputation - Rolling Average"
)

def_cols <- setNames(scales::hue_pal()(length(def_labels)), names(def_labels))


get_filtered_data_node <- function(definition, data_node) {
  if (definition == "D2T_RA_sens3") {
    data_node <- data_node %>%
      filter(Year_diagnosis >= 2006)}
  return(data_node)
}

data_inc_pooled <- map2_dfr(
  final_list,
  seq_along(final_list),
  function(data_node, i) {
    map_dfr(defs_all, function(def) {
      data_node_filtered <- get_filtered_data_node(def, data_node)
      get_incidence_data(def, data_node_filtered) %>%
        mutate(imputation = i) })
  })

data_km <- data_inc_pooled %>%
  dplyr::filter(status %in% c(0, 1),
                definition %in% names(def_labels)) %>%
  dplyr::mutate(definition = factor(definition, levels = names(def_labels))) # ensures data matches label/color sets and fix factor order


# Split data by imputation
km_data_split <- group_split(data_km, imputation)
time_breaks <- seq(0, 216, by = 24)  # up to 12 years every 2 years 
time_labels <- as.character(time_breaks / 12)

# Fit and extract summaries, converting to data frames
risk_tables <- map(km_data_split, function(df) {
  fit <- survfit(Surv(time, status) ~ definition, data = df)
  # use same breaks for the risk table summaries
  summary_fit <- summary(fit, times = time_breaks)
  # Convert to data.frame and keep needed components
  tibble(
    time = summary_fit$time,
    strata = summary_fit$strata,
    n.risk = summary_fit$n.risk,
    n.event = summary_fit$n.event )})

# Combine and average
# --- risk_df ---
risk_df <- bind_rows(risk_tables, .id = "imputation") %>%
  group_by(time, strata) %>%
  summarise(
    n.risk  = round(mean(n.risk,  na.rm = TRUE)),
    n.event = round(mean(n.event, na.rm = TRUE)), .groups = "drop"
  ) %>%
  mutate(
    Definition_id = sub("definition=", "", strata),
    time_years    = factor(time / 12)) %>%
  filter(Definition_id %in% names(def_labels)) %>%
  mutate(
    RowOrder = fct_rev(factor(Definition_id, levels = names(def_labels))))


# --- risk_table_plot ---
risk_table_plot <- ggplot(risk_df, aes(x = time_years, y = RowOrder, color = Definition_id)) +
  geom_text(aes(label = paste0(n.risk, " / ", n.event)), size = 3.5, show.legend = FALSE) +
  scale_color_manual(values = def_cols, guide = "none") +
  labs(x = "Years", y = NULL, title = "At Risk / Events (every 2 years)") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y  = element_blank(),
    axis.title.y = element_blank(),
    panel.grid   = element_blank(),
    plot.title   = element_text(hjust = 0.5),
    axis.ticks   = element_blank()
  )


# Fit survival models
km_models <- map(km_data_split, ~ survfit(Surv(time, status) ~ definition, data = .x))

km_df_list <- map2(km_models, km_data_split, function(fit, data) {
  surv_summary(fit, data = data) %>%
    mutate(
      cuminc = 1 - surv,
      lower = 1 - upper,
      upper = 1 - lower )})

# Combine and average across imputations
km_pooled <- bind_rows(km_df_list, .id = "imputation") %>%
  group_by(time, strata) %>%
  summarise(
    M = n(),
    cuminc_bar = mean(cuminc, na.rm = TRUE),
    W = mean((std.err)^2, na.rm = TRUE),     # within-imputation variance
    B = var(cuminc, na.rm = TRUE),           # between-imputation variance
    T_var = W + (1 + 1/M) * B,               # total variance
    SE = sqrt(T_var),
    LL = pmax(0, cuminc_bar - 1.96 * SE),
    UL = pmin(1, cuminc_bar + 1.96 * SE),
    .groups = "drop"
  ) %>%
  mutate(
    Definition = sub("definition=", "", strata) )

fit_km <- survfit2(Surv(time, status) ~ definition, data = data_km)

# Show up to 15 years = 180 months, and label in years
time_max    <- 180                        # months (15 years)
time_breaks <- seq(0, time_max, by = 24)  # every 2 years
time_labels <- as.character(time_breaks / 12)

# Manually convert survival to cumulative incidence
fit_km$surv  <- 1 - fit_km$surv
fit_km$lower <- 1 - fit_km$upper
fit_km$upper <- 1 - fit_km$lower

# Then plot as usual (skip `transform_surv_to_cuminc()`)
km_sens <- ggsurvfit(fit_km, linewidth = 1.2) +
  scale_x_continuous(
    breaks = time_breaks, 
    labels = time_labels,
    limits = c(0, time_max),     # <-- crop at 15 years
    expand = c(0, 0)
  ) +
  scale_y_continuous(limits = c(0, 0.35), expand = c(0, 0)) +
  labs(
    title = "Cumulative Incidence of D2T RA over Time with All Sensitivity Definitions",
    x = "Years from Risk Start",
    y = "Cumulative Incidence of D2T",
    color = "Definition"
  ) +
  scale_color_manual(
    values = def_cols,
    labels = def_labels,
    breaks = names(def_labels)
  ) +
  guides(fill = "none") +
  theme_ggsurvfit_default() +
  theme(legend.position = "bottom")

sens_plot <- km_sens / risk_table_plot + plot_layout(heights = c(3, 1))
print(sens_plot)


################################################################
# Stratified Cumulative Incidence Curves 
################################################################
# Select imputed dataset to use (e.g., first)
data_node <- final_list[[1]]

# Generate the `data_inc` for D2T_step3 only
data_inc <- get_incidence_data("D2T_step3", data_node)

get_stratified_plot <- function(data_node, def = "D2T_step3", strat_var_expr, title = "Stratified Plot") {
  # Generate incidence data for the selected definition
  data_inc <- get_incidence_data(def, data_node) %>%
    filter(status %in% c(0, 1)) %>%
    mutate(strata = {{ strat_var_expr }})
  
  # Fit stratified Kaplan-Meier
  fit <- survfit2(Surv(time, status) ~ strata, data = data_inc)
  
  # Extract summary and calculate cumulative incidence
  df <- surv_summary(fit) %>%
    mutate(
      cuminc = 1 - surv,
      lower = 1 - upper,
      upper = 1 - lower,
      strata = sub("^strata=", "", strata)
    )
  
  # Plot
  ggplot(df, aes(x = time / 12, y = cuminc, color = strata)) +
    geom_step(linewidth = 1.2) +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = strata), alpha = 0.1, color = NA) +
    labs(
      title = title,
      x = "Years from Diagnosis",
      y = "Cumulative Incidence",
      color = NULL,
      fill = NULL
    ) +
    scale_x_continuous(breaks = seq(0, 20, by = 2), expand = c(0, 0)) +  # ← Added this line
    coord_cartesian(xlim = c(0, 18)) + 
    scale_y_continuous(limits = c(0, 0.3), expand = c(0, 0)) +
    theme_minimal(base_size = 13) +
    theme(legend.position = "bottom")
}

# Figure 3 - Diagnosis Year 
i_diag_year <- get_stratified_plot(
  data_node = data_node,
  def = "D2T_step3",
  strat_var_expr = case_when(
    Year_diagnosis < 2006 ~ "<2006", Year_diagnosis <= 2010 ~ "2006–2010",
    Year_diagnosis <= 2015 ~ "2011–2015", Year_diagnosis <= 2020 ~ "2016–2020",
    Year_diagnosis > 2020 ~ ">2020",
    TRUE ~ "Missing"
  ),
  title = "Cumulative Incidence of D2T RA by Year of Diagnosis"
)

# Figure 4 - Age at Diagnosis
i_age <-get_stratified_plot(
  data_node = data_node,
  def = "D2T_step3",
  strat_var_expr = case_when(
    Age_diagnosis < 40 ~ "<40", Age_diagnosis < 50 ~ "40–50",
    Age_diagnosis < 60 ~ "50–60", Age_diagnosis < 70 ~ "60–70",
    Age_diagnosis < 80 ~ "70–80", Age_diagnosis >= 80 ~ ">80",
    TRUE ~ "Missing"
  ),
  title = "Cumulative Incidence of D2T RA by Age at Diagnosis"
)

# Figure 5 - RF-ACCP positivity 
i_markers <-get_stratified_plot(
  data_node = data_node,
  def = "D2T_step3",
  strat_var_expr = case_when(
    RF_positivity == 1 & anti_CCP == 1 ~ "Both+",
    RF_positivity == 1 & is.na(anti_CCP) ~ "All RF+",
    anti_CCP == 1 & is.na(RF_positivity) ~ "All CCP+",
    (RF_positivity == 1 & anti_CCP == 0) |
      (RF_positivity == 0 & anti_CCP == 1) ~ "Either+",
    RF_positivity == 0 & anti_CCP == 0 ~ "Seronegative",
    TRUE ~ "Missing info"
  ),
  title = "Cumulative Incidence of D2T RA by Serological Status"
)

print(i_diag_year)
print(i_age)
print(i_markers)


###############################################################################
#  Prevalence  overall and per gender/age category
###############################################################################
# Base prevalence over time for all 8 definitions
ydefs <- paste0("D2T_step", 0:6)
def_labels <- c(
  "D2T_step0" = "1. DMARDs ≥2", "D2T_step1" = "2. + DAS28 ≥ 3.2",
  "D2T_step2" = "3. + VAS >50", "D2T_step3" = "4. Rolling DAS28",
  "D2T_step4" = "5. VAS >75", "D2T_step5" = "6. ≤10y after 2000",
  "D2T_step6" = "7. ≤5y after 2000")

# Collect prevalence data across imputations
prev_all <- lapply(seq_along(final_list), function(i) { data_node <- final_list[[i]]
data_node %>%
  filter(current_year > 2005) %>%
  select(pat_ID, current_year, all_of(defs)) %>%
  pivot_longer(cols = defs, names_to = "Definition", values_to = "D2T") %>%
  group_by(pat_ID, current_year, Definition) %>%
  summarise(D2T = ifelse(all(is.na(D2T)), NA, max(D2T, na.rm = TRUE)), .groups = "drop") %>%
  group_by(current_year, Definition) %>%
  summarise(
    N = n(),
    p = mean(D2T, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(imputation = i)
}) %>%
  bind_rows()

prevalence_pooled <- prev_all %>%
  group_by(current_year, Definition) %>%
  summarise(
    M = n(),  # should be 10 imputations
    p_bar = mean(p),  u_bar = mean(p * (1 - p) / N),
    b = var(p), T_var = u_bar + (1 + 1/M) * b,
    SE = sqrt(T_var), LL = pmax(0, p_bar - 1.96 * SE),   UL = pmin(1, p_bar + 1.96 * SE),
    .groups = "drop"
  ) %>%
  mutate(
    Prevalence = round(p_bar * 100, 2), LL = round(LL * 100, 2),  UL = round(UL * 100, 2),
    Label = paste0(Prevalence, " (", LL, ", ", UL, ")"),
    Definition = factor(Definition, levels = defs, labels = def_labels)
  )

ggplot(prevalence_pooled, aes(x = current_year, y = Prevalence, color = Definition)) +
  geom_line(linewidth = 1.2) +
  labs(
    title = "Stepwise D2T RA Prevalence Over Calendar Years in the UMCU",
    x = "Calendar Year",
    y = "Prevalence per 100 Patients",
    color = "D2T Definition"
  ) +
  scale_x_continuous(
    limits = c(min(prevalence_pooled$current_year), 2024),
    breaks = seq(min(prevalence_pooled$current_year), 2024, by = 2)
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "right")


Table4 <- prevalence_pooled %>%
  select(current_year, Definition, Label) %>%
  pivot_wider(names_from = Definition, values_from = Label) %>%
  arrange(current_year)


##### Stratified Prevalence 
# Use first imputed dataset
data_node <- final_list[[1]]

# Summarise by patient and year
data_node_prev <- data_node %>%
  filter(current_year > 2005) %>%
  group_by(pat_ID, current_year) %>%
  summarise(
    D2T_step3 = max(D2T_step3, na.rm = TRUE),
    Sex = max(Sex),
    RF = first(RF_positivity, na.rm = TRUE),
    anti_CCP = first(anti_CCP, na.rm = TRUE),
    current_age_cat = first(current_age_cat),
    .groups = "drop"
  )

# Recode age categories (assuming they are 1 to 6)
age_cuts <- c(0, 40, 50, 60, 70, 80, Inf)
age_labels <- c("1" = "< 40", "2" = "40–50","3" = "50–60","4" = "60–70", "5" = "70–80", "6" = "> 80")

# Collect stratified prevalence across imputations
prev_age_list <- lapply(seq_along(final_list), function(i) {
  data_node <- final_list[[i]]
  data_node %>%
    filter(current_year > 2005) %>%
    select(pat_ID, current_year, current_age, D2T_step2) %>%
    mutate(age_group = cut(current_age, breaks = age_cuts, labels = age_labels, right = FALSE)) %>%
    group_by(pat_ID, current_year, age_group) %>%
    summarise(D2T = ifelse(all(is.na(D2T_step2)), NA, max(D2T_step2, na.rm = TRUE)), .groups = "drop") %>%
    group_by(current_year, age_group) %>%
    summarise(  N = n(),  p = mean(D2T, na.rm = TRUE),  .groups = "drop"
    ) %>%
    mutate(imputation = i)
}) %>%
  bind_rows()

# Pool prevalence across imputations
prev_age_pooled <- prev_age_list %>%
  group_by(current_year, age_group) %>%
  summarise( M = n(),  p_bar = mean(p),  u_bar = mean(p * (1 - p) / N),
             b = var(p),   T_var = u_bar + (1 + 1/M) * b,   SE = sqrt(T_var),  LL = pmax(0, p_bar - 1.96 * SE),   
             UL = pmin(1, p_bar + 1.96 * SE),  .groups = "drop"
  ) %>%
  mutate(  Prevalence = round(p_bar * 100, 2),   LL = round(LL * 100, 2),    UL = round(UL * 100, 2),
           Label = paste0(Prevalence, " (", LL, ", ", UL, ")")
  )
p_age <- ggplot(prev_age_pooled, aes(x = current_year, y = Prevalence, color = age_group)) +
  geom_line(linewidth = 1.2) +
  labs(
    title = "Prevalence of Base Definiton of D2T by Age Group Over Time",
    x = "Calendar Year",
    y = "Prevalence per 100 Patients",
    color = "Age Group"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "right")

p_age

# Collect stratified prevalence across imputations
prev_gender_list <- lapply(seq_along(final_list), function(i) {
  data_node <- final_list[[i]]
  
  data_node %>%
    filter(current_year > 2005) %>%
    select(pat_ID, current_year, Sex, D2T_step2) %>%
    group_by(pat_ID, current_year, Sex) %>%
    summarise(D2T = ifelse(all(is.na(D2T_step2)), NA, max(D2T_step2, na.rm = TRUE)), .groups = "drop") %>%
    group_by(current_year, Sex) %>%
    summarise(  N = n(),   p = mean(D2T, na.rm = TRUE),   .groups = "drop"
    ) %>%
    mutate(imputation = i)
}) %>%
  bind_rows()

# Pool prevalence across imputations
prev_gender_pooled <- prev_gender_list %>%
  group_by(current_year, Sex) %>%
  summarise(
    M = n(), p_bar = mean(p),  u_bar = mean(p * (1 - p) / N),  b = var(p),  T_var = u_bar + (1 + 1/M) * b,
    SE = sqrt(T_var), LL = pmax(0, p_bar - 1.96 * SE),    UL = pmin(1, p_bar + 1.96 * SE),
    .groups = "drop"
  ) %>%
  mutate(
    Prevalence = round(p_bar * 100, 2),  LL = round(LL * 100, 2),   UL = round(UL * 100, 2),
    Label = paste0(Prevalence, " (", LL, ", ", UL, ")")
  )

p_sex <- ggplot(prev_gender_pooled, aes(x = current_year, y = Prevalence, color = factor(Sex))) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(
    values = c("0" = "#1f77b4", "1" = "#ff7f0e"),  # choose your preferred colors
    labels = c("0" = "Male", "1" = "Female"),
    name = "Gender"
  ) +
  labs(
    title = "Prevalence of the base definition of D2T RA by Gender Over Time",
    x = "Calendar Year",
    y = "Prevalence per 100 Patients",
    color = "Gender"
  ) +
  scale_x_continuous(
    limits = c(min(prevalence_pooled$current_year), 2024),
    breaks = seq(min(prevalence_pooled$current_year), 2024, by = 2)
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "right")

p_sex 


prev_sero_list <- lapply(seq_along(final_list), function(i) {
  data_node <- final_list[[i]]
  
  data_node %>%
    filter(current_year > 2005) %>%
    select(pat_ID, current_year, D2T_step2, RF_positivity, anti_CCP) %>%
    mutate(
      RF_aCCP_group = case_when(
        RF_positivity == 1 & anti_CCP == 1 ~ "Both+",
        RF_positivity == 1 ~ "All RF+",
        anti_CCP == 1 ~ "All CCP+",
        (RF_positivity == 1 & is.na(anti_CCP)) |
          (anti_CCP == 1 & is.na(RF_positivity)) |
          (RF_positivity == 1 & anti_CCP == 0) |
          (RF_positivity == 0 & anti_CCP == 1) ~ "Either+",
        RF_positivity == 0 & anti_CCP == 0 ~ "Seronegative",
        is.na(RF_positivity) & is.na(anti_CCP) ~ "Missing info",
        TRUE ~ "Missing info"
      )
    ) %>%
    group_by(pat_ID, current_year, RF_aCCP_group) %>%
    summarise(D2T = ifelse(all(is.na(D2T_step2)), NA, max(D2T_step2, na.rm = TRUE)), .groups = "drop") %>%
    group_by(current_year, RF_aCCP_group) %>%
    summarise(
      N = n(),
      p = mean(D2T, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(imputation = i)
}) %>%
  bind_rows()

prev_sero_pooled <- prev_sero_list %>%
  group_by(current_year, RF_aCCP_group) %>%
  summarise(
    M = n(),
    p_bar = mean(p), 
    u_bar = mean(p * (1 - p) / N),
    b = var(p), 
    T_var = u_bar + (1 + 1/M) * b,
    SE = sqrt(T_var),
    LL = pmax(0, p_bar - 1.96 * SE),   
    UL = pmin(1, p_bar + 1.96 * SE),
    .groups = "drop"
  ) %>%
  mutate(
    Prevalence = round(p_bar * 100, 2), 
    LL = round(LL * 100, 2),  
    UL = round(UL * 100, 2),
    Label = paste0(Prevalence, " (", LL, ", ", UL, ")")
  )
p_markers <- ggplot(prev_sero_pooled, aes(x = current_year, y = Prevalence, color = RF_aCCP_group)) +
  geom_line(linewidth = 1.2) +
  labs(
    title = "Prevalence of the base definition of D2T RA by RF/aCCP Positivity",
    x = "Calendar Year",
    y = "Prevalence per 100 Patients",
    color = "Autoantibody Status"
  ) +
  scale_x_continuous(
    limits = c(min(prevalence_pooled$current_year), 2024),
    breaks = seq(min(prevalence_pooled$current_year), 2024, by = 2)
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom")

p_markers 




################################################################
# Cox regression
################################################################
# Step 1 : Build survival data across imputations 
imputed_surv_data <- lapply(final_list, function(data_node) {
  data_node %>%
    filter(!is.na(Visit_months_from_diagnosis)) %>%
    group_by(pat_ID) %>%
    summarise(
      event = as.integer(any(D2T_step2 == 1, na.rm = TRUE)),
      
      time_months = if (any(D2T_step2 == 1, na.rm = TRUE)) {
        # first event time
        min(Visit_months_from_diagnosis[D2T_step2 == 1], na.rm = TRUE)
      } else {
        # censor at last follow-up
        max(Visit_months_from_diagnosis, na.rm = TRUE)
      },
      # covariates
      age_diag = first(Age_diagnosis),
      gender   = as.numeric(first(Sex)),
      RF       = first(RF_positivity),
      CCP      = first(anti_CCP),
      Year_diagnosis = first(Year_diagnosis),
      .groups = "drop"
    ) %>%
    mutate(time = time_months /12,
           RF_CCP_group = case_when(
             RF == 1 & CCP == 1 ~ "Both Positive",
             (RF == 1 & CCP != 1) | (CCP == 1 & RF != 1) ~ "Either Positive",
             RF == 0 & CCP == 0 ~ "Seronegative",
             TRUE ~ "Missing"
           ),
           RF_CCP_group = factor(RF_CCP_group, levels = c("Seronegative", "Either Positive", "Both Positive", "Missing")),
           year_cat = case_when(
             Year_diagnosis < 2006 ~ "<2006",
             Year_diagnosis <= 2010 ~ "2006–2010",
             Year_diagnosis <= 2015 ~ "2011–2015",
             Year_diagnosis > 2015 ~ "2016–2024",
             TRUE ~ NA_character_
           ),
           year_cat = factor(year_cat, levels = c("<2006", "2006–2010", "2011–2015", "2016–2024"))
    )
})
# step 2: Fit Cox models across imputations

models <- lapply(imputed_surv_data, function(data) {
  coxph(Surv(time, event) ~ age_diag + gender + RF_CCP_group + year_cat, data = data)
})
# step 3: Pool results using mice:: pool ()
# Create mira object
mira_obj <- as.mira(models)

# Pool results
pooled_results <- pool(mira_obj)

# Tidy and format
multi_df <- summary(pooled_results, conf.int = TRUE) %>%
  mutate(
    Variable = recode(term,
                      "RF_CCP_groupEither Positive" = "Either vs Seronegative",
                      "RF_CCP_groupBoth Positive" = "Both vs Seronegative",
                      "RF_CCP_groupMissing" = "Missing vs Seronegative",
                      "gender" = "Female gender",
                      "age_diag" = "Age at Diagnosis",
                      "year_cat2006–2010" = "2006–2010 vs <2006",
                      "year_cat2011–2015" = "2011–2015 vs <2006",
                      "year_cat2016–2024" = "2016–2024 vs <2006"
    ),
    Multivariate = paste0(round(estimate, 2), " (", round(`2.5 %`, 2), " – ", round(`97.5 %`, 2), ")")
  ) %>%
  select(Variable, Multivariate)

get_univ <- function(var) {
  models_univ <- lapply(imputed_surv_data, function(data) {
    f <- as.formula(paste("Surv(time, event) ~", var))
    coxph(f, data = data)
  })
  
  pooled <- pool(as.mira(models_univ))
  pooled_summary <- summary(pooled, conf.int = TRUE)
  
  # Identify confidence interval columns
  conf_low_col <- if ("2.5 %" %in% names(pooled_summary)) "2.5 %" else "conf.low"
  conf_high_col <- if ("97.5 %" %in% names(pooled_summary)) "97.5 %" else "conf.high"
  
  pooled_summary %>%
    mutate(
      conf.low = .[[conf_low_col]],
      conf.high = .[[conf_high_col]],
      Variable = recode(term,
                        "gender" = "Female gender",
                        "age_diag" = "Age at Diagnosis",
                        "RF_CCP_groupEither Positive" = "Either vs Seronegative",
                        "RF_CCP_groupBoth Positive" = "Both vs Seronegative",
                        "RF_CCP_groupMissing" = "Missing vs Seronegative",
                        "year_cat2006–2010" = "2006–2010 vs <2006",
                        "year_cat2011–2015" = "2011–2015 vs <2006",
                        "year_cat2016–2024" = "2016–2024 vs <2006"
      ),
      Univariate = paste0(
        round(estimate, 2), " (",
        round(conf.low, 2), " – ",
        round(conf.high, 2), ")"
      )
    ) %>%
    select(Variable, Univariate)
}



univ_vars <- c("age_diag", "gender", "RF_CCP_group", "year_cat")
univ_all <- bind_rows(lapply(univ_vars, get_univ))

Table5 <- full_join(univ_all, multi_df, by = "Variable") %>%
  select(Variable, Univariate, Multivariate)
View(Table5)
