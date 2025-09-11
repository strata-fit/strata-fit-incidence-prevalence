###############################################################################
# Incidence and Prevalence of D2T RA 
# VERSION 3.0
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
data_node<- read_csv("Documents/data_node_1507.csv") # Upload your own data 

data_node <- data_node %>%
  filter(!is.na(Visit_months_from_diagnosis))

# Count JAKi as 1 class in tsDMARD
data_node <- data_node %>%
  mutate(tsDMARD = case_when(
    is.na(tsDMARD) ~ NA_real_,
    tsDMARD >= 1   ~ 1,
    TRUE           ~ 0
  ))

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
    VAS_physician_mean = mean(Ph_global, na.rm = TRUE)
  ) %>%
  mutate(
    visits_per_year = ifelse(is.na(FU_years) | FU_years == 0, NA, n_visits / FU_years)
  )

# Build summary table
Table1 <- list("Female, n (%)" = {
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
  "Initial calendar year of follow-up" = {
    m <- mean(patient_summary$Year_diagnosis, na.rm = TRUE)
    s <- sd(patient_summary$Year_diagnosis, na.rm = TRUE)
    paste0(round(m, 1), " (", round(s, 1), ")")},
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
    paste0(round(m, 1), " (", round(s, 1), ")")  })

# Convert to data frame
summary_table1 <- data.frame(Variable = names(Table1), Value = unlist(Table1), row.names = NULL)
print(summary_table1)

data_node <- data_node[order(data_node$pat_ID, data_node$Visit_months_from_diagnosis),]
rm(patient_summary, summary_table1)
n_distinct(data_node$pat_ID)

###############################################################################
#  Create derived variables needed for D2T RA definition and analyses 
###############################################################################

#  cumulative DMARD treatment variables 
data_node <- data_node %>%
  group_by(pat_ID) %>%
  mutate(
    cum_csDMARD1 = cumsum(!duplicated(csDMARD1) & !is.na(csDMARD1)),
    cum_csDMARD2 = cumsum(!duplicated(csDMARD2) & !is.na(csDMARD2)),
    cum_csDMARD3 = cumsum(!duplicated(csDMARD3) & !is.na(csDMARD3)),
    cum_bDMARD = cumsum(!duplicated(bDMARD) & !is.na(bDMARD)),
    cum_tsDMARD = cumsum(!duplicated(tsDMARD) & !is.na(tsDMARD)),
    cum_btsDMARD = cum_bDMARD + cum_tsDMARD + N_prev_bDMARD + N_prev_tsDMARD,
    cum_GC = cumsum(!duplicated(GC) & !is.na(GC)),
    cum_csDMARD = max(cum_csDMARD1, cum_csDMARD2, cum_csDMARD3, N_prev_csDMARD, na.rm= TRUE),
    cum_btsDMARDmin = cummin(cum_btsDMARD)) %>% # Define cum_btsDMARDmin (cum_btsDMARDs at start follow-up) for left/period censoring 
  ungroup()

#define current calender year variable
data_node$current_year <- round(data_node$Year_diagnosis + (data_node$month_diagnosis/12) +  (data_node$Visit_months_from_diagnosis)/12)

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
      rol_av_Pat_global = rollmean(Pat_global_imp, k = 2, fill = NA, align = "right"),
      rol_av_Ph_global = rollmean(Ph_global_imp, k = 2, fill = NA, align = "right")
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
# Apply D2T logic to all 10 datasets
final_list <- lapply(processed_list, function(data_node) {
  data_node <- data_node %>%
    mutate(
      # Criterion 1
      D2T_crit1 = ifelse(cum_csDMARD > 0 & cum_btsDMARD > 1, 1, 0),
      D2T_crit1a = ifelse(cum_csDMARD > 0 & cum_btsDMARD > 2, 1, 0),
      D2T_crit1b = ifelse(cum_csDMARD > 0 & cum_btsDMARD > 2 & FU_2000 < 120, 1, 0),
      
      # Criterion 2
      D2T_crit2 = ifelse(DAS28_imp > 3.2 | (changed_MOA == 1 & cum_btsDMARD > 1), 1, 0),
      D2T_crit2a = ifelse(rol_av_DAS28 > 3.2 | (changed_MOA == 1 & cum_btsDMARD > 1), 1, 0),
      D2T_crit2b = ifelse(DAS28 > 3.2 | (changed_MOA == 1 & cum_btsDMARD > 1), 1, 0),
      D2T_crit2sens1 = ifelse(DAS28_imp > 3.2, 1, 0),
      D2T_crit2sens2 = ifelse(DAS28 > 3.2 | (changed_MOA == 1 & cum_btsDMARD > 1), 1, 0),
      
      # Criterion 3
      D2T_crit3 = ifelse(Pat_global_imp > 50 | Ph_global_imp > 50, 1, 0),
      D2T_crit3a = ifelse(Pat_global_imp > 75 | Ph_global_imp > 75, 1, 0),
      D2T_crit3b = ifelse(Pat_global > 50 | Ph_global > 50, 1, 0),
      
      # D2T steps
      D2T_step0 = ifelse(D2T_crit1 == 1, 1, 0),
      D2T_step1 = ifelse(D2T_crit1 == 1 & D2T_crit2 == 1, 1, 0),
      D2T_step2 = ifelse(D2T_crit1 == 1 & D2T_crit2 == 1 & D2T_crit3 == 1, 1, 0),
      D2T_step3 = ifelse(D2T_crit1 == 1 & D2T_crit2a & D2T_crit3 == 1, 1, 0),
      D2T_step4 = ifelse(D2T_crit1 == 1 & D2T_crit2a & D2T_crit3a == 1, 1, 0),
      D2T_step5 = ifelse(D2T_crit1 == 1 & D2T_crit2a & D2T_crit3a == 1 & FU_2000 <= 120, 1, 0),
      D2T_step6 = ifelse(D2T_crit1 == 1 & D2T_crit2a & D2T_crit3a == 1 & FU_2000 <= 60, 1, 0),
      D2T_step7 = ifelse(D2T_crit1a == 1 & D2T_crit2a & D2T_crit3a == 1 & FU_2000 <= 60, 1, 0),   # Should we use 3rd MOA prior to step 5 ? 
      
      # Sensitivity analyses
      D2T_RA_sens1 = ifelse(D2T_crit1 == 1 & D2T_crit2sens1 == 1 & D2T_crit3 == 1, 1, 0),
      D2T_RA_sens2 = ifelse(D2T_crit1 == 1 & D2T_crit2sens2 == 1 & D2T_crit3b == 1, 1, 0),
      D2T_RA_sens3 = ifelse(D2T_crit1 == 1 & D2T_crit2a  & D2T_crit3 == 1, 1, 0)
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
calculate_incidence_per_dataset <- function(data_node, defs = paste0("D2T_step", 0:7), time_point = 60) {
  
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
        cum_btsDMARDmin = max(cum_btsDMARDmin, na.rm = TRUE),
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
Table2_pooled_results  <- incidence_all %>%
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
      cum_btsDMARDmin = max(cum_btsDMARDmin, na.rm = TRUE),
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

# ---- Step 2: Run survival for each dataset and definition ----
defs <- paste0("D2T_step", 0:7)

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

# Replace final_list with your actual list of imputed datasets (with D2T applied)


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
  "D2T_step6" = "7. ≤5y after 2000",
  "D2T_step7" = "8. MOA ≥3"
)

Inc_plot <- ggplot(pooled_curve, aes(x = time_rounded / 12, y = mean_inc, color = definition)) +
  geom_line(size = 1.2) +
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
  scale_y_continuous(limits = c(0, 0.40), expand = c(0, 0)) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "right")

Inc_plot

################################################################
# Baseline Characteristics at different D2T Definitions 
################################################################
# ---- Loop across imputed datasets and definitions ----
baseline_list <- lapply(final_list, function(data_node) {
  purrr::map_dfr(defs, function(def) {
    data_sub <- data_node %>%
      group_by(pat_ID) %>%
      summarise(
        definition = def,
        D2T_status = if (all(is.na(.data[[def]]))) NA else max(.data[[def]], na.rm = TRUE),
        Sex = first(Sex),
        RF = if (all(is.na(RF_positivity))) NA else max(RF_positivity, na.rm = TRUE),
        aCCP = if (all(is.na(anti_CCP))) NA else max(anti_CCP, na.rm = TRUE),
        Age_diag = first(Age_diagnosis),
        .groups = "drop"
      ) %>%
      filter(D2T_status == 1)
    
    tibble(
      Definition = def,
      N = nrow(data_sub),
      Female_pct = mean(data_sub$Sex == 1, na.rm = TRUE) * 100,
      RF_pos_pct = mean(data_sub$RF == 1, na.rm = TRUE) * 100,
      aCCP_pos_pct = mean(data_sub$aCCP == 1, na.rm = TRUE) * 100,
      Age_mean = mean(data_sub$Age_diag, na.rm = TRUE),
      Age_sd = sd(data_sub$Age_diag, na.rm = TRUE)
    )
  })
})

# ---- Combine into one dataframe ----
baseline_all <- bind_rows(baseline_list, .id = "Imputation")

# ---- Pooled summary by Rubin’s Rules (simple averaging) ----
Table3 <- baseline_all %>%
  group_by(Definition) %>%
  summarise(
    N_mean = mean(N),
    Female_pct_mean = mean(Female_pct),
    RF_pos_pct_mean = mean(RF_pos_pct),
    aCCP_pos_pct_mean = mean(aCCP_pos_pct),
    Age_mean = mean(Age_mean),
    Age_sd_pooled = sqrt(mean(Age_sd^2)),  # Pooled SD across datasets
    .groups = "drop"
  )

# ---- Print nicely ----
print(Table3)


################################################################
# Average Disease AFTER D2T RA 
################################################################
# Step 1: Combine the 10 imputed + processed datasets
pooled_data <- bind_rows(final_list, .id = "imputation_id")

# Definitions (D2T steps)
defs <- paste0("D2T_step", 0:7)

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

################################################################
# Persistance
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
Table4 <- persistence_summary %>%
  left_join(followup_summary, by = "definition")

# View the result
print(Table4)

################################################################
# Sensitvity analysis 
################################################################
defs_all <- c("D2T_step3", "D2T_RA_sens1", "D2T_RA_sens2", "D2T_RA_sens3")

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
  filter(status %in% c(0, 1)) %>%
  mutate(definition = factor(definition, levels = defs_all))


# Split data by imputation
km_data_split <- group_split(data_km, imputation)

# Fit and extract summaries, converting to data frames
risk_tables <- map(km_data_split, function(df) {
  fit <- survfit(Surv(time, status) ~ definition, data = df)
  summary_fit <- summary(fit, times = time_breaks)
  
  # Convert to data.frame and keep needed components
  tibble(
    time = summary_fit$time,
    strata = summary_fit$strata,
    n.risk = summary_fit$n.risk,
    n.event = summary_fit$n.event )})

    # Combine and average
    risk_df <- bind_rows(risk_tables, .id = "imputation") %>%
      group_by(time, strata) %>%
      summarise(
        n.risk = round(mean(n.risk, na.rm = TRUE)),
        n.event = round(mean(n.event, na.rm = TRUE)),
        .groups = "drop"
      ) %>%
      mutate(
        Definition = sub("definition=", "", strata),
        time_years = time / 12)
    
    risk_table_plot <- risk_df %>%
      mutate(time_years = as.factor(time_years)) %>%
      ggplot(aes(x = time_years, y = Definition)) +
      geom_text(aes(label = paste0(n.risk, " / ", n.event)), size = 3.5) +
      labs(x = "Years", y = NULL, title = "At Risk / Events") +
      theme_minimal(base_size = 12) +
      theme(
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.text.y = element_text(hjust = 1),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_text(margin = margin(t = 5)),
        axis.ticks = element_blank())

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

# Manually convert survival to cumulative incidence
fit_km$surv <- 1 - fit_km$surv
fit_km$lower <- 1 - fit_km$upper
fit_km$upper <- 1 - fit_km$lower

# Then plot as usual (skip `transform_surv_to_cuminc()`)
km_sens <- ggsurvfit(fit_km, linewidth = 1.2) +
  scale_x_continuous(breaks = time_breaks, labels = time_labels) +
  coord_cartesian(xlim = c(0, 15)) +  # Limit at 15 years of follow-up
  scale_y_continuous(limits = c(0, 0.35),expand = c(0, 0)) +
  labs(
    title = "Cumulative Incidence of D2T RA over Time with All Sensitivity Definitions",
    x = "Years from Risk Start",
    y = "Cumulative Incidence of D2T",
    color = "Definition") +
  scale_color_manual(
    values = c(
      "D2T_step3" = "#984ea3",
      "D2T_RA_sens1" = "#e41a1c",
      "D2T_RA_sens2" = "#4daf4a",
      "D2T_RA_sens3" = "#377eb8"),
    labels = c(
      "D2T_step3" = "Base Definition (Rolling Average)",
      "D2T_RA_sens1" = "Sens. 1: No 3rd MOA",
      "D2T_RA_sens2" = "Sens. 2: No Imputation",
      "D2T_RA_sens3" = "Sens. 3: Diagnosed ≥ 2006" )
  ) +
  guides(fill = "none") +
  theme_ggsurvfit_default() +
  theme(legend.position = "bottom")

Sens_plot <- km_sens / risk_table_plot + plot_layout(heights = c(3, 1))
print(Sens_plot)

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
    scale_x_continuous(breaks = seq(0, 20, by = 2), expand = c(0, 0)) +  
    coord_cartesian(xlim = c(0, 15)) +  # Limit at 15 years of follow-up
    scale_y_continuous(limits = c(0, 0.3), expand = c(0, 0)) +
    theme_minimal(base_size = 13) +
    theme(legend.position = "bottom")
}

# Figure 3 - Diagnosis Year 
Inc_diag_year <- get_stratified_plot(
  data_node = data_node,
  def = "D2T_step3",
  strat_var_expr = case_when(
    Year_diagnosis < 2006 ~ "<2006", Year_diagnosis <= 2010 ~ "2006–2010",
    Year_diagnosis <= 2015 ~ "2011–2015", Year_diagnosis <= 2020 ~ "2016–2020",
    Year_diagnosis > 2020 ~ ">2020",
    TRUE ~ "Missing"
  ),
  levels = c("<2006", "2006–2010", "2011–2015", "2016–2020", ">2020", "Missing"),
  title = "Cumulative Incidence by Year of Diagnosis"
)
Inc_diag_year

# Figure 4 - Ag at Diagnosis
Inc_diag_age <- get_stratified_plot(
  data_node = data_node,
  def = "D2T_step3",
  strat_var_expr = case_when(
    Age_diagnosis < 40 ~ "<40", Age_diagnosis < 50 ~ "40–50",
    Age_diagnosis < 60 ~ "50–60", Age_diagnosis < 70 ~ "60–70",
    Age_diagnosis < 80 ~ "70–80", Age_diagnosis >= 80 ~ ">80",
    TRUE ~ "Missing"
  ),
  levels = c("<40", "40–50", "50–60", "60–70", "70–80", ">80", "Missing"),
  title = "Cumulative Incidence by Age at Diagnosis"
)

Inc_diag_age

# Figure 5 - RF-ACCP positivity 
Inc_sero <- get_stratified_plot(
  data_node = data_node,
  def = "D2T_step3",
  strat_var_expr = case_when(
    RF_positivity == 1 & anti_CCP == 1 ~ "Both+",
    RF_positivity == 1  ~ "All RF+",
    anti_CCP == 1  ~ "All CCP+",
    (RF_positivity == 1 & anti_CCP == 0) |
      (RF_positivity == 0 & anti_CCP == 1) ~ "Either+",
    RF_positivity == 0 & anti_CCP == 0 ~ "Seronegative",
    TRUE ~ "Missing info"
  ),
  levels = c("Both+", "All RF+", "All CCP+", "Either+", "Seronegative", "Missing info"),
  title = "Cumulative Incidence by Serological Status"
)

Inc_sero
###############################################################################
#  Prevalence  overall and per gender/age category
###############################################################################
# Base prevalence over time for all 8 definitions
defs <- paste0("D2T_step", 0:7)
def_labels <- c(
  "D2T_step0" = "1. DMARDs ≥2", "D2T_step1" = "2. + DAS28 ≥ 3.2",
  "D2T_step2" = "3. + VAS >50", "D2T_step3" = "4. Rolling DAS28",
  "D2T_step4" = "5. VAS >75", "D2T_step5" = "6. ≤10y after 2000",
  "D2T_step6" = "7. ≤5y after 2000","D2T_step7" = "8. MOA ≥3")

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


Prev_plot <- ggplot(prevalence_pooled, aes(x = current_year, y = Prevalence, color = Definition)) +
  geom_line(linewidth = 1.2) +
  labs(
    title = "Stepwise D2T RA Prevalence Over Calendar Years",
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

Prev_plot

# Create Table 4 with prevalence labels
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
    RF = max(RF_positivity, na.rm = TRUE),
    anti_CCP = max(anti_CCP, na.rm = TRUE),
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
      time = max(Visit_months_from_diagnosis, na.rm = TRUE) / 12,
      event = max(D2T_step2, na.rm = TRUE),
      age_diag = first(Age_diagnosis),
      gender = as.numeric(first(Sex)),
      RF = first(RF_positivity),
      CCP = first(anti_CCP),
      Year_diagnosis = first(Year_diagnosis),
      .groups = "drop"
    ) %>%
    filter(!is.na(time) & !is.na(event)) %>%
    mutate(
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

