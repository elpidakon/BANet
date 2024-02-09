# Kontsioti, Maskell, Anderson & Pirmohamed, Identifying drug-drug interactions in spontaneous reports utilizing signal detection and biological plausibility aspects (2024)
# This script tests multiple sets of hyperparameters applied to the signal detection algorithm and evaluates their relative performance using ROC analysis

library(extraDistr)
library(matrixStats)
library(dplyr)
library(ROCR)
library(TailRank)
library(VGAM)
library(readxl)
library(DBI)
library(collections)
library(tidyverse)
library(qdapTools)
library(insurancerating)
library(ExtDist)
library(PerformanceAnalytics)
library(ggoutlier)
library(ggplot2)
library(reshape2)

# Set working directory (replace with path)
workdir <- "path/to/working/directory"
setwd(workdir)
# Set seed value
set.seed(13579)

### PART A: Data pre-processing (FAERS counts)

## Positive controls

# Read spreadsheet with FAERS counts for positive controls
pos_counts <- readxl::read_excel("data/evaluation/DR1_FAERS_COUNTS.xlsx")

# Keep instances with non-zero FAERS counts for the drug-drug-event triplet
pos_counts <- pos_counts[pos_counts$n_111 != 0,]

# Calculate the remaining counts
pos_counts$n_110 <- pos_counts$d1_d2_counter - pos_counts$n_111
pos_counts$n_100 <- pos_counts$d1_not_d2_counter - pos_counts$n_101
pos_counts$n_010 <- pos_counts$not_d1_d2_counter - pos_counts$n_011
pos_counts$n_000 <- pos_counts$not_d1_not_d2_counter - pos_counts$n_001

# Keep selected columns and convert to uppercase
pos_counts <- pos_counts %>%
  select(dde_tuple, n_111, n_101, n_011, n_001,
         n_110, n_100, n_010, n_000) %>%
  rename_with(toupper)

# Add class label
pos_counts$CLASS <- 1

## Negative controls

# Read spreadsheet with FAERS counts for negative controls
neg_counts <- readxl::read_excel("data/evaluation/DR2_FAERS_COUNTS.xlsx")

# Keep instances with non-zero FAERS counts for the drug-drug-event triplet
neg_counts <- neg_counts[neg_counts$n_111 != 0,]

# Calculate the remaining counts
neg_counts$n_110 <- neg_counts$d1_d2_counter - neg_counts$n_111
neg_counts$n_100 <- neg_counts$d1_not_d2_counter - neg_counts$n_101
neg_counts$n_010 <- neg_counts$not_d1_d2_counter - neg_counts$n_011
neg_counts$n_000 <- neg_counts$not_d1_not_d2_counter - neg_counts$n_001

# Keep selected columns and convert to uppercase
neg_counts <- neg_counts %>%
  select(dde_tuple, n_111, n_101, n_011, n_001,
         n_110, n_100, n_010, n_000) %>%
  rename_with(toupper)

# Add class label
neg_counts$CLASS <- 0

# Concatenate the FAERS counts of positive and negative controls to a single df
controls <- rbind(pos_counts, neg_counts)

controls$DDE_TUPLE_COPY <- controls$DDE_TUPLE
# Remove the first and last character and split concept_ids to separate columns
controls$DDE_TUPLE_COPY <- substr(
  controls$DDE_TUPLE_COPY, 2,
  nchar(controls$DDE_TUPLE_COPY) - 1
)
controls <- controls %>%
  separate(DDE_TUPLE_COPY, c(
    "DRUG_1_CONCEPT_ID", "DRUG_2_CONCEPT_ID",
    "EVENT_CONCEPT_ID"
  ), sep = ", ") %>% mutate(EVENT_CONCEPT_ID = as.numeric(EVENT_CONCEPT_ID))

### PART B: Estimate AE rate distribution

# Load the drug-event pair FAERS count tables from SQL
# Connect to PostgreSQL
# Replace with your PostgreSQL credentials
con <- DBI::dbConnect(
  odbc::odbc(),
  Driver   = "PostgreSQL Unicode(x64)",
  Server   = host,
  Database = database,
  UID = username,
  PWD = password,
  Port = port
)

# Find total report count (both ISR and primaryid codes)
standard_case_drug_ingr <-
  dbSendQuery(con,
              "SELECT * FROM faers.standard_case_drug_ingr_level_unique")
standard_case_drug_ingr_tbl <- dbFetch(standard_case_drug_ingr)
dbClearResult(standard_case_drug_ingr)

standard_case_outcome <-
  dbSendQuery(con, "SELECT * FROM faers.standard_case_outcome")
standard_case_outcome_tbl <- dbFetch(standard_case_outcome)
dbClearResult(standard_case_outcome)

# Find reports (both ISR and primaryid codes) that contain at least one drug and one PT (outcome)
intersect_primaryid <-
  intersect(unique(na.omit(standard_case_drug_ingr_tbl$primaryid)),
            unique(na.omit(standard_case_outcome_tbl$primaryid)))
intersect_isr <-
  intersect(unique(na.omit(standard_case_drug_ingr_tbl$isr)),
            unique(na.omit(standard_case_outcome_tbl$isr)))

# Keep only rows from reports that contain both drug(s) and PT(s)
standard_case_drug_ingr_tbl_new <-
  standard_case_drug_ingr_tbl[(standard_case_drug_ingr_tbl$isr %in% intersect_isr) |
                                (standard_case_drug_ingr_tbl$primaryid %in% intersect_primaryid),]
standard_case_outcome_tbl_new <-
  standard_case_outcome_tbl[(standard_case_outcome_tbl$isr %in% intersect_isr) |
                              (standard_case_outcome_tbl$primaryid %in% intersect_primaryid),]

total_report_count <-
  length(intersect_primaryid) + length(intersect_isr)

# For N_11 counts
standard_drug_ingr_outcome_count <-
  dbSendQuery(con, "SELECT * FROM faers.standard_drug_ingr_outcome_count")
standard_drug_ingr_outcome_count_tbl <-
  dbFetch(standard_drug_ingr_outcome_count)
dbClearResult(standard_drug_ingr_outcome_count)

# For N_1. counts
standard_drug_ingr_count_tbl <- standard_case_drug_ingr_tbl_new %>%
  group_by(standard_ingr_concept_id) %>%
  summarise(drug_count = n())

# For N_.1 counts
standard_outcome_count_tbl <- standard_case_outcome_tbl_new %>%
  group_by(outcome_concept_id) %>%
  summarise(outcome_count = n())

# Join drug and outcome counts
standard_drug_ingr_outcome_count_tbl <-
  standard_drug_ingr_outcome_count_tbl %>%
  left_join(standard_drug_ingr_count_tbl,
            by = c("drug_concept_id" = "standard_ingr_concept_id")) %>%
  left_join(standard_outcome_count_tbl, by = "outcome_concept_id") %>%
  rename(N_11 = drug_ingr_outcome_pair_count,
         N_1. = drug_count,
         N_.1 = outcome_count)

# Calculate N_10, N_01, and N_00
standard_drug_ingr_outcome_count_tbl$N_10 <-
  standard_drug_ingr_outcome_count_tbl$N_1. - standard_drug_ingr_outcome_count_tbl$N_11

standard_drug_ingr_outcome_count_tbl$N_01 <-
  standard_drug_ingr_outcome_count_tbl$N_.1 - standard_drug_ingr_outcome_count_tbl$N_11

standard_drug_ingr_outcome_count_tbl$N_00 <- total_report_count -
  standard_drug_ingr_outcome_count_tbl$N_11 -
  standard_drug_ingr_outcome_count_tbl$N_10 -
  standard_drug_ingr_outcome_count_tbl$N_01

standard_drug_ingr_outcome_count_tbl$ae_rate <-
  (
    standard_drug_ingr_outcome_count_tbl$N_11 +
      standard_drug_ingr_outcome_count_tbl$N_01
  ) / total_report_count

# Get a new tbl with all AEs reported to FAERS and their observed rates
ae_rate_unique_tbl <-
  standard_drug_ingr_outcome_count_tbl %>%
  #filter(outcome_concept_id %in% intersect_ctls_tbl$PT_CONCEPT_ID) %>%
  select(outcome_concept_id, N_.1, ae_rate) %>%
  distinct() %>%
  left_join(controls[, c("EVENT_CONCEPT_ID")],
            by = c("outcome_concept_id" = "EVENT_CONCEPT_ID")) %>%
  distinct()

# Function to estimate Beta parameters
estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}

hist(ae_rate_unique_tbl$ae_rate[ae_rate_unique_tbl$ae_rate < 0.00001],
     freq = FALSE, breaks = 50)
curve(dbeta(x, a, b),
      add = TRUE, col = "red", 
      lwd = 2)

# Histogram of AE observed rates
ggoutlier_hist(ae_rate_unique_tbl, 
               "ae_rate", 
               cut_off_ceiling = 0.002,
               binwidth = 0.0001) +
  stat_function(fun = function(x) dbeta(x, 0.05, 0.5), color = "red",
                size = 1)

curve(dbeta(x, a, b), col = "red", lwd = 2) +
curve(dbeta(x, 0.5*a, 0.5*b), col = "blue", lwd = 2) +
curve(dbeta(x, 1, 1), col = "black", lwd = 2)

rates_mu <- mean(ae_rate_unique_tbl[,"ae_rate"])
rates_var <- var(ae_rate_unique_tbl[,"ae_rate"])
# Estimate Beta hyperparameters
params <- estBetaParams(rates_mu, rates_var)

a <- params$alpha
b <- params$beta


### PART C: Modelling

# Function to calculate log-likelihoods in eq.63-70
# Input values: - the eight counts from FAERS;
#               - observed positive examples (i.e. event is present)
#                 for the different cases of presence/absence of the
#                 two drugs under consideration
#               - observed negative examples (i.e. event is not present)
#                 for the different cases of presence/absence of the
#                 two drugs under consideration
loglik <- function(n_001,n_000,n_101,n_100,n_011,n_010,n_111,n_110,
                   a_00,a_10,a_01,a_11,
                   b_00,b_10,b_01,b_11) {
  extraDistr::dbbinom(n_001, n_001 + n_000,
                      a_00,
                      b_00, log = TRUE) +
    extraDistr::dbbinom(n_101, n_101 + n_100,
                        a_10,
                        b_10, log = TRUE) +
    extraDistr::dbbinom(n_011, n_011 + n_010,
                        a_01,
                        b_01, log = TRUE) +
    extraDistr::dbbinom(n_111, n_111 + n_110,
                        a_11,
                        b_11, log = TRUE)
}

a <- 0.01714869
b <- 112.5411


hyperparam_1 <- list(a_0 = 1, a_1 = 1, a_2 = 1, a_12_1 = 0.5, a_12_2 = 0.5,
                     b_0 = 1, b_1 = 1, b_2 = 1, b_12_1 = 0.5, b_12_2 = 0.5)
hyperparam_2 <- list(a_0 = 1, a_1 = 1, a_2 = 1, a_12_1 = 0.5, a_12_2 = 0.5,
                     b_0 = 0.5, b_1 = 1, b_2 = 1, b_12_1 = 0.5, b_12_2 = 0.5)
hyperparam_3 <- list(a_0 = 3, a_1 = 1, a_2 = 1, a_12_1 = 0.5, a_12_2 = 0.5,
                     b_0 = 0.5, b_1 = 1, b_2 = 1, b_12_1 = 0.5, b_12_2 = 0.5)
hyperparam_4 <- list(a_0 = a, a_1 = 1, a_2 = 1, a_12_1 = 0.5, a_12_2 = 0.5,
                     b_0 = b, b_1 = 1, b_2 = 1, b_12_1 = 0.5, b_12_2 = 0.5)
hyperparam_5 <- list(a_0 = a, a_1 = 1, a_2 = 1, a_12_1 = 0.5, a_12_2 = 0.5,
                     b_0 = b, b_1 = 0.01, b_2 = 0.01, b_12_1 = 0.005, b_12_2 = 0.005)
hyperparam_6 <- list(a_0 = a*0.95, a_1 = 1, a_2 = 1, a_12_1 = 0.5, a_12_2 = 0.5,
                     b_0 = b*0.95, b_1 = 1, b_2 = 1, b_12_1 = 0.5, b_12_2 = 0.5)
hyperparameter_sets <- list(hyperparam_1, hyperparam_2, hyperparam_3,
                            hyperparam_4, hyperparam_5, hyperparam_6)

# Create an empty list to store the individual matrices
list_of_matrices <- vector("list", length(hyperparameter_sets))

# This loop calculates the log-likelihoods for each hypothesis and
# each set of hyperparameters.
for (h in 1:length(hyperparameter_sets)) {
  parameters <- hyperparameter_sets[[h]]
  alpha_0 <- parameters$a_0
  alpha_1 <- parameters$a_1
  alpha_2 <- parameters$a_2
  alpha_12_1 <- parameters$a_12_1
  alpha_12_2 <- parameters$a_12_1
  beta_0 <- parameters$b_0
  beta_1 <- parameters$b_1
  beta_2 <- parameters$b_2
  beta_12_1 <- parameters$b_12_1
  beta_12_2 <- parameters$b_12_2
  
  # Create an empty list to store the log-likelihoods
  loglik_values <- list()
  # Create an empty list to store the log posterior probabilities
  log_post_prob <- list()
  # Create an empty list to store the posterior probabilities for increasing rates
  log_post_prob_pos <- list()
    
  for (k in 1:nrow(controls)) {
    d <- as.list(controls[k,])
    
    # H_000
    # Hyperparameter setting (Eq.21-22)
    a0_H000 <- alpha_0 + alpha_1 + alpha_2 + alpha_12_1 + alpha_12_2
    b0_H000 <- beta_0 + beta_1 + beta_2 + beta_12_1 + beta_12_2
    # Log-likelihood calculation (Eq.68)
    loglik_H000 <- loglik(
      d$N_001,d$N_000,d$N_101,d$N_100,d$N_011,d$N_010,d$N_111,d$N_110,
      a0_H000,
      a0_H000 + d$N_001,
      a0_H000 + d$N_001 + d$N_101,
      a0_H000 + d$N_001 + d$N_101 + d$N_011,
      b0_H000,
      b0_H000 + d$N_000,
      b0_H000 + d$N_000 + d$N_100,
      b0_H000 + d$N_000 + d$N_100 + d$N_010
    )
    
    # H_001
    # Hyperparameter setting (Eq.25-28)
    a0_H001 <- alpha_0 + alpha_1 + alpha_2
    b0_H001 <- beta_0 + beta_1 + beta_2
    a12_H001 <- alpha_12_1 + alpha_12_2
    b12_H001 <- beta_12_1 + beta_12_2
    # Log-likelihood calculation (Eq.69)
    loglik_H001 <- loglik(
      d$N_001,d$N_000,d$N_101,d$N_100,d$N_011,d$N_010,d$N_111,d$N_110,
      a0_H001,
      a0_H001 + d$N_001,
      a0_H001 + d$N_001 + d$N_101,
      a12_H001,
      b0_H001,
      b0_H001 + d$N_000,
      b0_H001 + d$N_000 + d$N_100,
      b12_H001
    )
    
    # H_100
    # Hyperparameter setting (Eq.46-49)
    a0_H100 <- alpha_0 + alpha_2 + alpha_12_2
    b0_H100 <- beta_0 + beta_2 + beta_12_2
    a1_H100 <- alpha_1 + alpha_12_1
    b1_H100 <- beta_1 + beta_12_1
    # Log-likelihood calculation (Eq.72)
    loglik_H100 <- loglik(
      d$N_001,d$N_000,d$N_101,d$N_100,d$N_011,d$N_010,d$N_111,d$N_110,
      a0_H100,
      a1_H100,
      a0_H100 + d$N_001,
      a1_H100 + d$N_101,
      b0_H100,
      b1_H100,
      b0_H100 + d$N_000,
      b1_H100 + d$N_100
    )
    
    # H_101
    # Hyperparameter setting (Eq.53-58)
    a0_H101 <- alpha_0 + alpha_2
    b0_H101 <- beta_0 + beta_2
    a1_H101 <- alpha_1
    b1_H101 <- beta_1
    a12_H101 <- alpha_12_1 + alpha_12_2
    b12_H101 <- beta_12_1 + beta_12_2
    # Log-likelihood calculation (Eq.73)
    loglik_H101 <- loglik(
      d$N_001,d$N_000,d$N_101,d$N_100,d$N_011,d$N_010,d$N_111,d$N_110,
      a0_H101,
      a1_H101,
      a0_H101 + d$N_001,
      a12_H101,
      b0_H101,
      b1_H101,
      b0_H101 + d$N_000,
      b12_H101
    )
    
    # H_010
    # Hyperparameter setting (Eq.31-34)
    a0_H010 <- alpha_0 + alpha_1 + alpha_12_1
    b0_H010 <- beta_0 + beta_1 + beta_12_1
    a2_H010 <- alpha_2 + alpha_12_2
    b2_H010 <- beta_2 + beta_12_2
    # Log-likelihood calculation (Eq.70)
    loglik_H010 <- loglik(
      d$N_001,d$N_000,d$N_101,d$N_100,d$N_011,d$N_010,d$N_111,d$N_110,
      a0_H010,
      a0_H010 + d$N_001,
      a2_H010,
      a2_H010 + d$N_011,
      b0_H010,
      b0_H010 + d$N_000,
      b2_H010,
      b2_H010 + d$N_010
    )
    
    # H_011
    # Hyperparameter setting (Eq.38-43)
    a0_H011 <- alpha_0 + alpha_1
    b0_H011 <- beta_0 + beta_1
    a2_H011 <- alpha_2
    b2_H011 <- beta_2
    a12_H011 <- alpha_12_1 + alpha_12_2
    b12_H011 <- beta_12_1 + beta_12_2
    # Log-likelihood calculation (Eq.71)
    loglik_H011 <- loglik(
      d$N_001,d$N_000,d$N_101,d$N_100,d$N_011,d$N_010,d$N_111,d$N_110,
      a0_H011,
      a0_H011 + d$N_001,
      a2_H011,
      a12_H011,
      b0_H011,
      b0_H011 + d$N_000,
      b2_H011,
      b12_H011
    )
    
    # H_111
    # Hyperparater setting
    a0_H111 <- alpha_0
    b0_H111 <- beta_0
    a1_H111 <- alpha_1
    b1_H111 <- beta_1
    a2_H111 <- alpha_2
    b2_H111 <- beta_2
    a12_H111 <- alpha_12_1 + alpha_12_2
    b12_H111 <- beta_12_1 + beta_12_2
    # Log-likelihood calculation (Eq.75)
    loglik_H111 <- loglik(
      d$N_001,d$N_000,d$N_101,d$N_100,d$N_011,d$N_010,d$N_111,d$N_110,
      a0_H111,
      a1_H111,
      a2_H111,
      a12_H111,
      b0_H111,
      b1_H111,
      b2_H111,
      b12_H111
    )
    
    # H_110
    # Hyperparameter setting (Eq.62-67)
    a0_H110 <- alpha_0
    b0_H110 <- beta_0
    a1_H110 <- alpha_1 + alpha_12_1
    b1_H110 <- beta_1 + beta_12_1
    a2_H110 <- alpha_2 + alpha_12_2
    b2_H110 <- beta_2 + beta_12_2
    
    # Numerical approximation of integral in Eq.74
    n.sim <- 1000
    # Beta approximation for theta_1
    T1 <- rbeta(n.sim, a1_H110 + d$N_101, b1_H110 + d$N_100)
    # Beta approximation for theta_2
    T2 <- rbeta(n.sim, a2_H110 + d$N_011, b2_H110 + d$N_010)
    
    loglik.sum <- 0
    for (j in 1:n.sim) {
      loglik_j <- dbinom(
        x = d$N_111,
        size = d$N_111 + d$N_110,
        prob = T1[j] + (1 - T1[j]) * T2[j],
        log = TRUE
      )
      loglik.sum = loglik.sum + loglik_j
    }
    betabinom_11 <- loglik.sum / n.sim
    
    # Log-likelihood calculation (Eq.74)
    loglik_H110 <-
      extraDistr::dbbinom(d$N_001, d$N_001 + d$N_000,
                          a0_H110,
                          b0_H110, log = TRUE) +
      extraDistr::dbbinom(d$N_101, d$N_101 + d$N_100,
                          a1_H110,
                          b1_H110, log = TRUE) +
      extraDistr::dbbinom(d$N_011, d$N_011 + d$N_010,
                          a2_H110,
                          b2_H110, log = TRUE) +
      betabinom_11
    
    # Aggregate the log-likelihoods to a single vector
    loglik_vector <- c(
      'loglik_H_000' = loglik_H000,
      'loglik_H_001' = loglik_H001,
      'loglik_H_100' = loglik_H100,
      'loglik_H_101' = loglik_H101,
      'loglik_H_010' = loglik_H010,
      'loglik_H_011' = loglik_H011,
      'loglik_H_110' = loglik_H110,
      'loglik_H_111' = loglik_H111
    )
    
    # Store the log-likelihood vector back to the list
    loglik_values[[k]] <- loglik_vector
    
    # (a) Posteriors for "surprising" rates
    # Calculate and store the log posterior probabilities to the list
    # (we have assumed that the priors on all models are equal)
    # Eq.78
    log_post_prob_vector <- loglik_vector - logSumExp(loglik_vector)
    # Assign names to the vector elements
    names(log_post_prob_vector) <- c('log_post_H_000', 'log_post_H_001',
                                     'log_post_H_100', 'log_post_H_101',
                                     'log_post_H_010', 'log_post_H_011',
                                     'log_post_H_110', 'log_post_H_111')
    # Store the vector to the list of log posterior values
    log_post_prob[[k]] <- log_post_prob_vector
    
    # (b) Posteriors for increasing rates
    
    # H_001
    # Sample from betas
    # theta_0 (Eq.81)
    T0_H001 <- rbeta(1000,
                     d$N_101 + d$N_011 + d$N_001 + a0_H001,
                     d$N_100 + d$N_010 + d$N_000 + b0_H001)
    # theta_12 (Eq.82)
    T12_H001 <- rbeta(1000, d$N_111 + a12_H001, d$N_110 + b12_H001)
    
    # Use the distribution with the lower variance
    if (d$N_111 + a12_H001 + d$N_110 + b12_H001 < 
        d$N_101 + d$N_011 + d$N_001 + a0_H001 + 
        d$N_100 + d$N_010 + d$N_000 + b0_H001) {
      # Use samples from theta_0
      # Eq.79
      log_post_pos_H001 <-
        log(mean(1 - pbeta(
          T0_H001,
          shape1 = d$N_111 + a12_H001,
          shape2 = d$N_110 + b12_H001,
          log.p = FALSE
        )))
    } else {
      # Use samples from theta_12
      # Eq.80
      log_post_pos_H001 <-
        mean(
          pbeta(
            T12_H001,
            shape1 = d$N_101 + d$N_011 + d$N_001 + a0_H001,
            shape2 = d$N_100 + d$N_010 + d$N_000 + b0_H001,
            log.p = TRUE
          )
        )
    }
    
    # H_101
    # Sample from betas
    # theta_0 (Eq.88)
    T0_H101 <- rbeta(1000,
                     d$N_011 + d$N_001 + a0_H101,
                     d$N_010 + d$N_000 + b0_H101)
    # theta_1 (Eq.89)
    T1_H101 <- rbeta(1000,
                     d$N_101 + a1_H101,
                     d$N_100 + b1_H101)
    # theta_12 (Eq.90)
    T12_H101 <- rbeta(1000, 
                      d$N_111 + a12_H101, 
                      d$N_110 + b12_H101)
    
    # Use the distribution with the lower variance
    # First component of the product
    if (d$N_111 + a12_H101 + d$N_110 + b12_H101 > 
        d$N_011 + d$N_001 + a0_H101 + d$N_010 + d$N_000 + b0_H101) {
      # Eq.84
      log_post_pos_H101_first <-
        mean(
          pbeta(
            T12_H101,
            shape1 = d$N_011 + d$N_001 + a0_H101,
            shape2 = d$N_010 + d$N_000 + b0_H101,
            log.p = TRUE
          )
        )
    } else {
      # Eq.85
      log_post_pos_H101_first <-
        log(mean(1 - pbeta(
          T0_H101,
          shape1 = d$N_111 + a12_H101,
          shape2 = d$N_110 + b12_H101,
          log.p = FALSE
        ))) 
    }
    # Second component of the product
    if (d$N_111 + a12_H101 + d$N_110 + b12_H101 > 
        d$N_101 + a1_H101 + d$N_100 + b1_H101) {
      # Eq.86
      log_post_pos_H101_second <-
        mean(
          pbeta(
            T12_H101,
            shape1 = d$N_101 + a1_H101,
            shape2 = d$N_100 + b1_H101,
            log.p = TRUE
          )
        )
    } else {
      # Eq.87
      log_post_pos_H101_second <-
        log(mean(1 - pbeta(
          T1_H101,
          shape1 = d$N_111 + a12_H101,
          shape2 = d$N_110 + b12_H101,
          log.p = FALSE
        ))) 
    }
    
    #Eq.83
    log_post_pos_H101 <- log_post_pos_H101_first + log_post_pos_H101_second
    
    # H_011
    # Sample from betas
    # theta_0 (Eq.96)
    T0_H011 <- rbeta(1000, 
                     d$N_101 + d$N_001 + a0_H011,
                     d$N_100 + d$N_000 + b0_H011)
    # theta_2 (Eq.97)
    T2_H011 <- rbeta(1000, 
                     d$N_011 + a2_H011,
                     d$N_010 + b2_H011)
    # theta_12 (Eq.98)
    T12_H011 <- rbeta(1000, 
                      d$N_111 + a12_H011,
                      d$N_110 + b12_H011)
    
    # Use the distribution with the lower variance
    # First component of the product
    if (d$N_111 + a12_H011 + d$N_110 + b12_H011 > 
        d$N_101 + d$N_001 + a0_H011 + d$N_100 + d$N_000 + b0_H011) {
      # Eq.92
      log_post_pos_H011_first <-
        mean(
          pbeta(
            T12_H011,
            shape1 = d$N_101 + d$N_001 + a0_H011,
            shape2 = d$N_100 + d$N_000 + b0_H011,
            log.p = TRUE
          )
        )
    } else {
      # Eq.93
      log_post_pos_H011_first <-
        log(mean(1 - pbeta(
          T0_H011,
          shape1 = d$N_111 + a12_H011,
          shape2 = d$N_110 + b12_H011,
          log.p = FALSE)
        )
        ) 
    }
    
    # Second component of the product
    if (d$N_111 + a12_H011 + d$N_110 + b12_H011 > 
        d$N_011 + a2_H011 + d$N_010 + b2_H011 
    ) {
      # Eq.94
      log_post_pos_H011_second <-
        mean(
          pbeta(
            T12_H011,
            shape1 = d$N_011 + a2_H011,
            shape2 = d$N_010 + b2_H011,
            log.p = TRUE
          )
        )
    } else {
      # Eq.95
      log_post_pos_H011_second <-
        log(mean(1 - pbeta(
          T2_H011,
          shape1 = d$N_111 + a12_H011,
          shape2 = d$N_110 + b12_H011,
          log.p = FALSE
        ))
        ) 
    }
    # Eq.91
    log_post_pos_H011 <- log_post_pos_H011_first + log_post_pos_H011_second
    
    # H_111
    # Sample from betas
    # theta_0 (Eq.104)
    T0_H111 <- rbeta(1000,
                     d$N_001 + a0_H111,
                     d$N_000 + b0_H111)
    # theta_1 (Eq.105)
    T1_H111 <- rbeta(1000, 
                     d$N_101 + a1_H111,
                     d$N_100 + b1_H111)
    # theta_2 (Eq.106)
    T2_H111 <- rbeta(1000, 
                     d$N_011 + a2_H111,
                     d$N_010 + b2_H111)
    # theta_12 (Eq.107)
    T12_H111 <- rbeta(1000, 
                      d$N_111 + a12_H111,
                      d$N_110 + b12_H111)
    
    # Use the distribution with the lower variance
    # First component of the product
    if (d$N_111 + a12_H111 + d$N_110 + b12_H111 > 
        d$N_001 + a0_H111 + d$N_000 + b0_H111) {
      # Eq.100
      log_post_pos_H111_first <-
        mean(pbeta(
          T12_H111,
          shape1 = d$N_001 + a0_H111,
          shape2 = d$N_000 + b0_H111,
          log.p = TRUE
        )) 
    } else {
      # Eq.101
      log_post_pos_H111_first <- 
        log(mean(1 - pbeta(
          T0_H111,
          shape1 = d$N_111 + a12_H111,
          shape2 = d$N_110 + b12_H111,
          log.p = FALSE
        ))) }
    
    # Second component of the product
    if (d$N_111 + a12_H111 + d$N_110 + b12_H111 > 
        d$N_101 + a1_H111 + d$N_100 + b1_H111) {
      # Eq.102
      log_post_pos_H111_second <- 
        mean(pbeta(T12_H111,
                   shape1 = d$N_101 + a1_H111,
                   shape2 = d$N_100 + b1_H111,
                   log.p = TRUE))
    } else {
      # Eq.103
      log_post_pos_H111_second <- 
        log(mean(1 - pbeta(
          T1_H111,
          shape1 = d$N_111 + a12_H111,
          shape2 = d$N_110 + b12_H111,
          log.p = FALSE
        )))
    }
    # Third component of the product
    if (d$N_111 + a12_H111 + d$N_110 + b12_H111 > 
        d$N_011 + a2_H111 + d$N_010 + b2_H111) {
      # Eq.104
      log_post_pos_H111_third <- 
        mean(pbeta(
          T12_H111,
          shape1 = d$N_011 + a2_H111,
          shape2 = d$N_010 + b2_H111,
          log.p = TRUE
        ))
    } else {
      # Eq.105
      log_post_pos_H111_third <- 
        log(mean(1 - pbeta(
          T2_H111,
          shape1 = d$N_111 + a12_H111,
          shape2 = d$N_110 + b12_H111,
          log.p = FALSE
        )))
    }
    
    # Eq.99    
    log_post_pos_H111 <- log_post_pos_H111_first + log_post_pos_H111_second + 
      log_post_pos_H111_third
    
    # Eq.108
    postModel_pos_prob_001 <- log_post_pos_H001 + unname(log_post_prob[[k]][2])
    # Eq.109
    postModel_pos_prob_101 <- log_post_pos_H101 + unname(log_post_prob[[k]][4])
    # Eq.110
    postModel_pos_prob_011 <- log_post_pos_H011 + unname(log_post_prob[[k]][6])
    # Eq.111
    postModel_pos_prob_111 <- log_post_pos_H111 + unname(log_post_prob[[k]][8])
    
    # Aggregate the posterior probabilities for increasing rates
    # to a single vector
    log_post_pos_vector <- c(
      'log_post_pos_H_001' = postModel_pos_prob_001,
      'log_post_pos_H_101' = postModel_pos_prob_101,
      'log_post_pos_H_011' = postModel_pos_prob_011,
      'log_post_pos_H_111' = postModel_pos_prob_111
    )
    
    # Store the posterior probabilities for increasing rates to the list
    log_post_prob_pos[[k]] <- log_post_pos_vector
  }

  # Convert the list of loglikelihoods to a matrix
  log_lik <- do.call(rbind, loglik_values)
  # Convert the log posterior probabilities list to a matrix
  log_post_prob <- do.call(rbind, log_post_prob)
  # Convert the log posterior probabilities list for increasing rates to a matrix
  log_post_prob_pos <- do.call(rbind, log_post_prob_pos)
  # Append the matrices to the original df
  controls_post_prob <- cbind(controls, log_lik, log_post_prob, log_post_prob_pos)
  
  # Logsumexp of log-likelihoods
  # Positive outcome
  # LogSumExp of loglik_H_001, loglik_H_101, loglik_H_011, loglik_H_111 
  controls_post_prob$logsumexp_H_..1 <- rowLogSumExps(
    matrix(unlist(controls_post_prob[, c(15, 17, 19, 21)]), ncol=4))
  
  # Negative outcome
  # LogSumExp of loglik_H_000, loglik_H_100, loglik_H_010, loglik_H_110
  controls_post_prob$logsumexp_H_..0 <- rowLogSumExps(
    matrix(unlist(controls_post_prob[,c(14, 16, 18, 20)]), ncol=4))
  
  # Difference of LogSumExp values
  controls_post_prob$ratio <- controls_post_prob$logsumexp_H_..1 -
    controls_post_prob$logsumexp_H_..0
  
  # Add posterior probabilities
  controls_post_prob$log_post_H_..1 <- log(exp(controls_post_prob$log_post_H_001) + 
                                             exp(controls_post_prob$log_post_H_101) +
                                             exp(controls_post_prob$log_post_H_011) + 
                                             exp(controls_post_prob$log_post_H_111))
  controls_post_prob$log_post_H_..0 <- log(exp(controls_post_prob$log_post_H_000) + 
                                             exp(controls_post_prob$log_post_H_100) +
                                             exp(controls_post_prob$log_post_H_010) + 
                                             exp(controls_post_prob$log_post_H_110))
  
  controls_post_prob$sum_posterior_pos <- log(exp(controls_post_prob$log_post_pos_H_001) + 
                                                exp(controls_post_prob$log_post_pos_H_101) + 
                                                exp(controls_post_prob$log_post_pos_H_011) + 
                                                exp(controls_post_prob$log_post_pos_H_111))
  
  controls_post_prob$post_odds_ratio <- controls_post_prob$sum_posterior_pos - 
    controls_post_prob$log_post_H_..0
  
  # Append the matrices to the corresponding lists
  list_of_matrices[[h]] <- controls_post_prob
  
}

# Use ratio to plot the ROC curve
pred_hyper_set1 <- pROC::roc(list_of_matrices[[1]]$CLASS,
                  list_of_matrices[[1]]$ratio, 
                  direction = "<")
pred_hyper_set2 <- pROC::roc(list_of_matrices[[2]]$CLASS,
                             list_of_matrices[[2]]$ratio, 
                             direction = "<")
pred_hyper_set3 <- pROC::roc(list_of_matrices[[3]]$CLASS,
                             list_of_matrices[[3]]$ratio, 
                             direction = "<")
pred_hyper_set4 <- pROC::roc(list_of_matrices[[4]]$CLASS,
                             list_of_matrices[[4]]$ratio, 
                             direction = "<")
pred_hyper_set5 <- pROC::roc(list_of_matrices[[5]]$CLASS,
                             list_of_matrices[[5]]$ratio, 
                             direction = "<")
pred_hyper_set6 <- pROC::roc(list_of_matrices[[6]]$CLASS,
                             list_of_matrices[[6]]$ratio, 
                             direction = "<")


roc.list <- list('Hyperparameter set 1' = pred_hyper_set1,
                 'Hyperparameter set 2' = pred_hyper_set2,
                 'Hyperparameter set 3' = pred_hyper_set3,
                 'Hyperparameter set 4' = pred_hyper_set4,
                 'Hyperparameter set 5' = pred_hyper_set5,
                 'Hyperparameter set 6' = pred_hyper_set6)

# extract auc
data.auc <- roc.list %>% 
  map(~tibble(AUC = .x$auc)) %>% 
  bind_rows(.id = "name")

# generate labels labels
data.labels <- data.auc %>% 
  mutate(label_long=paste0(name," , AUC = ",paste(round(AUC,3))),
         label_AUC=paste0("AUC = ",paste(round(AUC,3))))

png("output/figures/sensitivity_analysis_hyperparameters.png",
    units = "in",
    width = 12,
    height = 5,
    res = 300)
pROC::ggroc(roc.list, legacy.axes=TRUE) +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1),
               color="darkgrey", linetype="dashed") +
  scale_color_discrete(labels=data.labels$label_long) +
  theme_minimal() + theme(legend.title=element_blank())
dev.off()

# Use posterior odds ratio to plot the ROC curve
pred <- pROC::roc(list_of_matrices[[1]][is.finite(list_of_matrices[[1]]$post_odds_ratio),]$CLASS, 
                  list_of_matrices[[1]][is.finite(list_of_matrices[[1]]$post_odds_ratio),]$post_odds_ratio,
                  direction = "<")
plot(pred, print.auc=TRUE)
title("Two drugs - ROC curve (posterior log odds ratio)", line = 3)

## Quality check
controls_post_prob$sums <- exp(controls_post_prob$log_post_H_000) +
  exp(controls_post_prob$log_post_H_001) + exp(controls_post_prob$log_post_H_100) +
  exp(controls_post_prob$log_post_H_101) + exp(controls_post_prob$log_post_H_010) +
  exp(controls_post_prob$log_post_H_011) + exp(controls_post_prob$log_post_H_110) +
  exp(controls_post_prob$log_post_H_111)

# Write output to a .csv file
write.csv(
  subset(controls_post_prob, select = -c(2:9,11:13, 41)),
  "output/bayesian_approach_two_drugs_scores.csv"
)
