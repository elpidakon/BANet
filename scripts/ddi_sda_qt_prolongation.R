# Kontsioti, Maskell, Anderson & Pirmohamed, Identifying drug-drug interactions in spontaneous reports utilizing signal detection and biological plausibility aspects (2024)
# This script applies the DDI signal detection algorithm for QT interval prolongation using drug pairs being co-reported with QT interval prolongation in at least 5 reports in FAERS.
# The output is a file with the respective algorithm scores for these drug pairs.

library(extraDistr)
library(matrixStats)
library(dplyr)
library(ROCR)

# Set working directory
#workdir <- "/Users/elpida/OneDrive - The University of Liverpool/BIOLOGICAL PLAUSIBILITY PAPER"
workdir <- "C:/Users/Elpida/OneDrive - The University of Liverpool"
setwd(workdir)
# Set seed value
set.seed(13579)

# Read data file with FAERS counts for TdP
#tdp_counts <- read.csv("tdp_faers_counts_min_5.csv")
counts <- read.csv("qt_prolongation_faers_counts_min_5.csv")
# Generate df with all counts
qt_prolongation_cases <- data.frame(dde_tuple = counts$dde_tuple,
                        N_111 = counts$n_111,
                        N_110 = counts$d1_d2_counter - counts$n_111,
                        N_101 = counts$n_101,
                        N_100 = counts$d1_not_d2_counter - counts$n_101,
                        N_011 = counts$n_011,
                        N_010 = counts$not_d1_d2_counter - counts$n_011,
                        N_001 = counts$n_001,
                        N_000 = counts$not_d1_not_d2_counter - counts$n_001)


# Set hyperparameter values
alpha_0 <- a
alpha_1 <- 1
alpha_2 <- 1
alpha_12_1 <- 0.5
alpha_12_2 <- 0.5
beta_0 <- b
beta_1 <- 1
beta_2 <- 1
beta_12_1 <- 0.5
beta_12_2 <- 0.5

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

# Create an empty list to store the log-likehoods
loglik_values <- list()
# Create an empty list to store the log posterior probabilities
log_post_prob <- list()
# Create an empty list to store the posterior probabilities for increasing rates
log_post_prob_pos <- list()

# This loop first calculates the log-likelihoods for each hypothesis
# and the posteriors in the case of: (a) changing (i.e. "surprising") rates
# and (b) increasing rates.
# This loop calculates the log-likelihoods for each hypothesis.
for (k in 1:nrow(qt_prolongation_cases)) {
  d <- as.list(qt_prolongation_cases[k,])
  
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
cases_log_lik <- do.call(rbind, loglik_values)
# Convert the log posterior probabilities list to a matrix
cases_log_post_prob <- do.call(rbind, log_post_prob)
# Convert the log posterior probabilities list for increasing rates to a matrix
cases_log_post_prob_pos <- do.call(rbind, log_post_prob_pos)
# Append the matrices to the original df
qt_prolongation_post_prob <- cbind(qt_prolongation_cases, cases_log_lik, 
                                   cases_log_post_prob, cases_log_post_prob_pos)

# Logsumexp of log-likelihoods
# Positive outcome
# LogSumExp of loglik_H_001, loglik_H_101, loglik_H_011, loglik_H_111 
qt_prolongation_post_prob$logsumexp_H_..1 <- rowLogSumExps(
  matrix(unlist(qt_prolongation_post_prob[, c(11, 13, 15, 17)]), ncol=4))

# Negative outcome
# LogSumExp of loglik_H_000, loglik_H_100, loglik_H_010, loglik_H_110
qt_prolongation_post_prob$logsumexp_H_..0 <- rowLogSumExps(
  matrix(unlist(qt_prolongation_post_prob[,c(10, 12, 14, 16)]), ncol=4))

# Difference of LogSumExp values
qt_prolongation_post_prob$ratio <- qt_prolongation_post_prob$logsumexp_H_..1 -
  qt_prolongation_post_prob$logsumexp_H_..0

# Add posterior probabilities
qt_prolongation_post_prob$log_post_H_..1 <- log(exp(qt_prolongation_post_prob$log_post_H_001) + 
                                           exp(qt_prolongation_post_prob$log_post_H_101) +
                                           exp(qt_prolongation_post_prob$log_post_H_011) + 
                                           exp(qt_prolongation_post_prob$log_post_H_111))
qt_prolongation_post_prob$log_post_H_..0 <- log(exp(qt_prolongation_post_prob$log_post_H_000) + 
                                           exp(qt_prolongation_post_prob$log_post_H_100) +
                                           exp(qt_prolongation_post_prob$log_post_H_010) + 
                                           exp(qt_prolongation_post_prob$log_post_H_110))

qt_prolongation_post_prob$sum_posterior_pos <- log(exp(qt_prolongation_post_prob$log_post_pos_H_001) + 
                                              exp(qt_prolongation_post_prob$log_post_pos_H_101) + 
                                              exp(qt_prolongation_post_prob$log_post_pos_H_011) + 
                                              exp(qt_prolongation_post_prob$log_post_pos_H_111))

qt_prolongation_post_prob$post_odds_ratio <- qt_prolongation_post_prob$sum_posterior_pos - 
  qt_prolongation_post_prob$log_post_H_..0

# Write output to a .csv file
write.csv(
  qt_prolongation_post_prob,
  "output/qt_prolongation_sda_scores.csv"
)
