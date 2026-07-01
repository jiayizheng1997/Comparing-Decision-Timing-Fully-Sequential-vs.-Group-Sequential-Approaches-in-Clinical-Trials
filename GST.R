# Function to generate arrival times for n patients
# - T0: end time for recruiting (all arrivals <= T0)
# - n: total number of patients
# Output: vector of length n with sorted arrival times tau_1 < tau_2 < ... < tau_n, all in [0, T0]
arrival_times <- function(T0, n) {
  # Generate n uniform random times in [0, T0]
  unsorted_times <- runif(n, min = 0, max = T0)
  
  # Sort them to get the ordered arrival times
  tau <- sort(unsorted_times)
  
  return(tau)
}

# Function: to attribute patients into 2 sites
# Input: n = total number of patients, lambda1 = Poisson parameter of site 1, lambda2 for site 2
# Output: vector of length n with ~n/2 ones and ~n/2 twos, randomly shuffled
site_attribution <- function(n, lambda1, lambda2) {
  # Number in each group (balanced as possible)
  n1 <- floor(n*(lambda1/(lambda1+lambda2)))
  n2 <- n - n1
  
  # Create two sequences
  group1 <- rep(1, n1)
  group2 <- rep(2, n2)
  
  # Combine and shuffle
  assignments <- c(group1, group2)
  assignments <- sample(assignments)  # random permutation
  
  return(assignments)
}

# Function: to attribute patients into 2 treatments.
# The assumed probability is 0.5, a fair coin distribution.
# p is the probability for treatment.
# but using balanced shuffling to guarantee (almost) equal group sizes.
treatment_attribution <- function(n, p =0.5) {
  # Number in each treatment (balanced as possible)
  n1 <- floor(n*p)          # treatment 1
  n2 <- n - n1                # treatment 2
  
  # Create two sequences
  treat1 <- rep(1, n1)
  treat2 <- rep(2, n2)
  
  # Combine and shuffle (random permutation)
  assignments <- c(treat1, treat2)
  assignments <- sample(assignments)
  
  return(assignments)
}

# Function to generate survival times (Y_k) for each patient
# Inputs:
# arrival_times: vector of arrival times (length n, tau_1 to tau_n)
# sites: vector of site attributions (1 or 2, length n)
# treats: vector of treatment attributions (1 or 2, length n)
# omega1, omega2: site effects (default examples)
# theta1, theta2: treatment effects (or use delta for reparameterization)
# p: treatment probability
# theta, delta: our hypotheses
# n: total patient number
# Output: data frame with patient_id, arrival_time, site, treatment, survival_time (Y)
generate_patient_data <- function(lambda1 = 2, lambda2 = 3, omega1 = 0.9, omega2 = 1.1, theta = 1, delta = 1.1, T0 = 1000, Tf = 500, n = 10000, p = 0.5) {
  
  # Step 1: Arrival time of patients
  arrival_time = arrival_times(T0,n)
  
  # Step 2: Site and treatment assignment (balanced shuffle)
  sites <- site_attribution(n, lambda1, lambda2)       
  treats <- treatment_attribution(n, p) 
  
  # Step 3: Generate true survival times Y_k
  thetas <- ifelse(treats == 1, 1/delta, delta) # can we assume theta1 to be one?
  omega <- c(omega1, omega2)
  Y <- numeric(n)
  for (k in 1:n) {
    rate <- omega[sites[k]] * thetas[treats[k]]
    # print(rate)
    Y[k] <- rexp(1, rate)
  }
  
  # Step 4: Create data frame (including possible censoring at Tf)
  df <- data.frame(
    patient_id     = 1:n,
    arrival_time   = arrival_time,
    site           = sites,
    treatment      = treats,
    true_survival  = Y,                     # uncensored Y_k
    max_followup   = Tf + T0 - arrival_time,    # max possible follow-up 
    censored_time  = pmin(Y, Tf + T0 - arrival_time),  # observed time if censored at Tf
    status = as.integer(Y <= Tf + T0 - arrival_time)  # 1 if event by Tf
  )
  
  return(df)
}

# Calculation of K11, K12 ...T21, T22
# start_row, end_row: divide data with GST; allowance of full data as first and last row

calc_KT <- function(data,start_row,end_row){
  
  data = data[start_row:end_row,] 
  
  K11 <- sum(data$status[data$site == 1 & data$treatment == 1])
  K12 <- sum(data$status[data$site == 1 & data$treatment == 2])
  K21 <- sum(data$status[data$site == 2 & data$treatment == 1])
  K22 <- sum(data$status[data$site == 2 & data$treatment == 2])
  
  T11 <- sum(data$censored_time[data$site == 1 & data$treatment == 1])
  T12 <- sum(data$censored_time[data$site == 1 & data$treatment == 2])
  T21 <- sum(data$censored_time[data$site == 2 & data$treatment == 1])
  T22 <- sum(data$censored_time[data$site == 2 & data$treatment == 2])
  
  return(c(K11,K12,K21,K22,T11,T12,T21,T22))
}

# Log-likelihood function (eq. 3.1 in chapter, with theta1=1/delta, theta2=delta)
# KTrow = c(K11,K12,K21,K22,T11,T12,T21,T22)
loglik <- function(delta,KTrow) {
  
  theta1 <- 1 / delta
  theta2 <- delta
  
  list(K11, K12, K21, K22, T11, T12, T21, T22) <- KTrow
  
  # Conditional MLE for omega_i (eq. 3.2)
  omega1 <- (K11 + K12) / (theta1 * T11 + theta2 * T12)
  omega2 <- (K21 + K22) / (theta1 * T21 + theta2 * T22)
  
  # Log-likelihood
  ll <- K11 * log(omega1 * theta1) - omega1 * theta1 * T11 +
    K12 * log(omega1 * theta2) - omega1 * theta2 * T12 +
    K21 * log(omega2 * theta1) - omega2 * theta1 * T21 +
    K22 * log(omega2 * theta2) - omega2 * theta2 * T22
  
  return(ll)
}

opt <- optim(par = 1, fn = loglik, method = "L-BFGS-B", lower = 0.01, upper = 10, control = list(fnscale = -1))



# Use formula 3.3 from Dr Kim's dissertation
# Function for parametric estimation (maximize log-likelihood for delta)

# data: the full survival dataset
# sub_n: the number of entries when considering for the current group in GST 

parametric_fit <- function(data) {
  # Compute K_ij and T_ij (failures and total time on test)
  K11 <- sum(data$status[data$site == 1 & data$treatment == 1])
  K12 <- sum(data$status[data$site == 1 & data$treatment == 2])
  K21 <- sum(data$status[data$site == 2 & data$treatment == 1])
  K22 <- sum(data$status[data$site == 2 & data$treatment == 2])
  
  T11 <- sum(data$censored_time[data$site == 1 & data$treatment == 1])
  T12 <- sum(data$censored_time[data$site == 1 & data$treatment == 2])
  T21 <- sum(data$censored_time[data$site == 2 & data$treatment == 1])
  T22 <- sum(data$censored_time[data$site == 2 & data$treatment == 2])
  
  # Log-likelihood function (eq. 3.1 in chapter, with theta1=1/delta, theta2=delta)
  loglik <- function(delta) {
    theta1 <- 1 / delta
    theta2 <- delta
    
    # Conditional MLE for omega_i (eq. 3.2)
    omega1 <- (K11 + K12) / (theta1 * T11 + theta2 * T12)
    omega2 <- (K21 + K22) / (theta1 * T21 + theta2 * T22)
    
    # Log-likelihood
    ll <- K11 * log(omega1 * theta1) - omega1 * theta1 * T11 +
      K12 * log(omega1 * theta2) - omega1 * theta2 * T12 +
      K21 * log(omega2 * theta1) - omega2 * theta1 * T21 +
      K22 * log(omega2 * theta2) - omega2 * theta2 * T22
    
    return(ll)
  }
  
  # Optimize (maximize loglik)
  opt <- optim(par = 1, fn = loglik, method = "L-BFGS-B", lower = 0.01, upper = 10, control = list(fnscale = -1))
  
  # Estimated delta and HR = delta^2
  est_delta <- opt$par
  
  return(est_delta)
}

# Then we need a function to calculate a and b, based on alpha and beta




set.seed(1234)
lambda1 = 2
lambda2 = 3
omega1 = 0.9
omega2 = 1.1
theta = 0.01
delta = 1.1
T0 = 10
Tf = 5
n = 10000
p = 0.5

# data generation
df = generate_patient_data(lambda1,lambda2,omega1,omega2,theta,delta,T0,Tf,n,p)

# write.csv # maybe we want to save a csv file. but it need only 1 second to run everything

parametric_fit(df)













