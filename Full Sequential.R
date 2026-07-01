# Corrected Parametric Full Sequential Simulation for Table 3.2: E[exp(-R_a)]

# lambda1, lambda2: poisson parameter for 2 sites.
# omega1, omega2: site effect for 2 sites
# delta: 
# a: 
# t0: max recruiting time
# tf: experiment end time


# Function (corrected indexing, full Lambda, edge handling)
simulate_parametric_exp_R <- function(lambda1 = 2, lambda2 = 3, omega1 = 0.9, omega2 = 1.1, delta = 0.5, a = 2, t0 = 1000, tf = 500) {
  lambda_total <- lambda1 + lambda2
  
  # Generate arrival times until we exceed T0
  
  arrival_times <- c(0)
  
  while (TRUE) {
    next_interarrival <- rexp(1, rate = lambda_total)
    next_arrival <- tail(arrival_times, 1) + next_interarrival
    
    if (next_arrival > t0) break
    
    arrival_times <- c(arrival_times, next_arrival)
  }
  
  max_n = length(arrival_times)
  
  # Keep in mind that p in rbinom is the probability for 1.
  site <- rbinom(max_n, 1, lambda2 / lambda_total) + 1 #lambda2 is 2 and lambda1 is 1.
  treat <- rbinom(max_n, 1, 0.5) + 1 # assume that patients get each treatment with p=0.5
  theta <- c(1/delta, delta)
  Y <- rep(0, max_n) #Y is the true life time for each patient.
  omega <- c(omega1, omega2)
  for (k in 1:max_n) {
    rate <- omega[site[k]] * theta[treat[k]]
    Y[k] <- rexp(1, rate)
  }
  
  Lambda <- rep(0, max_n)
  n_stop <- NA
  for (n in 1:max_n) {
    K <- matrix(0, 2, 2)
    T <- matrix(0, 2, 2)
    for (k in 1:n) {
      i <- site[k]
      j <- treat[k]
      obs_time <- min(Y[k], arrival_times[n+1] - arrival_times[k])  
      delta_k <- ifelse(Y[k] <= arrival_times[n+1] - arrival_times[k], 1, 0)
      K[i, j] <- K[i, j] + delta_k
      T[i, j] <- T[i, j] + obs_time
    }
    
    # MLE under H0 (theta=1)
    omega_tilde <- rep(0, 2)
    for (i in 1:2) {
      denom <- T[i,1] + T[i,2]
      omega_tilde[i] <- if (denom > 0) (K[i,1] + K[i,2]) / denom else 0
    }
    
    # MLE under H1 (theta=delta)
    omega_hat <- rep(0, 2)
    for (i in 1:2) {
      denom <- (1/delta) * T[i,1] + delta * T[i,2]
      omega_hat[i] <- if (denom > 0) (K[i,1] + K[i,2]) / denom else 0
    }
    
    # Full Lambda_n = ell_H1 - ell_H0 (from eq. 3.1)
    ell_H1 <- 0
    ell_H0 <- 0
    
    for (i in 1:2) { # We need to avoid log(0) here.
      # Under H0
      denom0 <- T[i,1] + T[i,2]
      if (denom0 > 1e-10 && omega_tilde[i] > 1e-10) {  # small threshold for floating point safety
        log_omega_tilde <- log(omega_tilde[i])
        ell_H0 <- ell_H0 + (K[i,1] + K[i,2]) * log_omega_tilde - omega_tilde[i] * denom0
      }
      # else: no contribution (K should also be 0, term = 0)
      
      # Under H1 — more careful with small omega_hat
      denom1 <- (1/delta) * T[i,1] + delta * T[i,2]
      if (denom1 > 1e-10 && omega_hat[i] > 1e-10) {
        ell_H1 <- ell_H1 + 
          K[i,1] * (log(omega_hat[i]) - log(delta)) - omega_hat[i] * T[i,1] / delta +
          K[i,2] * (log(omega_hat[i]) + log(delta)) - omega_hat[i] * delta * T[i,2]
      }
      # else: no contribution
    }
    Lambda[n] <- ell_H1 - ell_H0
    # print(Lambda[n])
    if (Lambda[n] > a) {
      n_stop <- n
      break
    }
  }
  
  if (!is.na(n_stop)) {
    R_a <- Lambda[n_stop] - a
    exp_neg_R <- exp(-R_a)
  } else {
    exp_neg_R <- NA
  }
  return(exp_neg_R)
}

# Simulation loop
set.seed(123)
reps <- 1000  # Increase to 4000
deltas <- c(0.5, 0.75, 1.25, 1.5)
as <- c(0.5, 1, 2, 4)
results <- matrix(0, nrow = length(deltas), ncol = length(as))
for (d in 1:length(deltas)) {
  for (aa in 1:length(as)) {
    exp_neg <- replicate(reps, simulate_parametric_exp_R(delta = deltas[d], a = as[aa]))
    results[d, aa] <- mean(exp_neg, na.rm = TRUE)
  }
}
colnames(results) <- paste("a=", as)
rownames(results) <- paste("delta=", deltas)
print(results)  # Empirical values; expect convergence like Table 3.2 (e.g., delta=0.5: 0.6 -> 0.3)















