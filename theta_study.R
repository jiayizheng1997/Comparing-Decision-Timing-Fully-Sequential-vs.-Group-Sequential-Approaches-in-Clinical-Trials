# Reminder: theta0 and theta1 are used in this code.

# *********************************************************************.
# *								    *.
# *								    *.
# *	m6data.fun(); generating two-site, two-treatment data       *.
# *								    *.
# *
# *********************************************************************.

# We used a different version of data functions here.
# t1, t2: theta1 and theta2

m6data.fun <- function(rate1, rate2, w1, w2, t1, t2, 
                       row.num, 
                       p = 0.5,   # prob of treatment 1 (same for both sites)
                       r = 0.5)   # prob of site 1
{
  # ================================================================
  # m6data.fun()  —  Corrected & completed version (March 2026)
  # Generates exactly 'row.num' patients for a two-site, two-treatment
  # survival trial with exponential arrivals and exponential lifetimes.
  #
  # Each patient is independently assigned:
  #   - Site 1 with probability r, Site 2 with 1-r
  #   - Treatment 1 with probability p, Treatment 2 with 1-p
  #
  # Arrival processes are independent Poisson processes per site.
  # Lifetimes depend on site × treatment × delta as in the rest of your code.
  # Returns a single data.frame sorted by arrival time.
  # ================================================================
  
  # 1. Assign site and treatment to each of the row.num patients
  site  <- sample(c(1, 2), size = row.num, replace = TRUE, 
                  prob = c(r, 1 - r))
  treat <- sample(c(1, 2), size = row.num, replace = TRUE, 
                  prob = c(p, 1 - p))
  
  # Status codes: 11, 12, 21, 22  (exactly as used in g6.fun, Ra6.fun, etc.)
  status <- 10 * site + treat
  
  # 2. Number of patients per site
  n1 <- sum(site == 1)
  n2 <- sum(site == 2)
  
  # 3. Generate arrival times separately for each site (Poisson process)
  arrivals <- numeric(row.num)
  
  if (n1 > 0) {
    arrivals[site == 1] <- cumsum(rexp(n1, rate = rate1))
  }
  if (n2 > 0) {
    arrivals[site == 2] <- cumsum(rexp(n2, rate = rate2))
  }
  
  # 4. Define the four hazard rates (exactly as in your original code)
  w11 <- w1 * t1
  w12 <- w1 * t2
  w21 <- w2 * t1
  w22 <- w2 * t2
  
  # 5. Generate lifetimes (survival times)
  life <- numeric(row.num)
  life[status == 11] <- rexp(sum(status == 11), rate = w11)
  life[status == 12] <- rexp(sum(status == 12), rate = w12)
  life[status == 21] <- rexp(sum(status == 21), rate = w21)
  life[status == 22] <- rexp(sum(status == 22), rate = w22)
  
  # 6. Build the data.frame
  df <- data.frame(
    life     = life,
    arrival  = arrivals,
    status   = status
  )
  
  # 7. Sort by arrival time (required by Ra6.fun, m6test.fun, etc.)
  ii <- order(df$arrival)
  df <- df[ii, ]
  
  # Optional: add row names like the old versions
  rownames(df) <- NULL
  
  return(df)
}

# *********************************************************************.
# *								    *.
# *								    *.
# *		m6test.fun()					    *.
# *								    *.
# *								    *.
# *********************************************************************.

# ================================================================
# UPDATED m6test.fun()  —  Now includes the 2nd-order slow term
# Uses the new g6.grad2.fun() you just got
# Allow different inputs of data delta and estimation delta
# delta_alter cannot be 1, or g, grad would always be zero
# to obtain gamma1, put alternative delta under delta_alter and delta_null
# to obtain gamma0, put alternative delta under delta_alter and 1 under delta_null
# ================================================================

m6test.fun <- function(rate1 = 1, rate2 = 1, 
                       w1 = 1, w2 = 1, 
                       theta1 = 1.25,          # H1: theta = delta
                       theta2 = 1,
                       cut.num = 3000,        # number of steps for min(Zk)
                       mov.num = 30)          # window for stationary correction
{
  # 1. Generate one LONG dataset
  total_n <- cut.num + 2 * mov.num + 10
  data1 <- m6data.fun(rate1 = rate1, rate2 = rate2, 
                      w1 = w1, w2 = w2, t1 = theta1, t2 = theta2, 
                      row.num = total_n)
  
  # 2. Compute all three ingredients exactly as in the theory
  delta_null = 1
  delta_alter = sqrt(theta2/theta1)
  w11 <- w1 * theta1;  w12 <- w1 * theta2
  w21 <- w2 * theta1;  w22 <- w2 * theta2
  
  p1 <- rate1 / (rate1 + rate2)
  p2 <- rate2 / (rate1 + rate2)
  
  nu      <- c(p1/2, p1/2, p2/2, p2/2,
               p1/(2*w11), p1/(2*w12),
               p2/(2*w21), p2/(2*w22))
  
  g       <- g6.fun(delta_alter, nu)
  g.grad  <- g6.grad.fun(delta_alter, nu)
  g.grad2 <- g6.grad2.fun(delta_alter, nu)     # ← NEW: the 8x8 Hessian
  
  # 3. Extract data
  xk     <- data1$life
  tou    <- c(0, data1$arrival)
  status <- data1$status
  B11 <- as.numeric(status == 11)
  B12 <- as.numeric(status == 12)
  B21 <- as.numeric(status == 21)
  B22 <- as.numeric(status == 22)
  
  # 4. Compute ξ̃₀ (initial stationary correction) — same as before
  xi0 <- numeric(8)
  posK <- (xk[1:mov.num] > (tou[2:(mov.num+1)] - tou[1]))
  xi0[1] <- -sum(posK * B11[1:mov.num])
  xi0[2] <- -sum(posK * B12[1:mov.num])
  xi0[3] <- -sum(posK * B21[1:mov.num])
  xi0[4] <- -sum(posK * B22[1:mov.num])
  xi0[5] <- -sum((xk[1:mov.num] - (tou[2:(mov.num+1)] - tou[1])) * posK * B11[1:mov.num])
  xi0[6] <- -sum((xk[1:mov.num] - (tou[2:(mov.num+1)] - tou[1])) * posK * B12[1:mov.num])
  xi0[7] <- -sum((xk[1:mov.num] - (tou[2:(mov.num+1)] - tou[1])) * posK * B21[1:mov.num])
  xi0[8] <- -sum((xk[1:mov.num] - (tou[2:(mov.num+1)] - tou[1])) * posK * B22[1:mov.num])
  
  Xi0 <- sum(g.grad * xi0)
  
  # 5. Main loop — now with full approximation including slow term
  Sn   <- numeric(cut.num)
  Xin  <- numeric(cut.num)
  slow <- numeric(cut.num)
  Zk   <- numeric(cut.num)
  
  for(n in 1:cut.num) {
    # Current sliding window for xin
    idx <- (n + 1):(n + mov.num)
    pos <- (xk[idx] > (tou[idx + 1] - tou[n + 1]))
    
    xin <- numeric(8)
    xin[1] <- -sum(pos * B11[idx])
    xin[2] <- -sum(pos * B12[idx])
    xin[3] <- -sum(pos * B21[idx])
    xin[4] <- -sum(pos * B22[idx])
    xin[5] <- -sum((xk[idx] - (tou[idx + 1] - tou[n + 1])) * pos * B11[idx])
    xin[6] <- -sum((xk[idx] - (tou[idx + 1] - tou[n + 1])) * pos * B12[idx])
    xin[7] <- -sum((xk[idx] - (tou[idx + 1] - tou[n + 1])) * pos * B21[idx])
    xin[8] <- -sum((xk[idx] - (tou[idx + 1] - tou[n + 1])) * pos * B22[idx])
    
    # Cumulative V_n^* up to n
    vvn <- c(sum(B11[1:n]), sum(B12[1:n]), sum(B21[1:n]), sum(B22[1:n]),
             sum(xk[1:n] * B11[1:n]), sum(xk[1:n] * B12[1:n]),
             sum(xk[1:n] * B21[1:n]), sum(xk[1:n] * B22[1:n]))
    
    wn <- vvn - n * nu
    
    # Random walk part
    Sn[n]   <- n * g + sum(g.grad * (vvn - n * nu))
    Xin[n]  <- sum(g.grad * xin)
    
    # NEW: slowly changing term ψ_n (eq. 3.12)
    slow[n] <- (1 / (2 * n)) * (t(wn) %*% g.grad2 %*% wn)
    
    # Full centered process Z_n (eq. 3.13 / 3.15b + slow term)
    Zk[n]   <- Sn[n] + Xi0 - Xin[n] + slow[n]
  }
  
  if (g < 0) Zk <- -Zk                   # flip sign so min becomes the correct M
  # (If you are computing γ₁ instead, do NOT flip)
  
  M <- min(Zk)
  
  cat("m6test.fun completed — M =", round(M, 4), 
      " (with slow term included)\n")
  
  return(list(M = M, 
              Xi0 = Xi0, 
              Sn = Sn, 
              Xin = Xin, 
              slow = slow, 
              Zk = Zk))
}

Ra6.two.fun <- function(rate1, rate2, w1, w2, data_t1, data_t2, grad_t1, grad_t2, a_upper, b_lower, m = m_emp) {
  
  cross_upper <- numeric(m)
  cross_lower <- numeric(m)
  
  for (i in 1:m) {
    
    # Generate one full dataset under the true delta (same as your original call)
    data1 <- m6data.fun(rate1 = rate1, rate2 = rate2,
                        w1 = w1, w2 = w2, t1 = data_t1, t2 = data_t2,
                        row.num = 500)   # you can increase this number if needed
    
    xk   <- data1$life
    tou  <- c(0, data1$arrival)
    status <- data1$status
    
    B11 <- as.numeric(status == 11)
    B12 <- as.numeric(status == 12)
    B21 <- as.numeric(status == 21)
    B22 <- as.numeric(status == 22)
    
    site.num <- length(xk)
    
    # Sequential monitoring — two-sided
    n <- 0
    log.ratio <- 0
    
    while (n < site.num) {
      n <- n + 1
      
      # Compute Kij and Tij with the stationary correction (exactly as in your original code)
      posK11 <- (xk[1:n] > (tou[n+1] - tou[1:n])) * B11[1:n]
      K11 <- sum(B11[1:n]) - sum(posK11)
      
      posK12 <- (xk[1:n] > (tou[n+1] - tou[1:n])) * B12[1:n]
      K12 <- sum(B12[1:n]) - sum(posK12)
      
      posK21 <- (xk[1:n] > (tou[n+1] - tou[1:n])) * B21[1:n]
      K21 <- sum(B21[1:n]) - sum(posK21)
      
      posK22 <- (xk[1:n] > (tou[n+1] - tou[1:n])) * B22[1:n]
      K22 <- sum(B22[1:n]) - sum(posK22)
      
      posT11 <- (xk[1:n] > (tou[n+1] - tou[1:n])) * B11[1:n]
      T11 <- sum(xk[1:n] * B11[1:n]) - sum((xk[1:n] - (tou[n+1] - tou[1:n])) * posT11)
      
      posT12 <- (xk[1:n] > (tou[n+1] - tou[1:n])) * B12[1:n]
      T12 <- sum(xk[1:n] * B12[1:n]) - sum((xk[1:n] - (tou[n+1] - tou[1:n])) * posT12)
      
      posT21 <- (xk[1:n] > (tou[n+1] - tou[1:n])) * B21[1:n]
      T21 <- sum(xk[1:n] * B21[1:n]) - sum((xk[1:n] - (tou[n+1] - tou[1:n])) * posT21)
      
      posT22 <- (xk[1:n] > (tou[n+1] - tou[1:n])) * B22[1:n]
      T22 <- sum(xk[1:n] * B22[1:n]) - sum((xk[1:n] - (tou[n+1] - tou[1:n])) * posT22)
      
      # Skip if any total time-on-test is zero (same safeguard you had)
      if ((T11 + T12) * (T21 + T22) == 0) next
      
      # CORRECT log-likelihood ratio statistic Λ_n (equation 3.3 in the dissertation)
      vn <- c(K11, K12, K21, K22, T11, T12, T21, T22)
      log.ratio <- n * g6.fun(sqrt(grad_t2/grad_t1), vn / n)
      
      # TWO-SIDED stopping rule
      if (log.ratio >= a_upper || log.ratio <= b_lower) {
        break
      }
    }
    
    # Record which boundary was crossed
    cross_upper[i] <- as.numeric(log.ratio >= a_upper)
    cross_lower[i] <- as.numeric(log.ratio <= b_lower)
  }
  
  list(emp_alpha = mean(cross_upper),
       emp_beta  = mean(cross_lower))
}


t1 = Sys.time()

set.seed(123)

library(parallel)

mc.cores = 40

# ================================================================
# FULL REPRODUCTION OF TABLE 3.3 
# ================================================================

# ====================== PARAMETERS  ======================

rate1 <- 1
rate2_list <- c(1,2)
w1    <- 1
w2    <- 1
alpha_target <- 0.05
beta_target  <- 0.1
theta1s <- c(8/9, 8/25, 9/16, 9/25)
theta2s <- c(1/2, 1/2, 1/3, 1/3)

m_gamma      <- 50000      # replications for gamma
m_emp        <- 1000      # replications for empirical error rates
cut.num      <- 3000
mov.num      <- 30




n_iter <- length(theta1s)*length(rate2_list)
boundary_table <- data.frame( # we store the results in this table
  delta       = numeric(n_iter),
  alpha       = numeric(n_iter), 
  beta        = numeric(n_iter),
  theta1      = numeric(n_iter),
  theta2      = numeric(n_iter),
  gamma0      = numeric(n_iter),
  gamma1      = numeric(n_iter),
  lower_b     = numeric(n_iter),
  upper_a     = numeric(n_iter),
  emp_alpha   = numeric(n_iter),
  emp_beta    = numeric(n_iter)
)



# ====================== simulation steps ==========================

cat("delta | nominal (lower, upper) | empirical α | empirical β\n")
cat("----------------------------------------------------------\n")

for (j in seq_along(rate2_list)) {
  for (i in seq_along(theta1s)) {
    theta1 <- theta1s[i]; theta2 <- theta2s[i]
    rate2 = rate2_list[j]
    d = sqrt(theta2/theta1)
    # --- gamma1 from m6test.fun (H1) ---
    M_alt <- unlist(mclapply(  X = 1:m_gamma,               # number of replications 
                               FUN = function(i) {m6test.fun(rate1=rate1, rate2=rate2, w1=w1, w2=w2, theta1 = theta1, theta2 = theta2, cut.num=cut.num, mov.num=mov.num)$M},mc.cores = mc.cores))
    print('Alter M simulation completed')
    
    # Compute nu, g, g.grad, g.grad2
    
    # w11 <- w1*theta1; w12 <- w1*theta2
    # w21 <- w2*theta1; w22 <- w2*theta2
    # p1 <- rate1/(rate1+rate2); p2 <- rate2/(rate1+rate2)
    # nu <- c(p1/2, p1/2, p2/2, p2/2, p1/(2*w11), p1/(2*w12), p2/(2*w21), p2/(2*w22))
    # g       <- g6.fun(d, nu)
    # g.grad  <- g6.grad.fun(d, nu)
    # g.grad2 <- g6.grad2.fun(d, nu)
    # 
    mu1     <- mean(pmax(M_alt,0))
    gamma1  <- mean(1 - exp(-pmax(M_alt, 0))) / mu1
    
    # --- gamma0 from m6test.fun (H0) ---
    M_null <- unlist(mclapply(  X = 1:m_gamma,               # number of replications 
                                FUN = function(i) {m6test.fun(rate1=rate1, rate2=rate2, w1=w1, w2=w2, theta1 = theta1, theta2 = theta2, cut.num=cut.num, mov.num=mov.num)$M},mc.cores = mc.cores))
    # Compute nu, g, g.grad, g.grad2
    # theta1 <- 1; theta2 <- 1
    # 
    # print('Null M simulation completed')
    # w11 <- w1*theta1; w12 <- w1*theta2
    # w21 <- w2*theta1; w22 <- w2*theta2
    # p1 <- rate1/(rate1+rate2); p2 <- rate2/(rate1+rate2)
    # nu <- c(p1/2, p1/2, p2/2, p2/2, p1/(2*w11), p1/(2*w12), p2/(2*w21), p2/(2*w22))
    # g       <- g6.fun(d, nu)
    # g.grad  <- g6.grad.fun(d, nu)
    # g.grad2 <- g6.grad2.fun(d, nu)
    
    mu0     <- mean(pmax(M_null,0))  #mu = g(v) #generate the whole M+, then take average
    gamma0  <- mean(1 - exp(-pmax(M_null, 0))) / mu0
    
    # --- Solve boundaries ---
    ab <- ab.fun(gamma0 = gamma0, gamma1 = gamma1, 
                 alpha = alpha_target, beta = beta_target)
    lower <- ab[1]
    upper <- ab[2]
    
    print('upper and lower bounds calculated')
    
    # --- Empirical errors: TWO separate simulations ---
    emp_alpha <- Ra6.two.fun(rate1=rate1, rate2=rate2, w1=w1, w2=w2, data_t1 = 1, data_t2 = 1, grad_t1 = theta1s[i], grad_t2 = theta2s[i], a_upper=upper, b_lower=lower, m=m_emp)$emp_alpha
    
    emp_beta  <- Ra6.two.fun(rate1=rate1, rate2=rate2, w1=w1, w2=w2, data_t1 = theta1s[i], data_t2 = theta2s[i], grad_t1 = theta1s[i], grad_t2 = theta2s[i], a_upper=upper, b_lower=lower, m=m_emp)$emp_beta
    
    boundary_table[length(theta1s)*(j-1)+i, ] <- c(d, alpha_target, beta_target, theta1s[i], theta2s[i], gamma0, gamma1, lower, upper, emp_alpha, emp_beta)
  }
  
}

# ====================== Display results ======================
print(boundary_table)
t2 = Sys.time()
write.csv(boundary_table, paste0("theta_study.csv"), append = FALSE, row.names = FALSE)









