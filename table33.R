t1 = Sys.time()

set.seed(123)

library(parallel)

mc.cores = 40

# ================================================================
# FULL REPRODUCTION OF TABLE 3.3 
# ================================================================

# ====================== PARAMETERS  ======================
rate1 <- 1
rate2 <- 1
w1    <- 1
w2    <- 1
alpha_target <- 0.05
beta_targets  <- c(0.1)
deltas <- seq(0.2,2.0,by=0.02)
deltas <- deltas[deltas != 1]


m_gamma      <- 50000      # replications for gamma
m_emp        <- 1000      # replications for empirical error rates
cut.num      <- 3000
mov.num      <- 30




n_iter <- length(deltas)*length(beta_targets)
boundary_table <- data.frame( # we store the results in this table
  delta       = numeric(n_iter),
  alpha       = numeric(n_iter), 
  beta        = numeric(n_iter),
  gamma0      = numeric(n_iter),
  gamma1      = numeric(n_iter),
  lower_b     = numeric(n_iter),
  upper_a     = numeric(n_iter),
  emp_alpha   = numeric(n_iter),
  emp_beta    = numeric(n_iter),
  M_null_med = numeric(n_iter),
  M_alt_med  = numeric(n_iter)
)



# ====================== simulation steps ==========================

cat("delta | nominal (lower, upper) | empirical α | empirical β\n")
cat("----------------------------------------------------------\n")

for (j in seq_along(beta_targets)) {
  for (i in seq_along(deltas)) {
    d <- deltas[i]
    beta_target = beta_targets[j]
    # --- gamma1 from m6test.fun (H1) ---
    M_alt <- unlist(mclapply(  X = 1:m_gamma,               # number of replications 
                               FUN = function(i) {m6test.fun(rate1=rate1, rate2=rate2, w1=w1, w2=w2,delta_alter=d, delta_null=d, cut.num=cut.num, mov.num=mov.num)$M},mc.cores = mc.cores))
    print('Alter M simulation completed')
    
    # Compute nu, g, g.grad, g.grad2
    theta1 <- 1/d; theta2 <- d
    w11 <- w1*theta1; w12 <- w1*theta2
    w21 <- w2*theta1; w22 <- w2*theta2
    p1 <- rate1/(rate1+rate2); p2 <- rate2/(rate1+rate2)
    nu <- c(p1/2, p1/2, p2/2, p2/2, p1/(2*w11), p1/(2*w12), p2/(2*w21), p2/(2*w22))
    g       <- g6.fun(d, nu)
    g.grad  <- g6.grad.fun(d, nu)
    g.grad2 <- g6.grad2.fun(d, nu)
    
    mu1     <- mean(pmax(M_alt,0))
    gamma1  <- mean(1 - exp(-pmax(M_alt, 0))) / mu1
    
    # --- gamma0 from m6test.fun (H0) ---
    M_null <- unlist(mclapply(  X = 1:m_gamma,               # number of replications 
                                FUN = function(i) {m6test.fun(rate1=rate1, rate2=rate2, w1=w1, w2=w2,delta_alter=d, delta_null=1, cut.num=cut.num, mov.num=mov.num)$M},mc.cores = mc.cores))
    print('Null M simulation completed')
    # Compute nu, g, g.grad, g.grad2
    theta1 <- 1; theta2 <- 1
    w11 <- w1*theta1; w12 <- w1*theta2
    w21 <- w2*theta1; w22 <- w2*theta2
    p1 <- rate1/(rate1+rate2); p2 <- rate2/(rate1+rate2)
    nu <- c(p1/2, p1/2, p2/2, p2/2, p1/(2*w11), p1/(2*w12), p2/(2*w21), p2/(2*w22))
    g       <- g6.fun(d, nu)
    g.grad  <- g6.grad.fun(d, nu)
    g.grad2 <- g6.grad2.fun(d, nu)
    
    mu0     <- mean(pmax(M_null,0))  #mu = g(v) #generate the whole M+, then take average
    gamma0  <- mean(1 - exp(-pmax(M_null, 0))) / mu0
    
    # --- Solve boundaries ---
    ab <- ab.fun(gamma0 = gamma0, gamma1 = gamma1, 
                 alpha = alpha_target, beta = beta_target)
    lower <- ab[1]
    upper <- ab[2]
    
    print('upper and lower bounds calculated')
    
    # --- Empirical errors: TWO separate simulations ---
    emp_alpha <- Ra6.two.fun(rate1=rate1, rate2=rate2, w1=w1, w2=w2,grad_delta = d,
                             data_delta=1,          # ← NULL hypothesis
                             a_upper=upper, b_lower=lower, m=m_emp)$emp_alpha
    
    emp_beta  <- Ra6.two.fun(rate1=rate1, rate2=rate2, w1=w1, w2=w2,grad_delta = d,
                             data_delta=d,          # ← ALTERNATIVE hypothesis
                             a_upper=upper, b_lower=lower, m=m_emp)$emp_beta
    
    boundary_table[length(deltas)*(j-1)+i, ] <- c(d, alpha_target, beta_target, gamma0, gamma1, lower, upper, emp_alpha, emp_beta, median(M_null), median(M_alt))
  }
  
}

# ====================== Display results ======================
cat("\n=== ALL CALCULATED BOUNDARIES STORED IN boundary_table ===\n")
print(boundary_table)

cat("\n=== Table 3.3 reproduced (boundaries stored) ===\n")

t2 = Sys.time()
# Optional: save to CSV
write.csv(boundary_table, paste0("Table_3.3_boundaries_rate",rate1,rate2,"_w",w1,w2,".csv"), append = FALSE, row.names = FALSE)









