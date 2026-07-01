# this is our source file.
# The upper part before line xxx are used codes
# The part below line xxx contains unused codes

# *********************************************************************.
# *								    *.
# *								    *.
# *		g6.fun()					    *.
# *		calling Splus function: cmp6.vec()		    *.
# *								    *.
# *********************************************************************.
g6.fun<-function(delta, vn)
{
  # name : g6.fun which is used in the function Ra6.fun.
  # to calculate log-likelihood ratio as a function of Kij and Tij.
  # H0 : theta=1 vs. H1 : theta = delta (delta is not equal to one)
  x11 <- vn[1]
  x12 <- vn[2]
  x21 <- vn[3]
  x22 <- vn[4]
  y11 <- vn[5]
  y12 <- vn[6]
  y21 <- vn[7]
  y22 <- vn[8]
  g <- (x11 + x12) * log((y11 + y12)/(y11 + delta^2 * y12)) + (x21 + x22) * 
    log((y21 + y22)/(y21 + delta^2 * y22)) + 2 * (x12 + x22) * log(
      delta)
  return(g)
}

# *********************************************************************.
# *								    *.
# *								    *.
# *		g6.grad.fun()					    *.
# *		calling Splus function: cmp6.vec()		    *.
# *								    *.
# *********************************************************************.
g6.grad.fun<-function(delta, nu)
{
  # name : g6.grad.fun
  # H0:theta=1 vs. H1:theta=delta.
  # g6.theta : log-likelihood function for simple vs. simple hypothesis.
  # to calculate 
  # g.grad : gradiant of g6
  # date :11/04/01
  # partial derivative is calculate with maple
  # reference : document/thesis/maple-simple.txt.
  #	g <- (x11 + x12) * log((y11 + y12)/(y11 + delta^2 * y12)) + (x21 + x22
  #		) * log((y21 + y22)/(y21 + delta^2 * y22)) + 2 * (x12 + x22) * 
  #		log(delta)
  theta <- delta
  x11 <- nu[1]
  x12 <- nu[2]
  x21 <- nu[3]
  x22 <- nu[4]
  y11 <- nu[5]
  y12 <- nu[6]
  y21 <- nu[7]
  y22 <- nu[8]
  dx11 <- log((y11 + y12)/(y11 + theta^2 * y12))
  dx12 <- log((y11 + y12)/(y11 + theta^2 * y12)) + 2 * log(theta)
  dx21 <- log((y21 + y22)/(y21 + theta^2 * y22))
  dx22 <- log((y21 + y22)/(y21 + theta^2 * y22)) + 2 * log(theta)
  dy11 <- ((x11 + x12) * (1/(y11 + theta^2 * y12) - (y11 + y12)/(y11 + 
                                                                   theta^2 * y12)^2))/(y11 + y12) * (y11 + theta^2 * y12)
  dy12 <- ((x11 + x12) * (1/(y11 + theta^2 * y12) - (y11 + y12)/(y11 + 
                                                                   theta^2 * y12)^2 * theta^2))/(y11 + y12) * (y11 + theta^2 * y12)
  dy21 <- ((x21 + x22) * (1/(y21 + theta^2 * y22) - (y21 + y22)/(y21 + 
                                                                   theta^2 * y22)^2))/(y21 + y22) * (y21 + theta^2 * y22)
  dy22 <- ((x21 + x22) * (1/(y21 + theta^2 * y22) - (y21 + y22)/(y21 + 
                                                                   theta^2 * y22)^2 * theta^2))/(y21 + y22) * (y21 + theta^2 * y22)
  g.grad <- as.vector(c(dx11, dx12, dx21, dx22, dy11, dy12, dy21, dy22))
  return(g.grad)
}

# *********************************************************************.
# *								    *.
# *								    *.
# *		g6.grad2.fun()					    *.
# *		          	    *.
# *								    *.
# *********************************************************************.

g6.grad2.fun <- function(delta, nu)
{
  # ================================================================
  # g6.grad2.fun(delta, nu)
  # Returns the 8x8 Hessian matrix of second partial derivatives
  # of g(delta, v) with respect to the 8 variables:
  #   v = c(K11, K12, K21, K22, T11, T12, T21, T22)
  # Evaluated at the point nu (exactly as needed for the slow term
  # in Ra6.fun and the dissertation equation (3.12))
  #
  # This completes the approximation:
  #   log-ratio ≈ Sn + Xin + (1/(2n)) * wn' * g.grad2 * wn
  # ================================================================
  
  x11 <- nu[1]; x12 <- nu[2]
  x21 <- nu[3]; x22 <- nu[4]
  y11 <- nu[5]; y12 <- nu[6]
  y21 <- nu[7]; y22 <- nu[8]
  
  H <- matrix(0, nrow = 8, ncol = 8)
  
  # ====================== SITE 1 (variables 1,2,5,6) ======================
  S1 <- y11 + y12
  C1 <- y11 + delta^2 * y12
  A1 <- x11 + x12
  
  if (S1 > 0 && C1 > 0) {
    # Cross terms: ∂²g / ∂K ∂T  (same for K11 and K12)
    dlog_dy11 <- 1/S1 - 1/C1          # = (C1 - S1)/(S1 * C1)
    dlog_dy12 <- 1/S1 - delta^2 / C1
    
    H[1,5] <- dlog_dy11;  H[5,1] <- dlog_dy11
    H[2,5] <- dlog_dy11;  H[5,2] <- dlog_dy11
    H[1,6] <- dlog_dy12;  H[6,1] <- dlog_dy12
    H[2,6] <- dlog_dy12;  H[6,2] <- dlog_dy12
    
    # Pure T second derivatives: A1 * Hessian of log(S1/C1) wrt (y11,y12)
    h11 <- A1 * (1/C1^2 - 1/S1^2)
    h12 <- A1 * (delta^2 / C1^2 - 1/S1^2)
    h22 <- A1 * (delta^4 / C1^2 - 1/S1^2)
    
    H[5,5] <- h11
    H[5,6] <- h12;   H[6,5] <- h12
    H[6,6] <- h22
  }
  
  # ====================== SITE 2 (variables 3,4,7,8) ======================
  S2 <- y21 + y22
  C2 <- y21 + delta^2 * y22
  A2 <- x21 + x22
  
  if (S2 > 0 && C2 > 0) {
    dlog_dy21 <- 1/S2 - 1/C2
    dlog_dy22 <- 1/S2 - delta^2 / C2
    
    H[3,7] <- dlog_dy21;  H[7,3] <- dlog_dy21
    H[4,7] <- dlog_dy21;  H[7,4] <- dlog_dy21
    H[3,8] <- dlog_dy22;  H[8,3] <- dlog_dy22
    H[4,8] <- dlog_dy22;  H[8,4] <- dlog_dy22
    
    h33 <- A2 * (1/C2^2 - 1/S2^2)
    h34 <- A2 * (delta^2 / C2^2 - 1/S2^2)
    h44 <- A2 * (delta^4 / C2^2 - 1/S2^2)
    
    H[7,7] <- h33
    H[7,8] <- h34;   H[8,7] <- h34
    H[8,8] <- h44
  }
  
  # All other entries (cross-site, K-K, etc.) remain 0
  return(H)
}




# *********************************************************************.
# *								    *.
# *								    *.
# *	m6data.fun(); generating two-site, two-treatment data       *.
# *								    *.
# *
# *********************************************************************.

m6data.fun <- function(rate1, rate2, w1, w2, delta, 
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
  t1  <- 1 / delta
  t2  <- delta
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
                       delta_alter = 1.25,          # H1: theta = delta
                       delta_null = 1,
                       cut.num = 3000,        # number of steps for min(Zk)
                       mov.num = 30)          # window for stationary correction
{
  # 1. Generate one LONG dataset
  total_n <- cut.num + 2 * mov.num + 10
  data1 <- m6data.fun(rate1 = rate1, rate2 = rate2, 
                      w1 = w1, w2 = w2, delta = delta_null, 
                      row.num = total_n)
  
  # 2. Compute all three ingredients exactly as in the theory
  theta1 <- 1 / delta_null
  theta2 <- delta_null
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

# ========================================================
# NEW STANDALONE FUNCTION — gamma.fun()
# Computes gamma1 or gamma0 exactly as described in the article
# Input : M = numeric vector of simulated M values 
#         (one vector per hypothesis)
# ========================================================

gamma.fun <- function(M) {
  
  # M must be a vector from many replications of m6test.fun(), under either null or alter hypo
  
  Mplus <- pmax(M, 0)                     # M⁺ = max(M, 0)
  
  # Compute the expectation that appears in (3.17)
  expect_term <- mean(1 - exp(-Mplus))
  
  return(expect_term)
}

ab.fun <- function(gamma0,gamma1,alpha=0.05,beta=0.05){
  eb <- beta / (gamma0 * (1 - alpha))          
  b  <- -log(eb)
  denom = (gamma0*gamma1*eb*(alpha-1))+gamma1
  a = -log(alpha/denom)
  return(c(-b,a))
}

Ra6.two.fun <- function(rate1, rate2, w1, w2, grad_delta, data_delta, a_upper, b_lower, m = m_emp) {
  
  cross_upper <- numeric(m)
  cross_lower <- numeric(m)
  
  for (i in 1:m) {
    
    # Generate one full dataset under the true delta (same as your original call)
    data1 <- m6data.fun(rate1 = rate1, rate2 = rate2,
                        w1 = w1, w2 = w2, delta = data_delta,
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
      log.ratio <- n * g6.fun(grad_delta, vn / n)
      
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


# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

# separation of used and unused codes.

# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


Ra6.fun0<-function(a, rate1, rate2, w1, w2, delta, row.num, nu, g, g.grad, 
                   g.grad2)
{
  # name : Ra6.fun
  # date : 11/1/01
  # to calculate Ra=log.likelihood ratio - a
  # H0:theta=1 vs. H1:theta=delta
  #	theta1 <- 1/delta
  #	theta2 <- delta
  #	w11 <- w1 * theta1
  #	w12 <- w1 * theta2
  #	w21 <- w2 * theta1
  #	w22 <- w2 * theta2
  #	p1 <- rate1/(rate1 + rate2)
  #	p2 <- rate2/(rate1 + rate2)
  log.ratio <- 0
  n <- 0	#	nu <- c(p1/2, p1/2, p2/2, p2/2, p1/(2 * w11), p1/(2 * w12), p2/(2 * w21), p2/(2 * w22))
  #	g <- g6.fun(delta, nu)
  #	g.grad <- g6.grad.fun(delta, nu)
  #	g.grad2 <- g6.grad2.fun(delta, nu)
  site1 <- data.fun(rate1, rate2, w1, w2, delta, row.num)$site1
  site2 <- data.fun(rate1, rate2, w1, w2, delta, row.num)$site2
  site1.max <- site1[row.num, 1]
  site2.max <- site2[row.num, 1]
  site.max.min <- min(site1.max, site2.max)
  site <- rbind(site1, site2)
  ii <- order(site[, 1])	# to sort by arrival time
  site <- site[ii,  ]
  site.num <- sum(site[, 1] < site.max.min) + 1
  tou <- c(0, site[, 1])
  xk <- site[, 2]
  B11 <- site[, 4]
  B12 <- site[, 5]
  B21 <- site[, 6]
  B22 <- site[, 7]
  while((log.ratio < a) || (n == site.num)) {
    n <- n + 1
    xin.K11 <- (-1) * sum((xk[1:n] > (tou[n + 1] - tou[1:
                                                         n])) * B11[1:n])
    vvn.K11 <- sum(B11[1:n])
    K11 <- vvn.K11 + xin.K11
    xin.K12 <- (-1) * sum((xk[1:n] > (tou[n + 1] - tou[1:
                                                         n])) * B12[1:n])
    vvn.K12 <- sum(B12[1:n])
    K12 <- vvn.K12 + xin.K12
    xin.K21 <- (-1) * sum((xk[1:n] > (tou[n + 1] - tou[1:
                                                         n])) * B21[1:n])
    vvn.K21 <- sum(B21[1:n])
    K21 <- vvn.K21 + xin.K21
    xin.K22 <- (-1) * sum((xk[1:n] > (tou[n + 1] - tou[1:
                                                         n])) * B22[1:n])
    vvn.K22 <- sum(B22[1:n])
    K22 <- vvn.K22 + xin.K22
    pos.T11 <- (xk[1:n] > (tou[n + 1] - tou[1:n])) * B11[
      1:n]
    xin.T11 <- (-1) * sum((xk[1:n] - (tou[n + 1] - tou[1:
                                                         n])) * pos.T11)
    vvn.T11 <- sum(xk[1:n] * B11[1:n])
    T11 <- vvn.T11 + xin.T11
    pos.T12 <- (xk[1:n] > (tou[n + 1] - tou[1:n])) * B12[
      1:n]
    xin.T12 <- (-1) * sum((xk[1:n] - (tou[n + 1] - tou[1:
                                                         n])) * pos.T12)
    vvn.T12 <- sum(xk[1:n] * B12[1:n])
    T12 <- vvn.T12 + xin.T12
    pos.T21 <- (xk[1:n] > (tou[n + 1] - tou[1:n])) * B21[
      1:n]
    xin.T21 <- (-1) * sum((xk[1:n] - (tou[n + 1] - tou[1:
                                                         n])) * pos.T21)
    vvn.T21 <- sum(xk[1:n] * B21[1:n])
    T21 <- vvn.T21 + xin.T21
    pos.T22 <- (xk[1:n] > (tou[n + 1] - tou[1:n])) * B22[
      1:n]
    xin.T22 <- (-1) * sum((xk[1:n] - (tou[n + 1] - tou[1:
                                                         n])) * pos.T22)
    vvn.T22 <- sum(xk[1:n] * B22[1:n])
    T22 <- vvn.T22 + xin.T22
    vn <- c(K11, K12, K21, K22, T11, T12, T21, T22)
    vvn <- c(vvn.K11, vvn.K12, vvn.K21, vvn.K22, vvn.T11, 
             vvn.T12, vvn.T21, vvn.T22)
    xin <- c(xin.K11, xin.K12, xin.K21, xin.K22, xin.T11, 
             xin.T12, xin.T21, xin.T22)
    wn <- vn - n * nu
    if((T11 + T12) * (T21 + T22) != 0) {
      log.ratio <- g6.fun(delta, vn)
      Sn <- n * g + sum(g.grad * (vvn - n * nu))
      Xin <- sum(g.grad * xin)
      slow <- (1/(2 * n)) * (t(wn) %*% g.grad2 %*% 
                               wn)
      log.ratio.app <- Sn + Xin + slow
    }
  }
  Ra <- log.ratio - a
  Ra.g <- c(Ra, n, Sn, Xin, slow)
  return(Ra.g, vn)
}



# *********************************************************************.
# *								    *.
# *								    *.
# *		Ra6.vec()					    *.
# *								    *.
# *								    *.
# *********************************************************************.
Ra6.vec<-function(a, rate1, rate2, w1, w2, delta, row.num, m)
{
  #function name:Ra6.vec
  #revised:11/04/01 
  #To make vectors of N = inf{n:log.ratio > a}, Ra. log.ratio, 
  #approximation part including random walk, stationary, slowlyching, and Kij, Tij
  #underlying distribution : exponential (xi_1*theta_j)
  #Patients arrive at site 1 with rate1, at site 2 with rate2 which are known.
  #
  cat(" ********", "\n")
  cat("Starting time:", date(), "\n")
  Ra.vec <- matrix(0, nrow = m, ncol = 5)
  Vn <- matrix(0, nrow = m, ncol = 8)
  theta1 <- 1/delta
  theta2 <- delta
  w11 <- w1 * theta1
  w12 <- w1 * theta2
  w21 <- w2 * theta1
  w22 <- w2 * theta2
  p1 <- rate1/(rate1 + rate2)
  p2 <- rate2/(rate1 + rate2)
  nu <- c(p1/2, p1/2, p2/2, p2/2, p1/(2 * w11), p1/(2 * w12), 
          p2/(2 * w21), p2/(2 * w22))
  g <- g6.fun(delta, nu)
  g.grad <- g6.grad.fun(delta, nu)
  g.grad2 <- g6.grad2.fun(delta, nu)
  for(i in 1:m) {
    Ra.vn <- Ra6.fun(a, rate1, rate2, w1, w2, delta, 
                     row.num, nu, g, g.grad, g.grad2)
    Ra.vec[i,  ] <- Ra.vn$Ra.g
    Vn[i,  ] <- Ra.vn$vn
    if(i/100 == floor(i/100)) {
      cat(i, " ")
    }
  }
  Ra <- Ra.vec[, 1]
  N <- Ra.vec[, 2]
  Sn <- Ra.vec[, 3]
  Xin <- Ra.vec[, 4]
  slow <- Ra.vec[, 5]
  cat("\n", "Ending time:", date(), "\n")
  cat("****** delta = ", delta, " ****** a = ", a, 
      "************", "\n")
  # cat("Mean of Ra :", round(mean(Ra), 4), "  Std.error of Ra :",
  #     round(se.fun(Ra), 4), "\n")
  # cat("Mean of e^(-Ra) :", round(mean(exp( - Ra)), 4), 
  #     "  Std.error of Ra :", round(se.fun(exp( - Ra)), 4), 
  #     "\n")
  # cat("mean of N : ", round(mean(N), 4), "  Std.error of N : ", 
  #     round(se.fun(N), 4), "\n")
  # cat("Mean of Sn :", round(mean(Sn), 4), 
  #     "  Std.error of Sn  :", round(se.fun(Sn), 4), "\n")
  # cat("Mean of Xin", round(mean(Xin), 4), 
  #     "  Std.error of Xin :", round(se.fun(Xin), 4), "\n")
  # cat("Mean of slow :", round(mean(slow), 4), 
  #     "  Std.error of slow :", round(se.fun(slow), 4), "\n"
  # )
  return(list(Ra,N,Sn,Xin,slow,Vn))
}


# *********************************************************************.
# *								    *.
# *								    *.
# *		mRa1.fun()					    *.
# *		#cutnum, movnum!				    *.
#   *								    *.
# *********************************************************************.
mRa1.fun<-function(delta = 1.25, lamda = 1, 
                   cut.num = 5000, mov.num = 40)
{
  # mRa1.fun
  # date : 9/10/01
  # To generate M for limiting distribution of Ra when arrival time is exponential with rate lamda
  yk <- rexp(cut.num + mov.num, 
             delta)
  inter.arrival <- c(0, rexp(
    cut.num + mov.num, lamda
  ))
  tou <- cumsum(inter.arrival)	
  #** tou is patient's arrival time.
  xk <- log(delta) - (delta - 1) * 
    yk[1:cut.num]
  yp0 <- yk[1:mov.num] + tou[2:(
    mov.num + 1)] - rep(tou[
      1], mov.num)
  yp0.pos <- yp0 > 0
  xi0 <- (delta - 1) * sum(yp0[
    yp0.pos]) - log(delta) * 
    sum(yp0.pos)
  zk <- cumsum(xk) + xi0	
  ###	cat("\n", "yk :", round(yk, 4))
  ###	cat("\n", "xk :", round(xk, 4))
  ###	cat("\n", "tou :", round(tou, 4))
  ###	cat("\n", "yp0 :", round(yp0, 4))
  ###	cat("\n", "xi0 :", round(xi0, 4))
  ###	cat("\n", "zk :", round(zk, 4))
  xi <- rep(0, cut.num)
  for(i in 1:cut.num) {
    yp <- yk[(i + 1):(
      mov.num + i)] + 
      tou[(i + 2):(
        mov.num + i + 1)
      ] - rep(tou[i + 
                    1], mov.num)
    yp.pos <- yp > 0
    xi[i] <- (delta - 1) * 
      sum(yp[yp.pos]) - 
      log(delta) * sum(
        yp.pos)	
    ###		cat("\n", "yp :", round(yp, 4))
    ###		cat("\n", "xi :", round(xi[i], 4))
  }
  zk <- zk - xi
  return(c(min(zk), xi0))
}

# *********************************************************************.
# *								    *.
# *								    *.
# *	m6data.fun(); generating two-site, two-treatment data       *.
# *								    *.
# *	STILL UNDER CONSTRUCTION				    *.
# * IT IS AN OLDER VERSION
# *********************************************************************.
m6data1.fun<-function(rate1, rate2, w1, w2, delta, row.num, p=0.5, r=0.5)
{
  #m6data.fun()
  #2/9/02
  #generates data for 2 site- 2 treatment situation;
  #calling function: m6test.func();
  #this is a new version; see m6data0.fun() for old one.
  #p=probability that treatment 1 is chosen (regardless of site)
  #r=probability that site 1 is chosen
  
  t1<-1/delta
  t2<-delta
  w11<-w1*t1 
  w12<-w1*t2
  w21<-w2*t1
  w22<-w2*t2
  
  ### framework of site data with row.num cases ;
  site<-data.frame(matrix(0,nrow=row.num, ncol=3))
  dimnames(site)[[2]]<-c("life","arrival","status")
  site$status[1:row.num]<-10 #site information; 1x=site 1; 2x=site 2;
  site2.ind<-as.numeric(runif(n=row.num)>r)
  site$status[1:row.num]<-site$status[1:row.num]+
    
    ### p is the treatment 1 probability  (same for both sites)
    # assignment of treatment for site 1
    treat2<-as.numeric(runif(n=row.num)>p) #random assignment of treatment2
  site$status[1:row.num]<-site$status[1:row.num]+treat2
  
  # assignment of treatment for site 2
  treat2<-as.numeric(runif(n=row.num)>p)
  site$status[(row.num+1):(2*row.num)]<-site$status[(row.num+1):(2*row.num)]+treat2
  
  ### arrival time assignment for each site
  arr<-cumsum(rexp(row.num,rate1))
  site$arrival[1:row.num]<-arr
  arr<-cumsum(rexp(row.num,rate2))
  site$arrival[(row.num+1):(2*row.num)]<-arr
  
  ### lifetime assignment for each site/treatment combination with different rates
  c11<-site$status==11
  site$life[c11]<-rexp(sum(c11),w11)
  c12<-site$status==12
  site$life[c12]<-rexp(sum(c12),w12)
  c21<-site$status==21
  site$life[c21]<-rexp(sum(c21),w21)
  c22<-site$status==22
  site$life[c22]<-rexp(sum(c22),w22)
  
  ### sorting the data in the arrival order
  ii<-order(site$arrival)
  site<-site[ii,]
  return(site)
}


# *********************************************************************.
# *								    *.
# *								    *.
# *	m6data0.fun(); generating 2 site, 2 treatment data	    *.
# *								    *.
# *								    *.
# *********************************************************************.
m6data0.fun<-function(rate1, rate2, w1, w2, delta, row.num, p=0.5)
{
  #m6data0.fun()
  #2/7/02
  #generates data for 2 site- 2 treatment situation;
  #calling function: m6test.func();
  #modification of m6data.fun()
  t1<-1/delta
  t2<-delta
  w11<-w1*t1 
  w12<-w1*t2
  w21<-w2*t1
  w22<-w2*t2
  
  ### framework of site data with first row.num cases from site 1 and the next from site 2;
  site<-data.frame(matrix(0,nrow=2*row.num, ncol=3))
  dimnames(site)[[2]]<-c("life","arrival","status")
  site$status[1:row.num]<-11 #site information; 1x=site 1; 2x=site 2;
  site$status[(row.num+1):(2*row.num)]<-21
  
  ### p is the treatment 1 probability  (same for both sites)
  # assignment of treatment for site 1
  treat2<-as.numeric(runif(n=row.num)>p) #random assignment of treatment2
  site$status[1:row.num]<-site$status[1:row.num]+treat2
  
  # assignment of treatment for site 2
  treat2<-as.numeric(runif(n=row.num)>p)
  site$status[(row.num+1):(2*row.num)]<-site$status[(row.num+1):(2*row.num)]+treat2
  
  ### arrival time assignment for each site
  arr<-cumsum(rexp(row.num,rate1))
  site$arrival[1:row.num]<-arr
  arr<-cumsum(rexp(row.num,rate2))
  site$arrival[(row.num+1):(2*row.num)]<-arr
  
  ### lifetime assignment for each site/treatment combination with different rates
  c11<-site$status==11
  site$life[c11]<-rexp(sum(c11),w11)
  c12<-site$status==12
  site$life[c12]<-rexp(sum(c12),w12)
  c21<-site$status==21
  site$life[c21]<-rexp(sum(c21),w21)
  c22<-site$status==22
  site$life[c22]<-rexp(sum(c22),w22)
  
  ### sorting the data in the arrival order
  ii<-order(site$arrival)
  site<-site[ii,]
  return(site[1:row.num,])
}

m6test0.fun <- function(rate1 = 1, rate2 = 1, 
                        w1 = 1, w2 = 1, 
                        delta = 1.25,          # H1: theta = delta
                        cut.num = 5000,        # number of steps for min(Zk)
                        mov.num = 30)          # window for stationary correction
{
  # 1. Generate one LONG dataset
  total_n <- cut.num + 2 * mov.num + 10
  data1 <- m6data.fun(rate1 = rate1, rate2 = rate2, 
                      w1 = w1, w2 = w2, delta = delta, 
                      row.num = total_n)
  
  # 2. Compute all three ingredients exactly as in the theory
  theta1 <- 1 / delta
  theta2 <- delta
  w11 <- w1 * theta1;  w12 <- w1 * theta2
  w21 <- w2 * theta1;  w22 <- w2 * theta2
  
  p1 <- rate1 / (rate1 + rate2)
  p2 <- rate2 / (rate1 + rate2)
  
  nu      <- c(p1/2, p1/2, p2/2, p2/2,
               p1/(2*w11), p1/(2*w12),
               p2/(2*w21), p2/(2*w22))
  
  g       <- g6.fun(delta, nu)
  g.grad  <- g6.grad.fun(delta, nu)
  g.grad2 <- g6.grad2.fun(delta, nu)     # ← NEW: the 8x8 Hessian
  
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

# *********************************************************************.
# *								    *.
# *								    *.
# *		m6test.g1.fun()					    *.
# *								    *.
# *								    *.
# *********************************************************************.

m6test.g1.fun <- function(rate1 = 1, rate2 = 1, 
                          w1 = 1, w2 = 1, 
                          delta = 1.25,          # H1: theta = delta
                          cut.num = 5000,        # number of steps for the min
                          mov.num = 30)          # window size for stationary approximation (dissertation uses ~30)
{
  # ================================================================
  # This function only involves first order gradient.
  # m6test.g1.fun() — Correct implementation of the limiting M 
  # (Chapter III, Section 3.3 of the dissertation)
  # Calls: m6data.fun(), g6.fun(), g6.grad.fun()
  # Returns M = inf_k Z_k  (used to approximate the limiting excess distribution)
  # ================================================================
  
  # 1. Generate one LONG dataset (enough patients for cut.num + windows)
  total_n <- cut.num + 2 * mov.num + 10   # safety margin
  data1 <- m6data.fun(rate1 = rate1, rate2 = rate2, 
                      w1 = w1, w2 = w2, delta = delta, 
                      row.num = total_n)
  
  # 2. Compute nu, g, g.grad exactly as in the theory (eq. 3.9 and 3.4)
  theta1 <- 1 / delta
  theta2 <- delta
  w11 <- w1 * theta1;  w12 <- w1 * theta2
  w21 <- w2 * theta1;  w22 <- w2 * theta2
  
  p1 <- rate1 / (rate1 + rate2)
  p2 <- rate2 / (rate1 + rate2)
  
  nu <- c(p1/2, p1/2, p2/2, p2/2, # check nu definition, nu is a specific parameter
          p1/(2*w11), p1/(2*w12),
          p2/(2*w21), p2/(2*w22))
  
  g     <- g6.fun(delta, nu)
  g.grad <- g6.grad.fun(delta, nu)   # 8-element gradient
  
  # 3. Extract data (same format used everywhere else)
  xk     <- data1$life
  tou    <- c(0, data1$arrival)      # tou[1] = 0
  status <- data1$status
  B11 <- as.numeric(status == 11) # B_{ij k} = indicator that kth patient is site i, treatment j
  B12 <- as.numeric(status == 12)
  B21 <- as.numeric(status == 21)
  B22 <- as.numeric(status == 22)
  
  # 4. Compute xi0 (stationary correction at "time 0") — using first mov.num patients
  #    (this is xi0 in the dissertation notation 3.8)
  xi0 <- numeric(8)
  # K components
  posK <- (xk[1:mov.num] > (tou[2:(mov.num+1)] - tou[1])) 
  xi0[1] <- -sum(posK * B11[1:mov.num])
  xi0[2] <- -sum(posK * B12[1:mov.num])
  xi0[3] <- -sum(posK * B21[1:mov.num])
  xi0[4] <- -sum(posK * B22[1:mov.num])
  # T components
  posT11 <- posK * B11[1:mov.num]
  xi0[5] <- -sum( (xk[1:mov.num] - (tou[2:(mov.num+1)] - tou[1])) * posT11 )
  posT12 <- posK * B12[1:mov.num]
  xi0[6] <- -sum( (xk[1:mov.num] - (tou[2:(mov.num+1)] - tou[1])) * posT12 )
  posT21 <- posK * B21[1:mov.num]
  xi0[7] <- -sum( (xk[1:mov.num] - (tou[2:(mov.num+1)] - tou[1])) * posT21 )
  posT22 <- posK * B22[1:mov.num]
  xi0[8] <- -sum( (xk[1:mov.num] - (tou[2:(mov.num+1)] - tou[1])) * posT22 )
  
  Xi0 <- sum(g.grad * xi0)
  
  # 5. Loop over cut.num steps — compute Sn, Xin, Zk exactly as in dissertation
  Sn  <- numeric(cut.num)
  Xin <- numeric(cut.num)
  Zk  <- numeric(cut.num)
  
  for(n in 1:cut.num) {
    # Current window for xin (patients n+1 to n+mov.num)
    idx <- (n+1):(n + mov.num)
    pos <- (xk[idx] > (tou[idx+1] - tou[n+1])) 
    
    xin <- numeric(8)
    # K components
    xin[1] <- -sum(pos * B11[idx])
    xin[2] <- -sum(pos * B12[idx])
    xin[3] <- -sum(pos * B21[idx])
    xin[4] <- -sum(pos * B22[idx])
    # T components
    xin[5] <- -sum( (xk[idx] - (tou[idx+1] - tou[n+1])) * (pos * B11[idx]) )
    xin[6] <- -sum( (xk[idx] - (tou[idx+1] - tou[n+1])) * (pos * B12[idx]) )
    xin[7] <- -sum( (xk[idx] - (tou[idx+1] - tou[n+1])) * (pos * B21[idx]) )
    xin[8] <- -sum( (xk[idx] - (tou[idx+1] - tou[n+1])) * (pos * B22[idx]) )
    
    # Cumulative counts up to n (V_n^* in dissertation 3.6)
    vvn <- c(sum(B11[1:n]), sum(B12[1:n]), sum(B21[1:n]), sum(B22[1:n]),
             sum(xk[1:n] * B11[1:n]), sum(xk[1:n] * B12[1:n]),
             sum(xk[1:n] * B21[1:n]), sum(xk[1:n] * B22[1:n]))
    
    # Random walk part
    Sn[n]  <- n * g + sum(g.grad * (vvn - n * nu)) # dissertation 3.11
    Xin[n] <- sum(g.grad * xin)
    
    # Centered process Z_n (dissertation 3.13)
    Zk[n]  <- Sn[n] + Xi0 - Xin[n]
  }
  
  M <- min(Zk)
  
  cat("m6test.fun finished — M =", round(M, 4), "\n")
  return(list(M = M, Xi0 = Xi0, Sn = Sn, Xin = Xin, Zk = Zk))
}

m6test0.fun<-function(data = data1, rate1 = 1, rate2 = 1, w1 = 1, w2 = 1, delta = 3, cut.num = 15, mov.num
                      = 5, nu = 2, g = 4, g.grad = 2.5)
{
  # name : m6test.fun
  # date : 2/7/02
  # to calculate limiting distribution of Ra
  # preliminary testing function for developing m6test.func(), a C implementation of m6.fun().
  # H0:theta=1 vs. H1:theta=delta
  #	theta1 <- 1/delta
  #	theta2 <- delta
  #	w11 <- w1 * theta1
  #	w12 <- w1 * theta2
  #	w21 <- w2 * theta1
  #	w22 <- w2 * theta2
  #	p1 <- rate1/(rate1 + rate2)
  #	p2 <- rate2/(rate1 + rate2)
  cut.num <- 20
  mov.num <- 6
  Sn <- rep(0, cut.num)
  Xin <- rep(0, cut.num)
  data1 <- data
  xk <- data1$life
  tou <- c(0, data1$arrival)
  status <- data1$status
  B11 <- as.numeric(status == 11)
  B12 <- as.numeric(status == 12)
  B21 <- as.numeric(status == 21)
  B22 <- as.numeric(status == 22)
  nu <- 2
  delta <- 3
  g.grad <- c(1, 3, 2, 5, 4, 3, 3, 2)
  g <- 4
  xi0.K11 <- (-1) * sum((xk[1:mov.num] > (tou[2:(mov.num + 1)] - tou[1])) * B11[1:mov.num
  ])
  xi0.K12 <- (-1) * sum((xk[1:mov.num] > (tou[2:(mov.num + 1)] - tou[1])) * B12[1:mov.num
  ])
  xi0.K21 <- (-1) * sum((xk[1:mov.num] > (tou[2:(mov.num + 1)] - tou[1])) * B21[1:mov.num
  ])
  xi0.K22 <- (-1) * sum((xk[1:mov.num] > (tou[2:(mov.num + 1)] - tou[1])) * B22[1:mov.num
  ])
  pos.T11 <- (xk[1:mov.num] > (tou[2:(mov.num + 1)] - tou[1])) * B11[1:mov.num]
  xi0.T11 <- (-1) * sum((xk[1:mov.num] - (tou[2:(mov.num + 1)] - tou[1])) * pos.T11)
  pos.T12 <- (xk[1:mov.num] > (tou[2:(mov.num + 1)] - tou[1])) * B12[1:mov.num]
  xi0.T12 <- (-1) * sum((xk[1:mov.num] - (tou[2:(mov.num + 1)] - tou[1])) * pos.T12)
  pos.T21 <- (xk[1:mov.num] > (tou[2:(mov.num + 1)] - tou[1])) * B21[1:mov.num]
  xi0.T21 <- (-1) * sum((xk[1:mov.num] - (tou[2:(mov.num + 1)] - tou[1])) * pos.T21)
  pos.T22 <- (xk[1:mov.num] > (tou[2:(mov.num + 1)] - tou[1])) * B22[1:mov.num]
  xi0.T22 <- (-1) * sum((xk[1:mov.num] - (tou[2:(mov.num + 1)] - tou[1])) * pos.T22)
  xi0 <- c(xi0.K11, xi0.K12, xi0.K21, xi0.K22, xi0.T11, xi0.T12, xi0.T21, xi0.T22)
  Xi0 <- sum(g.grad * xi0)
  
  # Main loop: compute S_n, Xi_n, Z_n for n = 1,..., cut.num
  
  for(n in 1:cut.num) {
    xin.K11 <- (-1) * sum((xk[(n + 1):(mov.num + n)] > (tou[(n + 2):(mov.num + n + 
                                                                       1)] - tou[n + 1])) * B11[(n + 1):(mov.num + n)])
    vvn.K11 <- sum(B11[1:n])
    xin.K12 <- (-1) * sum((xk[(n + 1):(mov.num + n)] > (tou[(n + 2):(mov.num + n + 
                                                                       1)] - tou[n + 1])) * B12[(n + 1):(mov.num + n)])
    vvn.K12 <- sum(B12[1:n])
    xin.K21 <- (-1) * sum((xk[(n + 1):(mov.num + n)] > (tou[(n + 2):(mov.num + n + 
                                                                       1)] - tou[n + 1])) * B21[(n + 1):(mov.num + n)])
    vvn.K21 <- sum(B21[1:n])
    xin.K22 <- (-1) * sum((xk[(n + 1):(mov.num + n)] > (tou[(n + 2):(mov.num + n + 
                                                                       1)] - tou[n + 1])) * B22[(n + 1):(mov.num + n)])
    vvn.K22 <- sum(B22[1:n])
    pos.T11 <- (xk[(n + 1):(mov.num + n)] > (tou[(n + 2):(mov.num + n + 1)] - tou[n +
                                                                                    1])) * B11[(n + 1):(mov.num + n)]
    xin.T11 <- (-1) * sum((xk[(n + 1):(mov.num + n)] - (tou[(n + 2):(mov.num + n + 
                                                                       1)] - tou[n + 1])) * pos.T11)
    pos.T12 <- (xk[(n + 1):(mov.num + n)] > (tou[(n + 2):(mov.num + n + 1)] - tou[n +
                                                                                    1])) * B12[(n + 1):(mov.num + n)]
    xin.T12 <- (-1) * sum((xk[(n + 1):(mov.num + n)] - (tou[(n + 2):(mov.num + n + 
                                                                       1)] - tou[n + 1])) * pos.T12)
    pos.T21 <- (xk[(n + 1):(mov.num + n)] > (tou[(n + 2):(mov.num + n + 1)] - tou[n +
                                                                                    1])) * B21[(n + 1):(mov.num + n)]
    xin.T21 <- (-1) * sum((xk[(n + 1):(mov.num + n)] - (tou[(n + 2):(mov.num + n + 
                                                                       1)] - tou[n + 1])) * pos.T21)
    pos.T22 <- (xk[(n + 1):(mov.num + n)] > (tou[(n + 2):(mov.num + n + 1)] - tou[n +
                                                                                    1])) * B22[(n + 1):(mov.num + n)]
    xin.T22 <- (-1) * sum((xk[(n + 1):(mov.num + n)] - (tou[(n + 2):(mov.num + n + 
                                                                       1)] - tou[n + 1])) * pos.T22)
    vvn.T11 <- sum(xk[1:n] * B11[1:n])
    vvn.T12 <- sum(xk[1:n] * B12[1:n])
    vvn.T21 <- sum(xk[1:n] * B21[1:n])
    vvn.T22 <- sum(xk[1:n] * B22[1:n])
    cat("*******************************\n")
    xin <- c(xin.K11, xin.K12, xin.K21, xin.K22, xin.T11, xin.T12, xin.T21, xin.T22
    )
    vvn <- c(vvn.K11, vvn.K12, vvn.K21, vvn.K22, vvn.T11, vvn.T12, vvn.T21, vvn.T22
    )
    cat(xin, "\n")
    cat(vvn, "\n")
    Sn[n] <- n * g + sum(g.grad * (vvn - n * nu))
    Xin[n] <- sum(g.grad * xin)
  }
  Zk <- Sn + Xi0 - Xin
  M <- min(Zk)
  return(M, Xi0, Sn, Xin, Zk)
}

# *********************************************************************.
# *								    *.
# *								    *.
# *		Ra6.fun()					    *.
# *								    *.
# *								    *.
# *********************************************************************.

Ra6.fun <- function(a, rate1, rate2, w1, w2, delta, row.num, 
                    nu, g, g.grad, g.grad2)
{
  # ================================================================
  # Ra6.fun() — Fixed version using m6data.fun()
  # Generates data with m6data.fun() and runs the sequential monitoring
  # until log-ratio >= a (one-sided)
  # ================================================================
  
  # 1. Generate the data using the correct function
  data1 <- m6data.fun(rate1 = rate1, rate2 = rate2, 
                      w1 = w1, w2 = w2, delta = delta, 
                      row.num = row.num)
  
  # 2. Extract the necessary vectors (exactly as used in the rest of your code)
  xk     <- data1$life
  tou    <- c(0, data1$arrival)          # tou[1] = 0
  status <- data1$status
  
  B11 <- as.numeric(status == 11)
  B12 <- as.numeric(status == 12)
  B21 <- as.numeric(status == 21)
  B22 <- as.numeric(status == 22)
  
  site.num <- row.num                    # total number of patients
  
  # 3. Sequential monitoring loop
  log.ratio <- 0
  n <- 0
  
  while (log.ratio < a && n < site.num) {   # corrected condition (&& instead of ||)
    n <- n + 1
    
    # --- Compute Kij and Tij up to time n (same logic as original) ---
    xin.K11 <- (-1) * sum((xk[1:n] > (tou[n+1] - tou[1:n])) * B11[1:n])
    vvn.K11 <- sum(B11[1:n])
    K11 <- vvn.K11 + xin.K11
    
    xin.K12 <- (-1) * sum((xk[1:n] > (tou[n+1] - tou[1:n])) * B12[1:n])
    vvn.K12 <- sum(B12[1:n])
    K12 <- vvn.K12 + xin.K12
    
    xin.K21 <- (-1) * sum((xk[1:n] > (tou[n+1] - tou[1:n])) * B21[1:n])
    vvn.K21 <- sum(B21[1:n])
    K21 <- vvn.K21 + xin.K21
    
    xin.K22 <- (-1) * sum((xk[1:n] > (tou[n+1] - tou[1:n])) * B22[1:n])
    vvn.K22 <- sum(B22[1:n])
    K22 <- vvn.K22 + xin.K22
    
    # Tij
    pos.T11 <- (xk[1:n] > (tou[n+1] - tou[1:n])) * B11[1:n]
    xin.T11 <- (-1) * sum((xk[1:n] - (tou[n+1] - tou[1:n])) * pos.T11)
    vvn.T11 <- sum(xk[1:n] * B11[1:n])
    T11 <- vvn.T11 + xin.T11
    
    pos.T12 <- (xk[1:n] > (tou[n+1] - tou[1:n])) * B12[1:n]
    xin.T12 <- (-1) * sum((xk[1:n] - (tou[n+1] - tou[1:n])) * pos.T12)
    vvn.T12 <- sum(xk[1:n] * B12[1:n])
    T12 <- vvn.T12 + xin.T12
    
    pos.T21 <- (xk[1:n] > (tou[n+1] - tou[1:n])) * B21[1:n]
    xin.T21 <- (-1) * sum((xk[1:n] - (tou[n+1] - tou[1:n])) * pos.T21)
    vvn.T21 <- sum(xk[1:n] * B21[1:n])
    T21 <- vvn.T21 + xin.T21
    
    pos.T22 <- (xk[1:n] > (tou[n+1] - tou[1:n])) * B22[1:n]
    xin.T22 <- (-1) * sum((xk[1:n] - (tou[n+1] - tou[1:n])) * pos.T22)
    vvn.T22 <- sum(xk[1:n] * B22[1:n])
    T22 <- vvn.T22 + xin.T22
    
    vn  <- c(K11, K12, K21, K22, T11, T12, T21, T22)
    vvn <- c(vvn.K11, vvn.K12, vvn.K21, vvn.K22, vvn.T11, vvn.T12, vvn.T21, vvn.T22)
    xin <- c(xin.K11, xin.K12, xin.K21, xin.K22, xin.T11, xin.T12, xin.T21, xin.T22)
    
    wn <- vn - n * nu
    
    if ((T11 + T12) * (T21 + T22) != 0) {
      log.ratio <- g6.fun(delta, vn)
      Sn   <- n * g + sum(g.grad * (vvn - n * nu))
      Xin  <- sum(g.grad * xin)
      slow <- (1 / (2 * n)) * (t(wn) %*% g.grad2 %*% wn)
      # log.ratio.app <- Sn + Xin + slow   # optional, kept for debugging
    }
  }
  
  Ra <- log.ratio - a
  Ra.g <- c(Ra, n, Sn, Xin, slow)
  
  return(list(Ra.g = Ra.g, vn = vn))
}
