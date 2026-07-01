# Plot for table xxx

data <- read.csv("Table_3.3_boundaries_rate11_w11.csv")
data_005 <- data[data$beta == 0.05, ]
data_010 <- data[data$beta == 0.1, ]

# -----------------------------
# Plot for beta = 0.05
# -----------------------------
png("boundaries_beta_0.05.png", width = 800, height = 600)  # optional: save as PNG
plot(data_005$delta, data_005$upper_a,
     type = "o", pch = 16, col = "blue", lwd = 2,
     xlab = expression(Delta), 
     ylab = "Boundary value",
     main = expression(paste("Sequential boundaries for ", beta, " = 0.05")),
     ylim = range(c(data_005$lower_b, data_005$upper_a)),
     cex.main = 1.2)

lines(data_005$delta, data_005$lower_b,
      type = "o", pch = 16, col = "red", lwd = 2)

abline(h = 0, lty = 2, col = "gray")

legend("topright", 
       legend = c("upper_a", "lower_b"),
       col = c("blue", "red"), 
       pch = 16, lty = 1, lwd = 2,
       bty = "n")

dev.off()   # close the PNG device (remove this line if you want the plot in R console)

# -----------------------------
# Plot for beta = 0.1
# -----------------------------
png("boundaries_beta_0.1.png", width = 800, height = 600)  # optional: save as PNG
plot(data_010$delta, data_010$upper_a,
     type = "o", pch = 16, col = "blue", lwd = 2,
     xlab = expression(Delta), 
     ylab = "Boundary value",
     main = expression(paste("Sequential boundaries for ", beta, " = 0.1")),
     ylim = range(c(data_010$lower_b, data_010$upper_a)),
     cex.main = 1.2)

lines(data_010$delta, data_010$lower_b,
      type = "o", pch = 16, col = "red", lwd = 2)

abline(h = 0, lty = 2, col = "gray")

legend("topright", 
       legend = c("upper_a", "lower_b"),
       col = c("blue", "red"), 
       pch = 16, lty = 1, lwd = 2,
       bty = "n")

dev.off()   # close the PNG device




### Plot for table xxx

# Open a single PNG with 2x2 layout
png("emp_errors_1x2_scatter.png", width = 1200, height = 500, res = 120)

par(mfrow = c(1, 2),          # 2 rows × 2 columns
    mar = c(4, 4, 3, 1),      # tighter margins
    oma = c(0, 0, 2, 0))      # outer margin for overall title

# # =============================================
# # Row 1: Nominal beta = 0.05
# # =============================================
# 
# # 1. emp_alpha (beta = 0.05)
# plot(data_005$delta, data_005$emp_alpha,
#      type = "p", pch = 16, col = "blue", lwd = 2, cex = 1.2,
#      xlab = "delta", 
#      ylab = "emp_alpha",
#      main = expression(paste("emp_alpha | nominal ", beta, " = 0.05")),
#      ylim = c(0, 0.12),
#      cex.main = 1.1)
# abline(h = 0.05, lty = 2, col = "red", lwd = 2)
# text(max(data_005$delta) * 0.95, 0.055, "nominal α = 0.05", 
#      col = "red", adj = c(1, 0), cex = 0.9)
# 
# # 2. emp_beta (beta = 0.05)
# plot(data_005$delta, data_005$emp_beta,
#      type = "p", pch = 16, col = "red", lwd = 2, cex = 1.2,
#      xlab = "delta", 
#      ylab = "emp_beta",
#      main = expression(paste("emp_beta | nominal ", beta, " = 0.05")),
#      ylim = c(0, 0.12),
#      cex.main = 1.1)
# abline(h = 0.05, lty = 2, col = "red", lwd = 2)
# text(max(data_005$delta) * 0.95, 0.055, "nominal β = 0.05", 
#      col = "red", adj = c(1, 0), cex = 0.9)

# =============================================
# Row 2: Nominal beta = 0.1
# =============================================

# 3. emp_alpha (beta = 0.1)
plot(data_010$delta, data_010$emp_alpha,
     type = "p", pch = 16, col = "blue", lwd = 2, cex = 1.2,
     xlab = "delta", 
     ylab = "emp_alpha",
     main = expression(paste("emp_alpha | nominal ", beta, " = 0.1")),
     ylim = c(0, 0.15),
     cex.main = 1.1)
abline(h = 0.05, lty = 2, col = "red", lwd = 2)
text(max(data_010$delta) * 0.95, 0.055, "nominal α = 0.05", 
     col = "red", adj = c(1, 0), cex = 0.9)

# 4. emp_beta (beta = 0.1)
plot(data_010$delta, data_010$emp_beta,
     type = "p", pch = 16, col = "red", lwd = 2, cex = 1.2,
     xlab = "delta", 
     ylab = "emp_beta",
     main = expression(paste("emp_beta | nominal ", beta, " = 0.1")),
     ylim = c(0, 0.15),
     cex.main = 1.1)
abline(h = 0.10, lty = 2, col = "red", lwd = 2)
text(max(data_010$delta) * 0.95, 0.105, "nominal β = 0.1", 
     col = "red", adj = c(1, 0), cex = 0.9)

# Overall title
mtext("Empirical Type I & Type II Error Rates vs. delta", 
      outer = TRUE, cex = 1.3, font = 2, line = 0.5)

dev.off()

cat("Plot saved as 'emp_errors_1x2_scatter.png'\n")

### Plot for table xxx

data_010 <- read.csv("Table_3.3_boundaries_rate21_lambda11.csv")

# Open a single PNG with 2x2 layout
png("emp_errors_1*2_ratestudy.png", width = 1200, height = 900, res = 120)

par(mfrow = c(1, 2),          # 2 rows × 2 columns
    mar = c(4, 4, 3, 1),      # tighter margins
    oma = c(0, 0, 2, 0))      # outer margin for overall title

plot(data_010$delta, data_010$emp_alpha,
     type = "p", pch = 16, col = "blue", lwd = 2, cex = 1.2,
     xlab = "delta", 
     ylab = "emp_alpha",
     main = expression(paste("emp_alpha | nominal ", beta, " = 0.1")),
     ylim = c(0, 0.15),
     cex.main = 1.1)
abline(h = 0.05, lty = 2, col = "red", lwd = 2)
text(max(data_010$delta) * 0.95, 0.055, "nominal α = 0.05", 
     col = "red", adj = c(1, 0), cex = 0.9)

plot(data_010$delta, data_010$emp_beta,
     type = "p", pch = 16, col = "red", lwd = 2, cex = 1.2,
     xlab = "delta", 
     ylab = "emp_beta",
     main = expression(paste("emp_beta | nominal ", beta, " = 0.1")),
     ylim = c(0, 0.15),
     cex.main = 1.1)
abline(h = 0.10, lty = 2, col = "red", lwd = 2)
text(max(data_010$delta) * 0.95, 0.105, "nominal β = 0.1", 
     col = "red", adj = c(1, 0), cex = 0.9)

# Overall title
mtext("Empirical Error Rates vs. delta when rate1:rate2 = 2:1", 
      outer = TRUE, cex = 1.3, font = 2, line = 0.5)

dev.off()

cat("Plot saved as 'emp_errors_1*2_ratestudy.png'\n")



