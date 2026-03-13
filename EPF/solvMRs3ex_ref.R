solvMRs3ex_ref <- function(B.ols, X, Y, m2, m4, m6, tol = 1e-4, max.iter = 100) {
  P <- length(B.ols)
  B.pmm3 <- B.ols
  JZs <- matrix(NA, nrow = P, ncol = P)
  
  # Pre-computed constants
  A <- 1
  B_coef <- -3 * Y
  C_coef <- 3 * Y^2 - (m6 - 3 * m4 * m2) / (m4 - 3 * m2^2)
  D_coef <- Y * (m6 - 3 * m4 * m2) / (m4 - 3 * m2^2) - Y^3
  
  for (i in 1:max.iter) {
    Yx <- X %*% B.pmm3
    Z1 <- A * Yx^3 + B_coef * Yx^2 + C_coef * Yx + D_coef
    Z_matrix <- matrix(Z1, nrow = length(Z1), ncol = ncol(X), byrow = TRUE)
    Z <- colSums(Z_matrix * X)  # Element-wise multiplication
    
    # Compute the Jacobian
    JZ11 <- 3 * A * Yx^2 + 2 * B_coef * Yx + C_coef
    JZs[1, 1] <- sum(JZ11)
    for (ii in 2:P) {
      JZ1i <- JZ11 * X[, ii]
      JZs[1, ii] <- sum(JZ1i)
      JZs[ii, 1] <- JZs[1, ii]
      for (jj in ii:P) {
        JZs[ii, jj] <- sum(JZ1i * X[, jj])
        if (ii != jj) {
          JZs[jj, ii] <- JZs[ii, jj]  # Exploit the symmetry
        }
      }
    }
    
    # Solve the system of linear equations for the update
    Q <- solve(JZs, Z)
    B.pmm3 <- B.pmm3 - Q
    
    # Check for convergence
    if (norm(Q, type = "2") < tol) {
      break
    }
  }
  
  # Check for a significant difference between the estimates
  if (norm(B.ols - B.pmm3, "2") > 10) {
    B.pmm3 <- B.ols
  }
  
  return(B.pmm3)
}
