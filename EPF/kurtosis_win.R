kurtosis_win <- function(Y,N) 
{
  n <- length(Y)
  # g4 <- kurtosis(Y[1:N], na = TRUE) - 3
  m1 <- mean(Y[1:N], na = TRUE)
  m2 <- mean((Y[1:N]-m1)^2, na = TRUE)
  m4 <- mean((Y[1:N]-m1)^4, na = TRUE)
  g4 <- m4/m2^2 - 3
  for (j in 1:(n-N))
  {
    m1[j] <- mean(Y[(1+j):(1+j+N)], na = TRUE)
    m2[j] <- mean((Y[(1+j):(1+j+N)]-m1[j])^2, na = TRUE)
    m4[j] <- mean((Y[(1+j):(1+j+N)]-m1[j])^4, na = TRUE)
    g4[j] <- m4[j]/m2[j]^2 - 3
    # g4[j] <- kurtosis(Y[(1+j):(1+j+N)], na.rm = TRUE) - 3
  }
  return(g4)
}