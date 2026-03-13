skewness_win <- function(Y,N) 
{
  n <- length(Y)
  g3 <- skewness(Y[1:N], na = TRUE)
  for (j in 1:(n-N))
  {
    g3[j] <- skewness(Y[(1+j):(1+j+N)], na = TRUE)
  }
  return(g3)
}