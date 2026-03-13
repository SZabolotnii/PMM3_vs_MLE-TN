step_pol <- function(B,X) 
{
  L <- length(B)
  Y <- B[1]
  for (i in 2:L)
  {
    Y <- Y + B[i]*X^(i-1)
  }
  return(Y)
}