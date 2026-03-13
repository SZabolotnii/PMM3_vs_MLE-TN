a_pmm3ex <- function(X) 
{
  A.mean <- mean(X)
  res <- X - A.mean
  m2 <- mean(res^2); m4 <- mean(res^4); m6 <- mean(res^6)
  S1 <- mean(X); S2 <- mean(X^2); S3 <- mean(X^3); 
  A=1;   B=-S1*3;   C=3*S2-(m6-3*m4*m2)/(m4-3*m2^2);   D=S1*(m6-3*m4*m2)/(m4-3*m2^2)-S3;
  A.pmm3 <- A.mean
  del.A <- sqrt(m2)
  while (del.A > sqrt(m2)/10^4)
  {
    del.A <- (A*A.pmm3^3+B*A.pmm3^2+C*A.pmm3+D)/(3*A*A.pmm3^2+2*B*A.pmm3+C)
    A.pmm3 <- A.pmm3 - del.A
  }
  return(A.pmm3)
}