# функция нахождения решения системы уравнений ММПл при S=2
solvMRs2 <- function(B.ols,X,Y,m2,m3,m4) 
{
  P <- length(B.ols)
  # k2 <- m2; g3 <- m3/m2^(3/2); g4 <- (m4 - 3*m2^2)/m2^2
  # A <- g3; B <- -((2*g3*Y-k2^(1/2)*(g4+2))); C <- ((g3*Y^2-k2^(1/2)*(g4+2)*Y-k2*g3));
  A <- m3; B <- m4 - m2^2 - 2*m3*Y; C <- m3*Y^2 - (m4 - m2^2)*Y - m2*m3;
  B.pmm2 <- B.ols 
  JZs <- matrix(data = NA, nrow = P, ncol=P)

  for (i in 1:1000)
  {  
      Yx <- B.pmm2[1]*X[,1]  
      for (r in 2:P)
         {Yx <- Yx + B.pmm2[r]*X[,r]}
      Z1 <- A*Yx^2 + B*Yx + C
      Z <- sum(Z1)
      for (r in 2:P)
         {
        Zr <- Z1*X[,r]
        Z <- rbind(Z, sum(Zr))
         }
      JZ11 <- 2*A*Yx + B
      JZs[1,1] <- sum(JZ11)
      
      for (ii in 2:P)
      {
        JZ1i <- JZ11*X[,ii]
        JZs[1,ii] <- sum(JZ1i);
        JZs[ii,1] <- JZs[1,ii];
          
        for (jj in 2:P)
            {
             JZij=JZ11*X[,ii]*X[,jj]
             JZs[ii,jj]=sum(JZij)
            }
      }
     
      Q <- solve(JZs,Z)
      B.pmm2 <- B.pmm2 - Q
      
      if (norm(matrix(Q)) < 10^-5) break      
      
  }
  return(B.pmm2)
}